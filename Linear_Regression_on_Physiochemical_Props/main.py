# -*- coding: utf-8 -*-

# -------------------------------REQUIREMENTS---------------------------------#
# 1. This program requires separate csvs for strong binder set and full peptide 
#    set
# 2. This program requires only a single csv for the full set and strong 
#    binders (no separate sets)
# 3. This program requires peptide sequences are contained in single column as 
#    a string in each csv
# 4. This program requires that binding affinity data is present in each 
#    respective csv
# 5. If count discrimination is needed, this program requires columns titled 
#    'CP1', 'CP2', 'CP3', 'CE'
# 6. Requires pandas, numpy, scikit-learn
# -----------------------------------INFO-------------------------------------#
# 1. This program calculates TSS for a set of peptides (given binders) and 
#    outputs it
# 2. This program predicts binding affinities in several cross validations and 
#    outputs it
# 3. All outputted files are saved to same directory as working directory
# ----------------------------------IMPORTS-----------------------------------#
from tss_calculator import TSS_Calculator
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LinearRegression
import warnings

warnings.simplefilter(action='ignore', category=Warning)
# ---------------------------------CONSTANTS-----------------------------------#
pep = None
bind = None
count_discrim = None
highest_count = None
COUNTS = [0, 1, 2, 3, 4]
SEQ = "Sequence"
TARGET = "Binding_Affinity"
NUM_CV = 5
WEAK_UPPER = 8
MED_UPPER = 17
STRONG_UPPER = 20
FEATURES = list("abcehipr")
TEST = 0.2
PRED = "Predicted"
TSS = [None for count in COUNTS]
LIN_REG_MODELS = [[None for cv in range(NUM_CV)] for count in COUNTS]
DATA = [[None for cv in range(NUM_CV)] for count in COUNTS]
# -------------------------------LAMBDA FXNS-----------------------------------#
scaler = MinMaxScaler(feature_range=(0, 1), copy=True)
format_data = lambda data: data[[SEQ, TARGET]]
get_seq = lambda data: pd.DataFrame(data[SEQ].apply(lambda x: list(x)).tolist())
get_TSS = lambda p, b: TSS_Calculator(p, b).calculate_TSS()
out_TSS = lambda tss, c: tss.to_csv('TSS_C' + str(c) + '.csv')
discrim = lambda d, count: d[d['CE'] + d['CP3'] + d['CP2'] + d['CP1'] > count]
scale = lambda d: pd.DataFrame(scaler.fit_transform(d, [0, 1]), index=d.index, columns=d.columns)
scale_mult = lambda *a: (scale(d) for d in a)
only_left = lambda d1, d2: pd.merge(d1, d2, how='outer', indicator=True)


# --------------------------------FUNCTIONS------------------------------------#
def train_test(peps, bind, state):
    p_TRAIN, p_TEST = train_test_split(peps, test_size=TEST, random_state=state)
    bind_in_p_TRAIN = p_TRAIN.merge(bind)
    return p_TRAIN, p_TEST, bind_in_p_TRAIN


def get_model(data, state):
    train, test = train_test_split(data, test_size=TEST, random_state=state)
    X_test = scale(test[FEATURES])
    X_train = scale(train[FEATURES])
    X_full = scale(data[FEATURES])
    X_test, X_train, X_full = scale_mult(test[FEATURES], train[FEATURES], data[FEATURES])
    model = LinearRegression().fit(X_train, train[TARGET])
    return model


def get_linreg_models(peps, bind, count, seed):
    # p_n = format_data(discrim(peps, count))
    # p_n = format_data(peps, count)
    p_n = format_data(peps)
    b_n = format_data(bind)
    # print(p_n.head())
    print(b_n.head())
    p_TRAIN, p_TEST = train_test_split(p_n, test_size=0.2, train_size=0.8, random_state=seed)
    b_n, b_test = train_test_split(b_n, test_size=0.8, random_state=seed)
    TSS[count] = get_TSS(p_TRAIN, b_n)
    out_TSS(TSS[count], count)
    for j in range(5):
        print("Counts > " + str(count) + ": CV" + str(j + 1))
        data = TSS[count].merge(p_n, left_on=TSS[count].index, right_on=SEQ)
        data.set_index(SEQ)
        LIN_REG_MODELS[count][j] = get_model(data, j)


def predict_and_classify(count):
    tss = TSS[count]
    X = tss[FEATURES]
    scaler = MinMaxScaler(feature_range=(0, 1), copy=True)
    for j in range(5):
        model = LIN_REG_MODELS[count][j]
        XT = pd.DataFrame(scaler.fit_transform(X[FEATURES], [0, 1]), index=X.index, columns=FEATURES)
        pred = model.predict(XT[FEATURES])
        pred = np.interp(pred, (pred.min(), pred.max()), (0, 20))
        XT['pred'] = pred
        conditions = [(XT['pred'] >= 0) & (XT['pred'] < WEAK_UPPER),
                      (XT['pred'] >= 7) & (XT['pred'] < MED_UPPER),
                      (XT['pred'] <= STRONG_UPPER)]
        choices = ['Weak', 'Medium', 'Strong']
        XT['Class'] = np.select(conditions, choices, default='Weak')
        DATA[count][j] = XT
        XT.to_csv('model_count+' + str(count) + '_cv' + str(j) + '.csv')


# -------------------------------MAIN------------------------------------------#
if __name__ == "__main__":
    # step 1: input path of full peptide set:
    # pep = pd.read_csv(input("1. Path of peptides csv: "))
    pep = pd.read_csv("/Users/mabochen/Desktop/Research/IALR/" +
                      "Linear_Regression_on_Physiochemical_Props/chimpanzee_Patr-A*0101_9.csv")

    # step 2: input path of full strong binder set:
    # bind = pd.read_csv(input("2. Path of binders csv: "))
    bind = pd.read_csv("/Users/mabochen/Desktop/Research/IALR/" +
                       "Linear_Regression_on_Physiochemical_Props/StrongPeps.csv")

    # step 4: what is the name of the "Sequence" column in the csv?
    # SEQ = str(input("4. What is the name of the column containing peptide sequences?: "))
    SEQ = "sequence"

    # step 5: what is the name of the column containing binding affinity data?
    # TARGET = str(input("5. What is the name of column containing binding affinity data?: "))
    TARGET = "meas"

    # step 6: how many cross validations do you wish to perform?
    # NUM_CV = int(input("6. How many cross validations do you wish to perform?: "))
    NUM_CV = 5

    # step 7: what is upper binding affinity limit of "weak" binding peptides?
    # WEAK_UPPER = int(input("7: what is upper binding affinity limit of \"weak\" binding peptides?: "))
    WEAK_UPPER = -6

    # step 8: what is the upper binding affinity limit of "medium" binding peptides?
    # MED_UPPER = int(input("8: what is upper binding affinity limit of \"medium\" binding peptides?: "))
    MED_UPPER = -4

    # step 9: what is the upper binding affinity limit of strong binding peptides?
    # STRONG_UPPER = int(input("9: what is upper binding affinity limit of strong binding peptides?: "))
    STRONG_UPPER = -1

    TSS = [None for count in COUNTS]
    LIN_REG_MODELS = [[None for cv in range(NUM_CV)] for count in COUNTS]
    DATA = [[None for cv in range(NUM_CV)] for count in COUNTS]

    # perform calculations
    for count in COUNTS:
        get_linreg_models(pep, bind, count, count)
        predict_and_classify(count)
