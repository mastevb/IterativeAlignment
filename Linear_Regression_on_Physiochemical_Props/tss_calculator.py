import pandas as pd
import numpy as np
import random


class TSS_Calculator(object):

    def __init__(self, p, b, TSS_path=None, seq_as_string=True):
        self.peptides = []
        self.binders = []
        self.length = 0
        self.AA = list("ACDEFGHIKLMNPQRSTVWY")
        self.cluster_keys = list("abcehipr")
        self.matrices = dict.fromkeys(self.cluster_keys)
        self.total_scores = None
        self.dist = None
        for key in self.cluster_keys:
            self.matrices[key] = pd.read_csv("./dist_matrices/cluster_" + key + ".csv", index_col=0)
        if TSS_path is not None:
            self.import_TSS(TSS_path)
        self.__store_values(p, b, seq_as_string)
        if p is not None and b is not None:
            self.binder_scores = self.get_binder_scores()

    def __store_values(self, p=None, b=None, seq_as_string=True):
        print('storing')
        if p is not None:
            self.p = p
            if "Unnamed: 0" in list(self.p.columns):
                self.p = self.p.set_index("Unnamed: 0")
            if not seq_as_string: self.length = len(list(self.p.columns))
            if seq_as_string:
                self.peptides = list(self.p['AA_seq'])
                self.length = len(self.peptides[0])
            else:
                self.peptides = [''.join(list(self.p.iloc[m, :]))
                                 for m in range(len(list(self.p.index)))]
        self.pep_num = len(self.peptides)
        if b is not None:
            self.b = b
            if "Unnamed: 0" in list(self.b.columns):
                self.b = self.b.set_index("Unnamed: 0")
            self.length = len(list(self.b.columns))
            if seq_as_string:
                self.binders = list(self.b['AA_seq'])
                self.length = len(self.peptides[0])
            else:
                self.binders = [''.join(list(self.b.iloc[m, :]))
                                for m in range(len(list(self.b.index)))]
        self.bin_num = len(self.binders)

    def get_binder_scores(self):
        print('gettingbinderscores')
        self.dist = dict.fromkeys(self.cluster_keys)
        self.total_scores = dict.fromkeys(self.cluster_keys)
        for key in self.dist.keys():
            self.dist[key] = dict.fromkeys(self.AA)
            self.total_scores[key] = dict.fromkeys(self.AA)
            for AA1 in self.AA:
                self.dist[key][AA1] = dict.fromkeys(self.AA)
                self.total_scores[key][AA1] = dict.fromkeys(range(self.length))
                for AA2 in self.AA:
                    self.dist[key][AA1][AA2] = self.matrices[key][AA1][AA2]
                for l in range(self.length):
                    self.total_scores[key][AA1][l] = sum(self.dist[key][AA1] \
                                                             [self.binders[n][l]] for n in range(self.bin_num))
        return self.total_scores

    def calculate_TSS(self):
        print("Calculating TSS...")
        np_ss = np.zeros(shape=(self.pep_num, len(self.cluster_keys)))
        for m in range(self.pep_num):
            if self.peptides[m] in self.binders:
                binder_index = self.binders.index(self.peptides[m])
                for i, key in enumerate(self.cluster_keys):
                    np_ss[m][i] = -sum(self.matrices[key][self.peptides[m][l]]
                                       [self.binders[binder_index][l]]
                                       for l in range(self.length))
            for i, key in enumerate(self.cluster_keys):
                np_ss[m][i] += sum(self.binder_scores[key][self.peptides[m][l]][l] \
                                   for l in range(self.length))
        self.TSS_df = pd.DataFrame(np_ss, index=self.peptides, columns=self.cluster_keys)
        return self.TSS_df

    def generate_test_set(self, num_peps, num_binders, l):
        self.peptides = [''.join([random.choice(self.AA) for j in range(l)]) for j in range(num_peps)]
        self.binders = [''.join([random.choice(self.AA) for j in range(l)]) for j in range(num_binders)]
        self.length = l
        self.__store_values()


'''
COMB = pd.read_csv('data/combined.csv')
COMB = COMB[COMB['Sequence'] != 0]
COMB = COMB[COMB.Conserved_Affinity>0]
COMB = COMB.set_index(COMB.Sequence)

BIND = COMB[COMB['Conserved_Affinity'] > 17]
PEPS = COMB[COMB['Conserved_Affinity'] <= 17]

TAR = pd.DataFrame(COMB.Conserved_Affinity)
TAR = TAR.set_index(COMB.index)

p1 = PEPS
b1 = BIND
p1 = pd.DataFrame(p1.Sequence.apply(lambda x: list(x)).tolist())
b1 = pd.DataFrame(b1.Sequence.apply(lambda x: list(x)).tolist())
p1_train, p1_test = train_test_split(p1, test_size=0.2, random_state=0)
ac = AlignmentCalculator(p1_train, b1, seq_as_string=False)
tss = ac.calculate_TSS()
tss.to_csv("FULL_SET_TSS.csv")
'''
