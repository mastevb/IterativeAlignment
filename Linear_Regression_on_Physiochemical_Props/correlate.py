from scipy import stats
import pandas as pd

predicted = pd.read_csv("model_count+0_cv0.csv")
predicted_sorted = predicted.sort_values(by="sequence")
preds = predicted_sorted.pred.tolist()
original = pd.read_csv("chimpanzee_Patr-A*0101_9.csv")
original_sorted = original.sort_values(by="sequence")
meas = original_sorted.meas.tolist()
# calculate the Spearman correlation coefficient
print(len(meas))
# correlation, pvalue = stats.spearmanr(preds, meas)
