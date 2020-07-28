import csv
import math
import os
import pandas as pd

path = "/Users/mabochen/Desktop/Research/20SU/Additional_Models_for_IEDB_Data/Iterative_Alignment/"
root = "/Users/mabochen/Desktop/Research/20SU/Additional_Models_for_IEDB_Data/Iterative_Alignment/bdata.20130222.mhci.txt"
data = pd.read_csv(root, sep="	",
                   header=None, names=["species", "mhc", "peptide_length", "sequence", "inequality", "meas"])

# only consider "=" in the inequality column
# data before filtering
# HLA-A*02:01
# min_meas = float("inf")
# max_meas = float("-inf")
# for i, row in data.iterrows():
#     if data[data.columns[4]][i] == "=" and data[data.columns[2]][i] == 9 and data[data.columns[1]][i] == "HLA-A*02:01":
#         min_meas = min(min_meas, -1 * math.log10(data[data.columns[5]][i]))
#         max_meas = max(max_meas, -1 * math.log10(data[data.columns[5]][i]))
#
# lower_bound = min_meas + 0.2 * (max_meas - min_meas)
# upper_bound = min_meas + 0.9 * (max_meas - min_meas)
# low_result = []
# hi_result = []
result = pd.DataFrame([])
for i, row in data.iterrows():
    #if data[data.columns[4]][i] != "=" or data[data.columns[2]][i] != 9 or data[data.columns[1]][i] != "HLA-A*02:01":
    if data[data.columns[4]][i] != "=":
        continue
    result = result.append(pd.DataFrame({data.columns.values[3]: data[data.columns[3]][i],
                                         data.columns.values[5]: -1 * math.log10(data[data.columns[5]][i])},
                                        index=[0]), ignore_index=True)
    # s = -1 * math.log10(data[data.columns[5]][i])
    # seq = data[data.columns[3]][i]
    # if s >= upper_bound:
    #     hi_result.append([char for char in seq])
    # elif s <= lower_bound:
    #     low_result.append([char for char in seq])

file = open('filtered_bdata.20130222.csv', 'a+', newline='')
# writing the data into the file
with file:
    write = csv.writer(file)
    write.writerows(result)

result.set_index('sequence', inplace=True)
names = result.index.values.tolist()
sortedAffinities = result['meas']
temp = sortedAffinities.sort_values(ascending=False)
high = temp[0 : int(len(temp) * .05)].keys().tolist()
low = temp[int(len(temp) * 0.95) : int(len(temp))].keys().tolist()
low_result = []
for seq in low:
    low_result.append([char for char in seq])
high_result = []
for seq in high:
    high_result.append([char for char in seq])
# opening the csv file in 'a+' mode
file = open('peps/1/low.csv', 'a+', newline='')
# writing the data into the file
with file:
    write = csv.writer(file)
    write.writerows(low_result)
# opening the csv file in 'a+' mode
file = open('peps/1/high.csv', 'a+', newline='')
# writing the data into the file
with file:
    write = csv.writer(file)
    write.writerows(high_result)

# file = open('peps/1/low.csv', 'a+', newline='')
# with file:
#     write = csv.writer(file)
#     write.writerows(low_result)
#
# file = open('peps/1/high.csv', 'a+', newline='')
# with file:
#     write = csv.writer(file)
#     write.writerows(hi_result)

result.to_csv("whole.csv", sep='\t')