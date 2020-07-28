# This program filters bdata.20130222.mhci.txt
# 1. separate by allele
# 2. remove binding affinities that aren't exact

import csv
import pandas as pd
import math

root = '/Users/mabochen/Desktop/Research/IALR/Iterative_Alignment/bdata.20130222.mhci.txt'
data = pd.read_csv(root, sep='	',
                   header=None, names=['species', 'mhc', 'peptide_length', 'sequence', 'inequality', 'mean'])

# remove binding affinities that aren't exact
filtered = pd.DataFrame([])
for i, row in data.iterrows():
    print(i)
    if data[data.columns[4]][i] != "=":
        continue
    filtered = filtered.append(pd.DataFrame({
            data.columns.values[0]: data[data.columns[0]][i],
            data.columns.values[1]: data[data.columns[1]][i],
            data.columns.values[2]: data[data.columns[2]][i],
            data.columns.values[3]: data[data.columns[3]][i],
            data.columns.values[5]: -1 * math.log10(data[data.columns[5]][i])},
                                        index=[0]), ignore_index=True)

# temporarily save the filtered result
#filtered.to_csv("/Users/mabochen/Desktop/Research/IALR/filtering.csv")

filtered = pd.read_csv('/Users/mabochen/Desktop/Research/IALR/filtering.csv', sep=',',
                   header=None, names=['species', 'mhc', 'peptide_length', 'sequence', 'meas'], skiprows=1)

# separate by allele
gb = filtered.groupby(['species', 'mhc', 'peptide_length'])
for name, group in gb:
    fname = name[0] + '_' + name[1] + '_' + str(name[2]) + '.csv'  # name  for the file
    fname = fname.replace('/', '-')
    group.to_csv('/Users/mabochen/Desktop/Research/IALR/filtered_data/' + fname)