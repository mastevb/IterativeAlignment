# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:24:21 2018

@author: jtfl2
"""
import random
import os
import copy
import easygui as eg
import re
#aa = 'A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V'.split(',')

filenames = []
eds = []
mat_filepath = "/Users/mabochen/Desktop/Research/20SU/Additional_Models_for_IEDB_Data/Iterative_Alignment/matrix"
pepsfolder = "/Users/mabochen/Desktop/Research/20SU/Additional_Models_for_IEDB_Data/Iterative_Alignment/peps"
saveFolder = "/Users/mabochen/Desktop/Research/20SU/Additional_Models_for_IEDB_Data/Iterative_Alignment/saves"
[ite, cut, dataselec, per] = eg.multenterbox(fields = ["Iterations", "Cutoff", "Dataset Selection", "Perturb Decimal (values greater than or equal to 1 will be counted as integer changes)"])
for file in os.listdir(mat_filepath):
    if file.endswith('.csv'):
        filenames.append(file.replace('.csv',''))
        matrix = []
        with open(mat_filepath + '/' + file, 'r') as g:
            string = g.read()
            lines = string.replace('b','').replace('\xef\xbb\xbf','').replace('ï»¿','').split('\n')
            aa = lines[0].replace(' ','').split(',')[1:21]
            for i in lines[1:len(lines)]:
                row = i.split(',')
                if len(row) >= 20:
                    matrix.append(row[1:21])
        newDict = {}
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                newDict[(aa[i], aa[j])] = float(re.sub('[^0-9],.','',matrix[i][j]))
        eds.append(newDict)

highpeps = []
lowpeps = []
pepnames = []
for data in os.listdir(pepsfolder):
    h = []
    l = []
    pepnames.append(data)
    with open(pepsfolder+ '/' + data+ '/high.csv', 'r') as g:
        string = g.readline()
        while string:
            h.append(''.join(string.replace('\n','').replace(',','').replace(' ','').replace('b','').replace('\xef\xbb\xbf','').replace('ï»¿','')))
            string = g.readline()
    with open(pepsfolder+ '/' + data+ '/low.csv', 'r') as g:
        string = g.readline()
        while string:
            l.append(''.join(string.replace('\n','').replace(',','').replace(' ','').replace('b','').replace('\xef\xbb\xbf','').replace('ï»¿','')))
            string = g.readline()
    highpeps.append(h)
    lowpeps.append(l)

def preprocess(peplist):
    blankDict = {}
    for i in aa:
        blankDict[i] = 0
    Dict_for_Location = []
    for i in range(len(peplist)): # why?????
        Dict_for_Location.append(copy.deepcopy(blankDict))
    for pep in peplist:
        for ind,i in enumerate(pep):
            Dict_for_Location[ind][i] = Dict_for_Location[ind][i] + 1
    return Dict_for_Location

def findTSS(pep, highDict, lowDict, ed):
    h = 0
    l = 0
    for ind, i in enumerate(pep):
        for a in aa:
            if i != a:
                h = h + highDict[ind][a]*ed[i,a]
            else:
                h = h + (highDict[ind][a]-1)*ed[i,a]
            l = l + lowDict[ind][a]*ed[i,a]
    return [h,l]

for ind, dataset in enumerate(pepnames):
    if not os.path.isdir(saveFolder + '/' + dataset):
        os.mkdir(saveFolder + '/' + dataset)
    for index in range(len(filenames)):
        filename = filenames[index]
        iterations = int(ite)
        cutoff = int(cut)
        old_unchanged = []
        highvalue = 0
        lowvalue = 0
        ed = eds[index]
        if float(per) >= 1:
            perturb = int(per)
        else:
            perturb = float(per) * max(ed.values())
        maxV = max(ed.values())
        minV = min(ed.values())
        high = highpeps[ind][1:int(dataselec)+1]
        low = lowpeps[ind][-int(dataselec):len(lowpeps[ind])]
        highBinding = preprocess(high)
        lowBinding = preprocess(low)
#        for d in range(len(highBinding)):
#            currentpep = highBinding[d]
#            for f in range(len(highBinding)):
#                if d != f:
#                    highvalue = highvalue + pairwise2.align.globalds(currentpep,highBinding[f],ed,-int(gap),-0.8)[0][2]
#            for f in range(len(lowBinding)):
#                lowvalue = lowvalue + pairwise2.align.globalds(currentpep,lowBinding[f],ed,-int(gap),-0.8)[0][2]
        for d in high:
            [hi, lo] = findTSS(d,highBinding, lowBinding, ed)
            highvalue = highvalue + hi
            lowvalue = lowvalue + lo
        iteration_highvalues = [highvalue]
        iteration_lowvalues  = [lowvalue]
        unchanged = 0
        if not os.path.isdir(saveFolder + '/' + dataset +"/" + filename):
            os.mkdir(saveFolder + '/' + dataset + "/" + filename)
        while len(iteration_highvalues) < iterations and unchanged < cutoff:
            highvalue = 0
            lowvalue = 0
            current = copy.deepcopy(ed)
            changesMade = 0
            while (changesMade == 0):
                for i in range(len(aa)):
                    j = i
                    k = random.randint(0,19)
                    sign = random.randint(0,2)
                    if sign == 1:
                        current_perturb = -perturb
                    elif sign == 2:
                        current_perturb = perturb
                    else:
                        current_perturb = 0
                    l = current[aa[j],aa[k]] + current_perturb
                    while (l<minV) or (l>maxV):
                        j = i
                        k = random.randint(0,19)
                        sign = random.randint(0,2)
                        if sign == 1:
                          current_perturb = -perturb
                        elif sign == 2:
                            current_perturb = perturb
                        else:
                            current_perturb = 0
                        l = current[aa[j],aa[k]] + current_perturb
                    if (current_perturb != 0):
                        changesMade = changesMade + 1
                    current[aa[j],aa[k]] = l
                    current[aa[k],aa[j]] = l
#            for d in range(len(highBinding)):
#                currentpep = highBinding[d];
#                for f in range(len(highBinding)):
#                    if d != f:
#                        highvalue = highvalue + pairwise2.align.globalds(currentpep,highBinding[f],current,-int(gap),-0.8)[0][2]
#                for f in range(len(lowBinding)):
#                     lowvalue = lowvalue + pairwise2.align.globalds(currentpep,lowBinding[f],current,-int(gap),-0.8)[0][2]
            for d in high:
                [hi, lo] = findTSS(d,highBinding, lowBinding, current)
                highvalue = highvalue + hi
                lowvalue = lowvalue + lo
            if highvalue > iteration_highvalues[-1] and lowvalue < iteration_lowvalues[-1]:
                iteration_highvalues.append(highvalue)
                iteration_lowvalues.append(lowvalue)
                old_unchanged.append(unchanged)
                unchanged = 0
                print("Improvements: " + str(len(iteration_highvalues)) + "                              ", end = '\r')
                with open(saveFolder + '/' + dataset +'/' + filename + '/new_ed'+str(len(iteration_highvalues))+'.csv', 'w+') as g:
                    for i in aa:
                        string = str(current[i,aa[0]])
                        for j in aa[1:len(aa)]:
                            string = string + ',' + str(current[i,j])
                        g.write(string + '\n')
                ed = copy.deepcopy(current)
            else:
                unchanged = unchanged + 1
            print('Improvements: ' + str(len(iteration_highvalues)) + ',' + ' unchaged: ' + str(unchanged), end = '\r')
        values =  [iteration_highvalues, iteration_lowvalues]
        with open(saveFolder+ '/'+ dataset +'/' + filename + '/values.csv', 'w') as g:
            g.write('HighValues, LowValues \n')
            for i in range(len(values[0])):
                g.write(str(values[0][i])+ ',' + str(values[1][i]) + '\n')
        with open(saveFolder+ '/' + dataset +'/' + filename + '/unchangedData.csv', 'w') as g:
            for i in old_unchanged:
                g.write(str(i) + '\n')
