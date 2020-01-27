import os
import sys
import random
import numpy
import math

wkdir = '/Users/apple/Data/multilayerNetwork/data/regrOutput/scExprs/'

subtype = sys.argv[1]
l = sys.argv[2]

# subtype = 'Mesenchymal'


# For each gene in the target list, find its TF and miRNA regulators
# Trim away the objfunc value at the end of the regression result vector
# Any F value within +/- 1e-4 of 1 will be considered as 1 (non-regulatory)



tfpath = wkdir + "/" + subtype + "/TF_expr.txt"
f = open(tfpath,'r')
TFExpr = f.readlines()
f.close()
TFExpr = [line.strip().split('\t') for line in TFExpr]
TFList = [line[0] for line in TFExpr[1:len(TFExpr)]]

#mipath = wkdir + "/" + subtype + "/mi_expr.txt"
#f = open(mipath,'r')
#miExpr = f.readlines()
#f.close()
#miExpr = [line.strip().split('\t') for line in miExpr]
#miList = [line[0] for line in miExpr[1:len(miExpr)]]

tgpath = wkdir + "/" + "targetGeneList.txt"
f = open(tgpath,'r')
TGList = f.readlines()
f.close()
TGList = [line.strip() for line in TGList]

#tmpath = wkdir + "/" + "targetmiRNAList.txt"
#f = open(tmpath,'r')
#TMList = f.readlines()
#f.close()
#TMList = [line.strip() for line in TMList]

TF_to_gene = {}
TF_to_mi = {}
mi_to_gene = {}

tol = 1e-3

for gene in TGList:
    gpath = wkdir + "/" + subtype + "/genes/" + gene
    tpath = gpath + "/TFs.txt"
    f = open(tpath,'r')
    tfs = f.readlines()
    tfs = [int(line.strip()) for line in tfs]
    f.close()
    mpath = gpath + "/miRNAs.txt"
    f = open(mpath,'r')
    mis = f.readlines()
    mis = [int(line.strip().split('\t')[0]) for line in mis]
    f.close()
    ltf = len(tfs)
    #lmi = len(mis)
    lmi = 0
    rpath = gpath + "/res_ri_nomi_" + str(l) + ".txt"
    f = open(rpath,'r')
    res = f.readlines()
    res = [float(line.strip().split('\t')[0]) for line in res]
    f.close()
    lres = len(res)
    if lres > ltf + lmi + 1:
        res = res[0:(lres-1)]
        lres -= 1
    cweight = res[0]
    tweights = res[1:(ltf+1)]
    mweights = res[(ltf+1):lres]
    if 'const' not in TF_to_gene:
        TF_to_gene['const'] = {gene:cweight}
    else:
        TF_to_gene['const'][gene] = cweight
    for i in range(ltf):
        currtf = TFList[tfs[i]-1]
        if currtf not in TF_to_gene:
            if abs(tweights[i]-1) > tol:
                TF_to_gene[currtf] = {gene:tweights[i]}
            else:
                TF_to_gene[currtf] = {gene:-1}
        else:
            if abs(tweights[i]-1) > tol:
                TF_to_gene[currtf][gene] = tweights[i]
            else:
                TF_to_gene[currtf][gene] = -1
    if lmi > 0:
        for i in range(lmi):
            currmi = miList[mis[i]-1]
            if currmi not in mi_to_gene:
                if abs(mweights[i]-1) > tol:
                    mi_to_gene[currmi] = {gene:mweights[i]}
                else:
                    mi_to_gene[currmi] = {gene:-1}
            else:
                if abs(mweights[i]-1) > tol:
                    mi_to_gene[currmi][gene] = mweights[i]
                else:
                    mi_to_gene[currmi][gene] = -1


for miRNA in TMList:
    mpath = wkdir + "/" + subtype + "/miRNAs/" + miRNA
    tpath = mpath + "/TFs.txt"
    f = open(tpath,'r')
    tfs = f.readlines()
    tfs = [int(line.strip()) for line in tfs]
    f.close()
    ltf = len(tfs)
    rpath = mpath + "/res_ri_nomi_" + str(l) + ".txt"
    f = open(rpath,'r')
    res = f.readlines()
    res = [float(line.strip().split('\t')[0]) for line in res]
    f.close()
    lres = len(res)
    if lres > ltf + 1:
        res = res[0:(lres-1)]
        lres -= 1
    cweight = res[0]
    tweights = res[1:(ltf+1)]
    if 'const' not in TF_to_mi:
        TF_to_mi['const'] = {miRNA:cweight}
    else:
        TF_to_mi['const'][miRNA] = cweight
    for i in range(ltf):
        currtf = TFList[tfs[i]-1]
        if currtf not in TF_to_mi:
            if abs(tweights[i]-1) > tol:
                TF_to_mi[currtf] = {miRNA:tweights[i]}
            else:
                TF_to_mi[currtf] = {miRNA:-1}
        else:
            if abs(tweights[i]-1) > tol:
                TF_to_mi[currtf][miRNA] = tweights[i]
            else:
                TF_to_mi[currtf][miRNA] = -1

# Output edges with non-(-1) weights
allpath = wkdir + "/" + subtype + "/" + subtype + "_all_ri_nomi_" + str(l) + ".txt"
f_all = open(allpath,'w')

tgopath = wkdir + "/" + subtype + "/" + subtype + "_TF_gene_ri_nomi_" + str(l) + ".txt"
f = open(tgopath,'w')
for TF in TF_to_gene:
    for gene in TF_to_gene[TF]:
        w = TF_to_gene[TF][gene]
        if w != -1:
            line = TF + "_t" + "\t" + gene + "_g" + "\t" + str(w) + "\n"
            f.write(line)
            f_all.write(line)

f.close()

mgopath = wkdir + "/" + subtype + "/" + subtype + "_mi_gene_ri_nomi_" + str(l) + ".txt"
f = open(mgopath,'w')
flag = 0
for mi in mi_to_gene:
    for gene in mi_to_gene[mi]:
        w = mi_to_gene[mi][gene]
        if w != -1:
            flag = 1
            line = mi + "_m" + "\t" + gene + "_g" + "\t" + str(w) + "\n"
            f.write(line)
            f_all.write(line)

if flag == 0:
	f.write("none\tnone\tnone\n")

f.close()

tmopath = wkdir + "/" + subtype + "/" + subtype + "_TF_mi_ri_nomi_" + str(l) + ".txt"
f = open(tmopath,'w')
for TF in TF_to_mi:
    for mi in TF_to_mi[TF]:
        w = TF_to_mi[TF][mi]
        if w != -1:
            line = TF + "_t" + "\t" + mi + "_m" + "\t" + str(w) + "\n"
            f.write(line)
            f_all.write(line)

f.close()

f_all.close()



# Generate randomly shuffled network
# Random shuffling preserves total number of non-zero edges among potential
# regulators of a target???

#nrepeat = 500
#
#for i in range(nrepeat):
#    TF_to_gene_r = {}
#    TF_to_mi_r = {}
#    mi_to_gene_r = {}






















