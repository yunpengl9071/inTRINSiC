import os
import sys
import random
#import numpy
import math

wkdir = '/ahg/regevdata/projects/txnRegModeling/regression'

subtype = sys.argv[1]
l = sys.argv[2]
rmd = sys.argv[3]

if int(rmd) > 0:
    idx = sys.argv[4]

# subtype = 'Mesenchymal'


# For each gene in the target list, find its TF and miRNA regulators
# Trim away the objfunc value at the end of the regression result vector
# Any F value within +/- 1e-4 of 1 will be considered as 1 (non-regulatory)


if int(rmd) == 1:
    tfpath = wkdir + "/" + subtype + "/CV_" + idx + "/TF_expr.txt"
elif int(rmd) == 2:
    tfpath = wkdir + "/" + subtype + "/randShfl_" + idx + "/TF_expr.txt"
else:
    tfpath = wkdir + "/" + subtype + "/TF_expr.txt"

tgpath = wkdir + "/" + "targetGeneList.txt"

f = open(tfpath,'r')
TFExpr = f.readlines()
f.close()
TFExpr = [line.strip().split('\t') for line in TFExpr]
TFList = [line[0] for line in TFExpr[1:len(TFExpr)]]

f = open(tgpath,'r')
TGList = f.readlines()
f.close()
TGList = [line.strip() for line in TGList]

TF_to_gene = {}

tol = 1e-3

for gene in TGList:
    if int(rmd) == 1:
        gpath = wkdir + "/" + subtype + "/CV_" + idx + "/genes/" + gene
    elif int(rmd) == 2:
        gpath = wkdir + "/" + subtype + "/randShfl_" + idx + "/genes/" + gene
    else:
        gpath = wkdir + "/" + subtype + "/genes/" + gene
    tpath = gpath + "/TFs.txt"
    f = open(tpath,'r')
    tfs = f.readlines()
    tfs = [int(line.strip()) for line in tfs]
    f.close()
    ltf = len(tfs)
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

# Output edges with non-(-1) weights
if int(rmd) == 1:
    allpath = wkdir + "/" + subtype + "/CV_" + idx + "/" + subtype + "_all_ri_nomi_" + str(l) + ".txt"
elif int(rmd) == 2:
    allpath = wkdir + "/" + subtype + "/randShfl_" + idx + "/" + subtype + "_all_ri_nomi_" + str(l) + ".txt"
else:
    allpath = wkdir + "/" + subtype + "/" + subtype + "_all_ri_nomi_" + str(l) + ".txt"

f_all = open(allpath,'w')

if int(rmd) == 1:
    tgopath = wkdir + "/" + subtype + "/CV_" + idx + "/" + subtype + "_TF_gene_ri_nomi_" + str(l) + ".txt"
elif int(rmd) == 2:
    tgopath = wkdir + "/" + subtype + "/randShfl_" + idx + "/" + subtype + "_TF_gene_ri_nomi_" + str(l) + ".txt"
else:
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

f_all.close()





















