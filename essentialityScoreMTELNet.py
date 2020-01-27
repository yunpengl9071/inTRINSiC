#import pandas as pd
import sklearn.linear_model as sklm
import numpy as np
#from matplotlib import pyplot as plt
import sys
import os
import pandas as pd
import random

np.random.seed(77922766)

path = '/Users/yunpengl/data/multilayerNetwork/data/perturbation/'

files = []
for r, d, f in os.walk(path):
    for file in f:
        if 'avail_scaled_' in file:
            files.append(os.path.join(r, file))
    #res = np.array(res)

ccleScoresFile = '/Users/yunpengl/data/multilayerNetwork/data/CCLE_DepMap_cns_scores_avail_TV.txt'
ccleScores = pd.read_table(ccleScoresFile,index_col=0).transpose()
# f = open(ccleScoresFile,'r')
# res = f.readlines()
# f.close()
# res = res[1:]
# res = [line.strip().split('\t') for line in res]
allTFs = list(ccleScores.columns)
# res = [[(float(y)) for y in x[1:]] for x in res]
# res = np.array(res)
# ccleScores = res.transpose()

resList = []
nGenes = 0
#for file in files:
#    TF = file.split('perturbation/avail_scaled_')[1].split('_')[0]
#    if TF in allTFs:
##        res = f.readlines()
##        res = res[1:]
#        res = [line.strip().split('\t') for line in res]
#        res = [[(float(y)) for y in x[1:]] for x in res]
#        print TF
#        resList.append(res)
TFList = []

for TF in allTFs:
    file = '/Users/yunpengl/data/multilayerNetwork/data/perturbation/avail_unscaled_delta_' + TF + '_lambda_0_0.1_RBE_CCLE_deltaPRScores.perturbed_act_inh.adjAdj0.1.txt'
    res_pd = pd.read_table(file,index_col=0)
    TF_file = '/Users/yunpengl/data/multilayerNetwork/data/perturbation/' + TF + 'availCovariates.txt'
    f = open(TF_file,'r')
    res = f.readlines()
    f.close()
    TFvec = [line.strip().split('\t')[0] for line in res]
    TFList.append(TFvec)
    #print(TFvec)
    #f = open(file,'r')
    #res = f.readlines()
    #f.close()
    #res = res[1:]
    #res = [line.strip().split('\t') for line in res]
    #res = [[(float(y)) for y in x[1:]] for x in res]
    print TF
    resList.append(res_pd)

#print(len(TFList))
# Try building a separate model for each TF
regrResList = []
for i in range(0,len(allTFs)):
    TF = allTFs[i]
    print(TF)
    clf = sklm.ElasticNetCV(cv = 4,max_iter=10000,verbose=0)
    clf.fit((resList[i].transpose())[TFList[i]].values,ccleScores[TF])
    elnet = sklm.ElasticNet(max_iter=10000,alpha=clf.alpha_)
    elnet.fit((resList[i].transpose())[TFList[i]].values,ccleScores[TF])
    paramVec = [elnet.intercept_] + elnet.coef_.tolist()
    paramVec = [[TF] + [str(x) for x in paramVec]]
    regrResList.append(paramVec)

outfile = '/Users/yunpengl/data/multilayerNetwork/data/perturbation/deltaScores_elasticNetCVRes_unscaled_NZCovariates.txt'
f = open(outfile,'w')

for TFRes in regrResList:
    line = '\t'.join(TFRes[0]) + '\n'
    f.write(line)

f.close()
#subtype = sys.argv[1]
#l = sys.argv[2]
#fn = '/ahg/regevdata/projects/txnRegModeling/regression/computeTFCor_lin/' + subtype + #'_TFGeneCapacities_ri_nomi__lambda' + l + '.txt'
#f = open(fn,'r')
#res = f.readlines()
#f.close()
