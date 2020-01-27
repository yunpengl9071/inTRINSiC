import pandas as pd
import sklearn.linear_model as sklm
import numpy as np
from matplotlib import pyplot as plt
import sys

subtypes = ['Classical', 'Neural', 'Proneural', 'Mesenchymal']
l = sys.argv[1]
sigTF = sys.argv[2]
sigTF_idx = int(sys.argv[3])
fn = sys.argv[4]
hl = sys.argv[5]
f = open(fn,'r')
idx = f.readlines()
f.close()
sigTgt_idx = [int(x.strip()) for x in idx]
cormat = np.zeros((len(subtypes),518))
count = 0
for subtype in subtypes:
    fn = '/users/yunpengl/Data/multilayerNetwork/data/regrOutput/paramSweep/nonlinear/' + subtype + '_TFGeneCapacities_ri_nomi__lambda' + l + '.txt'
    f = open(fn,'r')
    res = f.readlines()
    f.close()
    res = res[1:]
    res = [line.strip().split('\t') for line in res]
    res = [[np.log2(float(y)) for y in x[1:]] for x in res if x[0] != 'const']
    res = np.array(res)
    nTFs = res.shape[0]
    nTgts = res.shape[1]
    #ransac = sklm.RANSACRegressor(min_samples=0.8)
    X = res[sigTF_idx].reshape(-1,1)[sigTgt_idx]
    for i in range(0,nTFs):
        if i != sigTF_idx:
            y = res[i][sigTgt_idx]
            if X.size >= 5 and y.size >= 5:
                try:
                    #temp = ransac.fit(X,y)
                    #inlier_mask = ransac.inlier_mask_
                    #if np.sum(inlier_mask) > 3:
                    #cormat[count,i] = np.corrcoef(X[inlier_mask].reshape(1,-1)[0],y[inlier_mask])[0,1]
                    cormat[count,i] = np.corrcoef(X.reshape(1,-1)[0],y)[0,1]
                except Exception:
                    pass
    sys.stdout.write('.')
    sys.stdout.flush()
    count += 1

print(' ')
fn = '/users/yunpengl/Data/multilayerNetwork/analyses/0_0.1/' + sigTF + '_sigTgtCormat_' + hl + '_' + l + '.txt'
np.savetxt(fn,cormat,delimiter='\t')
