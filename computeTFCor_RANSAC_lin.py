import pandas as pd
import sklearn.linear_model as sklm
import numpy as np
from matplotlib import pyplot as plt
import sys

subtype = sys.argv[1]
l = sys.argv[2]
fn = '/users/yunpengl/Data/multilayerNetwork/data/regrOutput/paramSweep/linear/' + subtype + '_TFGeneCapacities_ri_nomi__lambda' + l + '.txt'
f = open(fn,'r')
res = f.readlines()
f.close()
res = res[1:]
res = [line.strip().split('\t') for line in res]
res = [[np.log2(float(y)) for y in x[1:]] for x in res if x[0] != 'const']
res = np.array(res)
nTFs = res.shape[0]
nTgts = res.shape[1]
cormat = np.zeros((nTFs,nTFs))

# Per scikit-learn docs: when in doubt, use RANSAC (good scalability with number
# of samples for this case).

# Use RANSAC to detect inliers and outliers. We do not
# need the regression results here but just need the inliers
# for correlation calculation back in R.
ransac = sklm.RANSACRegressor(min_samples=0.8)
# idx = (res[i]!=1) & (res[j]!=1)
# X = np.log2(res[i].reshape(-1,1)[idx])
# y = np.log2(res[j][idx])
# ransac.fit(X,y)
# inlier_mask = ransac.inlier_mask_
# np.corrcoef(X[inlier_mask].reshape(1,-1)[0],y[inlier_mask])


for i in range(0,nTFs):
    if i % 10 == 0:
        print(".")
    for j in range(i,nTFs):
        if j > i:
            idx = (res[i]!=0) & (res[j]!=0)
            X = res[i].reshape(-1,1)[idx]
            y = res[j][idx]
            if X.size > 10 and y.size > 10:
                try:
                    temp = ransac.fit(X,y)
                    inlier_mask = ransac.inlier_mask_
                    cormat[i,j] = np.corrcoef(X[inlier_mask].reshape(1,-1)[0],y[inlier_mask])[0,1]
                except Exception:
                    pass

fn = '/users/yunpengl/Data/multilayerNetwork/data/' + subtype + '_cormat_' + l + '_lin.txt'
np.savetxt(fn,cormat,delimiter='\t')
