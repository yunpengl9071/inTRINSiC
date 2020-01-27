import pandas as pd
import sklearn.linear_model as sklm
import numpy as np
from matplotlib import pyplot as plt
import sys
import random

random.seed(9164079)
np.random.seed(9164079)


subtype = sys.argv[1]
print(subtype)
#l = sys.argv[2]
TF = sys.argv[2]
#survival = sys.argv[3]
fn = '/users/yunpengl/Data/multilayerNetwork/analyses/0_0.1/' + TF + '_' + subtype +  '_score_vs_survival.txt'
f = open(fn,'r')
res = f.readlines()
f.close()
res = [line.strip().split('\t') for line in res]
res = [[float(y) for y in x] for x in res]
res = np.array(res)
# cormat = np.zeros((nTFs,nTFs))

# Per scikit-learn docs: when in doubt, use RANSAC (good scalability with number
# of samples for this case).

x = res[:,0]
y = res[:,1]
idx = (x!=0) & (y!=0)
x = (x[idx])
y = (y[idx])
ransac = sklm.RANSACRegressor(min_samples=0.5)
X = x.reshape(-1,1)
temp = ransac.fit(X,y)
inlier_mask = ransac.inlier_mask_
cor = np.corrcoef(X[inlier_mask].reshape(1,-1)[0],y[inlier_mask])[0,1]
print(cor)
x_inlier = X[inlier_mask]
y_inlier = y[inlier_mask]
fn = '/users/yunpengl/Data/multilayerNetwork/analyses/0_0.1/' + TF + '_' + subtype + '_corWithSurvival_RANSAC.txt'
f = open(fn,'w')
for i in range((x_inlier.size)):
    line = str(x_inlier[i,0]) + '\t' + str(y_inlier[i]) + '\n'
    f.write(line)
f.close()

# for i in range(0,nTFs):
#     if i % 10 == 0:
#         sys.stdout.write('.')
#         sys.stdout.flush()
#     for j in range(i,nTFs):
#         if j > i:
#             idx = (res[i]!=0) & (res[j]!=0)
#             X = res[i].reshape(-1,1)[idx]
#             y = res[j][idx]
#             if X.size > 10 and y.size > 10:
#                 try:
#                     temp = ransac.fit(X,y)
#                     inlier_mask = ransac.inlier_mask_
#                     if np.sum(inlier_mask) > 5:
#                         cormat[i,j] = np.corrcoef(X[inlier_mask].reshape(1,-1)[0],y[inlier_mask])[0,1]
#                 except Exception:
#                     pass

# print(cormat[0][1])
# print(' ')

# fn = '/users/yunpengl/Data/multilayerNetwork/data/' + subtype + '_cormat_' + l + '.txt'
# np.savetxt(fn,cormat,delimiter='\t')
