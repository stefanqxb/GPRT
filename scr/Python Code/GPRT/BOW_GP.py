__author__ = 'Administrator'
import sys
import numpy as np
from sklearn import grid_search, neighbors, svm
import data_generator, feature_extraction, feature_extraction_single_aa
import matplotlib.pyplot as plt
import scipy.io
import math
import GPy
import calc_mtw
import pylab
from scipy.stats.stats import pearsonr
import time

# import data
f = "data.txt"
g = "hydrophobicity.txt"
data = data_generator.read_data(f)
hydrophobicity_len = data_generator.read_data(g)
peptide = data.data[0]
rt = data.data[2]
hydrophobicity = hydrophobicity_len.data[2]
length = hydrophobicity_len.data[1]

# feature selection
print "selecting features....."
#feature = feature_extraction.feature_extra(peptide, hydrophobicity,length ,1)
feature = feature_extraction_single_aa.feature_extra(peptide,1)

# train model
row = len(peptide[0:2000])

[train_set_svr, train_tag_svr, test_set_svr, test_tag_svr] = data_generator.processing_Data(feature, rt, row,'svr')

max_t = max(test_tag_svr)
min_t = min(test_tag_svr)

print "training model......."
#########  SVR model ##########
svr1 = svm.SVR(C=700, coef0=0.0, degree=3, epsilon=0.1, gamma= 50, kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
svr = grid_search.GridSearchCV(svr1, param_grid ={"C":np.linspace(100,700,num = 7),"gamma":np.linspace(10,100,num = 15)})
t0 = time.time()
svr.fit(train_set_svr,train_tag_svr)
pv_svr = svr.predict(test_set_svr)
t_svr = time.time()-t0
#########  Gaussian Process model ############
[train_set, train_tag, test_set, test_tag] = data_generator.processing_Data(feature, rt, row,'gp')
t0 = time.time()
kern = GPy.kern.RBF(train_set.shape[1])
m = GPy.models.GPRegression(train_set, train_tag, kern)
print "optimizing the hyperparameters....."
m.optimize('bfgs')
pv_gp ,ps = m.predict(test_set)
t_gp = time.time() - t0

# result evaluation
print "evaluating result...."
step = 100
diff_svr = abs(pv_svr - test_tag_svr)
diff_gp = abs(pv_gp -test_tag)
histo_svr = plt.hist(diff_svr,step)
histo_gp = plt.hist(diff_gp,step)

print "calculating the minimal time window and correlation coefficient......"

mtw_gp =calc_mtw.mini_time_win(histo_gp[0],max_t,min_t,step)
corrcoef_gp = pearsonr(pv_gp,test_tag)
mtw_svr = calc_mtw.mini_time_win(histo_svr[0],max_t,min_t,step)
corrcoef_svr = pearsonr(pv_svr,test_tag_svr)

print 'corrcoef of GP = ',corrcoef_gp[0], 'corrcoef of SVR = ',corrcoef_svr[0]

pylab.figure(1)
plt.hold('on')
plt.subplot(121)
plt.plot(pv_svr,test_tag, 'g*')
plt.plot([1,max_t],[1,max_t],'y-',linewidth = 2)
plt.title('SVR')
plt.subplot(122)
plt.plot(pv_gp,test_tag, 'r*',label = "correlation coefficient of GP")
plt.plot([1,max_t],[1,max_t],'b-',linewidth = 2)
plt.title('GP')

pylab.figure(2)
plt.hold('on')
plt.subplot(121)
line1 = plt.hist(diff_gp,bins = step,label = 'GP histogram',color = 'red')
plt.legend()
plt.subplot(122)
line2 = plt.hist(diff_svr,bins = step,label = 'SVR histogram',color = 'green')
plt.legend()
pylab.show()






