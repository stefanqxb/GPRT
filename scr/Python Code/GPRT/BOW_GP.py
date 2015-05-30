__author__ = 'Administrator'
import sys
sys.path.append('elude/')
import numpy as np
from sklearn import grid_search, neighbors, svm
import data_generator, feature_extraction
import matplotlib.pyplot as plt
import GPy
import calc_mtw
import pylab
import time
from scipy.stats.stats import pearsonr

# feature selection
print "selecting features....."
feature,rt = feature_extraction.feature_extra(1,2,3000) # num: 1--bow 2--elude num2 1--1gram word 2--2gram word  subset_num: subset

row = len(feature)

print "training model......."
#########  SVR model ##########
[train_set_svr, train_tag_svr, test_set_svr, test_tag_svr] = data_generator.processing_Data(feature, rt, row,'svr')
t0 = time.time()
svr1 = svm.SVR(C=600, gamma=0.1, coef0=0.0, degree=3, epsilon=0.1, kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
svr = grid_search.GridSearchCV(svr1, param_grid ={"C":np.linspace(100,1000,num = 10),"gamma":np.linspace(0.01,10,num = 100)})
svr.fit(train_set_svr,train_tag_svr)
pv_svr = svr.predict(test_set_svr)
t_svr = time.time()-t0
#########  Gaussian Process model ############
[train_set, train_tag, test_set, test_tag] = data_generator.processing_Data(feature, rt, row,'gp')
t0 = time.time()
kern = GPy.kern.RBF(train_set.shape[1])
m = GPy.models.GPRegression(train_set, train_tag, kern)
print "optimizing the hyperparameters....."
m.optimize()
pv_gp ,ps = m.predict(test_set)
t_gp = time.time() - t0

# result evaluation
print "evaluating result...."
step = 10
diff_svr = abs(pv_svr - test_tag_svr)
diff_gp = abs(pv_gp -test_tag)
histo_svr = plt.hist(diff_svr,step)
histo_gp = plt.hist(diff_gp,step)

print "calculating the minimal time window and correlation coefficient......"


max_total_gp = max(test_tag)
max_total_svr = max(test_tag)

mtw_gp =calc_mtw.mini_time_win(histo_gp[0],diff_gp,step,max_total_gp)
corrcoef_gp = pearsonr(pv_gp,test_tag)
mtw_svr = calc_mtw.mini_time_win(histo_svr[0],diff_svr,step,max_total_svr)
corrcoef_svr = pearsonr(pv_svr,test_tag_svr)

print 'corrcoef of GP = ',corrcoef_gp[0], '   corrcoef of SVR = ',corrcoef_svr[0]
print 'mtw of GP = ',mtw_gp, '   mtw of SVR = ', mtw_svr
print 'running time of GP = ',t_gp, '   running time of SVR = ',t_svr

pylab.figure(1)
plt.hold('on')
plt.subplot(121)
plt.plot(pv_svr,test_tag, 'g*')
plt.plot([1,max_total_svr],[1,max_total_svr],'y-',linewidth = 2)
plt.title('SVR')
plt.subplot(122)
plt.plot(pv_gp,test_tag, 'r*',label = "correlation coefficient of GP")
plt.plot([1,max_total_gp],[1,max_total_gp],'b-',linewidth = 2)
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






