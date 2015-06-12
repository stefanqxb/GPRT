__author__ = 'Administrator'
import sys
sys.path.append('elude/')
import numpy as np
from sklearn import grid_search, neighbors, svm
import data_generator, feature_extraction
import GPy

import time
import Help_file
import Pre_option
import Evaluation_and_Ploting


request1,request12,request2 = Pre_option.option()

# feature selection
print "selecting features....."

File = 'data/retention_time_peptide.csv'
feature,rt = feature_extraction.feature_extra(File,request1,request12,request2)

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
print "Evaluating result...."

diff_svr ,diff_gp,step, max_total_gp, max_total_svr = Evaluation_and_Ploting.evaluation(pv_svr,test_tag_svr,pv_gp,test_tag,t_gp,t_svr)
Evaluation_and_Ploting.ploting(pv_svr,pv_gp,test_tag,max_total_svr,max_total_gp,diff_svr,diff_gp,step)







