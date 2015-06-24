__author__ = 'Administrator'
import sys
sys.path.append('elude/')
import numpy as np
from sklearn import grid_search, svm
import data_generator, feature_extraction
import GPy
import time
from optparse import OptionParser
import Evaluation_and_Ploting

usage = "usage: %prog [-f 'input file path'][-w 'feature extraction method'][-s 'number of subset']"
parser = OptionParser(usage=usage)
parser.add_option("-w", "--request1",
              type="int", dest="request1",
              help="Choose feature extraction method: 1-Bag-of-Words , 2-Elude")
parser.add_option("-g", "--request12",
              type= "int", dest="request12",
              help="Choose the number of gram for BOW, if request1==2, ignore this term")
parser.add_option("-s", "--request2",
              type = "int", dest="request2",
              help="Number of subset, set 0 for the whole dataset")
parser.add_option("-f", "--filename",
              type = "string", dest="filename",
              help="import the target file")
parser.print_help()

(options, args) = parser.parse_args()
request1 = options.request1
request12 = options.request12
request2 = options.request2
File = options.filename


print "selecting features....."
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







