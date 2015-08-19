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
import data_spilt


# usage = "usage: %prog [-f 'input file path'][-w 'feature extraction method'][-s 'number of subset']"
# parser = OptionParser(usage=usage)
# parser.add_option("-t", "--trafilename",type = "string", dest="trainfile",help="import the train file")
# parser.add_option("-e", "--testfilename",type = "string", dest="testfile",help="import the test file")
# parser.add_option("-w", "--request1",type="int", dest="request1",help="Choose feature extraction method: 1-Bag-of-Words , 2-Elude")
# parser.add_option("-b", "--request12",type= "int", dest="request12",help="Choose the number of gram for BOW, if request1==2, set 0 to it")
# parser.add_option("-s", "--request2",type = "int", dest="request2",help="Number of subset, set 0 for the whole dataset")
# parser.add_option("-g", "--ef",type = "int", dest="ge",help="Number of iteration in CG")
#
# parser.print_help()
#
# (options, args) = parser.parse_args()
# request1 = options.request1
# request12 = options.request12
# request2 = options.request2
# trainFile = options.trainfile
# testFile = options.testfile
# evn = options.ge

File = 'data/retention_time_peptide.csv'
trainFile = 'data/set1/train_set5000.csv'
testFile = 'data/set1/test_set.csv'
request1 = 2
request12 = 0
request2 = 0
evn = 10
model = 'gp'

#data_spilt.split_data(File)

print "forming features....."

# peptide, feature, rt = feature_extraction.feature_extra(File, request1,request12,request2)
# row = len(feature)

peptide, train_feature, train_rt = feature_extraction.feature_extra(trainFile, request1, request12, request2)
peptide, test_feature, test_rt = feature_extraction.feature_extra(testFile, request1, request12, request2)


print "training model......."

print "training svr......"
########  SVR model ##########
## original [train_set_svr, train_tag_svr, test_set_svr, test_tag_svr] = data_generator.processing_Data(feature, rt, row,'svr')
[train_set_svr, train_tag_svr, test_set_svr, test_tag_svr] = data_generator.processing_Data_tt(train_feature, test_feature, train_rt, test_rt,'svr')
t0 = time.time()
svr1 = svm.SVR(C=600, gamma=0.1, coef0=0.0, degree=3, epsilon=0.1, kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
svr = grid_search.GridSearchCV(svr1, param_grid ={"C":np.linspace(100,1000,num = 10),"gamma":np.linspace(0.01,10,num = 100)})
svr.fit(train_set_svr,train_tag_svr)
pv_svr = svr.predict(test_set_svr)
t_svr = time.time()-t0

print "training gp......"
#########  Gaussian Process model ############
#[train_set, train_tag, test_set, test_tag] = data_generator.processing_Data(feature, rt, row,'gp')
[train_set, train_tag, test_set, test_tag] = data_generator.processing_Data_tt(train_feature, test_feature, train_rt, test_rt,'gp')
t0 = time.time()
kern = GPy.kern.Poly(train_set.shape[1],order = 5) + GPy.kern.RBF(train_set.shape[1], lengthscale = 3 )
m = GPy.models.GPRegression(train_set, train_tag, kern)

print "optimizing the hyperparameters....."
m.optimize('scg',max_iters = evn)
pv_gp ,ps = m.predict(test_set)
t_gp = time.time() - t0

print "Evaluating result...."

step ,max_total_gp,max_total_svr = Evaluation_and_Ploting.evaluation(pv_gp,pv_svr,test_tag,test_tag_svr,t_gp,t_svr)
Evaluation_and_Ploting.ploting(pv_gp,pv_svr, test_tag,test_tag_svr,max_total_gp,max_total_svr,step)
