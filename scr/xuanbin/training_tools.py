#!/usr/bin/python

import numpy as np
import GPy
from sklearn import svm, grid_search
import time
import sys
sys.path.append('elude/')
import data_manager as dm
import random
import feature_extraction as fe

class rt_trainer:
    def __init__(self, peptide, str):

        self.peptide = []
        self.rt = np.zeros((len(peptide),1))
        self.ef = np.zeros((len(peptide), len(peptide[0].elude_descriptor)))
        self.bf = np.zeros((len(peptide), len(peptide[0].bow_descriptor)))

        for i, pp in peptide:
            self.peptide.append(pp.sequence)
            self.rt[i] = pp.rt
            self.ef[i][:] = pp.elude_descriptor
            self.bf[i][:] = pp.bow_descriptor



        if str == "svr":
            self.trainSet = np.zeros((1,0))
            self.testSet = np.zeros((1,0))
        elif str =='gp':
            self.trainSet = np.zeros([1,1])
            self.testSet = np.zeros([1,1])

        self.preTag = []
        self.preUnc = []
        self.time = 0


    def dataSet_split(self, trainSetNum, str):
        total = range(len(self.peptide))
        random.shuffle(total)
        seed = total[0:trainSetNum]
        nonseed = total[trainSetNum:]
        self.trainSet.append(self.ef[seed])
        self.trainTag.append(self.rt[seed])
        self.testSet.append(self.ef[nonseed])
        self.testTag.append(self.rt[nonseed])

        if str =='svr':
            train_tag_temp = np.zeros([len(self.trainSet)])
            test_tag_temp = np.zeros([len(self.testSet)])
        elif str =='gp':
            train_tag_temp = np.zeros([len(self.trainSet),1])
            test_tag_temp = np.zeros([len(self.testSet),1])

        for i in range(len(self.trainSet)):
            self.trainTag[i] = self.rt[i]
        for i in range(len(self.testSet)):
            self.testTag[i] = self.testTag[i]

        return train_set,train_tag,test_set,test_tag

    def train_svr(self):
        print "Training svr"
        t0 = time.time()
        svr1 = svm.SVR(C=600, gamma=0.1, coef0=0.0, degree=3, epsilon=0.1, kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
        svr = grid_search.GridSearchCV(svr1, param_grid={"C": np.linspace(100, 1000, num=10), "gamma": np.linspace(0.01,10, num = 100)})
        svr.fit(self.trainSet, self.trainTag)
        self.preTag = svr.predict(self.testSet)
        self.time = time.time()-t0

    def train_gp(self):
        print "Training GP"
        t0 = time.time()
        kern = GPy.kern.Poly(self.trainSet.shape[1],order = 5) + GPy.kern.RBF(self.trainSet.shape[1], lengthscale=3 )
        m = GPy.models.GPRegression(self.trainSet, self.trainTag, kern)
        m.optimize('scg', max_iters = 10)
        self.preTag, self.preUnc = m.predict(self.testSet)
        self.time = time.time() - t0



