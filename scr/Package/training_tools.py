#!/usr/bin/python

import numpy as np
import GPy
from sklearn import svm, grid_search
import time
import sys
sys.path.append('elude/')
import data_manager as dm
import random

class rt_trainer:
    def __init__(self,peptide):
        self.peptide = peptide
        self.feature = []
        self.trainSet = []
        self.trainTag = []
        self.testSet = []
        self.testTag = []
        self.preTag = []
        self.preUnc = []
        self.time = 0

    def extract(inst,num):
        peptide =[]
        rt = []
        for item in inst:
            if isinstance(item,dm.PSMDescription):
                peptide.append(item.peptide)
                rt.append(item.retentionTime)
            if len(rt) >= num:
                break
        return peptide,rt


    def dataSet_split(self, subset_num, str):
        total = range(len(self.feature))
        random.shuffle(total)
        if subset_num != 0:
            seed = total[0:subset_num]
            feature = self.feature[seed]
            peptide,rt = self.extract(dm.psmDescriptions,len(self.feature))
        else:
            seed = total
            feature = self.feature[total]
            peptide,rt = self.extract(dm.psmDescriptions,len(self.feature))

        for item in seed:
            self.trainSet.append(peptide[item])
            self.trainTag.append(rt[item])

        if str =='svr':
            train_tag_temp = np.zeros([len(train_rt)])
            test_tag_temp = np.zeros([len(test_rt)])
        elif str =='gp':
            train_tag_temp = np.zeros([len(train_rt),1])
            test_tag_temp = np.zeros([len(test_rt),1])

        for i in range(len(train_rt)):
            train_tag_temp[i] = train_rt[i]
        for i in range(len(test_rt)):
            test_tag_temp[i] = test_rt[i]

        train_set = train_fea
        train_tag = train_tag_temp
        test_set = test_fea
        test_tag = test_tag_temp
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



