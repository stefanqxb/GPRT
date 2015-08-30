#!/usr/bin/python

import numpy as np
import random
import feature_extraction
import GPy
import matplotlib.pyplot as plt

from scipy.stats.stats import pearsonr
from sklearn import svm, grid_search

from joblib import Parallel, delayed  
import multiprocessing

def pcv_train(i,bench):
	model = bench.train_model(i)
	score = bench.eval_model(i,model)
	return score

def parallel_cross_validataion( bench ):
	num_cores = multiprocessing.cpu_count()
	results = Parallel(n_jobs=num_cores)( delayed(pcv_train)(i,bench) for i in range(bench.parts.nfolds) )
	return np.mean( results ), np.std( results )

class partitions:
    def __init__(self, ndata, nfolds):
        self.ndata = ndata
        self.nfolds = nfolds
        np.random.seed(8766)

    def n_train(self):
        return len(self.train_parts[0])

    def n_test(self):
        return len(self.test_parts[0])

    def gen_cross_val(self):
        perm = np.random.permutation(self.ndata)
        self.train_parts = [];
        self.test_parts = [];
        for i in range(self.nfolds):
            train = []
            test = []
            for j in range(self.ndata):
                if j % self.nfolds == i:
                    test.append(perm[j])
                else:
                    train.append(perm[j])
            self.train_parts.append(np.array(train))
            self.test_parts.append(np.array(test))

    def gen_rand_splits(self, ratio):
        self.train_parts = []
        self.test_parts = []
        for i in range(self.nfolds):
            perm = np.random.permutation(self.ndata)
            train = []
            test = []
            for j in range(self.ndata):
                r = float(j) / self.ndata
                if r < ratio:
                    train.append(perm[j])
                else:
                    test.append(perm[j])
            self.train_parts.append(np.array(train))
            self.test_parts.append(np.array(test));

    def get_train_part(self, ind):
        return self.train_parts[ind]

    def get_test_part(self, ind):
        return self.test_parts[ind]

class eval_tools:
    def mean_square_error(self, actual, predicted):
        return np.sum(np.power(actual - predicted, 2)) / len(actual)

    def mini_time_window(self, hist, diff, step, max_total):
        max_t = max(diff)
        min_t = min(diff)
        total = sum(hist)
        threshold = round(0.95 * total)
        count = 0
        counter = 0
        for i in range(len(hist)):
            count = hist[i] + count
            if count >= threshold:
                counter = i
                break
        time_interval = 2 * counter * (max_t - min_t) / (step * max_total)
        return time_interval

    def delta_t(self, actual, predicted):
		step = 10
		diff = abs(actual - predicted)	
		histo = np.histogram(diff,step)
		#histo = plt.hist(diff, step)
		max_total = max(actual) - min(actual)
		mtw = self.mini_time_window(histo[0], diff, step, max_total)
		corrcoef = pearsonr(actual, predicted)
		return mtw

class rt_model:
    def __init__(self, feature,model_type, model, norm, voc, em):
        self.feature = feature
        self.model_type = model_type
        self.model = model
        self.norm = norm
        self.voc = voc
        self.em = em

    def eval(self, p):
        res = [];
        if self.feature == "bow":
            vec = p.bow_descriptor(self.voc)
        if self.feature == "elude":
            vec = p.elude_descriptor(self.em)
        self.norm.normalize(vec)
        vec = np.matrix(vec)
        vals = self.model.predict(np.array(vec))
        if self.model_type == 'gp':
            res = ( vals[0][0, 0], vals[1][0, 0] )
        elif self.model_type == 'svr':
            res = ( int(vals), 0)
        return res


class rt_benchmark:
    def __init__(self, peptides, feature, model_type, ntrain=-1, nfolds = 5 ):
		self.peptides = peptides
		self.feature = feature
		self.model_type = model_type
		self.ntrain = ntrain
		self.parts = partitions(len(peptides), nfolds)
		self.parts.gen_rand_splits(0.70)
		#self.parts.gen_cross_val();
		
		if ntrain < 0 or ntrain > self.parts.n_train():
			ntrain = self.parts.n_train()

    def train_model(self, ind):
        train_peptides = self.peptides[self.parts.get_train_part(ind)]
        mg = feature_extraction.model_generator(train_peptides)
        voc = mg.get_bow_voc(2)
        em = mg.get_elude_model()

        Y = [];
        X = [];
        for p in train_peptides[0:self.ntrain]:
            Y.append(p.rt)
            if self.feature == 'bow':
                X.append(p.bow_descriptor(voc))
            elif self.feature == 'elude':
                X.append(p.elude_descriptor(em))
        X = np.matrix(X)
        if self.model_type == 'gp':
            Y = np.transpose(np.matrix(Y))
        elif self.model_type == 'svr':
            Y = np.transpose(Y)

        norm = feature_extraction.normalizer()
        if self.feature == "elude":
            norm.normalize_maxmin(X);
        if self.model_type == "gp":
            m = GPy.models.GPRegression(X,Y)
            m.optimize_restarts(num_restarts=10,verbose=False)
        elif self.model_type == "svr":
            m = svm.SVR(C=600, gamma=0.1, coef0=0.0, degree=3, epsilon=0.1, kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)
            m = grid_search.GridSearchCV(m, param_grid={"C": np.linspace(100, 1000, num=10), "gamma": np.linspace(0.01,10, num = 100)})
            m.fit(X, Y)
        return rt_model(self.feature,self.model_type, m, norm, voc, em)

    def eval_model(self, ind, model):
        test_peptides = self.peptides[self.parts.get_test_part(ind)]
        actual = []
        predicted = []
        for p in test_peptides:
            actual.append(p.rt)
            v, s = model.eval(p)
            predicted.append(v)
        actual = np.array(actual)
        predicted = np.array(predicted)
        et = eval_tools()
        return et.delta_t(actual, predicted)
