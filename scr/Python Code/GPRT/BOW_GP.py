__author__ = 'Administrator'
import sys
#import time,pdb, random
import numpy as np
from Bio.SubsMat import MatrixInfo
from sklearn import grid_search, neighbors, svm
import data_generator, feature_extraction
import matplotlib.pyplot as plt
import scipy.io
import scipy.sparse as sp
import math
import PSO_GP
import GPy


# import data
f = "data.txt"
g = "hydrophobicity.txt"

#
data = data_generator.read_data(f)
hydrophobicity_len = data_generator.read_data(g)
peptide = data.data[0]
rt = data.data[2]
hydrophobicity = hydrophobicity_len.data[2]
length = hydrophobicity_len.data[1]

# feature selection
# 1---Bag of Word     2---Alignment Analysis  normalization is already done
feature = feature_extraction.feature_extra(peptide, hydrophobicity,length ,1)


# training model
row = len(peptide[0:3000])
ratio = 0.8
train_row = int(round(ratio * row))
end = row #len(feature)

train_tag_temp = np.zeros([train_row,1])
test_tag_temp = np.zeros([end-train_row,1])

for i in range(train_row):
    train_tag_temp[i] = rt[i]
for i in range(end-train_row):
    test_tag_temp[i] = rt[train_row+i]


train_set = feature[0:train_row][:]
train_tag = train_tag_temp
test_set = feature[train_row:end][:]
test_tag = test_tag_temp

m = GPy.models.GPRegression(train_set,train_tag)
m.optimize('bfgs')
pv,ps = m.predict(test_set)
temp = pv - test_tag
t = np.zeros([len(temp),1])
for i in range(len(temp)):
    t[i] = temp[i]**2

mse = np.mean(t)
print mse


# optimization of the hyper parameter








