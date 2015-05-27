__author__ = 'Administrator'
import numpy as np

def feature_extra(string, num):
    if num == 1:
        feature = BOW(string)
    else:
        print('Warning : Unknown type of method !!')

    return feature

def BOW(string):
    #fea_mat = extracting_feature(string)
    feature = forming_feature(string)
    return feature

def forming_feature(string):
    M = len(string)
    N = 20
    weight =[0.4,0.5,0.1]
    feature = np.zeros([M,N])
    refer = grab_word(string)
    for i in range(M):
        for aa in string[i]:
            if aa.find('-') == -1:
               index = refer.index(aa)
               if index != -1:
                   feature[i][index] = feature[i][index] +1
        feature[i][0:N] = float(weight[0])*feature[i][0:N]
        feature[i] = feature[i][:]/float(sum(feature[i],0))
    return feature

