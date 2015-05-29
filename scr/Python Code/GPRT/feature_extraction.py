__author__ = 'Administrator'
import numpy as np
import data_generator
import sys
sys.path.append('elude/')
import elude_features as el
import math

def feature_extra( num,num2,subset_num):
    if num == 1:
        # Bag of Word model
        f = "data/data.txt"
        g = "hydrophobicity.txt"
        data = data_generator.read_data(f)
        hydrophobicity_len = data_generator.read_data(g)
        peptide = data.data[0]
        rt = data.data[2]
        hydrophobicity = hydrophobicity_len.data[2]
        length = hydrophobicity_len.data[1]
        feature_temp = BOW(peptide, hydrophobicity,length,num2)
        feature = feature_temp[0:subset_num]
        rt = rt[0:subset_num]
    elif num == 2:
        #execute elude feature
        psmDescriptions, feature_temp = el.getFeatures('data/modify.csv')
        feature = feature_temp[0:subset_num]
        peptide,rt = data_generator.extract(psmDescriptions,subset_num) # grab 1000 as subset to test
    else:
        print('Warning : Unknown type of method !!')

    return feature,rt

def BOW(string, hydrop, length,num2):
    if num2 == 1:
        feature = forming_feature_single(string)  # 1-grams
    elif num2 == 2:
        feature = forming_feature(string,hydrop,length) # 2-grams
    return feature


def forming_feature(string,hydrop,length):
    M = len(string)
    N = 20
    weight =[0.4,0.5,0.1]
    feature = np.zeros([M,N*N+3])
    refer = grab_word(string)
    for i in range(M):
        T= len(string[i])
        for j in range(T-1):
            a = string[i][j]+string[i][j+1]
            if a.find('-') == -1:
               index = refer.index(a)
               if index != -1:
                   feature[i][index] = feature[i][index] +1
        feature[i][0:N*N] = float(weight[0])*feature[i][0:N*N]
        feature[i][N*N+1] = feature[i][N*N+1] + float(weight[1])*float(hydrop[i])
        feature[i][N*N+2] = feature[i][N*N+2] + float(weight[2])*float(length[i])
        feature[i] = feature[i][:]/float(sum(feature[i],0))
    return feature

def forming_feature_single(string):
    M = len(string)
    N = 20
    weight =[0.4,0.5,0.1]
    feature = np.zeros([M,N])
    refer = grab_word_single(string)
    for i in range(M):
        for aa in string[i]:
            if aa.find('-') == -1:
               index = refer.index(aa)
               if index != -1:
                   feature[i][index] = feature[i][index] +1
        feature[i][0:N] = float(weight[0])*feature[i][0:N]
        feature[i] = feature[i][:]/float(sum(feature[i],0))
    return feature


def grab_word(string):
    letter ='A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V'
    str1 = letter.split(',')
    str2 = letter.split(',')
    sub_word = []
    for s in str1:
        for j in str2:
            temp = s+j
            sub_word.append(temp)
    return sub_word


def grab_word_single(string):
    letter ='A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V'
    str1 = letter.split(',')
    sub_word = []
    for s in str1:
        sub_word.append(s)
    return sub_word


