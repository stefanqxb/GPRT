__author__ = 'Administrator'
import numpy as np
import data_generator
import math

def feature_extra(string ,hydrop,length, num):
    if num == 1:
        feature = BOW(string,hydrop,length)
    else:
        print('Warning : Unknown type of method !!')

    return feature

def BOW(string, hydrop, length):
    #fea_mat = extracting_feature(string)
    feature = forming_feature(string,hydrop,length)
    return feature

def extracting_feature(string):
    fea_mat = np.zeros([20,20])
    M = len(string)
    # 2-gram word
    for i in range(M):
        N = len(string[i])-1
        for j in range(N):
            a = [string[i][j] , string[i][j+1]]
            idx = np.zeros([2,1])
            for k  in range(2):
                temp = identify(a[k])
                if temp != -1:
                   idx[k] = temp
            fea_mat[int(idx[0])][int(idx[1])] = fea_mat[int(idx[0])][int(idx[1])] + 1

    return fea_mat

def identify(a):
    letter ='ARNDCQEGHILKMFPSTWYVOU'
    ind = letter.find(a)
    return ind

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





