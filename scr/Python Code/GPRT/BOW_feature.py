__author__ = 'Administrator'
import numpy as np
import data_generator
import sys
sys.path.append('elude/')
import elude_features as el
import retention_model as rm
import data_manager as dm
import data_generator as dg
import operator
import math

def getBOWfeature(file,num):
    hydrophobicity =[]
    peptide_length = []
    psmDescriptions, aaAlphabet = el.processTrainData(file)
    Index = rm.buildRetentionIndex(aaAlphabet,psmDescriptions,True)
    customIndex = dict(zip(aaAlphabet, Index))
    for i, psmd in enumerate(psmDescriptions):
         aas = dm.getAminoAcidList(psmd.peptide)
         peptide_length.append(len(psmd.peptide))
         hydrophobicity.append(bowindexSum(aas, customIndex))
    peptide,rt = dg.extract(psmDescriptions,len(psmDescriptions))
    bow_feature = BOW(peptide,hydrophobicity,peptide_length,num)
    return psmDescriptions ,bow_feature

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


def bowindexSum(aas,index):
    sum_ind = 0
    for aa in aas:
        temp = index[aa]
        sum_ind += temp
    return sum_ind
