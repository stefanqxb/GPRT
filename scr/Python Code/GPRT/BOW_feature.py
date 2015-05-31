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
    bow_feature = BOW(peptide,hydrophobicity,peptide_length,aaAlphabet,num)
    return psmDescriptions ,bow_feature

def BOW(string, hydrop,length,aaAlphabet,num2):
    if num2 == 1:
        feature = forming_feature_single(string,aaAlphabet,1)  # 1-grams
    elif num2 == 2:
        feature = forming_feature(string,hydrop,aaAlphabet,length,2)  # 2-grams
    return feature

def forming_feature(string,hydrop,aaAlphabet,length,num2):
    M = len(string)
    N = len(aaAlphabet)
    weight =[0.4,0.5,0.1]
    feature = np.zeros([M,N*N+3])
    refer = grab_word(aaAlphabet,num2)
    for i in range(M):
        idx = 0
        T =len(string[i])
        while idx < T-1:
            if  idx < T - 2 and string[i][idx+2] == '[':
                nextIdx = string[i][idx+1:].find(']') + idx + 2
                a = string[i][idx:nextIdx]
                idx +=1
            elif  string[i][idx+1] == '[':
                nextIdx = string[i][idx+1:].find(']') + idx + 2
                a = string[i][idx:nextIdx+1]
                idx = nextIdx
            else:
                a = string[i][idx]+ string[i][idx+1]
                idx += 1
            if a.find('-') == -1:
               index = refer.index(a)
               if index != -1:
                   feature[i][index] = feature[i][index] +1

        feature[i][0:N*N] = float(weight[0])*feature[i][0:N*N]
        feature[i][N*N+1] = feature[i][N*N+1] + float(weight[1])*float(hydrop[i])
        feature[i][N*N+2] = feature[i][N*N+2] + float(weight[2])*float(length[i])
        feature[i] = feature[i][:]/float(sum(feature[i],0))
    return feature

def forming_feature_single(string,aaAlphabet,num2):
    M = len(string)
    N = len(aaAlphabet)
    weight =[0.4,0.5,0.1]
    feature = np.zeros([M,N])
    refer = grab_word(aaAlphabet,num2)
    for i in range(M):
        idx = 0
        T =len(string[i])
        while idx < T-1:
            if  string[i][idx+1] == '[':
                nextIdx = string[i][idx+1:].find(']') + idx + 2
                a = string[i][idx:nextIdx]
                idx = nextIdx
            else:
                a = string[i][idx]
                idx += 1
            if a.find('-') == -1:
               index = refer.index(a)
               if index != -1:
                   feature[i][index] = feature[i][index] +1
        feature[i][0:N] = float(weight[0])*feature[i][0:N]
        feature[i] = feature[i][:]/float(sum(feature[i],0))
    return feature


def grab_word(aaAlphabet,num2):
    s1 = aaAlphabet
    num = len(aaAlphabet)
    sub_word = []
    if num2 ==1:
        for i in range(num):
            temp = s1[i]
            sub_word.append(temp)
    elif num2 ==2:
        for i in range(num):
           for j in range(num):
               temp = s1[i] + s1[j]
               sub_word.append(temp)
    return sub_word



def bowindexSum(aas,index):
    sum_ind = 0
    for aa in aas:
        temp = index[aa]
        sum_ind += temp
    return sum_ind
