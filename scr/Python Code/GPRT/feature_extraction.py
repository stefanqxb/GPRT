__author__ = 'Administrator'
import numpy as np
import data_generator
import sys
sys.path.append('elude/')
import elude_features as el
import BOW_feature as bf
import math

def feature_extra( num,num2,subset_num):
    if num == 1:
        # Bag of Word model
        #f = "data/data.txt"
        #g = "hydrophobicity.txt"
        #data = data_generator.read_data(f)
        #hydrophobicity_len = data_generator.read_data(g)
        #peptide = data.data[0]
        #rt = data.data[2]
        #hydrophobicity = hydrophobicity_len.data[2]
        #length = hydrophobicity_len.data[1]
        #feature_temp = BOW(peptide, hydrophobicity,length,num2)
        #feature = feature_temp[0:subset_num]
        #rt = rt[0:subset_num]
        psmDescriptions,feature_temp = bf.getBOWfeature('data/retention_time_peptide.csv',num2)
        feature = feature_temp[0:subset_num]
        peptide,rt = data_generator.extract(psmDescriptions,subset_num) # grab 1000 as subset to test
    elif num == 2:
        #execute elude feature / modify data
        psmDescriptions, feature_temp = el.getFeatures('data/modify.csv')
        feature = feature_temp[0:subset_num]
        peptide,rt = data_generator.extract(psmDescriptions,subset_num) # grab 1000 as subset to test
    else:
        print('Warning : Unknown type of method !!')

    return feature,rt




