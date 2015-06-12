__author__ = 'Administrator'
import numpy as np
import data_generator
import sys
sys.path.append('elude/')
import elude_features as el
import BOW_feature as bf
import math

def feature_extra(File, num,num2,subset_num):
    num = int(num)
    num2 = int(num2)
    subset_num = int(subset_num)

    if num == 1:
        # Bag of Word model
        psmDescriptions,feature_temp = bf.getBOWfeature(File,num2)
        if subset_num != 0:
            feature = feature_temp[0:subset_num]
            peptide,rt = data_generator.extract(psmDescriptions,subset_num)
        else:
            subset_num = len(feature_temp)
            feature = feature_temp[0:subset_num]
            peptide,rt = data_generator.extract(psmDescriptions,subset_num)
    elif num == 2:
        #execute elude feature / modify data
        psmDescriptions, feature_temp = el.getFeatures(File)
        if subset_num != 0:
            feature = feature_temp[0:subset_num]
            peptide,rt = data_generator.extract(psmDescriptions,subset_num)
        else:
            subset_num = len(feature_temp)
            feature = feature_temp[0:subset_num]
            peptide,rt = data_generator.extract(psmDescriptions,subset_num)
    else:
        print('Warning : Unknown type of method !!')

    return feature,rt




