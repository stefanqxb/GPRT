__author__ = 'Administrator'
import numpy as np
import data_generator
import sys
sys.path.append('elude/')
import elude_features as el
import BOW_feature as bf
import random

def feature_extra(File, num,num2,subset_num):
    num = int(num)
    num2 = int(num2)
    subset_num = int(subset_num)
    peptide_final =[]
    rt_final =[]

    if num == 1:
        # Bag of Word model
        psmDescriptions,feature_temp = bf.getBOWfeature(File,num2)
    elif num == 2:
        #execute elude feature / modify data
        psmDescriptions, feature_temp = el.getFeatures(File)
    else:
        print('Warning : Unknown type of method !!')

    # shuffle the data

    total = range(len(feature_temp))
    random.shuffle(total)

    if subset_num != 0:
        seed = total[0:subset_num]
        feature = feature_temp[seed]
        peptide,rt = data_generator.extract(psmDescriptions,len(feature_temp))
    else:
        seed = total
        feature = feature_temp[total]
        peptide,rt = data_generator.extract(psmDescriptions,len(feature_temp))

    for item in seed:
        peptide_final.append(peptide[item])
        rt_final.append(rt[item])

    return peptide_final, feature,rt_final




