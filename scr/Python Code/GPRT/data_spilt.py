__author__ = 'Administrator'
import sys
sys.path.append('elude/')
import elude_features as el
import data_generator

def split_data(File):
    pst,aa  = el.processTrainData(File)
    peptide,rt  = data_generator.extract(pst,len(pst))
    test_pep = peptide[0:1000]
    test_rt = rt[0:1000]
    data_generator.seal_data(test_pep,test_rt,0)

    peptide_train_tot = peptide[1001:]
    rt_train_tot = rt[1001:]
    train_size = [100,200,500,1000,2000,5000,8000,10000,12000]

    for i in range(9):
        data_generator.seal_data(peptide_train_tot,rt_train_tot,train_size[i])
