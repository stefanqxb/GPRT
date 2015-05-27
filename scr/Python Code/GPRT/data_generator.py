__author__ = 'Administrator'
import numpy as np
import sys
sys.path.append('elude/')
import data_manager

class function_group:
        def __init__( self, keys):
                self.keys = keys
                self.keys_dict = dict()
                for i in range( len(keys) ):
                        self.keys_dict[ keys[i] ] = i

                self.data = []
                for i in range( len(keys) ):
                        self.data.append( [] )
                self.ndata = 0

        def add_data( self, values ):
                assert len( values ) == len( self.keys )
                self.ndata = self.ndata + 1
                for i in range( len(self.keys) ):
                        dump_ind = self.keys_dict[ self.keys[i] ]
                        self.data[dump_ind].append( values[i] )

        def get_keys( self ):
                return self.keys

        def get_n_data( self ):
                return self.ndata

        def get_values_for_key( self, key, unique=0 ):
                full_set = self.data[ self.keys_dict[ key ] ]
                if unique == 1:
                        uni_set = set()
                        for i in range( len(full_set) ):
                                uni_set.add( full_set[i] )
                        return list(uni_set)
                else:
                        return full_set

        def get_sub_dataset( self, key, key_value ):
                sub_data = grab( self.keys )
                key_ind = self.keys_dict[ key ]
                for i in range( self.ndata ):
                        if self.data[ key_ind ][i] == key_value:
                                row = []
                                for j in range( len(self.keys) ):
                                        row.append( self.data[j][i] )
                                sub_data.add_data( list(row) )
                return sub_data

        def split_by_key(self,  key):
            uni_key_values = self.get_values_for_key(key,  1)
            return [ self.get_sub_dataset(key,  key_value) for key_value in uni_key_values ]

        def get_indices(self,  key,  key_value):
            ind = []
            key_ind = self.keys_dict[key]
            for i in range( self.ndata ):
                if self.data[ key_ind ][i] == key_value:
                    ind.append(i)
            return ind

def read_data(path):
    with open(path,'r') as f:
        keys = f.readline().split('\t')
        for i, key in enumerate(keys):
            keys[i] = grab(key)
        data = function_group(keys)
        row = f.readline()
        while len(row)> 0:
            row =  row.split('\t')
            for i , string in enumerate(row):
                row[i] = grab(string)
            data.add_data(row)
            row = f.readline()
    return data

def grab(string):
    string_new = ''
    num = '0123456789'
    g = num.find(string[0])
    if g == -1:
       for s in string:
           if (s!='\r') and (s!='\n') and (s!='.'):
                    string_new += s
    else:
        for s in string:
           if (s!='\r') and (s!='\n'):
                    string_new += s

    return string_new

def processing_Data(feature,rt,row,str):
    ratio = 0.8
    train_row = int(round(ratio * row))
    end = row #len(feature)
    if str =='svr':
       train_tag_temp = np.zeros([train_row])
       test_tag_temp = np.zeros([end-train_row])
    elif str =='gp':
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
    return train_set,train_tag,test_set,test_tag

def extract(inst,num):
    peptide =[]
    rt = []
    for item in inst:
        if isinstance(item,data_manager.PSMDescription):
            peptide.append(item.peptide)
            rt.append(item.retentionTime)
            if len(rt) >= num:
                break
    return peptide,rt