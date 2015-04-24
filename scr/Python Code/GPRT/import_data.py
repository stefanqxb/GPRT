__author__ = 'Administrator'
import numpy as np

class table:
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
                sub_data = table( self.keys )
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

def read_data( data_path):
    with open( data_path, 'r' ) as f :
        keys = f.readline().split('\t');
        for i,  key in enumerate(keys):
            keys[i] = del_quote(key)
        data = table(keys)
        row = f.readline()
        while len(row)>0:
            row = row.split('\t')
            for i,  string in enumerate(row):
                row[i] = del_quote(string)
            data.add_data(row)
            row = f.readline()
    return data

def del_quote(string):
    string_new = ''
    for s in string:
        if ( s!='"') and (s!='\r') and (s!='\n'):
            string_new += s
    return string_new