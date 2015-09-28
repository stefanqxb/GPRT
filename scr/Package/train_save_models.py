#!/usr/bin/python

import platform
import numpy as np
import data_tools
import ml_tools
import pickle as pk

if __name__ == "__main__":
    if platform.system() == 'Darwin' :
        data_root = "/Users/heydar/Stuff/tmp/gprt"
    else :
        data_root = "/media/hdd/heydar/data/gprt"
    peptides = data_tools.read_data()
    models_name = 'zero_mean'
    ntrain = [100, 200, 300, 500, 1000, 1500, 2000, 2500, 3000]
    for n in ntrain : 
        benchmark = ml_tools.rt_benchmark(peptides, 'elude', 'gp', n, 10)
        models = ml_tools.single_train_gp( benchmark )
        save_path = "%s/%s/models_ntrain_%d.pk" % ( data_root, models_name, n )
        print save_path
        with open( save_path, 'w' ) as ff :
            pk.dump( [ models ], ff )
            ff.close()
        models = None 
