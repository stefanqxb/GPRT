#!/usr/bin/pyhon

import platform
import os

class parameters:
    def __init__( self ):
        if platform.system() == 'Darwin' :
            self.data_root = "/Users/heydar/Stuff/tmp/gprt"
        else :
            self.data_root = "/media/hdd/heydar/data/gprt"

        #self.data_path = './data/20121212_S25_3ug_300min_LCMS_PM3_Tryp_GRADIENT_15sec_MZ.rtimes.tsv'
        #self.data_name = 'no_modif'

        self.data_path = './data/20120126_EXQ5_KiSh_SA_LabelFree_HeLa_Phospho_Control_rep1_Fr3.rtimes.tsv'
        self.data_name = 'ptm_modif'

        self.models_name = 'original'
        self.models_tag = '%s_%s' % ( self.models_name, self.data_name )
        self.models_path = '%s/%s' % ( self.data_root, self.models_tag )

        if not os.path.isdir( self.models_path ) : 
            os.mkdir( self.models_path )

        self.save_tmp = '%s/%s/models_ntrain_%d.pk'
        self.ntrain = [ 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 ];
        self.nparts = 10
        self.train_ratio = 0.5