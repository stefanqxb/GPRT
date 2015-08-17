#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as pp
pp.ion()
from joblib import Parallel, delayed  
import multiprocessing

def get_ith_descriptor( i,p, voc ):
	print i
	return p.wbow_descriptor( voc )


class word:
	def __init__( self, sequence, seq_indices ):
		self.char_seq = sequence
		self.ind_seq = seq_indices
	def __eq__(self,other):
		return self.char_seq, self.char_seq == other.char_seq, other.char_seq
	def __hash__(self):
		return hash( self.char_seq )

class vocabulary:
	def __init__( self, words ):
		self.char_seqs = [];
		self.ind_seqs = [];
		for w in words:
			self.char_seqs.append( w.char_seq )
			self.ind_seqs.append( w.ind_seq )
		self.nwords = len( self.char_seqs )
		self.k = len(self.ind_seqs[0])
	def char_seq_index( self, seq ):
		return self.char_seqs.index( seq )
	def ind_seq_score( self, inds ):
		v = np.array( [0.0] * self.nwords )
		for i in range( len( self.ind_seqs ) ):
			count = 0.0
			for a1,a2 in zip( inds, self.ind_seqs[i] ):
				if a1 == a2 :
					count = count + 1.0
			v[i] = count / len( inds )
		return v

class peptide:
	def __init__( self, sequence, rt ):
		self.pre = sequence[0]
		self.post = sequence[-1]
		self.sequence = sequence[2:-2]
		self.rt = rt
		self.amino_acids = self.get_amino_acids()
		self.aa_indices = [];

	def build_amino_acid_indices( self, aa_list ):
		self.aa_indices = [];
		for aa in self.amino_acids :
	   		self.aa_indices.append( aa_list.index(aa) )	

	def get_amino_acids( self ):
		ll = []
		seq = self.sequence + '*' ;
		i = 0
		while i < len(seq)-1 :
			if seq[i+1] == '[' :
				v = seq[i:i+5]
				i = i+5
			else :
				v = seq[i]
				i = i+1
			ll.append( v )
		return ll

	def get_k_mers( self, k ):
		k_mers = [];
		for i in range( len(self.aa_indices)-k+1 ):
			indices = self.aa_indices[i:i+k]
			seq = ''.join( self.amino_acids[i:i+k] )
			w = word( seq, indices )
			k_mers.append(w)
		return k_mers

	def bow_descriptor( self, voc ):
		desc = np.array([0.0] * ( voc.nwords ))
		k_mers = self.get_k_mers( voc.k )
		for w in k_mers : 
			ind = voc.char_seq_index( w.char_seq )
			desc[ ind ] = desc[ ind ] + 1.0
		desc = desc / (np.linalg.norm( desc ) + 1e-12)
		return desc	

	def sim_score( self, str1, str2 ):
		score = 0.0
		for s1, s2 in zip( str1, str2 ):
			if s1 == s2 :
				score = score + 1.0	
		return score / len( str1 )

	def wbow_descriptor( self, voc ):
		desc = np.array([0.0] * ( voc.nwords ))
		k_mers = self.get_k_mers(voc.k)
		for w in k_mers :
			v = voc.ind_seq_score( w.ind_seq )
			desc = desc + v
		desc = desc / (np.linalg.norm( desc ) + 1e-12)
		return desc

	def elute_descriptor( self ):
			

class feature_extractor:
	def __init__( self, peptides ):
		self.peptides = peptides;
		self.amino_list = self.amino_acid_list()
		for p in peptides :
			p.build_amino_acid_indices( self.amino_list )

	def amino_acid_list( self ):
		amino_list = set()
		for p in self.peptides :
			for a in p.amino_acids :
				amino_list.add( a )
		return list(amino_list)

	def bow_voc( self, k ):
		words = set()
		for p in self.peptides : 
			k_mers = p.get_k_mers(k)
			for s in k_mers :
				words.add( s )
		words = list( words )
		voc = vocabulary( words )
		return voc

	def elute_data( self ):


	def gen_features_parallel( self, voc ):
		num_cores = multiprocessing.cpu_count();
		res = Parallel( n_jobs=num_cores )( delayed(get_ith_descriptor)(i,self.peptides[i],voc) for i in range(5000) )
