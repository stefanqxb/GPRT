#!/usr/bin/python

import numpy as np

class rt_trainer :
	def __init__( self, peptides ):
		self.peptides = peptides
		self.targes = [];
		self.inputs [];

		for p in self.peptides :
			self.targets.append( p.rt )
		self.targets = np.array( self.targets )

	def ext_bow( self, voc ):
		self.inputs = [];
		for p in self.peptides :
			self.inputs.append( p.bow_descriptor( voc ) )
		self.inputs = np.array( self.inputs )

	def ext_elute( self ):
		print "Extracting Elute"

	def train_gp( self ):
		print "Training GP"
