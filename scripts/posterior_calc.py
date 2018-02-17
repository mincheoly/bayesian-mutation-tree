"""
This class implements the calculation of the posteriors from the given data and priors. 
The algorithm is discussed in Kim and Simon's paper also discussed in our submission. 

As of now, it does not account for allelic dropout rate and the false discovery rate

Authors: Christian Choe, Min Cheol Kim

Create Date: 5/27/2016
"""

import numpy as np
import scipy
import collections
from math import log

class posterior_calculator:

	# Constructor for the posterior calculator.
	# Inputs:
	# 	data - OrderedDict where keys are gene names('PDE4DIP(A->G)') and values are a list of genotypes at all 58 cells. 
	#	priors - OrderedDict where keys are (i, j, 'x~y') and the value is the probability
	def __init__(self, data, priors, n):
		self.data = data
		self.priors = priors
		self.posteriors = {}
		self.n = n

	# Function for calculating the likelihood (equation #1 in paper)
	def calculate_likelihood(self):
		data_keys = self.data.keys()
		data_values = self.data.values()
		for idx1 in range(len(self.data)): # For each pairing
			for idx2 in range(idx1+1, len(self.data)):
				gene1, genotype1 = data_keys[idx1], data_values[idx1]
				gene2, genotype2 = data_keys[idx2], data_values[idx2]
				assert(len(genotype1) == self.n) # make sure the length is the same as # of cells
				self.posteriors[(gene1, gene2)] = [0, 0, 0]
				skipped_cells = 0

				# keep track of probabilities in logs because they're too small
				for cell in range(self.n):
					if genotype1[cell] == '-' or genotype2[cell] == '-':
						skipped_cells += 1
						continue
					allele1 = int(genotype1[cell]) if genotype1[cell] != 2 else 1
					allele2 = int(genotype2[cell]) if genotype2[cell] != 2 else 1
					self.posteriors[(gene1, gene2)][0] = self.posteriors[(gene1, gene2)][0] + log(self.priors[(int(allele1) , int(allele2), 'x->y')])
					self.posteriors[(gene1, gene2)][1] = self.posteriors[(gene1, gene2)][1] + log(self.priors[(int(allele1) , int(allele2), 'x<-y')])
					self.posteriors[(gene1, gene2)][2] = self.posteriors[(gene1, gene2)][2] + log(self.priors[(int(allele1) , int(allele2), 'x><y')])
				self.posteriors[(gene1, gene2)][0] = self.posteriors[(gene1, gene2)][0]/float(self.n-skipped_cells)
				self.posteriors[(gene1, gene2)][1] = self.posteriors[(gene1, gene2)][1]/float(self.n-skipped_cells)
				self.posteriors[(gene1, gene2)][2] = self.posteriors[(gene1, gene2)][2]/float(self.n-skipped_cells)
		

	# Function for finding the posterior (equation #3 in paper)
	def calculate_posteriors(self):
		for pair, likelihoods in self.posteriors.iteritems():
			self.posteriors[pair][0] = -1*(self.posteriors[pair][0] + log(self.priors['x->y']))
			self.posteriors[pair][1] = -1*(self.posteriors[pair][1] + log(self.priors['x<-y']))
			self.posteriors[pair][2] = -1*(self.posteriors[pair][2] + log(self.priors['x><y']))

	# Getter function for the posteriors
	def return_posteriors(self):
		return self.posteriors


					