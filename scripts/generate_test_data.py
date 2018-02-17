"""
This class implements the generation of "fake" test data, to see how well our algorithm performs under random conditions

This module supports introducing errors as of 06/05/2016.

Authors: Christian Choe, Min Cheol Kim

Create Date: 5/27/2016
"""

import networkx as nx
import random as rand
import numpy as np
import scipy
import collections
import time
from networkx.algorithms.traversal.depth_first_search import dfs_tree


class data_generator:

	# Constructor
	def __init__(self, num_cells, num_muts, alpha=0.92, withError=False, allelic_dropout=0.4309, false_discovery_rate=6.04e-5):

		# Simulation parameters
		self.num_cells = num_cells
		self.num_muts = num_muts
		self.withError = withError
		self.ad = allelic_dropout
		self.fd = false_discovery_rate

		# Tree variables
		self.tree = nx.DiGraph()
		self.Tk = []
		self.alpha = alpha
		self.num_nodes = sum(range(1, num_cells+1)) + num_cells

		# pre-cache the nodes that are in each level
		self.lvl_nodes_lookup = {0:[0]}
		start_idx = 1
		end_idx = 1

		for lvl in range(1, num_cells+2):
		    self.lvl_nodes_lookup[lvl] = range(start_idx, end_idx+1) if lvl != num_cells+1 else range(start_idx, end_idx)
		    start_idx = start_idx + lvl
		    end_idx = end_idx + lvl + 1

	# This method generates Tk's, using exponential distributions and rate of (k choose 2)
	def initialize_Tk(self):
		temp = [np.random.exponential(1/scipy.misc.comb(k, 2)) for k in range(2,self.num_cells+1)]
		temp.insert(0, sum(temp)*self.alpha/(1-self.alpha))
		self.Tk = temp
		self.cumsum = np.cumsum(temp)

	# Create the lineage tree that will host ALL the mutations (unlike the simuation, where it only hosted two mutations)
	def create_lineage_tree(self):
		self.tree.clear()
		self.tree.add_node(0) #this is a fake node
		self.tree.add_edge(0, 1)
		self.tree.add_node(1)
		next_node_number = 2
		current_nodes = [1]
		for level in range(1,self.num_cells+1):
			next_nodes = []
			node_to_split = -1 if level == self.num_cells else rand.choice(current_nodes)
			for node in current_nodes:
				for i in range(1 + int(node==node_to_split)):
					self.add_child(node, next_node_number, next_nodes)
					next_node_number += 1
			current_nodes = next_nodes
		self.leaves = self.tree.nodes()[-self.num_cells:]
		self.leave_genotypes = {}
		for leaf in self.leaves:
			self.leave_genotypes[leaf] = [0 for idx in range(self.num_muts)]
		# Pre-cache the leaves from each node
		self.node_leaves_lookup = {}
		for node in range(1, self.num_nodes + 1):
			subtree_nodes = dfs_tree(self.tree, node)
			self.node_leaves_lookup[node] = [leaf for leaf in subtree_nodes if self.tree.out_degree(leaf) == 0]

	# Helper method for putting in a child node in the tree
	def add_child(self, node, next_node_number, next_nodes):
		self.tree.add_node(next_node_number)
		self.tree.add_edge(node, next_node_number)
		next_nodes.append(next_node_number)

	# generate the random times for all mutations. If level n is returned, it means that
	# the mutation occurs at the edge beween level n and n + 1 (tree is labeled starting at 1 as root)
	def generate_mutation_levels(self):
		mut_times = [np.random.uniform(low=0.0, high=sum(self.Tk)) for mut in range(self.num_muts)]
		self.mut_times = mut_times
		mut_lvls = [np.argmax(self.cumsum > mut_times[mut_idx]) for mut_idx in range(self.num_muts)]
		return sorted(mut_lvls)

	# Apply and propagate the mutation to the leaves, where the information is stored in self.leave_genotypes
	def apply_mutations(self):
		mut_levels = self.generate_mutation_levels()
		nodes_at_lvls = [self.lvl_nodes_lookup[mut_levels[mut_idx]] for mut_idx in range(self.num_muts)]
		self.mutation_locations = []
		for mut_idx in range(len(nodes_at_lvls)):
			nodes = nodes_at_lvls[mut_idx]
			edges = []
			for node in nodes:
				edges.extend(self.tree.edges(node))
			choice = rand.choice(edges)
			self.mutation_locations.append(choice)
			mut_root = choice[1]
			leaves_affected = self.node_leaves_lookup[mut_root]
			for leaf in leaves_affected:
				self.leave_genotypes[leaf][mut_idx] = 1

	# returns the genotype data as a dictionary of lists. Now supports introducing errors in the data.
	def return_genotype_data(self):
		temp = collections.OrderedDict()
		for mut_idx in range(self.num_muts):
			sitename = 'gene' + str(mut_idx)
			site_genotype = []
			for leaf in self.leaves:
				genotype = self.leave_genotypes[leaf][mut_idx]
				if self.withError: # let's introduce some errors
					if genotype == 1:
						if rand.random() < self.ad:
							genotype = 0
					else: # genotype was one
						if rand.random() < self.fd:
							genotype = 1
				site_genotype.append( genotype )
			temp[sitename] = site_genotype
		return temp

	def construct_true_mut_tree(self):
		muttree = nx.DiGraph()
		for idx, edge in enumerate(self.mutation_locations):
			parent = edge[0]
			child = edge[1]
			if idx == 0:
				muttree.add_node('gene' + str(idx))
			else:
				muttree.add_node('gene' + str(idx))
				indices = range(idx)
				for jdx in reversed(indices):
					subtree = dfs_tree(self.tree, self.mutation_locations[jdx][1])
					if parent in subtree:
						muttree.add_edge('gene'+str(jdx), 'gene'+str(idx))
						break
					elif self.mutation_locations[jdx][1] == child:
						muttree.add_edge('gene'+str(jdx), 'gene'+str(idx))
						break
		return muttree

	def get_ambiguities(self):
		track_ambiguities = {}
		# initialize
		for gene_num, edge in enumerate(self.mutation_locations):
			track_ambiguities[edge] = []
		for gene_num, edge in enumerate(self.mutation_locations):
			track_ambiguities[edge].append('gene' + str(gene_num))
		return track_ambiguities.values()







