"""
This class implements the simulation to calculate priors. The algorithm is discussed in Kim and Simon's
paper also discussed in our submission. 

Authors: Christian Choe, Min Cheol Kim
Create Date: 5/19/2016
"""

import networkx as nx
import random as rand
import numpy as np
import scipy
import collections
from networkx.algorithms.traversal.depth_first_search import dfs_tree
import time

class prior_calculator:
	def __init__(self, n, alpha, Btree=1, Bmut=1, withError=False, allelic_dropout=0.4309, false_discovery_rate=6.04e-5):
		# measurement variables
		self.x_to_y = collections.Counter()
		self.y_to_x = collections.Counter()
		self.x_not_y = collections.Counter()
		self.Bmut = Bmut
		self.Btree = Btree

		# error table lookup
		self.errors = [[1-6.04e-5, 0.21545], [6.04e-5, 0.5691]]

		# tree structure variables
		self.n = n
		self.num_nodes = sum(range(1, n+1)) + n
		self.alpha = alpha
		self.b_tree = Btree
		self.b_mut = Bmut
		self.tree = nx.DiGraph()
		self.Tk = []
		self.withError = withError
		self.ad = allelic_dropout
		self.fd = false_discovery_rate

		# pre-cache the nodes that are in each level
		self.lvl_nodes_lookup = {0:[0]}
		start_idx = 1
		end_idx = 1

		for lvl in range(1, n+2):
		    self.lvl_nodes_lookup[lvl] = range(start_idx, end_idx+1) if lvl != n+1 else range(start_idx, end_idx)
		    start_idx = start_idx + lvl
		    end_idx = end_idx + lvl + 1

	# This method performs the simulation with previously given parameters
	def simulate(self):	
		#start = time.time()
		for tree_trial in range(self.Btree):
			#if tree_trial%10 == 0:
				#print time.time()-start
			self.initialize_Tk()
			self.generate_new_tree()
			for mut_trial in range(self.Bmut):
				leaves = self.mutate_tree()
				if (1, 1) not in leaves:
					self.x_not_y.update(leaves)
				elif (1, 0) in leaves:
					self.x_to_y.update(leaves)
				elif (0, 1) in leaves:
					self.y_to_x.update(leaves)
				else:
					choice = rand.choice([1, 2])
					if choice == 1:
						self.x_to_y.update(leaves)
					else:
						self.y_to_x.update(leaves)

	# This function returns the final priors
	def return_priors(self):
		return self.pre_error_probs

	# This function computes priors after the simulations are finished.
	def compute_final_probs(self):
		self.compute_pre_error_probs()
		self.compute_prior_stats()

	# This method computes pre-error probabilities
	def compute_pre_error_probs(self):
		# compute the probability of genotypes given mutation order
		self.pre_error_probs = {}
		for genotype, count in self.x_to_y.iteritems():
			self.pre_error_probs[(genotype[0], genotype[1],'x->y')] = float(count)/sum(self.x_to_y.values())
		for genotype, count in self.y_to_x.iteritems():
			self.pre_error_probs[(genotype[0], genotype[1],'x<-y')] = float(count)/sum(self.y_to_x.values())
		for genotype, count in self.x_not_y.iteritems():
			self.pre_error_probs[(genotype[0], genotype[1],'x><y')] = float(count)/sum(self.x_not_y.values())

		total = sum(self.x_to_y.values()) + sum(self.y_to_x.values()) + sum(self.x_not_y.values())

		# Compute the probability of mutation orders
		self.pre_error_probs['x->y'] = sum(self.x_to_y.values())/float(total)
		self.pre_error_probs['x<-y'] = sum(self.y_to_x.values())/float(total)
		self.pre_error_probs['x><y'] = sum(self.x_not_y.values())/float(total)

		# Compute in the 'impossible' probabilities
		self.pre_error_probs[(1, 0, 'x<-y')] = 1e-300
		self.pre_error_probs[(0, 1, 'x->y')] = 1e-300
		self.pre_error_probs[(1, 1, 'x><y')] = 1e-300

	# This method generates Tk's, using exponential distributions and rate of (k choose 2)
	def initialize_Tk(self):
		temp = [np.random.exponential(1/scipy.misc.comb(k, 2)) for k in range(2,self.n+1)]
		temp.insert(0, sum(temp)*self.alpha/(1-self.alpha))
		self.Tk = temp
		self.cumsum = np.cumsum(temp)

	# This method generates a new tree scaffold
	def generate_new_tree(self):
		self.tree.clear()
		self.tree.add_node(0, x=0, y=0) #this is a fake node
		self.tree.add_edge(0, 1, x=0, y=0)
		self.tree.add_node(1, x=0, y=0)
		next_node_number = 2
		current_nodes = [1]
		for level in range(1,self.n+1):
			next_nodes = []
			node_to_split = -1 if level == self.n else rand.choice(current_nodes)
			for node in current_nodes:
				for i in range(1 + int(node==node_to_split)):
					self.add_child(node, next_node_number, next_nodes)
					next_node_number += 1
			current_nodes = next_nodes
		self.leaves = self.tree.nodes()[-self.n:]

		# Pre-cache the leaves from each node
		self.node_leaves_lookup = {}
		for node in range(1, self.num_nodes + 1):
			subtree_nodes = dfs_tree(self.tree, node)
			self.node_leaves_lookup[node] = [leaf for leaf in subtree_nodes if self.tree.out_degree(leaf) == 0]

	# Helper method for putting in a child node in the tree
	def add_child(self, node, next_node_number, next_nodes):
		self.tree.add_node(next_node_number, x=0, y=0)
		self.tree.add_edge(node, next_node_number, x=0, y=0)
		next_nodes.append(next_node_number)

	# generate the random times for mutations x and y. If level n is returned, it means that
	# the mutation occurs at the edge beween level n and n + 1 (tree is labeled starting at 1 as root)
	def generate_mutation_levels(self):
		x = np.random.uniform(low=0.0, high=sum(self.Tk))
		y = np.random.uniform(low=0.0, high=sum(self.Tk))
		x_lvl = np.argmax(self.cumsum > x)
		y_lvl = np.argmax(self.cumsum > y)
		return (x_lvl, y_lvl)

	# This is the method for mutating the tree
	def mutate_tree(self):
		mutated_edges = self.apply_mutation()
		mutated_nodes = self.propagate_mutation(mutated_edges)
		leaves = self.get_leaves()
		self.clean_tree(mutated_edges, mutated_nodes)
		return leaves

	# This helper method introduces the mutations to the tree
	def apply_mutation(self):
		(x_lvl, y_lvl) = self.generate_mutation_levels()
		nodes_at_lvl_x = self.lvl_nodes_lookup[x_lvl]
		nodes_at_lvl_y = self.lvl_nodes_lookup[y_lvl]
		y_edges = []
		x_edges = []
		for node in nodes_at_lvl_x:
			x_edges.extend(self.tree.edges(node))
		for node in nodes_at_lvl_y:
			y_edges.extend(self.tree.edges(node))
		x_choice = rand.choice(x_edges);
		y_choice = rand.choice(y_edges);
		self.tree.edge[x_choice[0]][x_choice[1]]['x'] = 1
		self.tree.edge[y_choice[0]][y_choice[1]]['y'] = 1
		return [x_choice, y_choice]

	# This helper method propagates the mutated edge all the way to the leaves
	def propagate_mutation(self, mutated_edges):
		mutated_nodes = []
		mut_root_x = mutated_edges[0][1]
		mut_root_y = mutated_edges[1][1]
		for node in self.node_leaves_lookup[mut_root_x]:
			self.tree.node[node]['x'] = 1
			mutated_nodes.append(node)
		for node in self.node_leaves_lookup[mut_root_y]:
			self.tree.node[node]['y'] = 1
			mutated_nodes.append(node)
		return set(mutated_nodes)

	# Returns the leaf genotypes of the current tree, mutate the leaves if withError is True
	def get_leaves(self):
		if not self.withError:
			return [ (self.tree.node[leaf]['x'], self.tree.node[leaf]['y']) for leaf in self.node_leaves_lookup[1] ]
		else: #mutate the leaves before returning
			temp = []
			for leaf in self.node_leaves_lookup[1]:
				x_value = self.tree.node[leaf]['x']
				y_value = self.tree.node[leaf]['y']
				if x_value == 1:
					if rand.uniform(0,1) < self.ad:
						x_value = 0
				else:
					if rand.uniform(0,1) < self.fd:
						x_value = 1
				if y_value == 1:
					if rand.uniform(0,1) < self.ad:
						y_value = 0
				else:
					if rand.uniform(0,1) < self.fd:
						y_value = 1
				temp.append((x_value, y_value))
			return temp				

	# "Cleans" the tree; set all x,y values to 0
	def clean_tree(self, mutated_edges, mutated_nodes):
		for edge in mutated_edges:
			self.tree.edge[edge[0]][edge[1]]['x'] = 0
			self.tree.edge[edge[0]][edge[1]]['y'] = 0
		for node in mutated_nodes:
			self.tree.node[node]['x'] = 0
			self.tree.node[node]['y'] = 0

	# This method incorporates the errors into the probabilities:
	def compute_prior_stats(self):
		self.prior_stats = collections.OrderedDict()
		# compute x->y stats
		for x, y in [(0, 0), (1, 0), (1, 1)]:
			self.prior_stats[(x, y, 'x->y')] = self.errors[x][0]*self.errors[y][0]*self.pre_error_probs[(0, 0, 'x->y')] + \
											   self.errors[x][1]*self.errors[y][0]*self.pre_error_probs[(1, 0, 'x->y')] + \
											   self.errors[x][1]*self.errors[y][1]*self.pre_error_probs[(1, 1, 'x->y')]
		x_to_y_total = sum([value for key, value in self.prior_stats.iteritems() if key[2] == 'x->y'])
		for x, y in [(0, 0), (1, 0), (1, 1)]:
			self.prior_stats[(x, y, 'x->y')] = self.prior_stats[(x, y, 'x->y')]/x_to_y_total

		# compute x<-y stats
		for x, y in [(0, 0), (0, 1), (1, 1)]:
			self.prior_stats[(x, y, 'x<-y')] = self.errors[x][0]*self.errors[y][0]*self.pre_error_probs[(0, 0, 'x<-y')] + \
											   self.errors[x][0]*self.errors[y][1]*self.pre_error_probs[(0, 1, 'x<-y')] + \
											   self.errors[x][1]*self.errors[y][1]*self.pre_error_probs[(1, 1, 'x<-y')]
		y_to_x_total = sum([value for key, value in self.prior_stats.iteritems() if key[2] == 'x<-y'])
		for x, y in [(0, 0), (0, 1), (1, 1)]:
			self.prior_stats[(x, y, 'x<-y')] = self.prior_stats[(x, y, 'x<-y')]/y_to_x_total

		# compute x><y stats
		for x, y in [(0, 0), (0, 1), (1, 0)]:
			self.prior_stats[(x, y, 'x><y')] = self.errors[x][0]*self.errors[y][0]*self.pre_error_probs[(0, 0, 'x><y')] + \
											   self.errors[x][0]*self.errors[y][1]*self.pre_error_probs[(0, 1, 'x><y')] + \
											   self.errors[x][1]*self.errors[y][0]*self.pre_error_probs[(1, 0, 'x><y')]
		x_not_y_total = sum([value for key, value in self.prior_stats.iteritems() if key[2] == 'x><y'])
		for x, y in [(0, 0), (0, 1), (1, 0)]:
			self.prior_stats[(x, y, 'x><y')] = self.prior_stats[(x, y, 'x><y')]/x_not_y_total











