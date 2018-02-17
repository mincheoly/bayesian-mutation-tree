"""
This class will generate a tree when given the maximum likelihood data between
every mutation pair. Using the ML data it first creates a tree and then the
minimum spanning tree will be found.

Authors: Christian Choe, Min Cheol Kim
Create Date: 5/19/2016
"""

import collections
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.tree.branchings import Edmonds
import numpy as np
import random as rand

class tree_generation:
	def __init__(self, ML_dictionary):
		self.dict = ML_dictionary
		self.init_graph()
		self.find_min_span_tree()

	#Opens the excel file and saves the two worksheets
	def init_graph(self):
		self.tree = nx.DiGraph()

		#keys are gene pairs (x,y) and values are the likelihood of mutation relationships
		# in the order of (x->y), (y->x), (x><y)
		for key, value in self.dict.iteritems():
			x = key[0]
			y = key[1]
			self.tree.add_node(x, pos = (rand.uniform(0,5),rand.uniform(0,5)))
			self.tree.add_node(y, pos = (rand.uniform(0,5),rand.uniform(0,5)))

			#if (x><y) draw no edge because they are not related
			if (value[2] < value[1]) and (value[2] < value[0]):
				continue

			#take the -log of the probabiliy to get the weights
			self.tree.add_edge(x, y, weight = value[0])
			self.tree.add_edge(y, x, weight = value[1])

	def find_min_span_tree(self):
		edmond_tree = Edmonds(self.tree)
		self.min_tree = edmond_tree.find_optimum(attr = 'weight', default = 1, kind = 'min', style = 'arborescence')

	def get_tree(self):
		return self.tree

	def get_min_tree(self):
		return self.min_tree

	def print_edge_and_value(self, tree):
		for edge in tree.edges():
			print edge
			print tree[edge[0]][edge[1]]['weight']

	def print_tree(self, tree):
		pos = nx.get_node_attributes(tree, 'pos')
		nx.draw_networkx(tree)
		labels = nx.get_edge_attributes(tree, 'weight')
		print labels
		print pos
		nx.draw_networkx_edge_labels(tree, pos, edge_labels = labels)
		plt.show()	