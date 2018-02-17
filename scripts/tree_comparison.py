"""
This class will be used to compare two trees. It will give a metric ~[0,1]
where 1 means the two graphs are identical (with regard to ambiguities). Our
metric compares every pair of nodes and checks the relationships in both trees.
(e.g. checks to make sure that node_1 is in the same lineage as node_2 as well
as node_1 being above node_2, etc.)

#first check if parents and children order is good
#then check if seperate lineages are good

Authors: Christian Choe, Min Cheol Kim
Create Date: 6/4/2016
"""

import itertools
import networkx as nx
import numpy as np
from scipy.special import comb

#Given the node relationships it will return a list of tuples listing all the ambiguities
def list_all_ambiguity(relation): 
    ambi = []
    for relationship in relation:
        if len(relationship) > 0: #start choosing two
            for subset in itertools.combinations(relationship, 2):
                ambi.append(subset)
    return ambi

def type_of_node(G, node_1, node_2):
    successors = nx.dfs_successors(G)
    if (nx.has_path(G, source=node_1, target=node_2)): #node_1 comes first
        return 1
    elif (nx.has_path(G, source=node_2, target=node_1)): #node_2 comes first
        return 2
    else: #not same lineage
        return 0

#Min Cheol thought of a sick idea, so sick he coughed
#When comparing the derived tree and true tree I will also get a list of nodes that are considered interchangeable.
#When I check those two nodes then I know that they must be in the same lineage so '1' or '2' is fine.
#If I get a block of ambiguity ('a', 'b', 'c') then choose 2 and go compare those nodes
def compare_mutation_order(G1, G2, relation = []):
    ambi = list_all_ambiguity(relation)
    nodes = G1.nodes()
    count = 0
    for i, node_1 in enumerate(nodes[0:-1]):
        for node_2 in nodes[i+1:]:
            if ((node_1, node_2) in ambi) or ((node_2, node_1) in ambi): #the two nodes can go in either order
                count += ((type_of_node(G1, node_1, node_2) * type_of_node(G2, node_1, node_2)) > 0)*1
            else:
                count += (type_of_node(G1, node_1, node_2) == type_of_node(G2, node_1, node_2))*1
    return float(count)/comb(len(nodes),2)
