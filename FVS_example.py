'''
This file contains the example from thesis section 3.3
written by Daniel Plaugher 9/10/22
'''

import networkx as nx
import random
import FVS



#A fixed graph example
G=nx.DiGraph()
G.add_edges_from([('x1', 'x2'),
('x2', 'x3'),
('x3', 'x4'), 
('x4', 'x2')])


#calculate FVS
G_FVS=FVS.FVS(G)
#print(G.edges())
print(G_FVS)
