#!/usr/bin/regina-python

#  original (adjusted so that it works on my setup): #!/usr/local/bin/regina-python

###############################################################################
# 1. Usage for big (.bz2-zipped) files
# bzcat <file>.sig.bz2 | ./<pythonFile>.py
#
# 2. normal files
# cat <file>.sig | ./<pythonFile>.py
###############################################################################
#
# LOCATION of 1-vtx solid tori:
#
# /media/data/dim3census/torus-1vtx/n.sig 1 <= n <= 8
###############################################################################

# to deal with file issues
import os.path

# import the pickle module
import pickle

# to handle external commands
import subprocess

# for random choices
import random

# import all functions etc. from regina
from regina import *

# regular expressions
import re

# gcd for Smith normal form
from fractions import gcd

#import networkx as nx

import sys
sys.setrecursionlimit(100)

##########################################
##########################################
############ HELPER FUNCTIONS ############
##########################################
##########################################

### computes the boundary of all triangles
### and returns a list of 1-chains (1-cycles)
class FindCellFailure(Exception):
    pass

class FindLayerFailure(Exception):
    pass

class LogicalMistake(Exception):
    pass

class FacePoset:
    def __init__(self):
        self.layers ={}

    class PosetNode:
        def __init__(self, dim, name, cell,
                     parents = None, children = None):
            self.cell = cell
            self.name = name
            self.dim = dim

            if not parents:
                self.parents = []
            else:
                self.parents = parents
            if not children:
                self.children = []
            else:
                self.children = children

        def add_node_to_list(self, node, li):
            if not li:
                li = [node]
            else:
                li.append(node)
        
        def remove_node_from_list(self, node, li, error = False):
            if node in li:
                li.remove(node)
            elif error:
                raise CellFailure('No such cell found')

        def add_child(self, node):
            self.children.append(node)
        
        def add_parent(self, node):
            # figure out what to do with multi-edges
            self.parents.append(node)
        
        def remove_child(self, node):
            #self.remove_node_from_list(node, self.children)
            if node in self.children:
                self.children.remove(node)
            elif error:
                raise CellFailure('No such cell found')

        def remove_parent(self, node):
            self.remove_node_from_list(node, self.parents)
            if node in self.parents:
                self.parents.remove(node)
            elif error:
                raise CellFailure('No such cell found')
    
    def get_node(self, dim, name):
        try:
            return self.layers[dim][name]
        except:
            raise FindCellFailure('Could not find cell '+str(name)+' in layer '+str(dim)) 

    def add_node(self, dim, name, cell):
        if not dim in self.layers.keys():
            self.layers[dim] = {}
        self.layers[dim][name] = self.PosetNode(name = name, cell = cell, dim = dim)

    def remove_node(self, dim, node_label, suppress_error = False):
        if not dim in self.layers.keys():
            if suppress_error:
                return
            string = 'Could not find layer of dimension '+str(dim)+' when deleting node '+str(node_label)
            raise FindLayerFailure()
        else:
            node = self.layers[dim][node_label]
            for child in node.children:
                child.remove_parent(node)
            for parent in node.parents:
                parent.remove_child(node)
            _ = self.layers[dim].pop(node_label)
            del node

    def add_arc(self, n1_tup, n2_tup):
        dim1, name1 = n1_tup
        dim2, name2 = n2_tup
        
        for dim in dim1, dim2:
            if not dim in self.layers.keys():
                raise FindLayerFailure('Could not find layer '+str(dim)+' for creating arc between '+str(n1_tup, n2_tup))

        for dim, name in n1_tup, n2_tup:
            if not name in self.layers[dim].keys():
                raise FindCellFailure('Could not find node '+str(dim, name)+' for creating arc between '+str(n1_tup, n2_tup))
        
        if dim1 == dim2:
            raise LogicalMistake('cannot place arc between two nodes of the same dimension')

        if abs(dim1 - dim2) != 1:
            raise LogicalMistake('cannot place arc between two nodes which are not in adjacent layers')

        if dim1 < dim2:
            dim1, name1, dim2, name2 = dim2, name2, dim1, name1

        # so after this point, the (dim1, name1) is the node with the higher dimension

        n1 = self.layers[dim1][name1]
        n2=  self.layers[dim2][name2]

        #print n1.name, n1.cell, n1.dim, n1.children, n1.parents
        #print n2.name, n2.cell, n2.dim, n2.children, n2.parents

        n1.add_child(n2)
        n2.add_parent(n1)

        #print n1.name, n1.cell, n1.dim, n1.children, n1.parents
        #print n2.name, n2.cell, n2.dim, n2.children, n2.parents

    def remove_arc(self, tup1, tup2):
        dim1, name1 = n1_tup
        dim2, name2 = n2_tup
        
        for dim in dim1, dim2:
            if not dim in self.layers.keys():
                raise FindLayerFailure('Could not find layer '+str(dim)+' for creating arc between '+str(n1_tup, n2_tup))

        for dim, name in n1_tup, n2_tup:
            if not name in self.layers[dim].keys():
                raise FindCellFailure('Could not find node '+str(dim, name)+' for creating arc between '+str(n1_tup, n2_tup))
        
        if dim1 == dim2:
            raise LogicalMistake('arc cannot exist between two nodes of the same dimension')

        if abs(dim1 - dim2) != 1:
            raise LogicalMistake('arc cannot go between two nodes which are not in adjacent layers')

        if dim1 < dim2:
            dim1, name1, dim2, name2 = dim2, name2, dim1, name1

        # so after this point, the (dim1, name1) is the node with the higher dimension

        n1 = self.layers[dim1][name1]
        n2=  self.layers[dim2][name2]

        n1.remove_child(n2)
        
        n2.remove_parent(n1)

    def output(self):
        print "This depicts the face poset diagram from the 0th dimension cells upwards. Arcs are expressed downwards."
        for dim in sorted(self.layers.keys()):
            print 'Dim '+str(dim)+ ':'
            if dim == 0:
                for cell_name in self.layers[dim].keys():
                    print cell_name
            else:
                for cell_name in self.layers[dim].keys():
                    #print self.layers[dim][cell_name].children
                    string = ', '.join([str(child.name) for child in self.layers[dim][cell_name].children])
                    if string == '':
                        print str(cell_name)
                    else:
                        print str(cell_name) + ': ' + string


fp = FacePoset()
fp.add_node(0, 0, 'hello')
fp.add_node(0, 1, 'hi')
fp.add_node(1,0, 'higher')
fp.add_arc((0,0), (1,0))

#n1 = fp.get_node(1,0)
#print n1.cell
fp.output()




##########################################
##########################################
#########END HELPER FUNCTIONS ############
##########################################
##########################################

#for line in sys.stdin:
#  if ctr > 0: # max added to stop loop
#    break
#  ctr += 1
  #if ctr%1000 == 0:
  #  print ctr
  #sig = re.search('[a-zA-Z0-9]*' ,line)
  #t = NTriangulation.fromIsoSig(line)
#  t = Triangulation3.fromIsoSig(line)
  #print t.isoSig()
  #print t.detail()
  #print list([v for v in t.vertex_iterator()])


  #vertices = t.getVertices() OG
#  vertices = t.faces(0)

  

  #edges = t.getEdges() OG
#  edges = t.faces(1)
  #faces = t.getTriangles() OG
#  faces = t.faces(2)
  #tet = t.getTetrahedra() OG
#  tet = t.simplices()

  