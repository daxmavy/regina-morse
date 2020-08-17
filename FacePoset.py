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

# for duplicating Hasse diagram for randomised algorithm
import copy

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

class InputError(Exception):
    pass



class FacePoset:
    """
    This class represents the face poset diagram. Nodes in the diagram are given by the PosetNode class below.
    FacePoset can be instantiated by passing a TriangulationN (eg Triangulation3) in the Regina format, along with
    the dimension of the triangulation

    This is designed so that the user should only interact with the FacePoset object, with the PosetNode objects 
    used only for the backend of the FacePoset object.

    Example usage:
    
    tri = Triangulation3.fromIsoSig('fLAMcbcbdeehxjqhr')
    fp = FacePoset(triangulation = tri, dim = 3)
    fp.strip_multi_edges()
    fp.output_poset()

    User methods:

    get_node
    get_cell
    add_node
    add_arc: this requires the two end nodes to exist
    remove_node
    remove_arc
    """
    class PosetNode:
        """
        PosetNodes have three main attributes:
        cell: a reference to the Regina object for the cell in the triangulation
        name: an integer index uniquely identifying this cell amongst cells of the same dimension
        dim: an integer corresponding to the dimension of the cell

        PosetNodes are uniquely identified inside the FacePoset by the tuple (self.dim, self.name). 
        This is implemented by the PosetNode.key() method, which returns the above tuple.

        PosetNodes also carry 'pointers' to other cells:
        parents:    a list of references to PosetNodes of dimension self.dim+1, which have an arc between
                    each parent and the current PosetNode
        children:   a list of references to PosetNodes of dimension self.dim-1, which have an arc between
                    each child and the current PosetNode

        Note that during creation of the FacePoset, a parent may have multiple faces which are identified.
        As such, the PosetNode.parents list may well have repeat entries. This is intended behaviour.

        To strip these edges, apply FacePoset.strip_multi_edges() to the FacePoset object. 
        Parent / child relationships will be stored in a separate 'irregular_parents' and 'irregular_children'
        lists as attributes of each PosetNode, similarly to the ordinary 'parents' and 'children' attributes.
        """

        def __init__(self, dim, name, cell,
                     parents = None, children = None):
            self.cell = cell
            self.name = name
            self.dim = dim
            self.irregular_parents = []
            self.irregular_children = []
            self.morse_matched = None 

            if not parents:
                self.parents = []
            else:
                self.parents = parents
            if not children:
                self.children = []
            else:
                self.children = children
        
        def key(self):
            return (self.dim, self.name)

        def __key(self):
            return (self.dim, self.name)
        
        def __hash__(self):
            return hash(self.__key())

        def __eq__(self, other):
            return self.__key() == other.__key()

        def add_child(self, node):
            self.children.append(node)
        
        def add_parent(self, node):
            self.parents.append(node)
        
        def remove_child(self, node, error = False):
            if node in self.irregular_children:
                self.irregular_children.remove(node)
            if node in self.children:
                self.children.remove(node)
            elif error:
                raise FindCellFailure('No such cell '+str((node.dim, node.name))+' found')

        def remove_parent(self, node, error = False):
            if node in self.irregular_parents:
                self.irregular_parents.remove(node)
            if node in self.parents:
                self.parents.remove(node)
            elif error:
                raise FindCellFailure('No such cell '+str((node.dim, node.name))+' found')
    
    def get_node(self, dim, name):
        try:
            return self.layers[dim][name]
        except:
            raise FindCellFailure('Could not find cell '+str(name)+' in layer '+str(dim)) 

    def get_cell(self, dim, name):
        return self.layers[dim][name].cell

    def add_node(self, dim, name, cell):
        if not dim in self.layers.keys():
            self.layers[dim] = {}
        self.layers[dim][name] = self.PosetNode(name = name, cell = cell, dim = dim)

    def remove_node(self, dim, node_label, suppress_error = False):
        if not dim in self.layers.keys():
            if suppress_error:
                return
            string = 'Could not find layer of dimension '+str(dim)+' when deleting node '+str(node_label)
            raise FindLayerFailure(string)
        if not node_label in self.layers[dim].keys():
            if suppress_error:
                return
            string = 'Could not find node '+str((dim, node_label))+' when trying to delete it'
            raise FindCellFailure(string)
        else:
            #print 'before remove_node'
            #fp.output_poset()
            
            node = self.layers[dim][node_label]
            #print type(node)
            #print 'stripping the neighbours of ',node.key(), 'before:'
            #print [(el.key(), [l.key() for l in el.parents]) for el in node.children]
            #print [(el.key(), [l.key() for l in el.children]) for el in node.parents]
            for child in node.children:
                child.remove_parent(node)
            for irr_child in node.irregular_children:
                irr_child.remove_parent(node)
            for parent in node.parents:
                parent.remove_child(node)
            for irr_parent in node.irregular_parents:
                irr_parent.remove_child(node)
            #print 'after:'
            #print [(el.key(), [l.key() for l in el.parents]) for el in node.children]
            #print [(el.key(), [l.key() for l in el.children]) for el in node.parents]
            #_ = self.layers[dim].pop(node_label)
            del self.layers[dim][node_label]
            del node
            #print 'after remove_node'
            #fp.output_poset()
            #print type(node)
            

    def add_arc(self, n1_tup, n2_tup):
        dim1, name1 = n1_tup
        dim2, name2 = n2_tup
        
        for dim in dim1, dim2:
            if not dim in self.layers.keys():
                raise FindLayerFailure('Could not find layer '+str(dim)+' for creating arc between '+str((n1_tup, n2_tup)))

        for dim, name in n1_tup, n2_tup:
            if not name in self.layers[dim].keys():
                raise FindCellFailure('Could not find node '+str((dim, name))+' for creating arc between '+str((n1_tup, n2_tup)))
        
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

    def remove_arc(self, tup1, tup2, error = False):
        dim1, name1 = tup1
        dim2, name2 = tup2
        
        for dim in dim1, dim2:
            if not dim in self.layers.keys():
                raise FindLayerFailure('Could not find layer '+str(dim)+' for creating arc between '+str(n1_tup, n2_tup))

        for dim, name in tup1, tup2:
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
        #print 'Test'
        #print n1.dim, n1.name, n1.children
        #print n2.dim, n2.name, n2.parents

        if error:
            n1.remove_child(n2, error = True)
            n2.remove_parent(n1, error = True)
        n1.remove_child(n2)
        n2.remove_parent(n1)

    def output_poset(self):
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

    def __init__(self, triangulation = None, dim = None):
        
        if triangulation and not dim:
            raise InputError('No dim arg given, must be given accompanying a triangulation input')
        self.layers ={}
        if not triangulation and not dim:
            return
        self.dim = dim
        for dimension in range(dim):
            self.layers[dimension] = { name: self.PosetNode(dimension, name, cell) for name, cell in enumerate(triangulation.faces(dimension))}
        self.layers[dim] = {name: self.PosetNode(dim, name, cell) for name, cell in enumerate(triangulation.simplices())}

        for dimension in range(dim, 0, -1):
            for name, node in self.layers[dimension].items():
                cell = node.cell 
                for j in range(dimension+1):
                    face = cell.face(dimension-1, j)
                    for face_name, face_node in self.layers[dimension-1].items():
                        if face_node.cell == face:
                            self.add_arc((dimension, name), (dimension-1, face_name))
    
    def strip_multi_edges(self):
        for dimension in range(self.dim, 0, -1):
            for name, node in self.layers[dimension].items():
                node.children, node.irregular_children = self.separate_duplicates(node.children)
        for dimension in range(self.dim):
            for name, node in self.layers[dimension].items():
                node.parents, node.irregular_parents = self.separate_duplicates(node.parents)


    def separate_duplicates(self, array):
        seen = {}
        for x in array:
            if x not in seen:
                seen[x] = 1
            else:
                seen[x] += 1

        not_duplicates = [el for (el, count) in seen.items() if count == 1]
        duplicates = [el for (el, count) in seen.items() if count > 1]
        return not_duplicates, duplicates

    def match(self, tup1, tup2):
        dim1, name1 = tup1
        dim2, name2 = tup2

        for dim in dim1, dim2:
            if not dim in self.layers.keys():
                raise FindLayerFailure('Could not find layer '+str(dim)+' for creating arc between '+str((n1_tup, n2_tup)))

        for dim, name in n1_tup, n2_tup:
            if not name in self.layers[dim].keys():
                raise FindCellFailure('Could not find node '+str((dim, name))+' for creating arc between '+str((n1_tup, n2_tup)))
        
        if dim1 == dim2:
            raise LogicalMistake('cannot place arc between two nodes of the same dimension')

        if abs(dim1 - dim2) != 1:
            raise LogicalMistake('cannot place arc between two nodes which are not in adjacent layers')

        n1 = self.layers[dim1][name1]
        n2 = self.layers[dim2][name2]

        n1.morse_matched = n2
        n2.morse_matched = n1

    def randomised_morse_matching(self):
        # for now, implement this in a destructive way
        morse_pairs = []
        critical = []
        critical_candidate = None
        while True:
            found_unmatched_uncritical = False
            match_made = False
            for dimension in range(self.dim, -1, -1):
                for name, node in self.layers[dimension].items():
                    if not found_unmatched_uncritical:
                        critical_candidate = node
                        #print 'this is the candidate', node.dim, node.name 
                        found_unmatched_uncritical = True
                    if len(node.parents) == 1 and len(node.irregular_parents) == 0:
                        #print critical
                        #print [ el.key() for el in self.layers[dimension].values()]
                        morse_pairs.append((node.key(), node.parents[0].key()))
                        #print [(node.key(), el.key()) for el in node.parents]
                        self.remove_node(node.parents[0].dim, node.parents[0].name)
                        #print [(node.key(), el.key()) for el in node.parents]
                        self.remove_node(node.dim, node.name)
                        #except:
                        #    print critical 
                        #    print morse_pairs
                        #    return
                        match_made = True
            # need to save a candidate for making a cell critical so I don't need to repeat the search
            if found_unmatched_uncritical and not match_made:
                #make cell critical
                #print 'this is before deletion', critical_candidate.dim, critical_candidate.name
                #print 'We are about to make ', critical_candidate.key(), ' critical. Its parents are:'
                #print self.keys(critical_candidate.parents)
                #print 'Irregular parents:'
                #print self.keys(critical_candidate.irregular_parents)
                #print 'Children:'
                #print self.keys(critical_candidate.children)
                #print 'Irregular children:'
                #print self.keys(critical_candidate.irregular_children)
                critical.append((critical_candidate.dim, critical_candidate.name))
                self.remove_node(critical_candidate.dim, critical_candidate.name)
                #print 'this is the face poset after deletion'
                #fp.output_poset()
            if not found_unmatched_uncritical:
                break
        return morse_pairs, critical

    def keys(self, array):
        return [el.key() for el in array]

#fp = FacePoset()
tri = Triangulation3.fromIsoSig('fLAMcbcbdeehxjqhr')
fp = FacePoset(triangulation = tri, dim = 3)
fp.strip_multi_edges()
fp.output_poset()


morse, critical = fp.randomised_morse_matching()
print 'Morse pairs:'
print morse
print 'Critical cells:'
print critical



#fp.add_node(0, 0, 'A')
#fp.add_node(0, 1, 'B')
#fp.add_node(1,0, 'C')
#fp.add_node(2,0, 'D')

#fp.add_arc((0,0), (1,0))
#fp.add_arc((1,0), (2,0))
#fp.output()
#fp.remove_arc((0,0), (1,0), error = True)
#fp.output()
#fp.remove_arc((1,0), (2,0), error = True)
#fp.output()

#n1 = fp.get_node(1,0)
#print n1.cell





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

  