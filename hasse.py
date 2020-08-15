#!/usr/bin/regina-python
from regina import *

import numpy as np

print np.array([1,2,3,4])

import sys
sys.setrecursionlimit(100)

def compute_hasse(isosig, dimension):
    if dimension ==3:
        t = Triangulation3.fromIsoSig(isosig)
    if dimension == 4:
        t = Triangulation4.fromIsoSig(isosig)
    cells = []
    dim = t.dimension
    for i in range(dim):
        cells.append(t.faces(i))
    cells.append(t.simplices())

    upward = [[] for _ in range(dim+1)]
    downward = [[] for _ in range(dim+1)]

    for j in range(dim+1):
        if j > 0:
            for i in range(len(cells[j])):
                downward[j].append([])
        if j < dim:
            for i in range(len(cells[j])):
                upward[j].append([])

    print "initial"
    for k in range(1,dim+1):
        print "k = ", k
        print downward[k]

    for k in range(1,dim+1):                    #dim chg
        for i in range(len(cells[k])):
            for j in range(k+1):
                idx = cells[k-1].index(cells[k][i].face(k-1,j))
                downward[k][i].append(idx)
                upward[k-1][idx].append(i)

    print "middle"
    for k in range(1,dim+1):
        print "k = ", k
        print downward[k]


    for k in range(dim+1):
        if k > 0:
            for i in range(len(downward[k])):
                for j in range(t.countFaces(k-1)):
                    if downward[k][i].count(j) > 1:
                        downward[k][i]=[y for y in downward[k][i] if y != j]
        if k < dim:
            for i in range(len(upward[k])):
                if k+1 == dim:
                    for j in range(t.size()):
                        if upward[k][i].count(j) > 1:
                            upward[k][i]=[y for y in upward[k][i] if y != j]
                else:
                    for j in range(t.countFaces(k+1)):
                        if upward[k][i].count(j) > 1:
                            upward[k][i]=[y for y in upward[k][i] if y != j]
    print "final"
    for k in range(1,dim+1):
        print "k = ", k
        print downward[k]

for isosig in sys.stdin:
    compute_hasse(isosig, 3)
    break