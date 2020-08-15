#!/usr/local/bin/regina-python

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

import sys
sys.setrecursionlimit(100)

##########################################
##########################################
############ HELPER FUNCTIONS ############
##########################################
##########################################

### computes the boundary of all triangles
### and returns a list of 1-chains (1-cycles)
def SCBdry(t):
  edgeEmbeddings=[]
  for j in range(t.getNumberOfEdges()):
    edgeEmbeddings.append([])
    e = t.getEdge(j)
    for i in range(e.getDegree()):
      emb = e.getEmbedding(i)
      tetIdx = t.tetrahedronIndex(emb.getTetrahedron())
      vtcs = emb.getVertices()
      edgeEmbeddings[j].append([tetIdx,vtcs[0],vtcs[1]])
  bdrys = []
  for j in range(t.getNumberOfTriangles()):
    f = t.getTriangle(j)
    emb=f.getEmbedding(0)
    tetIdx = t.tetrahedronIndex(emb.getTetrahedron())
    vtcs = emb.getVertices()
    b = [[],[],[]]
    if [tetIdx,vtcs[1],vtcs[2]] in edgeEmbeddings[t.edgeIndex(f.getEdge(0))]:
      b[0] = [t.edgeIndex(f.getEdge(0)),1]
    else:
      b[0] = [t.edgeIndex(f.getEdge(0)),-1]

    if [tetIdx,vtcs[0],vtcs[2]] in edgeEmbeddings[t.edgeIndex(f.getEdge(1))]:
      b[1] = [t.edgeIndex(f.getEdge(1)),-1]
    else:
      b[1] = [t.edgeIndex(f.getEdge(1)),1]

    if [tetIdx,vtcs[0],vtcs[1]] in edgeEmbeddings[t.edgeIndex(f.getEdge(2))]:
      b[2] = [t.edgeIndex(f.getEdge(2)),1]
    else:
      b[2] = [t.edgeIndex(f.getEdge(2)),-1]
    bdrys.append(b)
  return bdrys

### collapses a knot complement and computes a Morse function
### of a knot complement from an oriented Hasse diagram
def collKnotCompl(upward,downward,t):
  f=[1,0,0,0]
  Morse=[]
  critical=[[0],[],[],[]]
  available=[[0],range(t.getNumberOfEdges()),range(t.getNumberOfTriangles()),range(t.getNumberOfTetrahedra())]

  for iii in [2,1]:
    free=[]
    for i in range(len(upward[iii])):
      if len(upward[iii][i]) == 1:
        free.append(i)
    while available[iii+1] <> []:
      if free==[]:
        f[iii+1]+=1
        r=random.choice(available[iii+1])
        # keep track of Morse function
        Morse.append([iii+1,r])
        critical[iii+1].append(r)
        available[iii+1].pop(available[iii+1].index(r))
        # update upward Hasse diagram
        R=downward[iii+1][r]
        for i in R:
          upward[iii][i].pop(upward[iii][i].index(r))
          if len(upward[iii][i]) == 1:
            free.append(i)
      else:
        r=random.choice(free)
        free.pop(free.index(r))
        pairedFace=upward[iii][r][0]
        # keep track of Morse function
        Morse.append([iii,r])
        Morse.append([iii+1,pairedFace])
        available[iii].pop(available[iii].index(r))
        available[iii+1].pop(available[iii+1].index(pairedFace))
        R=downward[iii+1][pairedFace]

        # update upward Hasse diagram
        for i in R:
          upward[iii][i].pop(upward[iii][i].index(pairedFace))
          if len(upward[iii][i]) == 1:
            free.append(i)
          elif upward[iii][i] == []:
            if i in free:
              free.pop(free.index(i))
        if iii > 1:
          R=downward[iii][r]
          for i in R:
            upward[iii-1][i].pop(upward[iii-1][i].index(r))

  for i in available[1]:
    Morse.append([1,i])
    critical[1].append(i)
    f[1]+=1
  Morse.append([0,0])
  return [f,critical,Morse]

def SCAddCrits(s1,s2):
  for i in range(len(s1)):
    s1[i][1]=s1[i][1]+s2[i][1];
  return s1

def SCAddCrit(s,chain):
  for i in range(len(s)):
    if s[i][0]==chain[0]:
      s[i][1]=s[i][1]+chain[1]
  return s


def SCGradient(Morse,face,bdrys):
  pos=Morse.index([1,face])
  matching=[]
  for i in range(pos,len(Morse)):
    if Morse[i][0] <> 2: continue
    if not face in [bdrys[Morse[i][1]][0][0],bdrys[Morse[i][1]][1][0],bdrys[Morse[i][1]][2][0]]: continue
    matching=Morse[i][1]
    if face == bdrys[Morse[i][1]][0][0]:
      sign=bdrys[Morse[i][1]][0][1]
    elif face == bdrys[Morse[i][1]][1][0]:
      sign=bdrys[Morse[i][1]][1][1]
    elif face == bdrys[Morse[i][1]][2][0]:
      sign=bdrys[Morse[i][1]][2][1]
    break
  if matching == []:
    return 0
  else:
    return [matching,(-1)*sign]


def SCFindGradientPaths(Morse,chain,crits,mult,lookup,bdrys):
  dict1=lookup
  # check if oriented edge 'chain' was already computed
  check=dict1.get(tuple(chain))
  if check<>None:
    return [check,dict1]

  # empty chain
  s=[[x,0] for x in crits]

  gradTrig=SCGradient(Morse,chain[0],bdrys)
  if gradTrig==0:
    # if no gradient triangle found, update dictionary with empty chain 's'
    # positive case
    dict1.update({tuple(chain): s})
    # negative case
    dict1.update({tuple([chain[0],-1*chain[1]]): [[x[0],-1*x[1]] for x in s]})
    return [s,dict1]
  
  #print 'edge', chain[0], '->', 'triangle', gradTrig
  # case: gradient triangle found
  # adjust multiplicity of gradient triangle
  gradTrig[1]=gradTrig[1]*chain[1]
  boundaryGradTrig=bdrys[gradTrig[0]]
  # orient boundary
  if gradTrig[1] == -1:
    boundaryGradTrig[0][1]*=-1
    boundaryGradTrig[1][1]*=-1
    boundaryGradTrig[2][1]*=-1
  outgoingEdges=[x for x in boundaryGradTrig if x[0]<>chain[0] and Morse.index([1,chain[0]]) < Morse.index([1,x[0]])]
  adjacent=[x[0] for x in outgoingEdges]

  # case: outgoing edges are critical edges
  intersection=[x for x in adjacent if x in crits];
  for i in intersection:
    pos=adjacent.index(i)
    if mult == outgoingEdges[pos][1]:
      s=SCAddCrit(s,[i,1])
    else:
      s=SCAddCrit(s,[i,-1])
    outgoingEdges.pop(pos)
    adjacent.pop(pos)




  #print len(outgoingEdges), 'recursion(s) ahead'
  #print 'triangle', gradTrig, '->', 'edge(s)', outgoingEdges
  # loop over remaining outgoing edges
  for i in outgoingEdges:
    #print 'SCFindGradientPaths:', x, ',', i, ',',  crits, ',', lookup
    #print 'SCFindGradientPaths:', Morse.index([1,i[0]]), Morse.index([1,crits[0]]), Morse.index([1,crits[1]]), dict1
    new=SCFindGradientPaths(Morse,i,crits,mult,dict1,bdrys)
    dict1=new[1]
    s=SCAddCrits(s,new[0])
  dict1.update({tuple(chain): s})
  dict1.update({tuple([chain[0],-1*chain[1]]): [[x[0],-1*x[1]] for x in s]});
  return [s,dict1]


# compute bdry operator between triangles and edges
# in a knot complement
def SCBdryOp(Morse,critsUp,critsDown,t,bdrys):
  lookup = {}
  tau=[]
  for i in critsUp:
    s=[]
    f = t.getTriangle(i)
    bd = bdrys[i]
    for x in critsDown:
      s.append([x,0])
    for ii in bd:
      if ii[0] in critsDown:
        s=SCAddCrit(s,ii)
      else:
        new=SCFindGradientPaths(Morse,ii,critsDown,1,lookup,bdrys)
        s1=new[0]
        s=SCAddCrits(s,s1)
    tau.append(s)
  tau2=[]
  for i in range(len(critsDown)):
    row=[]
    for ii in range(len(critsUp)):
      row.append(tau[ii][i][1])
    tau2.append(row)
  return tau2


##########################################
##########################################
#########END HELPER FUNCTIONS ############
##########################################
##########################################

ctr=0
for line in sys.stdin:
  ctr += 1
  #if ctr%1000 == 0:
  #  print ctr
  #sig = re.search('[a-zA-Z0-9]*' ,line)
  t = NTriangulation.fromIsoSig(line)
  #print t.isoSig()
  vertices = t.getVertices()
  edges = t.getEdges()
  faces = t.getTriangles()
  tet = t.getTetrahedra()

#  print t.getNumberOfVertices(),
#  print t.getNumberOfEdges(),
#  print t.getNumberOfTriangles(),
#  print t.getNumberOfTetrahedra()
#  for i in vertices:
#    print i
#  for i in edges:
#    print i
#  for i in faces:
#    print i
#  for i in tet:
#    print i
  upward=[[],[],[]]
  downward=[[],[],[],[]]
  downward[0]=[]
  downward[1]=[]
  upward[0]=[]
#  for i in range(len(vertices)):
#    downward[0].append([])
#    upward[0].append([])
  for i in range(len(edges)):
    upward[1].append([])
#    downward[1].append([])
  for i in range(len(faces)):
    upward[2].append([])
    downward[2].append([])
  for i in range(len(tet)):
    downward[3].append([])
#  for i in range(len(edges)):
#    for j in range(2):
#      idx = vertices.index(edges[i].getVertex(j))
#      downward[1][i].append(idx)
#      upward[0][idx].append(i)
  for i in range(len(faces)):
    for j in range(3):
      idx = edges.index(faces[i].getEdge(j))
      downward[2][i].append(idx)
      upward[1][idx].append(i)
  for i in range(len(tet)):
    for j in range(4):
      idx = faces.index(tet[i].getTriangle(j))
      downward[3][i].append(idx)
      upward[2][idx].append(i)
  # remove double edges
#  for i in range(len(downward[1])):
#    for j in range(t.getNumberOfVertices()):
#      if downward[1][i].count(j) > 1:
#        downward[1][i]=[y for y in downward[1][i] if y != j]
  for i in range(len(downward[2])):
    for j in range(t.getNumberOfEdges()):
      if downward[2][i].count(j) > 1:
        downward[2][i]=[y for y in downward[2][i] if y != j]
  for i in range(len(downward[3])):
    for j in range(t.getNumberOfTriangles()):
      if downward[3][i].count(j) > 1:
        downward[3][i]=[y for y in downward[3][i] if y != j]
#  for i in range(len(upward[0])):
#    for j in range(t.getNumberOfEdges()):
#      if upward[0][i].count(j) > 1:
#        upward[0][i]=[y for y in upward[0][i] if y != j]        
  for i in range(len(upward[1])):
    for j in range(t.getNumberOfTriangles()):
      if upward[1][i].count(j) > 1:
       upward[1][i]=[y for y in upward[1][i] if y != j] 
  for i in range(len(upward[2])):
    for j in range(t.getNumberOfTetrahedra()):
      if upward[2][i].count(j) > 1:
        upward[2][i]=[y for y in upward[2][i] if y != j] 

  tmp = collKnotCompl(upward,downward,t)
  f = tmp[0]
  critical = tmp[1]
  Morse = tmp[2]
  bdrys = SCBdry(t)

  #test=[0,0,0,0]
  #for i in Morse:
  #  test[i[0]]+=1
  #print [t.getNumberOfVertices(),t.getNumberOfEdges(),t.getNumberOfTriangles(),t.getNumberOfTetrahedra()]==test
  #print Morse
  #print critical
  #print t.toStringLong()

################ NEW STUFF ###############

  critsUp=critical[2]
  critsDown=critical[1]

  # get boundary operator
#  print t.isoSig()
  tmp=SCBdryOp(Morse,critsUp,critsDown,t,bdrys)


  # here you can say how many critical cells you want to have at least before you output something. 
  # At the moment everything is printed.
  if len(critsUp) < 0:
    continue
  else:
    # checking for critical dunce hats and similar trivial examples
    for j in critsUp:
      if bdrys[j][0][1]==bdrys[j][1][1] and bdrys[j][0][1]==bdrys[j][2][1]: continue
    print '# isomorphism signature:',
    print t.isoSig(), '\n\n'
    print '# downward Hasse diagram (with multiple edges removed)'
    print '['
    print '# tetrahedra to triangles'
    print downward[3]
    print '# triangles to edges'
    print downward[2]
    print ']\n\n'
    print 'Morse function ([i,j] means j-th face of dimension i):\n'
    print Morse, '\n\n'
    print 'critical triangle(s):\t', critsUp, '\tcritical edges:\t', critsDown, '\n\n'
    print 'oriented boundaries of triangles (k-th entry [[i_0,s_0],[i_1,s_1],[i_2,s_2]]: i_j = j-th edge of triangle k, s_j = orientation of i_j)\n'
    for i in bdrys:
      print i
    print '\n'
    print 'induced boundary operator between critical triangles (columns) and critical edges (rows)\n'
    for i in tmp:
      print i

#    print t.toStringLong(), "\n\n\n\n\n"
    print "\n\n\n\n\n"






