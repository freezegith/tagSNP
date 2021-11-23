# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:39:04 2020

@author: Awery
"""



import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import operator
from numpy import *
from numpy.linalg import *
import random

cc = 112                        #non-redundant SNPs
A=zeros((cc,cc),dtype=float)             
f=open('***')
lines = f.readlines()       
A_row = 0                 
for line in lines:         
    list0 = line.strip('\n').split(' ')    
    A[A_row:] = list0[0:cc]           
    A_row+=1                   
    #print(line)
f.close()


G=nx.Graph()                    #将矩阵转成图

temp=A

#print(temp)
print(len(temp))

#savetxt(r'***',temp,fmt='%d')

for i in range(len(temp)):
   for j in range(len(temp)):
      if(temp[i][j] == 1):
          G.add_edge(i, j)

nx.draw(G,with_labels=True,circular_layout=True)                  #绘制图    


'''
part = community.best_partition(A)
print (part)
'''


#按节点的度排序
def degrees():
    dc=G.degree()
    print("节点编号及其按度由大到小排序：")
    return(sorted(dc,key=lambda x:x[1],reverse=True))
print(degrees())
'''
#按聚类系数排序
def clustering(G,u):
    return sum([1 for x in G[u] for y in G[x] if y in G[u]]) / (G.degree[u] * (G.degree[u] - 1)) if G.degree[u] > 1 else 0
sortc=[(0,0)]*len(A)
for v in range(len(A)):
    sortc[v]=[v,clustering(G,v)]
print("节点编号及其按聚类系数由小到大排序：")
print(sorted(sortc,key=operator.itemgetter(1),reverse=False))
'''
'''
#按聚类系数排序new
def clustering():
    cc=nx.clustering(G)
    print("节点编号及其按聚类系数由小到大排序：")
    cc_list = sorted(cc.items(),key=operator.itemgetter(1),reverse=False)
    return (cc_list)
print (clustering())
'''
'''
#按特征向量中心性排序
def centrality():
    ec=nx.eigenvector_centrality(G)
    print("节点编号及其按特征向量中心性由大到小排序：")
    ec_list = sorted(ec.items(),key=operator.itemgetter(1),reverse=True)
    return (ec_list)
print (centrality())
'''
#按接近中心性排序
def closeness_centrality():
    clc=nx.closeness_centrality(G)
    print("节点编号及其按接近中心性由大到小排序：")
    clc_list = sorted(clc.items(),key=operator.itemgetter(1),reverse=True)
    return (clc_list)
print (closeness_centrality())
#按点介数排序
def betweenness():
  bc=nx.centrality.betweenness_centrality(G,normalized=False)
  print("节点编号及其按点介数由大到小排序：")
  bc_list = sorted(bc.items(),key=operator.itemgetter(1),reverse=True)
  return(bc_list)
print(betweenness())







