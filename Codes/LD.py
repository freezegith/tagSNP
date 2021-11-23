# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:35:53 2020

@author: Awery
"""



import networkx as nx
import matplotlib.pyplot as plt
import operator
from numpy import *
from numpy.linalg import *
import random


#一
#基于原始数据
set_printoptions(suppress=True)
A=zeros((176,2000),dtype=float)              
f=open('***')                       #输入文件位置
lines = f.readlines()       
A_row = 0                    
for line in lines:           
    list0 = line.strip('\n').split(' ')
    A[A_row:] = list0[0:2000]           
    A_row+=1                  
    #print(line)
f.close()
#print(lines)


n=len (A)
#删除满度行
for i in range (n-1,-1,-1):
    if sum(A[i]) == 2000:
        A=delete(A,i,axis=0)
        n=n-1
#删除零度行
for i in range (n-1,-1,-1):
    if sum(A[i]) == 0:
        A=delete(A,i,axis=0)
        n=n-1
        
#删除重复行
for i in range (n-1,-1,-1):
    for j in range (i-1,-1,-1):
        if (A[i] == A[j]).all():
            A=delete(A,i,axis=0)
            n=n-1

#删除互补行
for i in range (n-1,-1,-1):
    for j in range (i-1,-1,-1):
        if (A[i] != A[j]).all():
            A=delete(A,i,axis=0)
            n=n-1

savetxt(r'***', A,fmt='%d')

print (n)


#计算连锁不平衡
E=zeros(n,dtype=float)                    #E[i]表示第i行均值=P(a) 1-E[i]=P(A)
for i in range (n):
    E[i]=mean(A[i])

R=zeros((n,n),dtype=float)       #R为连锁不平衡系数
S=0
for j in range (n):
    for i in range (j+1,n):
        for k in range (2000):
            if A[i][k]==A[j][k]==0:
                S+=1
        R[i][j]=((S/2000-(1-E[i])*(1-E[j])))**2/(E[i]*(1-E[i])*E[j]*(1-E[j]))
        R[j][i]=R[i][j]
        S=0
        
for i in range(len(R)):            
   for j in range(len(R)):
       if R[i][j] > 0.8:            #阈值
           R[i][j] = 1
       else:
           R[i][j] = 0
           
#R = R.astype(int32)      #数据类型转换

savetxt(r'***', R,fmt='%d')

print (n)
