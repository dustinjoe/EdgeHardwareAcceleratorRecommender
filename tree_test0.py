#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:35:15 2020

@author: ubuntu
"""

from anytree import Node, RenderTree, AsciiStyle, PreOrderIter



def createTree0(inMu0,inStd0):
    # static topology test
    cloud = Node("cloud")
    
    fog0 = Node("fog0", parent=cloud)
    fog1 = Node("fog1", parent=cloud)
    fog2 = Node("fog2", parent=cloud)
    
    
    edge0_0 = Node("edge0_0", parent=fog0, inMu=inMu0, inStd=inStd0)
    edge0_1 = Node("edge0_1", parent=fog0, inMu=inMu0, inStd=inStd0)
    edge0_2 = Node("edge0_2", parent=fog0, inMu=inMu0, inStd=inStd0)
    edge0_3 = Node("edge0_3", parent=fog0, inMu=inMu0, inStd=inStd0)
    
    edge1_0 = Node("edge1_0", parent=fog1, inMu=inMu0, inStd=inStd0)
    edge1_1 = Node("edge1_1", parent=fog1, inMu=inMu0, inStd=inStd0)
    edge1_2 = Node("edge1_2", parent=fog1, inMu=inMu0, inStd=inStd0)
    edge1_3 = Node("edge1_3", parent=fog1, inMu=inMu0, inStd=inStd0)  
    
    edge2_0 = Node("edge2_0", parent=fog2, inMu=inMu0, inStd=inStd0)
    edge2_1 = Node("edge2_1", parent=fog2, inMu=inMu0, inStd=inStd0)
    edge2_2 = Node("edge2_2", parent=fog2, inMu=inMu0, inStd=inStd0)
    edge2_3 = Node("edge2_3", parent=fog2, inMu=inMu0, inStd=inStd0)
    
    
    #print(cloud)
    print(RenderTree(cloud, style=AsciiStyle()).by_attr())
    # 
    return cloud

def createTree1(inMu0,inStd0):
    # static topology test
    cloud = Node("cloud")
    
    fog0 = Node("fog0", parent=cloud)
    fog1 = Node("fog1", parent=cloud)
    fog2 = Node("fog2", parent=cloud)
    fog3 = Node("fog3", parent=cloud)
    
    edge0_0 = Node("edge0_0", parent=fog0, inMu=inMu0, inStd=inStd0)
    edge0_1 = Node("edge0_1", parent=fog0, inMu=inMu0, inStd=inStd0)
    edge0_2 = Node("edge0_2", parent=fog0, inMu=inMu0, inStd=inStd0)
    
    edge1_0 = Node("edge1_0", parent=fog1, inMu=inMu0, inStd=inStd0)
    edge1_1 = Node("edge1_1", parent=fog1, inMu=inMu0, inStd=inStd0)
    edge1_2 = Node("edge1_2", parent=fog1, inMu=inMu0, inStd=inStd0)
    
    edge2_0 = Node("edge2_0", parent=fog2, inMu=inMu0, inStd=inStd0)
    edge2_1 = Node("edge2_1", parent=fog2, inMu=inMu0, inStd=inStd0)
    edge2_2 = Node("edge2_2", parent=fog2, inMu=inMu0, inStd=inStd0)
    
    edge3_0 = Node("edge3_0", parent=fog3, inMu=inMu0, inStd=inStd0)
    edge3_1 = Node("edge3_1", parent=fog3, inMu=inMu0, inStd=inStd0)
    edge3_2 = Node("edge3_2", parent=fog3, inMu=inMu0, inStd=inStd0)
    
    #print(cloud)
    print(RenderTree(cloud, style=AsciiStyle()).by_attr())
    # 
    return cloud
'''
cloud = createTree0(inMu,inStd)
print(cloud)
print(RenderTree(cloud, style=AsciiStyle()).by_attr())
'''  

#%%
from sympy import Symbol
from sympy import solveset, S
from sympy.abc import x

from sympy.parsing.sympy_parser import parse_expr

from scipy.special import comb 

def twoINthreeRelia(x):
    return x**3 + 3*(x**2)*(1-x) #parse_expr("x**3 + 3*(x**2)*(1-x)", evaluate=False)

def twoINfourRelia(x):
    return x**4 + 4*(x**3)*(1-x) + 6*(x**2)*(1-x)**2



def k_IN_n_ReliaExpr(x,k,n):
    # generate reliability function for k out of n system 
    # x is the reliability for one node
    # Reliable: in n components, no less than k is working
    x = str(x) 
    expr = x+"**"+str(n)
    if k==0:
        return expr
    # from 1 fail(n-1 working) to (n-k) fail (k working)
    for i in range(n-k): 
        num_fail = i+1
        num_comb = int(comb(n,num_fail))
        num_work = int(n-num_fail)
        expr = expr + '+'+ str(num_comb) + '*('+ x + '**'+str(num_work)+')*(1-'+x+')**'+str(num_fail)
    return expr    

def k_IN_n_ReliaCompute(x,k,n):
    # compute reliability prob value for k out of n system 
    # x is the reliability for one node
    # Reliable: in n components, no less than k is working
    p_sys = x**n
    if k==0:
        return p_sys
    # from 1 fail(n-1 working) to (n-k) fail (k working)
    for i in range(n-k): 
        num_fail = i+1
        num_comb = int(comb(n,num_fail))
        num_work = int(n-num_fail)
        p_sys = p_sys + num_comb * (x**num_work) * ((1-x)**num_fail) 
        #expr = expr + '+'+ str(num_comb) + '*('+ x + '**'+str(num_work)+')*(1-'+x+')**'+str(num_fail)
    return p_sys    
    
#print(k_IN_n_Relia(x,2,4))    

def solveNodeRelia(rsys_func,Rsys):
    result = solveset(rsys_func - Rsys, x, domain=S.Reals)
    result = list(result)
    Rnode = 0
    for i in range(len(result)):
        if result[i]>=0 and result[i]<=1:
            Rnode = float(result[i])
    return Rnode

def k_IN_n_solver(Rsys,x,k,n):
    exprStr = k_IN_n_ReliaExpr(str(x),k,n)
    expr = parse_expr(exprStr, evaluate=False)
    pnode = solveNodeRelia(rsys_func=expr,Rsys=Rsys)
    #print(pnode)
    return pnode

    
'''
R_sys = 0.936 
k=2
n=3
p = k_IN_n_solver(R_sys,x,k,n)

r_func = twoINthreeRelia(x)
r_func2 = twoINfourRelia(x)
'''

#%% 
from math import sqrt,ceil
from scipy.stats import norm
#norm.ppf(0.95, loc=0, scale=1)
# compute number of a type of device needed given an input level and confidence
# theoretical computation power estimation regardless of placement location
def getDevNum(muDev,stdDev, muIn,cvIn, conf):
    Pr = conf
    stdIn = muIn*cvIn
    # ppf(q, loc=0, scale=1) Percent point function (inverse of cdf — percentiles).
    idxR = norm.ppf(Pr, loc=0, scale=1)
    std_G = sqrt( stdDev**2 + stdIn**2 )
    mu_G = idxR*std_G
    # mu_G=numDev*muDev-muIn
    numDev = (mu_G + muIn)/muDev
    numDev = int( ceil(numDev) )
    return numDev
    


def numDevTotal(muDevFreq,stdDevFreq,level,tree,conf0):
    # bandwidth assumption using Wifi
    # https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8016573
    #B_e2f = 135*1024/8  #135Mbps Wifi  IEEE802.11n convert to kb/s
    #B_f2c = 135*1024/8  #135Mbps Wifi  IEEE802.11n convert to kb/s	
    #t_e2f = S_dat/B_e2f
    #t_f2c = S_dat/B_f2c	
	
    #conf = 0.99
    numDev = 0
    
    muInTotal = 0
    stdInTotal= 0
    varIn = 0
    
    numFog = len(tree.children)
    
    Pr = conf0
    Pf = 1-conf0
    
    if level == 0:
        stdInTotal = 0
        for i in range(numFog):
            numEdge = len(tree.children[i].children)
            varIn = 0
            for j in range(numEdge):
                muInTotal = muInTotal+ tree.children[i].children[j].inMu
                stdij = tree.children[i].children[j].inStd
                varIn = varIn + stdij*stdij
                stdInTotal = stdInTotal + stdij
            #stdInTotal = stdInTotal + sqrt(varIn)
        print('muInTotal',muInTotal)
        
        #stdInTotal = sqrt(varIn)
        stdInTotal = stdInTotal*2 
        print('stdInTotal',stdInTotal)
        cvInTotal = stdInTotal/muInTotal
        conf = conf0
        numDev = getDevNum(muDevFreq,stdDevFreq, muInTotal,cvInTotal, conf)
    elif level == 1:        
        for i in range(numFog):
            numEdge = len(tree.children[i].children)
            muIn_fog = 0
            varIn_fog = 0
            for j in range(numEdge):
                muIn_fog = muIn_fog + tree.children[i].children[j].inMu
                stdij = tree.children[i].children[j].inStd
                varIn_fog = varIn + stdij*stdij
            #print(muIn_fog)            
        
            stdIn_fog = sqrt(varIn_fog)
            #print(stdIn_fog)
            cvIn_fog = stdIn_fog/muIn_fog
            numDev_fog = 0
            numDev_fog = getDevNum(muDevFreq,stdDevFreq, muIn_fog,cvIn_fog, conf)
        numDev = numDev+numDev_fog
    elif level == 2:
        Pr_fog = k_IN_n_solver(Pr,x,numFog,numFog)
        #Pr_fog = Pr**(1/numFog)
        Pf_fog = 1-Pr_fog
        print('Pf_fog',Pf_fog)
        for i in range(numFog):
            numEdge = len(tree.children[i].children)
            #Pf_fog = Pf**(1/numFog)
            #print(Pf_fog)           
            thresEdge = 1
            if numEdge>2:
                thresEdge = ceil(numEdge/2)
            Pr_edge = k_IN_n_solver(Pr_fog,x,thresEdge,numEdge)
            print('Pf_edge',1-Pr_edge)
            conf = Pr_edge
            for j in range(numEdge):
                muIn_ij = tree.children[i].children[j].inMu
                std_ij = tree.children[i].children[j].inStd
                cv_ij = std_ij/muIn_ij
                #Pf_edge = Pf_fog**(1/numEdge)                
                #conf = 1-Pf_edge                          
                
                
                #print(Pr_edge)
                numDev_ij = getDevNum(muDevFreq,stdDevFreq, muIn_ij,cv_ij, conf)
                #print(numDev_ij)
                numDev = numDev+numDev_ij
                

        
        
    return numDev
'''
inMu =10
inStd = 15

cloud = createTree0(inMu,inStd)
print(cloud)
print(RenderTree(cloud, style=AsciiStyle()).by_attr())

conf = 0.99

muRes = [ 0.479,7.621,34.5046,4.587,26.0323,3.733]
stdRes = [0.013,0.962,1.19,0.0607,3.37367,0.081]
muYolo = [0.3486,10.443,43.51,4.2035,16.9898,5.9358]
stdYolo = [0.0082,0.83738,1.848,0.051,0.59264,3.2745]

# 
muDev=muRes[4]
stdDev=stdRes[4]
n_cloud_1060 = numDevTotal(muDevFreq=muDev,stdDevFreq=stdDev,level=0,tree=cloud,conf0=conf)
print( 'n_cloud_1060',n_cloud_1060 )
      
muDev=muRes[1]
stdDev=stdRes[1]
n_edge_jetnano = numDevTotal(muDevFreq=muDev,stdDevFreq=stdDev,level=2,tree=cloud,conf0=conf) 
print( 'n_edge_jetnano',n_edge_jetnano  )

muDev=muRes[2]
stdDev=stdRes[2]
n_edge_u96 = numDevTotal(muDevFreq=muDev,stdDevFreq=stdDev,level=2,tree=cloud,conf0=conf) 
print('n_edge_u96',n_edge_u96)

muDev=muRes[2]
stdDev=stdRes[2]
n_cloud_u96 = numDevTotal(muDevFreq=muDev,stdDevFreq=stdDev,level=0,tree=cloud,conf0=conf) 
print('n_cloud_u96',n_cloud_u96)  

'''
