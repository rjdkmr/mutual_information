#!/usr/bin/env python
import re
import sys
import numpy as np
import mi
from ctypes import *

def read_file(FileName,n=2):
    infile = open(FileName,'r')
    x, data = [],[]
    for line in infile:
        if(re.match('#|@',line)==None):
            temp = line.split()
            data.append(np.array(temp))
    for j in range(0,n):
        x_temp =[]
        for i in range(len(data)):
            x_temp.append(np.float64((data[i][j])))
        x.append(x_temp)
    return x

x1 = read_file(sys.argv[1])
a1 = mi.new_doubleArray(len(x1[1]))

x2 = read_file(sys.argv[2]) 
a2 = mi.new_doubleArray(len(x2[1]))

t1 = np.random.normal(0,1,15000)
ta1 = mi.new_doubleArray(len(t1))
np.random.seed(2541822)
t2 = np.random.normal(0,1,15000)
ta2 = mi.new_doubleArray(len(t2))

#r = np.corrcoef(x1[0],x2[0])
#test = -(3/2)*np.log(1-(r[0][1]*r[0][1]))
#print r
#exit()


for i in range(len(t1)):
    mi.doubleArray_setitem(ta1,i,np.float64(t1[i]))
    mi.doubleArray_setitem(ta2,i,np.float64(t2[i]))

for i in range(len(x1[1])):
    mi.doubleArray_setitem(a1,i,x1[1][i])
    mi.doubleArray_setitem(a2,i,x2[1][i])

r1 = mi.mi(a1, a2, len(x1[1]),10)
r2 = mi.mi(a1, a1, len(x2[1]),10)
r3 = mi.mi(a2, a2, len(x2[1]),10)
print r1, r2, r3
#mi.mi(ta1, ta1, len(t1))
#mi.mi(ta1, ta2, len(t1))
#mi.mi(ta2, ta1, len(t1))
