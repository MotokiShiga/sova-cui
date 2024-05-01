# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 13:54:37 2022

@author: H. Morita
"""

import math
import numpy as np

length = [6.18704,  2.91581,  7.28313] # a, b, c
degree_angle = [81.6066,  99.4557,  76.4221] # alpha, beta, gamma
angle = [i * math.pi/180.0 for i in degree_angle] # degree --> rad

a_axis = [0.0] * 3
b_axis = [0.0] * 3
c_axis = [0.0] * 3

a_axis[0] = length[0]
a_axis[1] = 0.00000
a_axis[2] = 0.00000

b_axis[0] = length[1]*math.cos(angle[2])
b_axis[1] = length[1]*math.sin(angle[2])
b_axis[2] = 0.00000

print('angle[2]',angle[2],math.cos(angle[2]),math.acos(b_axis[0]/length[1])*180/math.pi)
print(b_axis[0])
c_axis[0] = length[2]*math.cos(angle[1])
A = (math.cos(angle[0])-math.cos(angle[1])*math.cos(angle[2]))/math.sin(angle[2])
c_axis[1] = length[2]*A
c_axis[2] = length[2]*math.sqrt(1-(math.cos(angle[1]))**2-A**2)

#print("TV",'{:>14.10f}'.format(a_axis[0]),'{:>14.10f}'.format(a_axis[1]),'{:>14.10f}'.format(a_axis[2]))
#print("TV",'{:>14.10f}'.format(b_axis[0]),'{:>14.10f}'.format(b_axis[1]),'{:>14.10f}'.format(b_axis[2]))
#print("TV",'{:>14.10f}'.format(c_axis[0]),'{:>14.10f}'.format(c_axis[1]),'{:>14.10f}'.format(c_axis[2]))

m = np.identity(3)
m[0,0]=a_axis[0]; m[0,1]=a_axis[1]; m[0,2]=a_axis[2]
m[1,0]=b_axis[0]; m[1,1]=b_axis[1]; m[1,2]=b_axis[2]
m[2,0]=c_axis[0]; m[2,1]=c_axis[1]; m[2,2]=c_axis[2]

a = m[0,0]
v = m[1,:2]
b = np.linalg.norm(v)
cosg = m[1,0]/b
sing = m[1,1]/b
v = m[2,:]
c = np.linalg.norm(v)
cosb = m[2,0]/c
cosa = (m[2,0]*cosg+m[2,1]*sing)/c

print('a : ', a)
print('b : ', b)
print('c : ', c)
print('alpha : ', math.acos(cosa)*180./math.pi)
print('beta : ', math.acos(cosb)*180./math.pi)
print('gamma : ', math.acos(cosg)*180./math.pi)

