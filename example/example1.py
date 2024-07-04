# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.structure_factor import (histogram,gr,total_gr,SQ,total_SQ,total_FQ,
                                               ncoeff,xcoeff,Gr,Tr,Nr)
import matplotlib.pyplot as plt

path = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()
dr = 0.05

# calculate histogram
#_, hist = histogram(atoms,dr,symbols=['Si','O'])
r, hist = histogram(atoms,dr)

fig = plt.figure(figsize=(12, 4)) 
for i in range(3):
    ax = fig.add_subplot(1, 3, i+1)
    ax.bar(r, hist.T[i], width=dr*0.8, label=atoms.pairs[i])
    ax.set_xlim(0.0,5.0)
    ax.set_ylim(0,500)
    ax.set_xlabel('r(Å)')
    ax.set_ylabel('Number')
    ax.legend()
plt.subplots_adjust(wspace=0.3)
plt.show()

# calculate g(r)
r, gr = gr(atoms,hist,dr)

# calculate Total g(r)
coeff = ncoeff(atoms.symbols,atoms.frac)
total_gr = total_gr(gr,coeff)

# calculate S(Q)
dq = 0.05
qmin = 0.3
qmax = 25.0
q, sq = SQ(atoms,gr,qmin,qmax,dr,dq)

# calculate Total S(Q)
total_sq = total_SQ(sq,coeff)

# calculate F(Q)
coeff = xcoeff(atoms.symbols,atoms.frac,q)
total_fq = total_FQ(sq,coeff)

# calculate Gr
rho = atoms.rho
_Gr = Gr(r,total_gr,rho)

# calculate Tr
_Tr = Tr(r,total_gr,rho)

# calculate Nr
_Nr = Nr(r,_Tr)

# show graph
fig = plt.figure(figsize=(18, 8)) 
ax = fig.add_subplot(2, 4, 1)
for i in range(3):    
    ax.plot(r, gr.T[i], label=atoms.pairs[i])
    ax.set_xlabel('r(Å)')
    ax.set_ylabel('Partial g(r)')
    ax.legend()

ax = fig.add_subplot(2, 4, 2)
ax.set_xlabel('r(Å)')
ax.set_ylabel('Total g(r)')
ax.plot(r,total_gr)

ax = fig.add_subplot(2, 4, 3)
for i in range(3):    
    ax.plot(q, sq.T[i], label=atoms.pairs[i])
    ax.set_xlabel('Q(Å^-1)')
    ax.set_ylabel('Partial S(Q)')
    ax.legend()

ax = fig.add_subplot(2, 4, 4)
ax.set_xlabel('Q(Å^-1)')
ax.set_ylabel('Total Neutron S(Q)')
ax.plot(q, total_sq)

ax = fig.add_subplot(2, 4, 5)
ax.set_xlabel('Q(Å^-1)')
ax.set_ylabel('Total X-ray S(Q)')
ax.plot(q, total_fq)

ax = fig.add_subplot(2, 4, 6)
ax.set_xlabel('r(Å)')
ax.set_ylabel('G(r)')
ax.plot(r, _Gr)

ax = fig.add_subplot(2, 4, 7)
ax.set_xlabel('r(Å)')
ax.set_ylabel('T(r)')
ax.plot(r, _Tr)

ax = fig.add_subplot(2, 4, 8)
ax.set_xlabel('r(Å)')
ax.set_title('N(r)')
ax.plot(r, _Nr)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.3)
plt.show()
