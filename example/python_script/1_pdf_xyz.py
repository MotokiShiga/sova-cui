from sovapy.core.file import File
from sovapy.computation.structure_factor import (histogram,gr,total_gr,SQ,total_SQ,total_FQ,
                                    ncoeff,xcoeff,Gr,Tr,Nr)
import matplotlib.pyplot as plt


# Load structural information from a xyz file
# The second line (CUB 24.713) in the xyz file indicates the shape of cell and its length.
# (CUB means cubics.)
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

print("Atom symbols:", atoms.symbols)

print("Is the periodicity information of the cell available?:")
print(atoms.volume.periodic)
# SOVA requires the periodicity information to compute PDF functions.
# Only histograms of distances between atom pairs can be computed without it.


# Histograms of distances between atom pairs
dr = 0.05   # bin width 
r, hist = histogram(atoms,dr) # calculate histograms
# Input symbols option to determine plot order of atoms 
#r, hist = histogram(atoms,dr,symbols=['Si','O'])

# Plot histograms of pair distance
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


# Calculate PDF functions
# Calculate g(r)
r, gr = gr(atoms,hist,dr)

# Calculate Total g(r)
coeff = ncoeff(atoms.symbols,atoms.frac)
total_gr = total_gr(gr,coeff)

# Calculate S(Q)
dq = 0.05
qmin = 0.3
qmax = 25.0
q, sq = SQ(atoms,gr,qmin,qmax,dr,dq)

# Calculate Total S(Q)
total_sq = total_SQ(sq,coeff)

# Calculate F(Q)
coeff = xcoeff(atoms.symbols,atoms.frac,q)
total_fq = total_FQ(sq,coeff)

# Calculate Gr
rho = atoms.rho
_Gr = Gr(r,total_gr,rho)

# Calculate Tr
_Tr = Tr(r,total_gr,rho)

# Calculate Nr
_Nr = Nr(r,_Tr)


# Plot functions g(r), total g(r), et al.
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
