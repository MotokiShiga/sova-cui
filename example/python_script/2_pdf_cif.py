from sovapy.core.file import File
from sovapy.computation.structure_factor import (atom_pair_hist, pair_dist_func, 
                                                 partial_structure_factor, atomc_pair_dist_func_neutron,
                                                 structure_factor_neutron, structure_factor_xray,
                                                 reduced_pair_dist_func, total_corr_fun, radial_dist_fun,
                                                 ncoeff, xcoeff)
import matplotlib.pyplot as plt

import sovapy
print('sovapy ver: ', sovapy.__version__)

# Load structural information from a cif file
structure_file = "../data/crystal/sio2_beta_cristobalite222.cif"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

print("Atom symbols:", atoms.symbols)

# CIF files have periodicity information.
print("Is the periodicity information of the cell available?:")
print(atoms.volume.periodic)

#Histogram of atom-pair distances
# (The number of atoms at a distance between r and r+dr from a given atom.)
dr = 0.05   # bin width 
r, hist = atom_pair_hist(atoms, dr) # calculate histograms
# Input symbols option to determine plot order of atoms 
#r, hist = histogram(atoms,dr,symbols=['Si','O'])

# Plot histograms of pair distance
fig = plt.figure(figsize=(12, 4)) 
for i in range(3):
    ax = fig.add_subplot(1, 3, i+1)
    ax.bar(r, hist.T[i], width=dr*0.8, label=atoms.pairs[i])
    ax.set_xlim(0.0,5.0)
    ax.set_ylim(0,500)
    ax.set_xlabel('r (Angstrom)')
    ax.set_ylabel('Number of atoms')
    ax.legend()
plt.subplots_adjust(wspace=0.3)
plt.show()


# Calculate Pair distribution function (PDF) g_{ab}(r) functions
r, partial_gr = pair_dist_func(atoms, hist, dr)

#Calculate related to neutron diffraction
coeff_neutron = ncoeff(atoms.symbols, atoms.frac)

# Calculate atomic pair distribution function for (neutron beam) g(r)
gr_neutron = atomc_pair_dist_func_neutron(partial_gr, coeff_neutron)

# Calculate Partial structure factors S{ab}(Q)
dq = 0.05
qmin = 0.3
qmax = 25.0
q, partial_sq = partial_structure_factor(atoms, partial_gr, qmin, qmax, dr, dq)

# Calculate structure factor by neutron beam diffraction S_N(Q)
sq_neutron = structure_factor_neutron(partial_sq, coeff_neutron)

#Calculate related to X-ray diffraction
coeff_xray = xcoeff(atoms.symbols, atoms.frac,q)

# Calculate structure factor by X-ray beam diffraction S_X(Q)
sq_xray = structure_factor_xray(partial_sq, coeff_xray)

# Atomic number density
rho = atoms.rho 

# Reduced atomic pair distribution function by neutron beam G(r)
Gr_neutron = reduced_pair_dist_func(r, gr_neutron, rho)

# Total correlation function by neutron beam T(r)
Tr_neutron = total_corr_fun(r, gr_neutron, rho)

# Calculate radial_dist_fun by neutron beam N(r)
Nr_neutron = radial_dist_fun(r, gr_neutron, rho)

# Plot functions g(r), total g(r), et al.
fig = plt.figure(figsize=(18, 8)) 
ax = fig.add_subplot(2, 4, 1)
for i in range(3):    
    ax.plot(r, partial_gr.T[i], label=atoms.pairs[i])
    ax.set_xlabel('r (Angstrom)')
    ax.set_ylabel('Partial g(r)')
    ax.legend()

ax = fig.add_subplot(2, 4, 2)
ax.set_xlabel('r (Angstrom)')
ax.set_ylabel('Atomic PDF (Neutron) g(r)')
ax.plot(r, gr_neutron)

ax = fig.add_subplot(2, 4, 3)
for i in range(3):    
    ax.plot(q, partial_sq.T[i], label=atoms.pairs[i])
    ax.set_xlabel('Q (Angstrom^(-1))')
    ax.set_ylabel('Partial structure factor S(Q)')
    ax.legend()

ax = fig.add_subplot(2, 4, 4)
ax.set_xlabel('Q (Angstrom^(-1))')
ax.set_ylabel('Structure factor by Neutron SN(Q)')
ax.plot(q, sq_xray)

ax = fig.add_subplot(2, 4, 5)
ax.set_xlabel('Q (Angstrom^(-1))')
ax.set_ylabel('Structure factor by X-ray SX(Q)')
ax.plot(q, sq_xray)

ax = fig.add_subplot(2, 4, 6)
ax.set_xlabel('r (Angstrom)')
ax.set_ylabel('Reduced atomic PDF G(r)')
ax.plot(r, Gr_neutron)

ax = fig.add_subplot(2, 4, 7)
ax.set_xlabel('r (Angstrom)')
ax.set_ylabel('Total correlation function T(r)')
ax.plot(r, Tr_neutron)

ax = fig.add_subplot(2, 4, 8)
ax.set_xlabel('r (Angstrom)')
ax.set_ylabel('Radial distribution function N(r)')
ax.plot(r, Nr_neutron)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.3)
plt.show()
