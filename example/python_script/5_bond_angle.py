from sovapy.core.file import File
from sovapy.computation.structure_factor import bond_angle_hist
import matplotlib.pyplot as plt

import sovapy
print('sovapy ver: ', sovapy.__version__)

# Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.get_atoms()

print("Atom symbols:", atoms.symbol_set)

# Build chemical bonds by setting maximum distances
bond_lengths = {('Si', 'O') : 1.77, ('Si', 'Si') : 1.32, ('O', 'O') : 2.22}
atoms.set_bond_lengths(bond_lengths) #build chemical bonds

# Calculate bond angle distributions
angle, hist = bond_angle_hist(atoms, num_bins=30, norm_sin=True, prob=True)

# Generate text labels of trios 
label_trios = list()
for trio in atoms.trios:
    label_trios.append(trio[0] + '-' + trio[1] + '-' + trio[2])

#Plot bond angle distributions
w = angle[1]-angle[0]
fig = plt.figure(figsize=(12, 5)) 
fig.suptitle("{} : TRIPLETS - bond angle correlations".format(structure_file))
for i, trio in enumerate(atoms.trios):
    ax = fig.add_subplot(2, 3, i+1)
    ax.bar(angle, hist[i], width=w*0.8, label=label_trios[i])
    ax.set_xlim(0.0,180.0)
    ax.set_ylim(0,1.0)
    ax.set_xlabel('Angle (Degree)')
    ax.set_ylabel('Probability')
    ax.legend()
plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.3)
plt.show()


# # Standard output 
# print('\nOutput from TRIPLETS - bond angle correlations')
# print()
# print('Maximum r values : ')
# s = ''
# for r in rmax:
#     s += "%17.8F" % r
# print(s)
# print('\n')
# for i, trio in enumerate(atoms.trios):
#     print('{0}\n'.format(trio))
#     print('Angle (degree), Probability ')
#     for x, y in zip(angle, hist[i]): 
#         print("%17.8F" % x + "%17.8F" % y)
#     print()
