import os, sys, h5py, re
import numpy as np
import matplotlib.pyplot as plt

from .data import ResultsFile

# from ..computation.structure_factor import (histogram,gr,total_gr,SQ,total_SQ,total_FQ,
#                                     ncoeff,xcoeff,Gr,Tr,Nr)
from sovapy.computation.structure_factor import (atom_pair_hist, pair_dist_func, 
                                                 partial_structure_factor, atomc_pair_dist_func_neutron,
                                                 structure_factor_neutron, structure_factor_xray,
                                                 reduced_pair_dist_func, total_corr_fun, radial_dist_fun,
                                                 ncoeff, xcoeff)
from ..computation.structure_factor import bond_angle_hist
from ..computation.structure_factor import tetrahedra
from ..computation.rings import RINGs, Ring
from ..computation.cavity import Cavity

# Get the current version number:
cur_dirname =  os.path.dirname(__file__)
file_version = os.path.join(cur_dirname, '..', '__init__.py')
with open(file_version) as fv:
    VERSION = re.search("__version__ = '(.*)'", fv.read()).group(1)

class Analysis(object):

    def __init__(self, atoms=None):
        self.atoms = atoms #Atoms object
        pass

    # Input settings and run all computations for analysis
    def run(self):
        pass

    # Plot all analysis results
    def plot(self):
        pass

    # Save analysis results to hdf5
    def save_atoms_to_hdf5(self, mode="w"):

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, mode) as f:
            f.create_dataset('name', data=["sovapy"])
            f.create_dataset('version', data=[VERSION])

            group = f.create_group("atoms")
            group['positions'] = self.atoms.positions
            group['radii'] = self.atoms.radii
            group['symbols'] = self.atoms.symbols  # group['elements'] = self.atoms.elements.tolist()
            group['volume'] = [str(self.atoms.volume)]
            group['volume_origin'] = self.atoms.volume.origin.tolist()

            #Original text data
            if (self.atoms.original_file_data is not None):
                group = f.create_group('original_data')
                dt = h5py.string_dtype()
                ds_name = group.create_dataset('file_name', (1, ), dtype=dt, compression="gzip")
                ds_name[0] = self.atoms.original_file_data.file_name
                ds_txt = group.create_dataset('structure_text', (1, ), dtype=dt, compression="gzip")
                ds_txt[0] = self.atoms.original_file_data.structure_text
                
            if self.atoms.bond_lengths is not None:
                list_bond_lengths = list()
                for a in self.atoms.bond_lengths.items():
                    list_bond_lengths.append(list(a[0])+[str(a[1])])
                if "bond_lengths" in f["atoms"]:
                    del f["atoms"]["bond_lengths"]
                f["atoms"]["bond_lengths"] = list_bond_lengths

    # Save analysis results to hdf5
    def save_to_hdf5(self):
        pass 

    # Load analysis results from hdf5
    def load_hdf5(self):
        pass


class PDFAnalysis(Analysis):

    def __init__(self, atoms, dr=0.05, dq=0.05, qmin=0.3, qmax=25.0):
        #Atoms object
        self.atoms = atoms

        if not self.atoms.volume.periodic:
            print('Error: Periodicity is required for PDF analysis!')
            return 

        # Computation settings
        self.dr = dr                       
        self.dq = dq
        self.qmin = qmin
        self.qmax = qmax

        # Variables to save results
        self.r = None
        self.q = None
        self.pdf_atom_pairs = None
        self.partial_gr = None
        self.gr_neutron = None
        self.partial_sq = None
        self.sq_neutron = None
        self.sq_xray = None
        self.Gr_neutron = None
        self.Tr_neutron = None
        self.Nr_neutron = None

    # Input settings and run all computations for analysis
    def run(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot run PDF analysis because of no periodicity!')
            return 

        #Histogram of atom-pair distances
        # (The number of atoms at a distance between r and r+dr from a given atom.)
        self.r, hist = atom_pair_hist(self.atoms, self.dr)

        # Calculate Pair distribution function (PDF) g_{ab}(r) functions
        self.partial_gr = pair_dist_func(self.atoms, self.r, hist)

        #Calculate related to neutron diffraction
        coeff_neutron = ncoeff(self.atoms)

        # Calculate atomic pair distribution function for (neutron beam) g(r)
        self.gr_neutron = atomc_pair_dist_func_neutron(self.partial_gr, coeff_neutron)

        # Calculate Partial structure factors S{ab}(Q)
        self.q, self.partial_sq = partial_structure_factor(self.atoms, self.r, self.partial_gr, self.qmin, self.qmax, self.dq)

        # Calculate structure factor by neutron beam diffraction S_N(Q)
        self.sq_neutron = structure_factor_neutron(self.partial_sq, coeff_neutron)

        #Calculate related to X-ray diffraction
        coeff_xray = xcoeff(self.atoms, self.q)

        # Calculate structure factor by X-ray beam diffraction S_X(Q)
        self.sq_xray = structure_factor_xray(self.partial_sq, coeff_xray)

        # Atomic number density
        rho = self.atoms.atom_number_density

        # Reduced atomic pair distribution function by neutron beam G(r)
        self.Gr_neutron = reduced_pair_dist_func(self.r, self.gr_neutron, rho)

        # Total correlation function by neutron beam T(r)
        self.Tr_neutron = total_corr_fun(self.r, self.gr_neutron, rho)

        # Calculate radial_dist_fun by neutron beam N(r)
        self.Nr_neutron = radial_dist_fun(self.r, self.gr_neutron, rho)

    # Plot all analysis results
    def plot(self, figsize=(18, 8)):
        if not self.atoms.volume.periodic:
            print('Error: Cannot plot because any computed results were NOT found!')
            return 
        
        num_pairs = len(self.atoms.pairs)

        # Labels of atom pairs
        label_pairs = []
        for i in range(num_pairs):
            txt = self.atoms.pairs[i][0]+'-'+self.atoms.pairs[i][1]
            label_pairs.append(txt)

        fig = plt.figure(figsize=figsize) 
        ax = fig.add_subplot(2, 4, 1)
        for i in range(num_pairs):    
            ax.plot(self.r, self.partial_gr.T[i], label=label_pairs[i])
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Partial PDF g(r)')
        ax.legend()

        ax = fig.add_subplot(2, 4, 2)
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Atomic PDF (Neutron) g(r)')
        ax.plot(self.r, self.gr_neutron)

        ax = fig.add_subplot(2, 4, 3)
        for i in range(num_pairs):    
            ax.plot(self.q, self.partial_sq.T[i], label=label_pairs[i])
        ax.set_xlabel('Q (Angstrom^(-1))')
        ax.set_ylabel('Partial structure factor S(Q)')
        ax.legend()

        ax = fig.add_subplot(2, 4, 4)
        ax.set_xlabel('Q (Angstrom^(-1))')
        ax.set_ylabel('Structure factor by Neutron SN(Q)')
        ax.plot(self.q, self.sq_neutron)

        ax = fig.add_subplot(2, 4, 5)
        ax.set_xlabel('Q (Angstrom^(-1))')
        ax.set_ylabel('Structure factor by X-ray SX(Q)')
        ax.plot(self.q, self.sq_xray)

        ax = fig.add_subplot(2, 4, 6)
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Reduced atomic PDF G(r)')
        ax.plot(self.r, self.Gr_neutron)

        ax = fig.add_subplot(2, 4, 7)
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Total correlation function T(r)')
        ax.plot(self.r, self.Tr_neutron)

        ax = fig.add_subplot(2, 4, 8)
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Radial distribution function N(r)')
        ax.plot(self.r, self.Nr_neutron)

        plt.subplots_adjust(wspace=0.3)
        plt.subplots_adjust(hspace=0.3)
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot save to hdf file because any computed results were NOT found!')
            return 

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, 'a') as f:
            print('saving pdf.....')
            if 'pdf' in f:
                del f['pdf']
            group               = f.create_group("pdf")
            group['dr']         = self.dr
            group['dq']         = self.dq
            group['qmin']       = self.qmin
            group['qmax']       = self.qmax
            group['r']          = self.r
            group['q']          = self.q
            group['partial_gr'] = self.partial_gr
            group['gr_neutron'] = self.gr_neutron
            group['partial_sq'] = self.partial_sq
            group['sq_neutron'] = self.sq_neutron
            group['sq_xray']    = self.sq_xray
            group['Gr_neutron'] = self.Gr_neutron
            group['Tr_neutron'] = self.Tr_neutron
            group['Nr_neutron'] = self.Nr_neutron

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:   
            if 'pdf' in f:
                self.dr         = float(np.array(f['pdf']['dr']))
                self.dq         = float(np.array(f['pdf']['dr']))
                self.qmin       = float(np.array(f['pdf']['qmin']))
                self.qmax       = float(np.array(f['pdf']['qmax']))
                self.r          = np.array(f['pdf']['r']) 
                self.q          = np.array(f['pdf']['q'])
                self.partial_gr = np.array(f['pdf']['partial_gr'])
                self.gr_neutron = np.array(f['pdf']['gr_neutron'])
                self.partial_sq = np.array(f['pdf']['partial_sq'])
                self.sq_neutron = np.array(f['pdf']['sq_neutron'])
                self.sq_xray    = np.array(f['pdf']['sq_xray'])
                self.Gr_neutron = np.array(f['pdf']['Gr_neutron'])
                self.Tr_neutron = np.array(f['pdf']['Tr_neutron'])
                self.Nr_neutron = np.array(f['pdf']['Nr_neutron'])


class CoordinationNumberAnalysis(Analysis):

    def __init__(self, atoms):
        #Atoms object
        self.atoms = atoms

        # Variables to save results 
        self.atom_symbol_set = atoms.symbol_set #The set of atoms (elements)) 
        self.coord_num       = None  # List of coordination numbers
        self.list_counts     = None  # Counts (elems, coord_num)

    # Input settings and run all computations for analysis
    def run(self):

        # Element list
        # self.elems = list(set(self.atoms.elements))
        # self.symbol_set

        # List of coordination numbers
        cnums = np.array([len(b) for b in self.atoms.bonds])
        # Range of coordination numbers
        cmin, cmax  = cnums.min(), cnums.max() 
        # List of coordination numbers
        self.coord_num = np.arange(cmin, cmax+1)


        # Count coordination numbers of each element
        num_elems = len(self.atom_symbol_set)  # The number of elements (atoms)
        num_crange = cmax-cmin+1
        self.list_counts = np.zeros((num_elems, num_crange))
        for cnt, elem in enumerate(self.atom_symbol_set):
            ids = np.array([i for i, s in enumerate(self.atoms.symbols) if s == elem])
            cnums_elem = cnums[ids]
            j=0
            for cn in self.coord_num:
                self.list_counts[cnt,j] = np.sum(cnums_elem==cn)
                j += 1


    # Plot all analysis results
    def plot(self, figsize=None):
        
        num_elems = len(self.atom_symbol_set) # The number of elements (atoms)
        row_num = (num_elems//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        fig = plt.figure(figsize=figsize) 
        for i, elem in enumerate(self.atom_symbol_set):
            ax = fig.add_subplot(row_num, 3, i+1)
            plt.bar(self.coord_num, self.list_counts[i,:])
            plt.xlabel('Coordination number of {:} atom'.format(elem))
            plt.ylabel('Counts')
            plt.xticks(self.coord_num)
        plt.tight_layout()
        plt.show()
        

    # Save analysis results to hdf5
    def save_to_hdf5(self):

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, "a") as f:
            if 'coordination_number' in f:
                del f['coordination_number']         
            group = f.create_group("coordination_number")

            num_elems = len(self.atom_symbol_set)
            str_dtype = h5py.special_dtype(vlen=str)
            dataset = group.create_dataset('atom_symbol_set', shape=(num_elems,), dtype=str_dtype)
            for i, a in enumerate(self.atom_symbol_set):
                dataset[i] = a

            group['coord_num'] = self.coord_num
            group['list_counts'] = self.list_counts

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:                
            if 'coordination_number' in f:
                atom_symbols = list(f['coordination_number']['atom_symbol_set'])
                self.atom_symbol_set = [a.decode() for a in atom_symbols]
                self.coord_num = np.array(f['coordination_number']['coord_num'])
                self.list_counts = np.array(f['coordination_number']['list_counts'])


class BondAngleAnalysis(Analysis):

    def __init__(self, atoms, bins=100):
        #Atoms object
        self.atoms = atoms
        
        # Computation settings
        self.num_bins   = bins  # The number of bins (Old name: nth)

        # Variables to save results
        
        self.trios       = None  # List of atom triplets
        self.angles      = None  # Angle values  (Old name: cth)
        self.hist_angles = None  # The number of counts  (Old name: ncth)

    # Input settings and run all computations for analysis
    def run(self):

        # Calculate bond angles
        # (Angle values, histograms for all triplets)
        # self.cth, self.ncth = bond_angle_hist(self.atoms, nth=self.nth) 
        self.angles, self.hist_angles = bond_angle_hist(self.atoms, num_bins=self.num_bins, norm_sin=True, prob=True)

        # Remove redundant triplets that does not include in the setting
        self.trios = np.array(self.atoms.trios)
        flag_triplet = list()
        for hist in self.hist_angles:
            _sum = np.sum(np.array(hist))
            flag_triplet.append(_sum>0)
        self.trios = self.trios[flag_triplet]
        self.hist_angles = np.array(self.hist_angles)
        self.hist_angles = self.hist_angles[flag_triplet] 

        log_str = 'Calculating bond angle histogram has completed.\n'
        if len(self.trios)>0:
            log_str += 'Extracted atom triplets: {:}'.format(self.trios)
        else:
            log_str += 'There is no triplet! Check the bond setting.'
        print(log_str)


    # Plot all analysis results
    def plot(self, figsize=None):

        w = self.angles[1]-self.angles[0] # Angle interval
        num_trio = len(self.trios)
        row_num = (num_trio//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        fig = plt.figure(figsize=figsize) 
        for i, trio in enumerate(self.trios):
            ax = fig.add_subplot(row_num, 3, i+1)
            txt = trio[0] + '-' + trio[1] + '-' + trio[2]
            ax.bar(self.angles, self.hist_angles[i], width=w*0.8, label=txt)
            ax.set_xlim(0.0,180.0)
            # ax.set_ylim(0,1.0)
            ax.set_xlabel('Angle (Degree)')
            ax.set_ylabel('Probability')
            ax.legend()
        plt.subplots_adjust(wspace=0.3)
        plt.subplots_adjust(hspace=0.3)
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, "a") as f:
            if 'bond_angle' in f:
                del f['bond_angle']         
            group = f.create_group("bond_angle")
            group['num_bins'] = self.num_bins
            group['angles'] = self.angles
            group['hist_angles'] = self.hist_angles

            list_trios = list()
            for trio in self.trios:
                # list_trios.append(trio.split('-'))
                list_trios.append([str(s) for s in trio])
            group['trios'] = list_trios

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:                
            if 'bond_angle' in f:
                self.num_bins   = int(np.array(f['bond_angle']['num_bins']))
                self.angles   = np.array(f['bond_angle']['angles'])
                self.hist_angles  = np.array(f['bond_angle']['hist_angles'])
                list_trios = list(f['bond_angle']['trios'])
                trios = list()
                for trio in list_trios:
                    trio = [t.decode() for t in trio]
                    # trios.append('-'.join(trio))
                    trios.append(trio)
                self.trios = np.array(trios)


class TetrahedralOrderAnalysis(Analysis):

    def __init__(self, atoms, num_bins=100, list_cc_dist=None):
        #Atoms object
        self.atoms = atoms

        if not self.atoms.volume.periodic:
            print('Error: Periodicity is required for tetrahedral order analysis!')
            return 

        # Computation settings
        self.num_bins = num_bins #Old name: self.num_bins
        if list_cc_dist is not None:
            self.list_cc_dist = list_cc_dist
        else:
            print("Example of parameters of Si-O4 tetrahedron: [['Si','O', 2.0]]")

        # Variables to save results
        self.list_q = None
        self.list_idx_center = None

    # Input settings and run all computations for analysis
    def run(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot run tetrahedral order analysis because of no periodicity!')
            return
        self.list_q = list()
        self.list_idx_center = list()
        if self.list_cc_dist is None:
            print('ERROR: Set parameter list_cc_dist!')
            return
        for elem1, elem2, rmax in self.list_cc_dist:
            polys = tetrahedra(self.atoms, center=elem1, around=elem2, rmax=float(rmax))
            q = []
            idx_center = []
            for poly in polys:
                if poly.q is not None:
                    q.append(poly.q)
                    idx_center.append(poly.center)
            print('Number of tetrahedra ({:}-{:}4): {:d}'.format(elem1, elem2, len(idx_center))) 
            self.list_q.append(q)
            self.list_idx_center.append(idx_center)

    # Plot all analysis results
    def plot(self, figsize=None):
        if not self.atoms.volume.periodic:
            print('Error: Cannot plot because any computed results were NOT found!')
            return 

        num_tet = len(self.list_cc_dist)
        row_num = (num_tet//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        plt.figure(figsize=figsize) 
        for i in range(num_tet):
            plt.subplot(row_num, 3, i+1)
            label = self.list_cc_dist[i][0] + '-' + self.list_cc_dist[i][1] +'4'
            plt.hist(self.list_q[i], bins=self.num_bins, range=(0, 1), label=label)
            plt.xlim(0.0,1.0)
            plt.xlabel('Tetrahedral order (q-value)')
            plt.ylabel('Counts')
            plt.legend()
        plt.tight_layout()
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot save to hdf file because any computed results were NOT found!')
            return 

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, "a") as f:                
            if 'tetrahedra' in f:                
                del f['tetrahedra']
            group = f.create_group("tetrahedra")
            max_len = 0
            list_idx_center = list()
            for idx_center in self.list_idx_center:
                list_idx_center.append(idx_center.copy())
            list_q = list()
            for q in self.list_q:
                list_q.append(q.copy())

            #Add unusual values to be same length of lists
            for idx_center in list_idx_center:
                max_len=max([len(idx_center), max_len])
            for idx_center in list_idx_center:
                for i in range(max_len-len(idx_center)):
                    idx_center.append(-1) 
            for q in list_q:
                for i in range(max_len-len(q)):
                    q.append(10)
            group['list_idx_center'] = list_idx_center
            group['list_q']          = list_q
            list_cc_dist = list()
            for center,corner,dist in self.list_cc_dist:
                list_cc_dist.append([center,corner,str(dist)])
            group['list_cc_dist']    = list_cc_dist

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):

        with h5py.File(hdf5_path, "r") as f:                
            if 'tetrahedra' in f: 
                self.list_idx_center  = list(np.array(f['tetrahedra']['list_idx_center']))
                list_idx_center = list()
                for idx_center in self.list_idx_center:
                    list_idx_center.append(idx_center[idx_center>=0])
                self.list_idx_center = list_idx_center
                
                self.list_q = list(np.array(f['tetrahedra']['list_q']))
                list_q = list()
                for q in self.list_q:
                    list_q.append(q[q<1.1])
                self.list_q = list_q

                list_tmp = list(f['tetrahedra']['list_cc_dist'])
                list_cc_dist = list()
                for cc_dist in list_tmp:
                    cc_dist = [s.decode() for s in cc_dist]
                    cc_dist[2] = float(cc_dist[2])
                    list_cc_dist.append(cc_dist)
                self.list_cc_dist     = list_cc_dist
            else:
                print('ERROR: "tetrahedral order data does NOT contain!"')


class RingAnalysis(Analysis):

    def __init__(self, atoms, guttman=True, king=False, primitive=False, cutoff_primitive=24, num_parallel=-1, close=True):
        #Atoms object
        self.atoms = atoms

        if not self.atoms.volume.periodic:
            print('Error: Periodicity is required for Ring analysis!')
            return 
        
        # Computation settings
        self.flag_guttman = guttman
        self.flag_king = king
        self.flag_primitive = primitive
        self.cutoff_primitive = cutoff_primitive   # Cut-off ring size of primitive rings
        self.num_parallel = num_parallel # The number of cores in parallel computation (-1: single core) 
        self.flag_close = close # Enumerate rings closed in the real space


        # Variables to save results
        self.guttman_ring   = None #Output of sovapy.Computation.Rings.calculate()
        self.king_ring      = None #Output of sovapy.Computation.Rings.calculate()
        self.primitive_ring = None #Output of sovapy.Computation.Rings.calculate()

    # Input settings and run all computations for analysis
    def run(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot run Ring analysis because of no periodicity!')
            return 

        if self.flag_guttman:
            self.guttman_ring = RINGs(self.atoms)
            print('Calculating Guttman rings....')
            self.guttman_ring.calculate(ring_type=RINGs.RingType.GUTTMAN, num_parallel=self.num_parallel)
            if self.flag_close:
                rings = list()
                for ring in self.guttman_ring.rings:
                    if ring.close:
                        rings.append(ring)
                self.guttman_ring.rings = rings

            print('Done.\n\n')
        if self.flag_king:
            self.king_ring = RINGs(self.atoms)
            print('Calculating King rings....')
            self.king_ring.calculate(ring_type=RINGs.RingType.KING, num_parallel=self.num_parallel)
            if self.flag_close:
                rings = list()
                for ring in self.king_ring.rings:
                    if ring.close:
                        rings.append(ring)
                self.king_ring.rings = rings
            print('Done.\n\n')
        if self.flag_primitive:
            print('Calculating Primitive rings....')
            self.primitive_ring = RINGs(self.atoms)
            self.primitive_ring.calculate(ring_type=RINGs.RingType.PRIMITIVE,
                                                 cutoff_size=self.cutoff_primitive, num_parallel=self.num_parallel)
            if self.flag_close:
                rings = list()
                for ring in self.primitive_ring.rings:
                    if ring.close:
                        rings.append(ring)
                self.primitive_ring.rings = rings
            print('Done.\n\n')

    # Plot all analysis results
    def plot(self, figsize=None):
        if not self.atoms.volume.periodic:
            print('Error: Cannot plot because any computed results were NOT found!')
            return 
        row_num =  int(self.guttman_ring is not None)
        row_num += int(self.king_ring is not None)
        row_num += int(self.primitive_ring is not None)

        if figsize is None:
            figsize=(12, row_num*2.5)

        plt.figure(figsize=figsize) 
        cnt = 1
        if self.guttman_ring is not None:
            s_num, hist_size, r_roundness, r_roughness = self.calculate_ring_stats(self.guttman_ring.rings)
            plt.subplot(row_num, 3, cnt)
            plt.bar(s_num, hist_size)
            plt.xlabel('The number of atoms')
            plt.ylabel('Counts')
            plt.xticks(s_num)
            cnt += 1

            ### Roundness and roughness distributions
            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roundness, bins=np.linspace(0,1,20))
            plt.xlabel('Roundness')
            plt.ylabel('Counts')
            cnt += 1

            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roughness, bins=np.linspace(0,1,20))
            plt.xlabel('Roughness')
            plt.ylabel('Counts')
            cnt += 1

        if self.king_ring is not None:
            s_num, hist_size, r_roundness, r_roughness = self.calculate_ring_stats(self.king_ring.rings)
            plt.subplot(row_num, 3, cnt)
            plt.bar(s_num, hist_size)
            plt.xlabel('The number of atoms')
            plt.ylabel('Counts')
            plt.xticks(s_num)
            cnt += 1

            ### Roundness and roughness distributions
            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roundness, bins=np.linspace(0,1,20))
            plt.xlabel('Roundness')
            plt.ylabel('Counts')
            cnt += 1

            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roughness, bins=np.linspace(0,1,20))
            plt.xlabel('Roughness')
            plt.ylabel('Counts')
            cnt += 1
        
        if self.primitive_ring is not None:
            s_num, hist_size, r_roundness, r_roughness = self.calculate_ring_stats(self.primitive_ring.rings)
            plt.subplot(row_num, 3, cnt)
            plt.bar(s_num, hist_size)
            plt.xlabel('The number of atoms')
            plt.ylabel('Counts')
            plt.xticks(s_num)
            cnt += 1

            ### Roundness and roughness distributions
            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roundness, bins=np.linspace(0,1,20))
            plt.xlabel('Roundness')
            plt.ylabel('Counts')
            cnt += 1

            plt.subplot(row_num, 3, cnt)
            plt.hist(r_roughness, bins=np.linspace(0,1,20))
            plt.xlabel('Roughness')
            plt.ylabel('Counts')
            cnt += 1

        plt.tight_layout()
        plt.show()

    def calculate_ring_stats(self, rings):
        r_size, r_roundness, r_roughness = list(), list(), list()
        for r in rings:
            r_size.append(r.number) # the number of atoms in a ring
            r_roundness.append(r.roundness)
            r_roughness.append(r.roughness)
            
        r_size = np.array(r_size)
        r_roundness = np.array(r_roundness)
        r_roughness = np.array(r_roughness)

        #### Ring size distribution
        # Maximum ring size
        s_max = r_size.max()

        # Calculate the histogram of ring size
        hist_size = np.zeros(s_max +1, dtype='int')
        for s in range(s_max+1):
            hist_size[s] = np.sum(r_size==s)

        s_num = np.arange(s_max+1)
        return s_num, hist_size, r_roundness, r_roughness

    # Save analysis results to hdf5
    def save_to_hdf5(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot save to hdf file because any computed results were NOT found!')
            return 

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, "a") as f_hdf5:
            if self.guttman_ring is not None:
                if "rings_guttman" in f_hdf5:
                    del f_hdf5["rings_guttman"]
                group = f_hdf5.create_group("rings_guttman")
                # get max size
                max_size = 0
                for ring in self.guttman_ring.rings:
                    max_size = max(max_size, len(ring.indexes))
                irings = []
                for ring in self.guttman_ring.rings:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings
                
            if self.king_ring is not None:
                if "rings_king" in f_hdf5:
                    del f_hdf5["rings_king"]
                group = f_hdf5.create_group("rings_king")
                # get max size
                max_size = 0
                for ring in self.king_ring.rings:
                    max_size = max(max_size, len(ring.indexes))
                irings = []
                for ring in self.king_ring.rings:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings

            if self.primitive_ring is not None:
                if "rings_primitive" in f_hdf5:
                    del f_hdf5["rings_primitive"]
                group = f_hdf5.create_group("rings_primitive")
                # get max size
                max_size = 0
                for ring in self.primitive_ring.rings:
                    max_size = max(max_size, len(ring.indexes))
                irings = []
                for ring in self.primitive_ring.rings:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings
                
    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        
        with h5py.File(hdf5_path, "r") as f:
            # # Load rings
            if 'rings_guttman' in f:
                self.guttman_ring = RINGs(self.atoms)
                group = f['rings_guttman']
                irings = np.array(group['indexes']).tolist()
                self.guttman_ring.rings = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.guttman_ring.rings.append(ring)
            if 'rings_king' in f:
                self.king_ring = RINGs(self.atoms)
                group = f['rings_king']
                irings = np.array(group['indexes']).tolist()
                self.king_ring.rings = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.king_ring.rings.append(ring)
            if 'rings_primitive' in f:
                self.primitive_ring = RINGs(self.atoms)
                group = f['rings_primitive']
                irings = np.array(group['indexes']).tolist()
                self.primitive_ring.rings = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.primitive_ring.rings.append(ring)


class CavityAnalysis(Analysis):

    def __init__(self, atoms, resolution=128, cutoff_radii=2.8 ):
        #Atoms object
        self.atoms = atoms

        if not self.atoms.volume.periodic:
            print('Error: Periodicity is required for Cavity analysis!')
            return 
        
        # Computation settings
        self.resolution       = resolution
        if isinstance(cutoff_radii, dict):
            self.cutoff_radii = cutoff_radii
        elif isinstance(cutoff_radii, float) or isinstance(cutoff_radii, int):
            cutoff_radii = float(cutoff_radii)
            self.cutoff_radii = dict()
            # for elem in atoms.elements_kind:
            for elem in atoms.symbol_set:
                self.cutoff_radii[str(elem)] = cutoff_radii
        else:
            print('ERROR: cutoff_radii should be float or dict!')

        # Variables to save results
        self.cavity          = None 

    # Input settings and run all computations for analysis
    def run(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot run Cavity analysis because of no periodicity!')
            return 
        self.cavity = Cavity(self.atoms)
        self.cavity.calculate(resolution=self.resolution, cutoff_radii=self.cutoff_radii, 
                              gyration_tensor_parameters=True)

    # Plot all analysis results
    def plot(self, figsize=(20,12)):
        if not self.atoms.volume.periodic:
            print('Error: Cannot plot because any computed results were NOT found!')
            return 
        cavities = [self.cavity.domains, self.cavity.center_cavities, self.cavity.surface_cavities]
        names = ["Cavity domain", "Center cavity", "Surface cavity"]
        attrs = ['volumes', 'surface_areas', 'squared_gyration_radii', 
         'asphericities', 'acylindricities', 'anisotropies']

        plt.figure(figsize=figsize)
        cnt = 1
        for cavity, name in zip(cavities, names):
            for n, attr in enumerate(attrs):
                plt.subplot(3,6,cnt)
                vals = getattr(cavity, attr)
                plt.hist(vals)
                plt.xlabel(attr)
                plt.ylabel('Counts')
                if n==0:
                    plt.title(name)
                cnt += 1
        plt.tight_layout()
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):
        if not self.atoms.volume.periodic:
            print('Error: Cannot save to hdf file because any computed results were NOT found!')
            return 

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'
        
        with h5py.File(hdf5_path, "a") as f_hdf5:
            overwrite = True
            if self.cavity.domains is not None:
                if "cavity_setting" in f_hdf5:
                    del f_hdf5["cavity_setting"]
                group_setting = f_hdf5.create_group("cavity_setting")
                group_setting["resolution"] = self.cavity.resolution
                if self.cavity.cutoff_radii is None:
                    group_setting["cutoff_radii"] = -1
                elif type(self.cavity.cutoff_radii) is dict:
                    list_radii = list()
                    for k, d in self.cavity.cutoff_radii.items():
                        list_radii.append([k,str(d)])
                    group_setting["cutoff_radii"] = list_radii
                else:
                    group_setting["cutoff_radii"] = self.cavity.cutoff_radii

                if "domains" in f_hdf5:
                    del f_hdf5["domains"]
                group = f_hdf5.create_group("domains")
                self.cavity.domains.tohdf(group, overwrite)

                    # Critical-domains are not saved automatically in sovapy.
                if "critical_domains" in f_hdf5:
                    del f_hdf5["critical_domains"]
                f_hdf5["critical_domains"] = self.cavity.domains.critical_domains

            if self.cavity.surface_cavities is not None:
                if "surface_cavities" in f_hdf5:
                    del f_hdf5["surface_cavities"]
                group = f_hdf5.create_group("surface_cavities")
                self.cavity.surface_cavities.tohdf(group, overwrite)
            if self.cavity.surface_cavities is not None:
                if "center_cavities" in f_hdf5:
                    del f_hdf5["center_cavities"]
                group = f_hdf5.create_group("center_cavities")
                self.cavity.center_cavities.tohdf(group, overwrite)


    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):

        with ResultsFile(hdf5_path, 'r') as f:
            if f.cavity is not None:
                self.cavity = f.cavity
                