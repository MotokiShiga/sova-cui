import os, h5py, re
import numpy as np
import matplotlib.pyplot as plt

from .data import ResultsFile

from ..computation.structure_factor import (histogram,gr,total_gr,SQ,total_SQ,total_FQ,
                                    ncoeff,xcoeff,Gr,Tr,Nr)
from ..computation.structure_factor import triplets
from ..computation.structure_factor import polyhedra
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
            group['elements'] = self.atoms.elements.tolist()
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

        # Computation settings
        self.dr = dr                       
        self.dq = dq
        self.qmin = qmin
        self.qmax = qmax

        # Variables to save results
        self.r = None
        self.q = None
        self.gr = None
        self.total_gr = None
        self.sq = None
        self.total_sq = None
        self.fq = None
        self.Gr = None
        self.Tr = None
        self.Nr = None

    # Input settings and run all computations for analysis
    def run(self):
        # Histograms of distances between atom pairs
        r, hist = histogram(self.atoms, self.dr)

        # Calculate PDF functions
        # Calculate g(r)
        self.r, self.gr = gr(self.atoms, hist, self.dr)

        # Calculate Total g(r)
        coeff = ncoeff(self.atoms.symbols, self.atoms.frac)
        self.total_gr = total_gr(self.gr, coeff)

        # Calculate S(Q)
        self.q, self.sq = SQ(self.atoms, self.gr, self.qmin, self.qmax, self.dr, self.dq)

        # Calculate Total S(Q)
        self.total_sq = total_SQ(self.sq, coeff)

        # Calculate F(Q)
        coeff = xcoeff(self.atoms.symbols, self.atoms.frac, self.q)
        self.fq = total_FQ(self.sq, coeff)

        # Calculate Gr
        rho = self.atoms.rho
        self.Gr = Gr(self.r, self.total_gr, rho)

        # Calculate Tr
        self.Tr = Tr(self.r, self.total_gr, rho)

        # Calculate Nr
        self.Nr = Nr(self.r, self.Tr)

    # Plot all analysis results
    def plot(self, figsize=(18, 8)):
        fig = plt.figure(figsize=figsize) 
        ax = fig.add_subplot(2, 4, 1)
        for i in range(3):    
            ax.plot(self.r, self.gr.T[i], label=self.atoms.pairs[i])
            ax.set_xlabel('r(Å)')
            ax.set_ylabel('Partial g(r)')
            ax.legend()

        ax = fig.add_subplot(2, 4, 2)
        ax.set_xlabel('r(Å)')
        ax.set_ylabel('Total g(r)')
        ax.plot(self.r, self.total_gr)

        ax = fig.add_subplot(2, 4, 3)
        for i in range(3):    
            ax.plot(self.q, self.sq.T[i], label=self.atoms.pairs[i])
            ax.set_xlabel('Q(Å^-1)')
            ax.set_ylabel('Partial S(Q)')
            ax.legend()

        ax = fig.add_subplot(2, 4, 4)
        ax.set_xlabel('Q(Å^-1)')
        ax.set_ylabel('Total Neutron S(Q)')
        ax.plot(self.q, self.total_sq)

        ax = fig.add_subplot(2, 4, 5)
        ax.set_xlabel('Q(Å^-1)')
        ax.set_ylabel('Total X-ray S(Q)')
        ax.plot(self.q, self.fq)

        ax = fig.add_subplot(2, 4, 6)
        ax.set_xlabel('r(Å)')
        ax.set_ylabel('G(r)')
        ax.plot(self.r, self.Gr)

        ax = fig.add_subplot(2, 4, 7)
        ax.set_xlabel('r(Å)')
        ax.set_ylabel('T(r)')
        ax.plot(self.r, self.Tr)

        ax = fig.add_subplot(2, 4, 8)
        ax.set_xlabel('r(Å)')
        ax.set_title('N(r)')
        ax.plot(self.r, self.Nr)

        plt.subplots_adjust(wspace=0.3)
        plt.subplots_adjust(hspace=0.3)
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):

        ori_fname = self.atoms.original_file_data.file_name
        hdf5_path = os.path.splitext(ori_fname)[0] + '.hdf5'

        with h5py.File(hdf5_path, 'a') as f:
            print('saving pdf.....')
            if 'pdf' in f:
                del f['pdf']
            group             = f.create_group("pdf")
            group['dr']       = self.dr
            group['dq']       = self.dq
            group['qmin']     = self.qmin
            group['qmax']     = self.qmax
            group['r']        = self.r
            group['q']        = self.q
            group['gr']       = self.gr
            group['total_gr'] = self.total_gr
            group['sq']       = self.sq
            group['total_sq'] = self.total_sq
            group['fq']       = self.fq
            group['Gr']       = self.Gr
            group['Tr']       = self.Tr
            group['Nr']       = self.Nr

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:   
            if 'pdf' in f:
                self.dr = float(np.array(f['pdf']['dr']))
                self.dq = float(np.array(f['pdf']['dr']))
                self.qmin = float(np.array(f['pdf']['qmin']))
                self.qmax = float(np.array(f['pdf']['qmax']))
                self.r = np.array(f['pdf']['r']) 
                self.q = np.array(f['pdf']['q'])
                self.gr = np.array(f['pdf']['gr'])
                self.total_gr = np.array(f['pdf']['total_gr'])
                self.sq = np.array(f['pdf']['sq'])
                self.total_sq = np.array(f['pdf']['total_sq'])
                self.fq = np.array(f['pdf']['fq'])
                self.Gr = np.array(f['pdf']['Gr'])
                self.Tr = np.array(f['pdf']['Tr'])
                self.Nr = np.array(f['pdf']['Nr'])

class CoordinationNumberAnalysis(Analysis):

    def __init__(self, atoms):
        #Atoms object
        self.atoms = atoms

        # Variables to save results
        
        self.elems        = None  # List of elements (atoms)
        self.coord_num    = None  # List of coordination numbers
        self.list_counts  = None  # Counts (elems, coord_num)

    # Input settings and run all computations for analysis
    def run(self):

        # Element list
        self.elems = list(set(self.atoms.elements))

        # List of coordination numbers
        cnums = np.array([len(b) for b in self.atoms.bonds])
        # Range of coordination numbers
        cmin, cmax  = cnums.min(), cnums.max() 
        # List of coordination numbers
        self.coord_num = np.arange(cmin, cmax+1)


        # Count coordination numbers of each element
        num_elems = len(self.elems)
        num_crange = cmax-cmin+1
        self.list_counts = np.zeros((num_elems,num_crange))
        for cnt, elem in enumerate(self.elems):
            ids = np.array([i for i, s in enumerate(self.atoms.elements) if s == elem])
            cnums_elem = cnums[ids]
            j=0
            for cn in self.coord_num:
                self.list_counts[cnt,j] = cn
                self.list_counts[cnt,j] = np.sum(cnums_elem==cn)
                j += 1


    # Plot all analysis results
    def plot(self, figsize=None):
        
        num_elems = len(self.elems)
        row_num = (num_elems//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        fig = plt.figure(figsize=figsize) 
        for i, elem in enumerate(self.elems):
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

            num_elems = len(self.elems)
            str_dtype = h5py.special_dtype(vlen=str)
            dataset = group.create_dataset('elements', shape=(num_elems,), dtype=str_dtype)
            for i, a in enumerate(self.elems):
                dataset[0] = a

            group['coord_num'] = self.coord_num
            group['list_counts'] = self.list_counts

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:                
            if 'coordination_number' in f:
                elems = list(f['coordination_number']['elements'])
                self.elems = [a.decode() for a in elems]
                self.coord_num = np.array(f['coordination_number']['coord_num'])
                self.list_counts = np.array(f['coordination_number']['list_counts'])


class BondAngleAnalysis(Analysis):

    def __init__(self, atoms, bins=100):
        #Atoms object
        self.atoms = atoms
        
        # Computation settings
        self.nth   = bins  # The number of bins

        # Variables to save results
        
        self.trios = None  # List of atom triplets
        self.cth   = None  # Angle values (nth +1)
        self.ncth  = None  # The number of counts

    # Input settings and run all computations for analysis
    def run(self):

        # Calculate bond angles
        # (Angle values, histograms for all triplets)
        self.cth, self.ncth = triplets(self.atoms, nth=self.nth) 
        
        # Remove redundant triplets that does not include in the setting
        self.trios = np.array(self.atoms.trios)
        flag_triplet = list()
        for hist in self.ncth:
            _sum = np.sum(np.array(hist))
            flag_triplet.append(_sum>0)
        self.trios = self.trios[flag_triplet]
        self.ncth = np.array(self.ncth)
        self.ncth = self.ncth[flag_triplet] 

        log_str = 'Calculating bond angle histogram has completed.\n'
        if len(self.trios)>0:
            log_str += 'Extracted atom triplets: {:}'.format(self.trios)
        else:
            log_str += 'There is no triplet! Check the bond setting.'
        print(log_str)

        # Save results in shared object
        self.nth   = self.nth
        self.trios = self.trios
        self.cth   = self.cth
        self.ncth  = self.ncth

    # Plot all analysis results
    def plot(self, figsize=None):

        w = self.cth[1]-self.cth[0] # Angle interval
        num_trio = len(self.trios)
        row_num = (num_trio//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        fig = plt.figure(figsize=figsize) 
        for i, trio in enumerate(self.trios):
            ax = fig.add_subplot(row_num, 3, i+1)
            ax.bar(self.cth, self.ncth[i], width=w*0.8, label=trio)
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
            group['nth'] = self.nth
            group['cth'] = self.cth
            group['ncth'] = self.ncth

            list_trios = list()
            for trio in self.trios:
                list_trios.append(trio.split('-'))
            group["trios"] = list_trios 

    # Load analysis results from hdf5
    def load_hdf5(self, hdf5_path):
        with h5py.File(hdf5_path, "r") as f:                
            if 'bond_angle' in f:
                self.nth   = int(np.array(f['bond_angle']['nth']))
                self.cth   = np.array(f['bond_angle']['cth'])
                self.ncth  = np.array(f['bond_angle']['ncth'])
                list_trios = list(f['bond_angle']['trios'])
                trios = list()
                for trio in list_trios:
                    trio = [t.decode() for t in trio]
                    trios.append('-'.join(trio))
                self.trios = np.array(trios)


class TetrahedralOrderAnalysis(Analysis):

    def __init__(self, atoms, bins=100, list_cc_dist=None):
        #Atoms object
        self.atoms = atoms

        # Computation settings
        self.bins = bins
        if list_cc_dist is not None:
            self.list_cc_dist = list_cc_dist
        else:
            print("Example of parameters of Si-O4 tetrahedron: [['Si','O', 2.0]]")

        # Variables to save results
        self.list_q = None
        self.list_idx_center = None

    # Input settings and run all computations for analysis
    def run(self):
        self.list_q = list()
        self.list_idx_center = list()
        if self.list_cc_dist is None:
            print('ERROR: Set parameter list_cc_dist!')
            return
        for elem1, elem2, rmax in self.list_cc_dist:
            polys = polyhedra(self.atoms, center=elem1, around=elem2, rmax=float(rmax))
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

        num_tet = len(self.list_cc_dist)
        row_num = (num_tet//3)+1

        if figsize is None:
            figsize=(12, row_num*2.5)

        plt.figure(figsize=figsize) 
        for i in range(num_tet):
            plt.subplot(row_num, 3, i+1)
            label = self.list_cc_dist[i][0] + '-' + self.list_cc_dist[i][1] +'4'
            plt.hist(self.list_q[i], bins=self.bins, range=(0, 1), label=label)
            plt.xlim(0.0,1.0)
            plt.xlabel('Tetrahedral order (q-value)')
            plt.ylabel('Counts')
            plt.legend()
        plt.tight_layout()
        plt.show()

    # Save analysis results to hdf5
    def save_to_hdf5(self):

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
        
        # Computation settings
        self.resolution       = resolution
        if isinstance(cutoff_radii, dict):
            self.cutoff_radii = cutoff_radii
        elif isinstance(cutoff_radii, float) or isinstance(cutoff_radii, int):
            cutoff_radii = float(cutoff_radii)
            self.cutoff_radii = dict()
            for elem in atoms.elements_kind:
                self.cutoff_radii[str(elem)] = cutoff_radii
        else:
            print('ERROR: cutoff_radii should be float or dict!')

        # Variables to save results
        self.cavity          = None 

    # Input settings and run all computations for analysis
    def run(self):
        self.cavity = Cavity(self.atoms)
        self.cavity.calculate(resolution=self.resolution, cutoff_radii=self.cutoff_radii, 
                              gyration_tensor_parameters=True)

    # Plot all analysis results
    def plot(self, figsize=(20,12)):
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
                