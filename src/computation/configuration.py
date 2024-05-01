import math
import numpy as np

try:
    from computation.histogram import histogram as hist
    has_hist = True
except ImportError:
    import computation.calc_histogram as hist
    has_hist = False
    print('Warning calculate histogram is slow')


class Configuration(object):

    def __init__(self,ni,vectors,positions):
        self.ni = ni
        self.vectors = vectors
        self.positions = positions
        self.ntypes = len(ni)
        self.nmol = np.sum(self.ni)
        self.npar = int(self.ntypes*(self.ntypes+1)/2)
        self.r = None
        self.q = None
        self._gr = None
        self.total_gr = None
        self.sq = None
        self.total_sq = None
        self._metric = None
        self.truncated = False
        self._volume = None
        
    @property
    def metric(self):
        if self._metric is None:
            self._metric = np.identity(3)
            for i in range(3):
                for j in range(3):
                    self._metric[i][j] = 0.0
                    for k in range(3):
                        self._metric[i][j] = self._metric[i][j] + \
                            self.vectors[k][i] * self.vectors[k][j]

        return self._metric

    @property
    def volume(self):
        if self._volume is None:    
            triprod = self.vectors[0][0]*self.vectors[1][1]*self.vectors[2][2] \
                    + self.vectors[1][0]*self.vectors[2][1]*self.vectors[0][2] \
                    + self.vectors[2][0]*self.vectors[0][1]*self.vectors[1][2] \
                    - self.vectors[2][0]*self.vectors[1][1]*self.vectors[0][2] \
                    - self.vectors[1][0]*self.vectors[0][1]*self.vectors[2][2] \
                    - self.vectors[0][0]*self.vectors[2][1]*self.vectors[1][2]
            self._volume = 8.0*abs(triprod)
            if(self.truncated == True):
                self._volume = self._volume/2.0

        return self._volume

    @property
    def d(self):
        triprod = self.vectors[0][0]*self.vectors[1][1]*self.vectors[2][2] \
                + self.vectors[1][0]*self.vectors[2][1]*self.vectors[0][2] \
                + self.vectors[2][0]*self.vectors[0][1]*self.vectors[1][2] \
                - self.vectors[2][0]*self.vectors[1][1]*self.vectors[0][2] \
                - self.vectors[1][0]*self.vectors[0][1]*self.vectors[2][2] \
                - self.vectors[0][0]*self.vectors[2][1]*self.vectors[1][2]

        axb1=self.vectors[1][0]*self.vectors[2][1]-self.vectors[2][0]*self.vectors[1][1]
        axb2=self.vectors[2][0]*self.vectors[0][1]-self.vectors[0][0]*self.vectors[2][1]
        axb3=self.vectors[0][0]*self.vectors[1][1]-self.vectors[1][0]*self.vectors[0][1]
        bxc1=self.vectors[1][1]*self.vectors[2][2]-self.vectors[2][1]*self.vectors[1][2]
        bxc2=self.vectors[2][1]*self.vectors[0][2]-self.vectors[0][1]*self.vectors[2][2]
        bxc3=self.vectors[0][1]*self.vectors[1][2]-self.vectors[1][1]*self.vectors[0][2]
        cxa1=self.vectors[1][2]*self.vectors[2][0]-self.vectors[2][2]*self.vectors[1][0]
        cxa2=self.vectors[2][2]*self.vectors[0][0]-self.vectors[0][2]*self.vectors[2][0]
        cxa3=self.vectors[0][2]*self.vectors[1][0]-self.vectors[1][2]*self.vectors[0][0]
        d1 = triprod/math.sqrt(axb1**2+axb2**2+axb3**2)
        d2 = triprod/math.sqrt(bxc1**2+bxc2**2+bxc3**2)
        d3 = triprod/math.sqrt(cxa1**2+cxa2**2+cxa3**2)
        _d = min(d1,d2,d3)
        if (self.truncated):
            
            d1 = 1.5*triprod/math.sqrt( \
                (axb1+bxc1+cxa1)**2+(axb2+bxc2+cxa2)**2+(axb3+bxc3+cxa3)**2)
            d2 = 1.5*triprod/math.sqrt( \
                (axb1-bxc1+cxa1)**2+(axb2-bxc2+cxa2)**2+(axb3-bxc3+cxa3)**2)
            d3 = 1.5*triprod/math.sqrt( \
                (axb1+bxc1-cxa1)**2+(axb2+bxc2-cxa2)**2+(axb3+bxc3-cxa3)**2)
            d4 = 1.5*triprod/math.sqrt( \
                (axb1-bxc1-cxa1)**2+(axb2-bxc2-cxa2)**2+(axb3-bxc3-cxa3)**2)
            _d = math.min(_d,d1,d2,d3,d4)
        
        return _d
    
    @property
    def rho(self):        
        return self.nmol/self.volume

    def gr(self,dr,coeff):
        """
        calculate g(r)

        Parameters
        ----------
        dr : float
            delta r.
        coeff : list
            coefficient to calcuatate g(r).

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """
        nxn = np.zeros(self.npar)
        
        ic = 0
        for itype in range(self.ntypes):
            for jtype in range(itype, self.ntypes):
                nxn[ic] = self.ni[itype]*self.ni[jtype]
                if(itype != jtype):
                    nxn[ic] = nxn[ic]*2
                ic = ic+1
        
        _d = self.d
        _volume = self.volume
        nr = int(_d/dr)+1
        self.r =np.array([(float(i)*dr) for i in range(nr) ])
        truncated=False
        
        if has_hist == True:
            atoms = []
            for pos in self.positions:
                atoms.append(list(pos))
            _metric = []
            for m in self.metric:
                _metric.append(list(m))
            
            histogram = hist.calc_histogram(atoms,_metric,self.ni,_d,dr,truncated)
        else:
            # TODO modify calc_histogram
            histogram = hist.calc_histogram(cfg, dr)
        
        self._gr = np.zeros_like(histogram, dtype=float)
        self.total_gr = np.zeros(nr)

        for ir in range(1, nr):   
            gnorm = (3.0*ir*ir+0.25)*dr*dr*dr*2.0*math.pi/(3.0*_volume)
            for ic in range(self.npar):
                if has_hist== True:
                    self._gr[ir,ic] = histogram[ir][ic]/(gnorm*nxn[ic])
                else:
                    self._gr[ir,ic] = histogram[ir,ic]/(gnorm*nxn[ic])

        for ir in range(nr):
            for ic in range(self.npar):
                self.total_gr[ir] += coeff[ic]*(self._gr[ir,ic]-1.0)

        return self.r, self._gr, self.total_gr
        
    def SQ(self,dr,dq,qmin,qmax,coeff):
        """
        calculate S(Q)

        Parameters
        ----------
        dr : float
            delta radius.
        dq : float
            delta q.
        qmin : float
            minimum Q.
        qmax : float
            maximum Q.
        coeff : list
            coefficient to calcuatate S(Q).

        Returns
        -------
        q : TYPE
            DESCRIPTION.
        sq : TYPE
            DESCRIPTION.
        sq_tot : TYPE
            DESCRIPTION.

        """
        nq = int(math.ceil((qmax-qmin)/dq))+1
        self.q =[ (qmin+float(i)*dq) for i in range(nq) ]
        nr = self._gr.shape[0]
        sqr = np.zeros((nr, nq+1), dtype=float)
        #s = np.zeros_like(sqr)

        self.sq = np.zeros((nq, self.npar), dtype=float)
        self.total_sq = np.zeros(nq)

        n = self.nmol
        volume = self.volume
        
        for iq in range(nq):
            s = np.zeros(self.npar)
            for ir in range(1, nr):
                r = float(ir)*dr
                sqr = 4.0*np.pi*float(n)/volume*r*np.sin(r*self.q[iq])/self.q[iq]*dr
                    
                for ic in range(self.npar):
                    s[ic] += (self._gr[ir, ic]-1.0)*sqr
            
            self.sq[iq] = s

        for iq in range(nq):
            for ic in range(self.npar):
                self.total_sq[iq] += coeff[ic]*self.sq[iq,ic]

        return self.q, self.sq, self.total_sq

class DistributionSettings(object):
    
    def __init__(self, partial_gr=False,partial_sq=False):
        self.partial_gr = partial_gr
        self.partial_sq = partial_sq
        self.neutron_sq = True
        self.neutron_gr = True
        self.xray_sq = True
        self.Gr = True
        self.Tr = True
        self.Nr = True
        self.dr = 0.05
        self.dq = 0.05
        self.qmin = 0.3
        self.qmax = 25.0
        