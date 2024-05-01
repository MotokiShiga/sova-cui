# -*- coding: utf-8 -*-

import numpy as np
import math
from . import element_data as elem
import sys

def calc_xcoeff(frac, symbol, q, option):
    """
    Options for ceofficients are: 
        1 - Use as calculated
        2 - Divide by <f**2)
        3 - Divide by <f>**2
    """
    
    nq = len(q)
    ntypes = len(symbol)
    npar = int(ntypes*(ntypes+1)/2)
    ff = [0.0]*npar
    xcoeff = [0.0]*nq
    
    for iq in range(nq):    
        q2 = q[iq] * q[iq] / (16.0 * math.pi * math.pi)
        fsq = 0.0
        f = 0.0
        local_xcoeff = [0.0]*npar
        
        for itype in range(ntypes):
            e = symbol[itype]
            a = elem.Xray[e].a
            b = elem.Xray[e].b
            c = elem.Xray[e].c
            
            _ff = 0.0
            for n in range(5):
                _ff = _ff + a[n]*math.exp(-b[n]*q2)
            _ff = _ff + c
            ff[itype] = _ff
            
            f = f + frac[itype]*_ff
            fsq = fsq + frac[itype]*_ff*_ff
         
        for itype in range(ntypes):
            for jtype in range(itype, ntypes):
                ic = int(itype*(2*ntypes-itype-1)/2) + jtype

                local_xcoeff[ic] = frac[itype]*frac[jtype]*ff[itype]*ff[jtype]
                if (itype != jtype):
                    local_xcoeff[ic] = local_xcoeff[ic]*2.0
                if (option == 2):
                    local_xcoeff[ic] = local_xcoeff[ic]/fsq
                if (option == 3):
                    local_xcoeff[ic] = local_xcoeff[ic]/(f*f)
                  
        xcoeff[iq] = local_xcoeff 
        
    return xcoeff
    

def coeff(frac, elements, in_path, out_path, option):
    """
    Options for ceofficients are: 
        1 - Use as calculated
        2 - Divide by <f**2)
        3 - Divide by <f>**2
    """
    q = []
    s = []
    try:
        with open(in_path.encode("utf-8"), "r") as f:
            try:
                nq = int(f.next())
                line = f.next()
                for i in range(nq):
                    item = f.next()
                    _q, _s = item.split()[:2]
                    q.append(float(_q))
                    s.append(float(_s))
            except StopIteration:
                pass
    except IOError:
        raise

    ntypes = len(frac)
    npar = ntypes*(ntypes+1)/2
    ff = [0.0]*npar
    coeff = [0.0]*npar

    fw = open(out_path, 'w')
    fw.write("%5d\n" % nq)
    fw.write("%s" % line)

    for iq in range(nq):

        q2 = q[iq] * q[iq] / (16.0 * math.pi * math.pi)
        fsq = 0.0
        f = 0.0
        for itype in range(ntypes):
            e = elements[itype]
            a = elem.Xray[e].a
            b = elem.Xray[e].b
            c = elem.Xray[e].c
            _ff = 0.0
            for n in range(5):
                _ff = _ff + a[n]*math.exp(-b[n]*q2)
            _ff = _ff + c
            ff[itype] = _ff
            f = f + frac[itype]*_ff
            fsq = fsq + frac[itype]*_ff*_ff
         
        for itype in range(ntypes):
            for jtype in range(itype, ntypes):
                ic = itype*(2*ntypes-itype-1)/2+jtype
                coeff[ic] = frac[itype]*frac[jtype]*ff[itype]*ff[jtype]
                if (itype != jtype):
                    coeff[ic] = coeff[ic]*2.0
                if (option == 2):
                    coeff[ic] = coeff[ic]/fsq
                if (option == 3):
                    coeff[ic] = coeff[ic]/(f*f)

        fw.write("%12.7F %12.7F" % (q[iq], s[iq]))
        for ic in range(npar):
            fw.write(" %12.7F" % coeff[ic])
        fw.write("\n")

    fw.close()

if __name__ == '__main__':

    print("  1: coeff \n  2:calc_xcoeff")
    call = input("Select a call function. :  ")
#    ntypes=3

#    for itype in range(ntypes):
#        for jtype in range(itype, ntypes):
#            ic = itype*(2*ntypes-itype-1)/2+jtype
#            print itype, jtype, ic
    
    #sys.exit()C:\python\pyModelling\rmc
      
    if call == 1: 
        in_path = './rmc/data/igzosqro.dat'
        out_path = './rmc/igzosqr.dat'
        frac = []
        elements = []
            
        frac.append(0.33)
        frac.append(0.67)
        elements.append('Si')
        elements.append('O')
           
        coeff(frac, elements, in_path, out_path, 1)


    #=====↓↓calc_xcoeff用=====
    if call == 2:
        in_path = './aaaaa.dat'
        # inputファイルの1列目がq、2列目がs
        out_path = './aaaaa.out'
        frac = []
        symbol = []
            
        frac.append(0.33)
        frac.append(0.67)
        symbol.append('Si')
        symbol.append('O')
        
        #↓qの読み込み, 確認用
        q = []
        try:
            with open(in_path.encode("utf-8"), "r") as f:
                try:
                    nq = int(f.next())
                    line = f.next()
                    for i in range(nq):
                        item = f.next()
                        _q, _s = item.split()[:2]
                        q.append(float(_q))
                except StopIteration:
                    pass
        except IOError:
            raise

        #↓関数の実行、書き出し
        result = calc_xcoeff(frac, symbol, q, 1)
        #print(result)

        fw = open(out_path, 'w')
        for i in range(len(result)):
            local_result = result[i]

            for j in range(len(local_result)):
                fw.write(" %12.7F" % (local_result[j]))

            fw.write("\n")
        fw.close

    else:
        print("Select a correct number.")


