#!/usr/bin/env python

__version__ = "0.3.0"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "04-10-2007"
__copyright__ = "Copyright (c) 2007 Pierre Legrand"
__license__ = "LGPL"

import sys
import time
import Numeric

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def reciprocal(cell):
    r2d = 180/Numeric.pi
    cell = Numeric.array(cell)
    sa, sb, sg = Numeric.sin(cell[3:6]/r2d)
    ca, cb, cg = Numeric.cos(cell[3:6]/r2d)
    v = cell[0]*cell[1]*cell[2]*(1-ca**2-cb**2-cg**2+2*ca*cb*cg)**0.5
    rc = Numeric.zeros(6, Numeric.Float)
    rc[0] = cell[1]*cell[2]*sa/v
    rc[1] = cell[2]*cell[0]*sb/v
    rc[2] = cell[0]*cell[1]*sg/v
    rc[3] = Numeric.arccos((cb*cg-ca)/(sb*sg)) * r2d
    rc[4] = Numeric.arccos((ca*cg-cb)/(sa*sg)) * r2d
    rc[5] = Numeric.arccos((ca*cb-cg)/(sa*sb)) * r2d
    return rc

def setmat(cell):
    r2d = 180/Numeric.pi
    cell = Numeric.array(cell,Numeric.Float)
    rcell = reciprocal(cell)
    sr = Numeric.sin(rcell[3:6]/r2d)
    cr = Numeric.cos(rcell[3:6]/r2d)
    a  = Numeric.zeros((3,3),Numeric.Float)
    a[0,0] = rcell[0]*sr[1]
    a[0,1] = rcell[1]*(cr[2]-cr[0]*cr[1])/sr[1]
    a[1,1] = ((rcell[1]*sr[0])**2-a[0,1]**2)**0.5
    a[2,0] = rcell[0]*cr[1]
    a[2,1] = rcell[1]*cr[0]
    a[2,2] = rcell[2]
    return a

def read_resolution(txt):
    idx1 = txt.index("!INCLUDE_RESOLUTION_RANGE=")
    idx2 = txt.index("\n", idx1)
    return map(float,txt[idx1+26:idx2].split())
    
def get_cell(txt):
    idx1 = txt.index("!UNIT_CELL_CONSTANTS=")
    idx2 = txt.index("\n", idx1)
    return map(float,txt[idx1+22:idx2].split())


def getHKLv2(inp_str):
    """ Return a list of hkls. Functional Programming style."""
    headerEnd_markup = "!END_OF_HEADER\n"
    dataEnd_markup = "!END_OF_DATA\n"
    id_endHead = inp_str.find(headerEnd_markup)+len(headerEnd_markup)
    hklOnly = lambda s: map(int, s.split()[:3])
    lines = inp_str[id_endHead:-len(dataEnd_markup)].splitlines()
    return map(hklOnly, lines)

def calc_reso(hkl,cell):
    sum = Numeric.sum
    a = setmat(cell)
    return  1/(sum(a[0,:]*hkl,1)**2 +
               sum(a[1,:]*hkl,1)**2 + 
               sum(a[2,:]*hkl,1)**2)**0.5

def get_info(hklFileName, t0=None):
    
    refl = opReadCl(hklFileName)
    cell = get_cell(refl[:2500])
    #res = read_resolution(refl[:2500])
    
    if t0:
        t1 = time.time()
        print "\n... Reading:          %.3fs" % (t1 - t0)
    hkls = getHKLv2(refl)
    if t0:
        t2 = time.time()
        print "... Extracting HKLs:  %.3fs" % (t2 - t1)
    
    try:
        1/0 
        res = read_resolution(refl[:2500])
        
    except:
	
        import Numeric
        hkls = Numeric.array(hkls)
        if t0:
            t3 = time.time()
            print "... Creating Array:   %.3fs" % (t3 - t2)
        res = calc_reso(hkls, cell)
        if t0:
            t4 = time.time()
            print "... Calculate Reso:   %.3fs" % (t4 - t3)

    res_low, res_high =  max(res), min(res)
    return {"res_low": res_low,
            "res_high": res_high,
            "refl_numb": len(hkls),
            "cell": cell,
            "hklFileName": hklFileName}

info_fmt = """
>>> Reflection File:          %(hklFileName)s
>>> Cell parameters:          %(cell_str)s
>>> Number of reflections:    %(refl_numb)d
>>> Low  resolution limit:    %(res_low).2f
>>> High resolution limit:     %(res_high).2f
"""

if __name__=='__main__':
    t0 = None
    if "-t" in sys.argv: 
       t0 = time.time()
       sys.argv.remove("-t")
    info = get_info(sys.argv[1], t0)
    info["cell_str"] = 6*"%.2f  " % tuple(info["cell"])
    print info_fmt % info
     
