#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.2"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "07-12-2015"
__copyright__ = "Copyright (c) 2010 Pierre Legrand"
__license__ = "LGPL"

import sys
import time
import numpy as N

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def reciprocal(cell):
    r2d = 180/N.pi
    cell = N.array(cell)
    sa, sb, sg = N.sin(cell[3:6]/r2d)
    ca, cb, cg = N.cos(cell[3:6]/r2d)
    v = cell[0]*cell[1]*cell[2]*(1-ca**2-cb**2-cg**2+2*ca*cb*cg)**0.5
    rc = N.zeros(6, N.float)
    rc[0] = cell[1]*cell[2]*sa/v
    rc[1] = cell[2]*cell[0]*sb/v
    rc[2] = cell[0]*cell[1]*sg/v
    rc[3] = N.arccos((cb*cg-ca)/(sb*sg)) * r2d
    rc[4] = N.arccos((ca*cg-cb)/(sa*sg)) * r2d
    rc[5] = N.arccos((ca*cb-cg)/(sa*sb)) * r2d
    return rc

def setmat(cell):
    r2d = 180/N.pi
    cell = N.array(cell,N.float)
    rcell = reciprocal(cell)
    sr = N.sin(rcell[3:6]/r2d)
    cr = N.cos(rcell[3:6]/r2d)
    a  = N.zeros((3,3),N.float)
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


def getHKLv2(inp_str, maxcolumn=3, dtype=int):
    """ Return a list of hkls. Functional Programming style."""
    headerEnd_markup = "!END_OF_HEADER\n"
    dataEnd_markup = "!END_OF_DATA\n"
    id_endHead = inp_str.find(headerEnd_markup)+len(headerEnd_markup)
    hklOnly = lambda s: map(dtype, s.split()[:maxcolumn])
    lines = inp_str[id_endHead:-len(dataEnd_markup)].splitlines()
    return map(hklOnly, lines)

def calc_reso(hkl, cell):
    sum = N.sum
    a = setmat(cell)
    return  1/(sum(a[0,:]*hkl[:,:3],1)**2 +
               sum(a[1,:]*hkl[:,:3],1)**2 +
               sum(a[2,:]*hkl[:,:3],1)**2)**0.5

def low_res_info(hklFileName, t0=None, filterR=False):
    "Read the header and return some details."
    refl = opReadCl(hklFileName)
    cell = get_cell(refl[:2500])
    if t0:
        t1 = time.time()
        print "\n... Reading:          %.3fs" % (t1 - t0)
    hkls = getHKLv2(refl, maxcolumn=8, dtype=float)
    if t0:
        t2 = time.time()
        print "... Extracting HKLs:  %.3fs" % (t2 - t1)
    #try:
    #    1/0
    #    res = read_resolution(refl[:2500])
    #except:
    if True:
        hkls = N.array(hkls)
        if t0:
            t3 = time.time()
            print "... Creating Array:   %.3fs" % (t3 - t2)
        res = calc_reso(hkls, cell)
        isort = res.argsort()
        res = res[isort]
        hkls = hkls[isort]
        if t0:
            t4 = time.time()
            print "... Calculate Reso:   %.3fs" % (t4 - t3)

    res_low, res_high =  max(res), min(res)
    #     "          45.18   -4   0   0   6.260e+02   2.613e+01      23.96
    print "          RESOL    H   K   L   INTENSITY    SIGMA    INTENSITY/SIGMA\n"
    for i in -1*N.array(range(1,N_LIST_LOW)):
        Iover6 = hkls[i,3]/abs(hkls[i,4])
        R = W = " "
        if hkls[i,4]<0:
            R = "R"
        if Iover6 < 3.:
            W = "W"
        print " %1s %1s %10.2f" % (R, W, res[i]), \
              "%5d%5d%5d%12.3e%12.3e" % tuple(hkls[i,:5]), \
              "%9.2f" % Iover6
    print "          RESOL    H   K   L   INTENSITY    SIGMA    INTENSITY/SIGMA\n"
    for i in -1*N.array(range(1,N_LIST_LOW)):
        Iover6 = hkls[i,3]/abs(hkls[i,4])
        R = W = " "
        if hkls[i,4]<0:
            R = "R"
        if Iover6 < 3.:
            W = "W"
        print " %1s %1s %10.2f %6.1f" % (R, W, res[i], Iover6), \
              "%5d%5d%5d" % tuple(hkls[i,:3]), \
              "%8.1f%8.1f%8.1f 3 3 3" % tuple(hkls[i,-3:])
              
             
    return {"res_low": res_low,
            "res_high": res_high,
            "refl_numb": len(hkls),
            "cell": cell,
            "hklFileName": hklFileName}


def get_info(hklFileName, t0=None):
    "Read the header and return some details."
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
    #try:
    #    1/0 
    #    res = read_resolution(refl[:2500])
    if 1:
        hkls = N.array(hkls)
        if t0:
            t3 = time.time()
            print "... Creating Array:   %.3fs" % (t3 - t2)
        res = calc_reso(hkls, cell)
        if t0:
            t4 = time.time()
            print "... Calculate Reso:   %.3fs" % (t4 - t3)

    res_low, res_high =  max(res), min(res)
    #print res[:10], res.shape
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
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t", "--test",
                  action="store_true", dest="test", default=False,
                  help="Performe and report timing tests.")
    parser.add_option("-f", "--file", metavar="FILE",
                  dest="file", default="XDS_ASCII.HKL" ,
                  help="Input reflection file.")
    parser.add_option("-r", "--lowres",
                  action="store_true", dest="lowres", default=False,
                  help="Report on quality of lower resolution reflections.")
    parser.add_option("-F","--filter",
                  action="store_true", dest="filterR", default=False,
                  help="Try to better filter lower resolution reflections.")
    parser.add_option("-n", "--num_list_refl",
                  dest="nlistlow", default=40, type="int", metavar="N_LIST_LOW",
                  help="Number of listed lower reflections.")

    (options, args) = parser.parse_args()
    N_LIST_LOW = options.nlistlow
    
    t0 = None
    if options.test: 
       t0 = time.time()
    if options.filterR:
        info = low_res_info(options.file, t0, filterR=True)
        info["cell_str"] = 6*"%.2f  " % tuple(info["cell"])
        print info_fmt % info
    elif options.lowres:
        info = low_res_info(options.file, t0)
        info["cell_str"] = 6*"%.2f  " % tuple(info["cell"])
        print info_fmt % info
    else:
        info = get_info(options.file, t0)
        info["cell_str"] = 6*"%.2f  " % tuple(info["cell"])
        print info_fmt % info
