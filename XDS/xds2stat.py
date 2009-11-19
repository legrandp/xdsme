#!/usr/bin/env python
"""
    09/07/03 legrand@embl-grenoble.fr
    
    TODO:
"""

import Numeric
import math, sys, os.path

_progname = os.path.split(sys.argv[0])[1]
_usage = """
  Calculate I and Sigma means vs resolution on XDS reflection files.
  It can read INTEGRATE.HKL, XDS_ASCII.HKL and XSCALE reflection files.

    Usage:   %s   reflection_file [-n10] [-h2.4]
    
  -nXX   specifies the number of resolution bins (default = 20).
  -hX.X  specifies a high resolution cutoff (default is no cutoff).
        
""" % _progname

                    

class XP: pass

FMT = XP()
FMT.dotx_mat = 3*"%15.8f"+3*"%10.6f"+"\n"
FMT.dotx_par = 4*"%12.5f"+4*"%10.4f"+6*"%2d"+"\n"
FMT.refl = ("%4d%4d%4d" +
            " 0%8.1f   100.0   1.00%6.1f %5.3f" +
            "%7.1f%7.1f 1.000    99.9\n")
                   
def spot_profile(boxe_size, background_radius, spot_radius):
    x0, y0 = boxe_size[0]/2.-0.5,boxe_size[1]/2.-0.5
    def cercle_2D(x, y):
        return Numeric.sqrt((x-x0)**2 + (y-y0)**2) <= radius
    radius = background_radius
    background_boxe = (Numeric.fromfunction(cercle_2D, boxe_size) == 0)
    radius = spot_radius
    spot_boxe = 2*Numeric.fromfunction(cercle_2D, boxe_size)
    boxe = spot_boxe+background_boxe
    boxe_str = ""
    line_fmt = " "+ boxe.shape[1]*"%1d"+ "\n"
    for l in boxe:
         boxe_str += line_fmt % tuple(l)
    return  boxe_str
    
def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r
    
def getXpar(_str,match,limit=80,func=float):
    start = _str.index(match)+len(match)
    tmp = _str[start:start+limit].split("\n")[0].split()
    return map(func,tmp)

def read_xdsascii_head(file_name_in):
    friedels_law = "UNKNOWN"
    if not os.path.exists(file_name_in):
        print "ERROR! Can't find file %s.\nSTOP.\n"
        sys.exit()
    raw = open(file_name_in)
    line = raw.readline()
    if line.count("OUTPUT_FILE=INTEGRATE.HKL"):
        ftype = "INTEGRATE"
    elif line.count("FORMAT=XDS_ASCII"):
        ftype = "XDS_ASCII"
    else: ftype = "UNKNOWN"
    while line[0] == "!":
        if line.count("UNIT_CELL_CONSTANTS=") == 1:
            cell = line[line.index("=")+1:-1].strip()
        if line.count("SPACE_GROUP_NUMBER=") == 1:
            sym = int(line[line.index("=")+1:-1].strip())
        if line.count("_LAW=") == 1:
            friedels_law = line[line.index("_LAW=")+5:-1].strip()
        if line.count("NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD="):
            ncolumn = int(line[line.index("RECORD=")+7:-1])
        line = raw.readline()
    return ftype, cell, sym, friedels_law, ncolumn

def read_xds_hkl(hkl_str, ncol):
    hkl_str = hkl_str[hkl_str.index("!END_OF_HEADER")+15:
                      hkl_str.index("!END_OF_DATA")]
    hkl_num = Numeric.array(map(float,hkl_str.split()))
    _len = len(hkl_num)
    hkl_num.shape = (_len/ncol,ncol)
    return hkl_num

def volum(cell):
    cell = Numeric.array(cell,Numeric.Float)
    ca = Numeric.cos(cell[3:6]*math.pi/180)
    va = (1+2*ca[0]*ca[1]*ca[2]-ca[0]**2-ca[1]**2-ca[2]**2)
    return cell[0]*cell[1]*cell[2]*math.sqrt(va)

def reciprocal(cell):
    cell = Numeric.array(cell,Numeric.Float)
    v = volum(cell)
    mod = [0,1,2,0,1] ; rcell = Numeric.zeros(6,Numeric.Float)
    sa = Numeric.sin(cell[3:6]*math.pi/180)
    ca = Numeric.cos(cell[3:6]*math.pi/180)
    for i in range(3):
       j=mod[i+1] ; k = mod[i+2]
       rcell[i]=cell[j]*cell[k]*sa[i]/v
       rcell[i+3]=180/math.pi*math.acos((ca[j]*ca[k]-ca[i])/(sa[j]*sa[k]))
    return rcell

def setmat(cell):
    cell = Numeric.array(cell, Numeric.Float)
    rcell = reciprocal(cell)
    sr = Numeric.sin(rcell[3:6]*math.pi/180)
    cr = Numeric.cos(rcell[3:6]*math.pi/180)
    a  = Numeric.zeros((3,3),Numeric.Float)
    a[0,0] = rcell[0]*sr[1]
    a[0,1] = rcell[1]*(cr[2]-cr[0]*cr[1])/sr[1]
    a[1,1] = math.sqrt((rcell[1]*sr[0])**2-a[0,1]**2)
    a[2,0] = rcell[0]*cr[1]
    a[2,1] = rcell[1]*cr[0]
    a[2,2] = rcell[2]
    return a
   
def calc_reso(hkl,cell):
    #hkl = Numeric.array(hkl)
    a = setmat(cell)
    return  1/Numeric.sqrt(Numeric.sum(a[0,:]*hkl,1)**2 +
        Numeric.sum(a[1,:]*hkl,1)**2 + Numeric.sum(a[2,:]*hkl,1)**2)

def res_bin(dmin,nbin) :
   resoRange = []
   for bin in range(1,nbin + 1):
       resoRange.append(1/(((1/dmin)**3*bin/nbin)**(1./3)))
   return resoRange

if __name__=='__main__':

    if len(sys.argv) == 1:
        print _usage
        sys.exit()
    refl_file_name = sys.argv[1]
    
    nbin = 20
    highest_res = None
    if len(sys.argv) >= 3:
        for arg in sys.argv:
            if arg[:2] == "-n":
                nbin = int(arg[2:])
            if arg[:2] == "-h":
                highest_res = float(arg[2:])
    
    print ">>> Selected Reflection file: ", refl_file_name
    ftype, cell, sym, friedels_law, ncol = read_xdsascii_head(refl_file_name)
    
    if ftype == "UNKNOWN":
        print ">>> ERROR: Unrecognize file type."
        sys.exit()
    
    print ">>>    File type:           ", ftype
    print ">>>    Cell:                ", cell
    print ">>>    Symmetry number:      %d\n" % sym
    
    cell = map(float, cell.split())
    refl_file = opReadCl(refl_file_name)
    hkletc = read_xds_hkl(refl_file, ncol)
    nrefl = hkletc.shape[0]
    print ">>> Total number of reflections:  ", nrefl
    
    S = calc_reso(hkletc[:,:3], cell)
    
    # Add a resolution column in the table
    S.shape = (nrefl,1)
    hkletc = Numeric.concatenate((hkletc[:,:5],S),1)
    
    # Sort the table by order of the last column= resolution
    hkletc = Numeric.take(hkletc, Numeric.argsort(hkletc[:,-1]))
    
    if not highest_res:
         highest_res = hkletc[0,-1]
    elif highest_res < hkletc[0,-1]:
         highest_res = hkletc[0,-1]
         
    reso_bins = res_bin(highest_res,nbin)
    reso_bins.insert(0, 999.99)
    s_bins = 1/(2*Numeric.array(reso_bins))
    print "\n  Resol.\n limit     S          N       <I>        <SIGMA>   <I/SIGMA>\n"
    for n in range(nbin):
         lowr, highr = reso_bins[n], reso_bins[n+1]
         select = (lowr > hkletc[:,-1]) * (hkletc[:,-1] >=  highr)
         nselect = Numeric.sum(select)
         Imean = Numeric.sum(hkletc[:,3]*select)/nselect
         SIGmean = Numeric.sum(hkletc[:,4]*select)/nselect
         IoSIGmean = Numeric.sum(hkletc[:,3]/hkletc[:,4]*select)/nselect
         fmt = "%6.2f %6.3f %9d %11.2f %11.2f %9.2f"
         print fmt % (highr, s_bins[n+1], nselect, Imean, SIGmean, IoSIGmean)
         
    select = (hkletc[:,-1] >=  highest_res)
    nselect = Numeric.sum(select)
    Imean = Numeric.sum(hkletc[:,3]*select)/nselect
    SIGmean = Numeric.sum(hkletc[:,4]*select)/nselect
    IoSIGmean = Numeric.sum(hkletc[:,3]/hkletc[:,4]*select)/nselect
    print "\n total        %9d %11.2f %11.2f %9.2f" % (nselect, Imean, SIGmean, IoSIGmean)
