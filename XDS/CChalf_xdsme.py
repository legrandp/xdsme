#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 19:42:16 2017

@author: lpecqueur
"""

import sys
import os

try:
    import numpy
except ImportError:
    raise ImportError('''%s
             Missing module numpy needed for option --CChalf
             %s'''%('='*47,'='*47))
try :
    import scipy
except :
    raise ImportError('''%s
             Missing module scipy needed for option --CChalf
             %s'''%('='*47,'='*47))
    
###############################################################################
class AimlessLogParser:
    """ A Parser for log file from AIMLESS.
    """
    def __init__(self, filename="", run_dir="",
                 verbose=False, raiseErrors=True):
        self.results = {}
        self.info = "AIMLESS Parser"
        self.fileType = "AIMLESS"
        self.verbose = verbose
        self.resolution = []
        self.CChalf= []
        #
        if not run_dir:
            run_dir = "./"
        self.run_dir = run_dir
        #
        full_filename = os.path.join(self.run_dir, filename)
        #
        if filename:
            try:
                fp = open(full_filename, "r")
                self.lp = fp.read()
                fp.close()
            except:
                raise IOError, "Can't read file: %s" % full_filename
        else:
            self.lp = ""

        if full_filename.count("aimless.log"):
            self.get_aimless_CChalf()
        else:
            if filename:
                raise IOError, "Don't know how to parse file: %s" % \
                               full_filename

    def get_aimless_CChalf(self):
        '''Parse AIMLESS log file for CChalf'''

        ### Select Diffraction range.
        sp1 = self.lp.index("  N  1/d^2    Dmid CCanom    Nanom   RCRanom   CC1/2   NCC1/2   Rsplit     CCfit CCanomfit   $$ $$")
        sp2 = self.lp.index("                   CCanom    Nanom   RCRanom   CC1/2   NCC1/2   Rsplit     CCfit CCanomfit", sp1)
        _table = self.lp[sp1:sp2].splitlines()[1:-2]
        _table = [ map(float, l.split()[1:]) for l in _table ]
        resolution, CChalf=[], []
        for l in _table:
            self.resolution.append(l[0])
            self.CChalf.append(l[5])
        return resolution, CChalf

def func(x,d0,r):
    '''x and y are lists of same size
       x is 1/d**2, y is CC1/2
       d0 is 1/d**2 at half decrease
       r is the steepness of the falloff
    '''
    from scipy import tanh
    return 0.5*(1 -tanh((x - d0)/r))
    
def tanh_fit(Ex, Ey):
    '''Ex and Ey are lists of same size
       x is 1/d**2, y is CC1/2
       d0 is 1/d**2 at half decrease
       r is the steepness of the falloff
    '''
    from scipy.optimize import curve_fit
#    from pylab import *
    
    #initializing parameters
    halffo=(max(Ey)-min(Ey))/2.
    DELTA=(max(Ey)-min(Ey))
    nearestpoint=int()
    for i in Ey:
        delta=abs(i - halffo)
        if delta<DELTA:
            DELTA=delta
            nearestpoint=Ey.index(i)
    d0=Ex[nearestpoint]
    r=0.1
    init_vals=[d0, r]
    p, pcov = curve_fit(func,Ex,Ey, p0=init_vals)
    d0_opt, r_opt=p[0], p[1]
#    plot(Ex, Ey, "wo", Ex, func(Ex,*p), "r-") # Plot of the data and the fit
    return [d0_opt, r_opt]
 
def EstimateCC_NewHighRes(cutoff, d0, r, Ey):
    '''Estimate a New High resolution
    cutoff based on CC1/2=value
    Ey is used to check if experimental data
    contain a CC1/2<Threshold, if not
    function is stopped
    '''
    from math import sqrt
    from scipy import tanh
    from numpy import linspace    

#Check if calculation is sensible or not
    if cutoff<=min(Ey):
        print "WARNING: data don't reach this cutoff value...skipping CC1/2 analysis"
        return
        
    cutoff=float(cutoff)
    d0=float(d0)
    r=float(r)
    HighRes=None
    CC_calc, delta=[], []
    #Searching for CC1/2 resolution by minimizing CC1/2-cutoff        
    x = linspace(0.,1.,1001).tolist()
    for val in x:
        CC_calc.append((0.5*(1 -tanh((val - d0)/r))))

    for val in CC_calc:
        delta.append(abs(val-cutoff))
    HighRes=sqrt(1/x[delta.index(min(delta))])
    CalculatedCutoff=CC_calc[delta.index(min(delta))]
#==============================================================================
#     tolerance= linspace(0.005,0.03,6).tolist()  
#     for val in CC_calc:
#         if True in [cutoff<=val<cutoff+i for i in tolerance]:
#                 HighRes=sqrt(1/x[CC_calc.index(val)])
#                 CalculatedCutoff=val
#         else: continue
#==============================================================================
          
    print '''   %s
   ->  Suggested New High Resolution Limit: %.2f A for CC1/2= %.2f <-
   %s
   ''' %('='*66,HighRes, CalculatedCutoff,'='*66)
#    del y    
    return HighRes

def CalculateAimlessHighRes(filename, run_dir="./", verbose=1, CChalf=0.3):
    aimlessdata=AimlessLogParser(filename, run_dir, verbose) 
    fit=tanh_fit(aimlessdata.resolution,aimlessdata.CChalf)
    Newh=EstimateCC_NewHighRes(CChalf, fit[0], fit[1], aimlessdata.CChalf)
    return Newh
    
if __name__ == "__main__":   
#    inputfile=sys.argv[len(sys.argv)-1]
    res = AimlessLogParser(filename="puck7_10_1_aimless.log", run_dir=".", verbose=1)
    res2=tanh_fit(res.resolution,res.CChalf)
    EstimateCC_NewHighRes(0.30, res2[0], res2[1], res.CChalf)