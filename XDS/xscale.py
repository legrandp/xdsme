#!/usr/bin/env python

__version__ = "0.6.6"
__author__ = "Pierre Legrand (legrand@emble-grenoble.fr)"
__date__ = "21-02-2005"
__copyright__ = "Copyright (c) 2003-2005 Pierre Legrand"
__license__ = "LGPL"

usage   = """

>>>   Usage : xscale.py [-h] [-u] [-a/n] [-N nbins] [reso_high ]
     
      -h         Print this message
      -u         Output unmerged reflections (defaulf is unmerge)
      -m         Output merged reflections (defaulf is unmerge)
      -a/n       -a: Set Friedel's law as "FALSE"
                 -n: Set Friedel's law as "TRUE"
                    (default is to read this information in the reflection
                     file header)
      -n nbins   Number of resolution bins for printing statistics.
                   (default is 12)
                   
>>>   Cell parameters and space group number are taken either from the
      "XDS_ASCII.HKL", "GXPARM.XDS" or "XPARM.XDS" file.
>>>   If not given, reso_high is read from either 
      "DEFPIX.LP" or "INTEGRATE.LP" if present or is set to 1.5A.\n
"""
import sys, os
from xupy import *


#### Default values
float_arg = []
merge = "FALSE"
#friedels_law = "TRUE"
friedels_law = 1
nbin = 12
hklout = "normal.hkl"
dmin = get_resmax_limit()

# Source file for spg, cell...
# Order XDS_ASCII.HKL, GXPARM, XPARM

xparmfile = ""
hklfile = "XDS_ASCII.HKL"
if  os.path.exists(hklfile):
    H = read_xdsascii_head(hklfile)
    if H["friedels_law"] == "FALSE":
        friedels_law = 0
        hklout = "anomal.hkl"
    cell = H["cell"]
    spaceGroupNum = H["sym"]
else:
    if os.path.exists("GXPARM.XDS"): xparmfile = "GXPARM.XDS"
    elif os.path.exists("XPARM.XDS"): xparmfile = "XPARM.XDS"
    if xparmfile:
        line = (open("XPARM.XDS",'r').readlines()[7]).split()
        spaceGroupNum, cell  = line[0],line[1:]

# To handle old XDS.DATA
if os.path.exists("XDS.HKL") and not os.path.exists("XDS_ASCII.HKL"):
    hklfile = "XDS.HKL"

if os.path.exists("XSCALE.INP") and not os.path.exists("XSCALE.INP.bck"):
    os.system("mv XSCALE.INP XSCALE.INP.bck")
if sys.argv.count("-h"):
    print usage
    sys.exit()
if sys.argv.count("-a"):
    sys.argv.remove("-a")
    file_name, friedels_law = "anomal.hkl", 0
if  sys.argv.count("-n"):
    sys.argv.remove("-n")
    file_name, friedels_law = "normal.hkl", 1
if sys.argv.count("-u"):
    sys.argv.remove("-u")
    merge = "FALSE"
if sys.argv.count("-m"):
    sys.argv.remove("-m")
    merge = "TRUE"
if sys.argv.count("-N"):
    arg_ind = sys.argv.index("-N")
    nbin = sys.argv[arg_ind+1]
    sys.argv.remove(nbin)
    sys.argv.remove("-N")
    try: nbin = int(nbin)
    except:
        print "\nERROR! Wrong value given for the number of bins:", nbin
        print "\tNumber of bins will be set to: 12\n"
        nbin = 12


dirname = os.path.split(os.getcwd())[1]
for arg in sys.argv[1:]:
   float_arg.append(float(arg))

if  len(sys.argv) == 2: dmin = float_arg[0]
elif len(sys.argv) >= 3:
        print usage
        sys.exit()

xlatt = Lattice(cell, "Unknown", symmetry=spaceGroupNum,
                      dmin=dmin, friedels_law=friedels_law)
run_xscale((hklfile,),hklout, xlatt, nbin=nbin, merge=merge)



