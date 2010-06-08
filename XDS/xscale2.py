#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.0"
__author__ = "Pierre Legrand (legrand@emble-grenoble.fr)"
__date__ = "22-10-2009"
__copyright__ = "Copyright (c) 2003-2009 Pierre Legrand"
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
import sys
import os
from xupy import exec_prog, resum_scaling
from XDSReflectionFile import XDSReflectionFile

template_str = """!
!      A single output file is generated from scaling 2 input
!      files making use of most of XSCALE's input parameters.
!      To activate an input parameter remove all "!" left of it.

!MAXIMUM_NUMBER_OF_PROCESSORS=16
!RESOLUTION_SHELLS= 10 6 4 3 2.5 2.0 1.8 1.7 1.6
!SPACE_GROUP_NUMBER= %(spg)s
!UNIT_CELL_CONSTANTS= %(cell)s
!REIDX=-1 0 0 0    0 -1 0 0    0 0 -1 0
!REFERENCE_DATA_SET= fae-rm.ahkl

!MINIMUM_I/SIGMA=3.0
!REFLECTIONS/CORRECTION_FACTOR=50   !minimum #reflections/correction_factor
!0-DOSE_SIGNIFICANCE_LEVEL=0.10

OUTPUT_FILE= %(hklout)s
STRICT_ABSORPTION_CORRECTION= TRUE 
"""

_fmt_final_stat = """
        Refined Parameters and Scaling Statistics
        =========================================\n

        Space group   number    %(spg_num)d
        symbol    %(spg_sym)s

        Cell parameters     %(cell)s

        Resolution           %(LowestReso)8.2f -%(reso)6.2f\
        (%(resoL).2f - %(reso).2f)

        Completeness                    %(compl)5.1f%%   (%(complL).1f%%)
        I/sigma(I)                     %(isig)6.1f    (%(isigL).1f)
        Rmeas                          %(rmeas)6.1f%%   (%(rmeasL).1f%%)
        Rsym                          %(rsym)7.2f%%   (%(rsymL).1f%%)
        Multiplicity               %(multiplicity)10.1f
        Compared                   %(compar)10d    (%(comparL)d)
        Measured                   %(total)10d
        Unique                     %(unique)10d
        Rejected misfits           %(misfit)10d
        Wilson scaling (B/Corr)    %(wilson_b)10.1f    (%(wilson_corr).2f)
"""

xdshome = os.getenv("XDSHOME")
if xdshome:
    xscale_exec = os.path.join(xdshome,"xscale_par")
else:
    xscale_exec = "xscale_par"

xds_input_files = []
print "File selected for scaling:\n"

for arg in sys.argv[1:]:
    try:
        if os.path.isfile(arg):
            print "for %s: " % arg
            f = open(arg)
            flines = f.readlines()
            f.close()
            if ("FORMAT=XDS_ASCII" in flines[0]) and \
               ("!END_OF_DATA" in flines[-1]):
                xds_input_files.append(arg)
                print arg
    except:
        pass

if xds_input_files == []:
    print "None!"
    print "\nError no valid reflection files given!"
    sys.exit()

# the first reflection file is taken as the reference.
hklref = XDSReflectionFile(xds_input_files[0])
refdic = {}
refdic['cell'] = hklref.header["UNIT_CELL_CONSTANTS"]
refdic['spg'] = hklref.header["SPACE_GROUP_NUMBER"]
refdic['hklout'] = "XSCALE.HKL"

f = open("XSCALE.INP","w")
f.write(template_str % refdic)
f.write("    MERGE=FALSE\n\n")

for file in xds_input_files:
    hklf = XDSReflectionFile(file)
    f.write("   INPUT_FILE= %s\n" % file)
    f.write("      INCLUDE_RESOLUTION_RANGE= %s\n" % hklf.header["INCLUDE_RESOLUTION_RANGE"])
    f.write("      FRIEDEL'S_LAW=            %s\n" % hklf.header["FRIEDEL'S_LAW"])
    f.write("      ! CORRECTIONS= DECAY MODULATION ABSORPTION\n")

f.close()


exec_prog("%s" % (xscale_exec), stdout=None, stderr="xupy.stderr")
s = resum_scaling(lpf="XSCALE.LP")
if not s:
    print "\nERROR while running XSCALE"
    sys.exit()
    #else:
    #    print s
print _fmt_final_stat % vars(s)
if s.absent:
   print _fmt_AbsIav % vars(s)
if ["FRIEDEL'S_LAW"] == "FALSE":
   print _fmt_anomal % vars(s)
            
            
sys.exit()


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

