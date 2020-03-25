#!/usr/bin/env python2
# -*- coding: utf-8 -*-

__version__ = "0.6.3"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "25-09-2018"
__copyright__ = "Copyright (c) 2003-2018 Pierre Legrand"
__license__ = "LGPL"

USAGE   = """

>>>   Usage : xscale2.py [-h] [-u] [-a/n] [-N nbins] [reso_high ]

      -h         Print this message
      -m         Output merged reflections (defaulf is unmerge)
      -a         Set Friedel's law to "FALSE"
      -n         Set Friedel's law to "TRUE"
                 Default is to read this information in the reflection
                     file header
      -r   FLOAT Set the high resolution limit to FLOAT value.
      -S0        Scaling schem: 0 = No correction.
      -z         Apply Zerro Dose extrapolation.
                 

>>>   Cell parameters and space group number are taken either from the
      "XDS_ASCII.HKL", "GXPARM.XDS" or "XPARM.XDS" file.
>>>   If not given, reso_high is read from either 
      "DEFPIX.LP" or "INTEGRATE.LP" if present or is set to 1.5A.\n
"""
import sys
import os
from xupy import exec_prog, resum_scaling
from XDSReflectionFile import XDSReflectionFile

INPUT_HEAD = """!
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

"""

OUTPUT_STR = """\nOUTPUT_FILE= %(hklout)s   ! For wavelength= %(wavelength).5f
STRICT_ABSORPTION_CORRECTION= TRUE
      MERGE=%(merge)s

"""

FMT_FINAL_STAT = """
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

RESOLUTION = None
SCALE_TYPE = None
ZERO_DOSE = None
MERGE = False
FRIEDELS_LAW = True
FRIEDELS_LAW_TRUE = False

if sys.argv.count("-h"):
    print USAGE
    sys.exit()
if sys.argv.count("-S0"):
    sys.argv.remove("-S0")
    SCALE_TYPE = 0
if sys.argv.count("-z"):
    sys.argv.remove("-z")
    ZERO_DOSE = "XTAL1"
if sys.argv.count("-m"):
    sys.argv.remove("-m")
    MERGE = True
if sys.argv.count("-a"):
    sys.argv.remove("-a")
    FRIEDELS_LAW = False
    FRIEDELS_LAW_TRUE = True
if sys.argv.count("-n"):
    sys.argv.remove("-n")
    FRIEDELS_LAW = True
    FRIEDELS_LAW_TRUE = True
if sys.argv.count("-r"):
    p = sys.argv.index("-r")
    sys.argv.remove("-r")
    RESOLUTION = float(sys.argv[p])
    sys.argv.remove(sys.argv[p])

for arg in sys.argv[1:]:
    try:
        if os.path.isfile(arg):
            f = open(arg)
            flines = f.readlines()
            f.close()
            if ("FORMAT=XDS_ASCII" in flines[0]) and \
               ("!END_OF_DATA" in flines[-1]):
                xds_input_files.append(arg)
    except:
        pass

if xds_input_files == []:
    print USAGE 
    sys.exit()

# the first reflection file is taken as the reference.
hklref = XDSReflectionFile(xds_input_files[0])
refdic = {}
refdic['cell'] = hklref.header["UNIT_CELL_CONSTANTS"]
refdic['spg'] = hklref.header["SPACE_GROUP_NUMBER"]
refdic['hklout'] = "XSCALE.HKL"
refdic['merge'] = MERGE

try:
    # work for XDS_ASCII from XDS
    refdic['wavelength'] = float(hklref.header["X-RAY_WAVELENGTH"])
except:
    refdic['wavelength'] = float(hklref.header["ISET_1"])
    
f = open("XSCALE.INP","w")
f.write(INPUT_HEAD % refdic)
if SCALE_TYPE == 0:
    f.write("SNRC=50\n")

hklf_files = []
wavelengths = []
muliwavelength = False

for _file in xds_input_files:
    hklf = XDSReflectionFile(_file)
    hklf_files.append(hklf)
    try:
        wavelength = float(hklf.header["X-RAY_WAVELENGTH"])
    except:
        wavelength = float(hklf.header["ISET_1"])
    wavelengths.append(wavelength)
    print "INPUT: %s    WAVELENGTH: %.4f" % (_file, wavelength)

w0 = wavelengths[0]
if len(wavelengths) > 1:
    for w in wavelengths[1:]:
        if abs(w - w0) > 0.0001:
            muliwavelength = True
            print 1
            break

if not muliwavelength:
    f.write((OUTPUT_STR % refdic).upper())

mad_fnames = []
print "MultiWavelength:", muliwavelength
for hklf in hklf_files:
    if muliwavelength:
        refdic['wavelength'] = float(hklf.header["X-RAY_WAVELENGTH"])
        fname = "XSCALE_%.4fA.HKL" % refdic['wavelength']
        refdic['hklout'] = (fname)
        print fname, mad_fnames
        if fname not in mad_fnames:
            f.write((OUTPUT_STR % refdic).upper())
            mad_fnames.append(fname)
    f.write("   INPUT_FILE= %s\n" % hklf.fileName)
    if ZERO_DOSE:
        f.write("CRYSTAL_NAME=%s\n" % "XTAL1")
    if RESOLUTION:
        res_range= "100 %.3f" % (RESOLUTION)
        f.write("      INCLUDE_RESOLUTION_RANGE= %s\n" % res_range)
    elif "INCLUDE_RESOLUTION_RANGE" in hklf.header:
        f.write("      INCLUDE_RESOLUTION_RANGE= %s\n" % \
                hklf.header["INCLUDE_RESOLUTION_RANGE"])
    if not FRIEDELS_LAW_TRUE:
        f.write("      FRIEDEL'S_LAW=        %s\n" % \
                                 hklf.header["FRIEDEL'S_LAW"])
    elif FRIEDELS_LAW == True:
        f.write("      FRIEDEL'S_LAW=        TRUE\n")
    elif FRIEDELS_LAW == False :
        f.write("      FRIEDEL'S_LAW=        FALSE\n")
    if SCALE_TYPE == 0:
        f.write("      CORRECTIONS= NONE\n")      
    else:
        f.write("      ! CORRECTIONS= DECAY MODULATION ABSORPTION\n")
f.close()


exec_prog("%s" % (xscale_exec), stdout=None, stderr="xupy.stderr")
s = resum_scaling(lpf="XSCALE.LP")
if not s:
    print "\nERROR while running XSCALE"
    sys.exit()
    #else:
    #    print s
print FMT_FINAL_STAT % vars(s)
#if ["FRIEDEL'S_LAW"] == "FALSE":
#   print _fmt_anomal % vars(s)
