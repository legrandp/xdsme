#!/usr/bin/env python2
"""    
    Uses some denzo convertion recipies taken from the rotgen CCP4 program
    by John W. Campbell.
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
    
    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "13-11-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "New BSD License"
__version__ = "0.4.4"

import sys
import os.path

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Converting Mosflm crystal orientation informations to Denzo format.

   A program to convert the orientation matix from Mosflm (.mat file)
   to Denzo format:

   USAGE:   %s  [OPTION]... FILE
    
      FILE is the Mosflm crystal orientation file (.mat).
          
   OPTIONS:
             
    -h
    --help
         Print this help message.

    -p
    --pg-permutations
         Print out the other equivalent crystal orientation
         informations based on the point group allowed permutations.
    
    -S VECTOR
    --spindle VECTOR
         Defines the spindle orientation in denzo frame.
         Default is  0,0,1 for spindle. For example: --spindle 0,1,0
    
    -V VECOTR
    --vertical VECTOR   
         Defines the vertical orientation in denzo frame.
         Default is 1,0,0 for vertical. For example: --vertical 1,0,0
         
    -v
    --verbose
         Turn on verbose output.     
""" % _progname

Qmos2dnz = mat3(ey, ez, ex).transpose() 

FMT_dotx_instructions = """
WAVELENGTH %(wavel).5f 
SPINDLE AXIS%(spindleAxisStr)s VERTICAL AXIS%(verticalAxisStr)s
UNIT CELL%(cellStr)s
CRYSTAL ROTX    %(rotx).3f ROTY  %(roty).3f ROTZ  %(rotz).3f 
"""

if __name__=='__main__':

    import getopt
    
    _debug = False
    _do_PG_permutations = False
    _verbose = False
    _template = "mos2dnz"
    _vertical = (1, 0, 0)
    _spindle =  (0, 0, 1)

    short_opt =  "dhpvS:V:"
    long_opt = ["debug",
                "help",
                "pg-permutations",
                "verbose",
                "spindle=",
                "vertical="]

    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)
        
    for o, a in opts:
        if o in ("-d", "--debug"):
            _debug = True
        if o in ("-v", "--verbose"):
            _verbose = True
        if o in ("-h", "--help"):
            print _usage
            sys.exit()
        if o in ("-p","--pg-permutations"):
            _do_PG_permutations = True
        if o in ("-S","--spindle"):
            try:
               _spindle = tuple(map(float,a.split(',')))
            except:
                print "Error! Can't parse -S or --spindle option."
                print _usage
                sys.exit(2)
        if o in ("-V","--vertical"):
            try:
               _vertical = tuple(map(float,a.split(',')))
            except:
                print "Error! Can't parse -V or --vertical option."
                print _usage
                sys.exit(2)

    MOSi = MosflmParser(inputf[0])

    DNZi = DenzoParser()
    DNZi.verticalAxis = vec3(_vertical)
    DNZi.spindleAxis =  vec3(_spindle)
    DNZi.wavel = MOSi.UB_to_wavelength()
    
    DNZi.cellStr = (6*"%9.3f") % tuple(MOSi.cell)
    DNZi.verticalAxisStr = "%4d%4d%4d" % tuple(DNZi.verticalAxis)
    DNZi.spindleAxisStr =  "%4d%4d%4d" % tuple(DNZi.spindleAxis)
        
    RmisSet = ThreeAxisRotation2(map_d2r(MOSi.missetingAngles),inversAxesOrder=1)
    # Here we combine both orientation informations,
    # from U and the misseting angles
    UBmos_recalc = RmisSet.tensor * MOSi.U * MOSi.B * DNZi.wavel
    assert rootSquareSum(UBmos_recalc - MOSi.UB) < 2e-6
    
    if _debug:
        Udecomp = MOSi.UB.decompose()[0]
        ErrMat =  MOSi.UB.decompose()[0] - MOSi.U
        fmt = 9*"%10.6f"
        print fmt % tuple(MOSi.U.mlist)
        print fmt % tuple(ErrMat.mlist)
        ErrDecomp = rootSquareSum(MOSi.UB.decompose()[0] - MOSi.U)
        print ">>> UB decomposition error: %.1e" % ErrDecomp
        assert ErrDecomp < 2e-5
        
    DNZi.UB = Qmos2dnz * MOSi.UB / DNZi.wavel
    DNZi.U = DNZi.Adnz_to_Udnz(DNZi.UB)
    DNZi.rotx, DNZi.roty, DNZi.rotz = DNZi.UB_to_Rotxyz()
    
    print "\n!--------------------- Adnz (UB) and U header matrices ----------------------\n"
    for i in range(3):
        ml = list(DNZi.UB.getRow(i)) + list(DNZi.U.getRow(i))
        print (3*"%15.8f"+ 3*"%10.6f") % tuple(ml)
    
    print "\n!----------------------------- Denzo instructions ---------------------------"
    print FMT_dotx_instructions % vars(DNZi)
