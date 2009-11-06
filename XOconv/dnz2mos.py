#!/usr/bin/env python
"""    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "21-10-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "GPL"
__version__ = "0.4.8"

import sys
import os.path

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Converting Denzo crystal orientation informations to Mosflm format.

   A program to convert the orientation matix from Denzo (.x file
   output files) to a Mosflm matrix file and a simple Mosflm input file:

   USAGE:   %s  [OPTION]... FILE
    
      FILE is the Denzo dot.x reflection file.
          
   OPTIONS:
             
    -h
    --help
         Print this help message.

    -p
    --pg-permutations
         Print out the other equivalent crystal orientation
         informations based on the point group allowed permutations.

    -s
    --start-mosflm
         Start "ipmosflm < dnz2mos.inp". Than clic on the "Predict" buton
         to verify predictions.
             
    -v
    --verbose
         Turn on verbose output.     
""" % _progname

DNZAxes = ey, -1*ex, -1*ez
Qdnz2mos = mat3(ez, ex, ey).transpose()

mosflmInpTemplate = """
TILE      Input created from  %(title)s

MATRIX            %(matrix_file)s
WAVELENGTH        %(wavelength)12.5f
DISTANCE          %(distance)12.3f
BEAM              %(beam_x)12.3f %(beam_y)12.3f
SYMMETRY          %(symmetry)s
MOSAICITY         %(mosaicity).3f

IMAGE             %(image_name)s
GO
"""

def PARS_dnz2mos(dnzPar):
    "Convert Denzo output parameters to Mosflm input parameters."
    
    mosPar = {}
    mosPar["title"] = "dnz2mos   version: %s" % (__version__)
    mosPar["distance"] = dnzPar.distance
    mosPar["wavelength"] = dnzPar.wavel
    mosPar["symmetry"] = dnzPar.spg.upper()
    
    #mosPar["omega"] = xdsPar["omega"]*r2d
    #mosPar["twotheta"] = xdsPar["twotheta"]*r2d
        
    mosPar["beam_x"] = dnzPar.beam_x
    mosPar["beam_y"] = dnzPar.beam_y
    mosPar["template"] = dnzPar.template
    mosPar["mosaicity"] = dnzPar.mosaic
        
    return mosPar

if __name__ == '__main__':

    import getopt
    
    short_opt =  "aho:pvs"
    long_opt = ["angles",
                "help",
                "output=",
                "pg-permutations",
                "start-mosflm",
                "verbose"]
    
    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)
    
    _angles = False
    _do_PG_permutations = False
    _start_mosflm = False
    _verbose = False
    
    if inputf:
        base = ".".join(inputf[0].split(".")[:-1])
        matName = os.path.basename(base + ".dnz2mos.umat")
        inpName = os.path.basename(base + ".dnz2mos.inp")

    for o, a in opts:
        if o == "-v":
            _verbose = True
        if o in ("-h", "--help"):
            print _usage
            sys.exit()
        if o in ("-a", "--angles"):
            _angles = True
            "  Writing Crystal setting angles in place of the U matrix."
        if o in ("-o", "--output"):
            matName = a
        if o in ("-p","--pg-permutations"):
            _do_PG_permutations = True
        if o in ("-s", "--start-mosflm"):
            _start_mosflm = True
        
    x = DenzoParser(inputf[0])
    B = BusingLevy(x.cell_r)
    
    MOSi = MosflmParser()
    MOSi.cell = x.cell
    MOSi.UB = Qdnz2mos * x.UB * x.wavel
    MOSi.U = MOSi.UB * B.inverse() / x.wavel
    
    print
    print  ">>> Space Group : %s" % (x.spg.upper())
    print  ">>> Unit Cell :   %s" % (6*"%.2f  " % (x.cell))
    
    printmat(MOSi.UB/x.wavel, '\n>>> UBmos/x.wavel')
    printmat(MOSi.U, '\n>>> Umos')
    printmat(B, '\n>>> Bmos')
    
    # There is two way to write the crystal orientation in mosflm:
    # using the U matrix or the phi1, phi2, phi3 misseting angles (in degree).
    if not _angles:
        MOSi.missetingAngles = 0, 0, 0
    else:
        solve = ThreeAxisRotation2(MOSi.U.mlist, inversAxesOrder=1).getAngles()
        MOSi.missetingAngles = map_r2d(solve[0])
        print "\n  Setting angles:  PHIX %8.2f  PHIY %8.2f PHIZ %8.2f\n" % \
                  tuple(MOSi.missetingAngles)
        Umos = mat3(1)

    MOSi.write_umat(matName)
    
    MOSpar = PARS_dnz2mos(x)
    MOSpar["matrix_file"] = matName
    MOSpar["image_name"] =  MOSpar["template"].replace("###","001")
    openWriteClose(inpName, mosflmInpTemplate %  MOSpar)

    if _start_mosflm:
        os.system("mosflm <  %s" % inpName)
