#!/usr/bin/env python
"""
    24/09/04 1st version pierre.legrand@synchrotron-soleil.fr
    
    Note: This program assum that the most commons XDS frame is choosen
          (see the distributed XDS input files):
          Beam is along the Z axis, Y axis point verticaly down (like gravity),
          and the X axis is defined to yield an orthonormal right handed
          laboratory coordinate system {X,Y,Z}.

    TODO:
        - detector conversion
    
    Uses the ScientificPython module by Konrad Hinsen
    http://starship.python.net/crew/hinsen/scientific.html
    
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "16-11-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "GPL"
__version__ = "0.2.2"

import math
import sys
import os.path

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Converting Mosflm crystal orientation informations to XDS format.

   The program convert the orientation matix from Mosflm (.mat file)
   to XDS pseudo XPARM.XDS file:

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

    -v
    --verbose
         Turn on verbose output.
""" % _progname

Qmos2xds = mat3(ez,-ey,ex).transpose() 

XPARM_fmt = """%(first_frame)6d%(phi_init)12.6f%(delta_phi)12.6f%(spindle)s
%(wavelength)15.6f       0.000000       0.000000%(inv_wavelength)15.6f
%(nx)10d%(ny)10d%(qy)10.5f%(qy)10.5f
%(distance)15.6f%(beam_x)15.6f%(beam_y)15.6f
%(detector_orientation)s
%(spg)10d%(a)10.3f%(b)10.3f%(c)10.3f%(alpha)10.3f%(beta)10.3f%(gamma)10.3f
%(xds_crystal_orientation)s
"""
XPARM_fmt ="""     1    PHI_INIT    DELT_PHI    1.000000    0.000000    0.000000
%(wavelength)15.6f       0.000000       0.000000%(inv_wavelength)15.6f
    NX_PIX    NY_PIX  X_PIXSIZ  Y_PIXSIZ
     DISTANCE      X_BEAM_CENT    Y_BEAM_CENT
       1.000000       0.000000       0.000000
       0.000000       1.000000       0.000000
       0.000000       0.000000       1.000000
         1%(cellStr)s
%(xds_crystal_orientation)s"""

        
def getCellParameters(UB):
    """Return an array containing the cell parameters with angles en degree"""
    Ar = Vector(UB[:,0])
    Br = Vector(UB[:,1])
    Cr = Vector(UB[:,2])
    return Numeric.array([Ar.length(), Br.length(), Cr.length(),
        Br.angle(Cr)*r2d, Cr.angle(Ar)*r2d, Ar.angle(Br)*r2d])

if __name__ == '__main__':

    import getopt
    
    _debug = False
    _do_PG_permutations = False
    _verbose = False

    short_opt =  "dhpv"
    long_opt = ["debug",
                "help",
                "pg-permutations",
                "verbose"]

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

    MOSi = MosflmParser(inputf[0])
    
    xdsPar = {}
    xdsPar["cellStr"] = 6*"%10.3f" % tuple(MOSi.cell)
    
    wavelength = MOSi.UB_to_wavelength()
    xdsPar["wavelength"] = wavelength
    xdsPar["inv_wavelength"] = 1/wavelength
    
    UBxds = Qmos2xds * MOSi.UB / wavelength
    Ar = vec3(UBxds.getColumn(0))
    Br = vec3(UBxds.getColumn(1))
    Cr = vec3(UBxds.getColumn(2))
    
    volumInv = Ar * Br.cross(Cr)
    A = Br.cross(Cr)/volumInv
    B = Cr.cross(Ar)/volumInv
    C = Ar.cross(Br)/volumInv
    
    #print A.length(), B.length(), C.length(), B.angle(C)*r2d, A.angle(C)*r2d, A.angle(B)*r2d
    
    fmtM = 3 * "%15.6f" + "\n"
    xdsPar["xds_crystal_orientation"] = fmtM % tuple(A) + \
                                        fmtM % tuple(B) + \
                                        fmtM % tuple(C)
    print XPARM_fmt % xdsPar,
