#!/usr/bin/env python
"""
    20/09/04 1st version pierre.legrand \at synchrotron-soleil.fr
    
    Note: This program assum that the most commons XDS frame is choosen
          (see the distributed XDS input files):
          Beam is along the Z axis, Y axis point verticaly down (like gravity),
          and the X axis is defined to yield an orthonormal right handed
          laboratory coordinate system {X,Y,Z}.

    TODO:
        - detector conversion
        - finding the good film rotation for each detector...
	
    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "22-10-2005"
__copyright__ = "Copyright (c) 2005-2009  Pierre Legrand"
__license__ = "New BSD License"
__version__ = "0.4.4"


import sys
import os.path

from XOconv import *


_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Convert Mosflm crystal orientation informations to XDS format.

   The program convert the orientation matix from Denzo (.x file)
   to XDS pseudo XPARM.XDS file:

   USAGE:   %s  [OPTION]... FILE

      FILE is the Denzo integration file (.x format).

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

Qdnz2xds = mat3(-ey, ex, -ez).transpose()

class XP: pass

XPARM_FMT = """%(first_frame)6d%(phi_init)12.6f%(delta_phi)12.6f%(spindle)s
%(wavelength)15.6f       0.000000       0.000000%(inv_wavelength)15.6f
%(nx)10d%(ny)10d%(qy)10.5f%(qy)10.5f
%(distance)15.6f%(beam_x)15.6f%(beam_y)15.6f
%(detector_orientation)s
%(spg)10d%(a)10.3f%(b)10.3f%(c)10.3f%(alpha)10.3f%(beta)10.3f%(gamma)10.3f
%(xds_crystal_orientation)s
"""

detector_orientation = str_mat(mat3(1),format="%15.6f")[:-1]
        
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
            
    x = DenzoParser(sys.argv[1])
    xdspar = XP()

    xdspar.first_frame = 1 # We can't realy guess that...
    xdspar.delta_phi = x.phi1 - x.phi0
    xdspar.phi_init = x.phi0 - (x.sector - xdspar.first_frame) * xdspar.delta_phi
    xdspar.spindle = "    1.000000    0.000000    0.000000"
    xdspar.wavelength = x.wavel
    xdspar.inv_wavelength = 1/x.wavel
    xdspar.distance = x.distance
    xdspar.nx = 3072
    xdspar.ny = 3072
    xdspar.qx = xdspar.qy = xdspar.distance/x.xtod
    xdspar.beam_x = x.beam_x/xdspar.qx
    xdspar.beam_y = x.beam_y/xdspar.qy
    xdspar.detector_orientation = detector_orientation
    xdspar.a, xdspar.b, xdspar.c, xdspar.alpha, xdspar.beta, xdspar.gamma = tuple(x.cell)
    xdspar.spg = SPGlib2[x.spg]
    
    UBxds = Qdnz2xds * x.UB

    Arxds = vec3(UBxds.getColumn(0))
    Brxds = vec3(UBxds.getColumn(1))
    Crxds = vec3(UBxds.getColumn(2))

    Axds = tuple(Brxds.cross(Crxds)*x.volum)
    Bxds = tuple(Crxds.cross(Arxds)*x.volum)
    Cxds = tuple(Arxds.cross(Brxds)*x.volum)

    fmt_mat = 3 * "%15.6f" + "\n"
    xdspar.xds_crystal_orientation = ""
    xdsXO = ""
    for v in (Axds, Bxds, Cxds): xdsXO += fmt_mat % v
    xdspar.xds_crystal_orientation = xdsXO[:-1] 
    print XPARM_FMT % vars(xdspar),
