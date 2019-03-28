#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    Compare orientation matrices

    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).

    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "18-06-2012"
__copyright__ = "Copyright (c) 2010  Pierre Legrand"
__license__ = "New BSD License"
__version__ = "0.0.3"


import sys
import os
import math

from XOconv import *
from AxisAndAngle import axis_and_angle, R2D

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   compare crystal orientation matrices.

   A program to convert the orientation matix, extract information
   from XDS output files and write a mosflm input file:

   USAGE:   %s  [OPTION]... FILE

      FILE can be one of these XDS output files:

         XPARM.XDS, GXPARM.XDS, IDXREF.LP, CORRECT.LP

   OPTIONS:

   -a
   --angles
         Writes out the crystal orientation in xds2mos.umat as
         setings angles (default is as a U matrix)

   -h
   --help
         Print this help message.

""" % _progname

FMT_AXE = 3*"%15.6f"

def kappaVector(alpha=49.64/r2d):
    return vec3([-1*math.cos(alpha), 0, math.sin(alpha)])

kappa = Rotation2(kappaVector(), 30./r2d)

def PARS_xds2mos(xdsPar):
    "Convert XDS output parameters to Mosflm input parameters."
    
    mosPar = {}
    mosPar["title"] = "xds2mos   version: %s" % (__version__)
    mosPar["distance"] = abs(xdsPar["distance"])
    mosPar["wavelength"] = 1/vec3(xdsPar["beam"]).length()
    mosPar["symmetry"] = spg_num2symb[xdsPar["symmetry"]]
    mosPar["omega"] = xdsPar["omega"]*r2d
    mosPar["twotheta"] = xdsPar["twotheta"]*r2d

    xc = xdsPar["origin"][0]
    yc = xdsPar["origin"][1]

    cosOmega = math.cos(xdsPar["omega"])
    sinOmega = math.sin(xdsPar["omega"])

    mosPar["beam_x"] = xc*cosOmega + yc*sinOmega
    mosPar["beam_y"] = xc*sinOmega + yc*cosOmega

    if "detector_type" in xdsPar.keys():
        mosPar["detector"] = detector2scanner[xdsPar["detector_type"]]

    mosPar["pixel_x"] =  xdsPar["pixel_size"][1]
    mosPar["pixel_y"] =  xdsPar["pixel_size"][0]
    mosPar["template"] = xdsPar["template"]
    mosPar["extention"] = xdsPar["template"].split(".")[-1]
    mosPar["image_numb"] = xdsPar["num_init"]
    mosPar["phi_i"] = xdsPar["phi_init"]
    mosPar["phi_f"] = xdsPar["phi_init"] + xdsPar["delta_phi"]

    if  "mosaicity" in xdsPar:
        mosPar["mosaicity"] = mosaicity_conversion_factor*xdsPar["mosaicity"]

    return mosPar


def GonioRotation(omega, kappa, phi):
    "input: gonio_angle in degree and return Goniometer Rotation operator."
    Rphi = Rotation2(ex, phi/r2d)
    Romega = Rotation2(ex, omega/r2d)
    Rkappa = Rotation2(kappaVector(), kappa/r2d)
    return Rphi*Rkappa*Romega

if __name__=='__main__':

    import getopt

    _debug = False
    _write_out_angles = False
    DO_PG_PERMUTATIONS = True
    _start_mosflm = False
    _verbose = False

    short_opt =  "ahpsv"
    long_opt = ["angles", "help", "pg-permutations", "start-mosflm", "verbose"]

    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)

    for o, a in opts:
        if o in ("-v", "--verbose"):
            _verbose = True
        if o in ("-h", "--help"):
            print _usage
            sys.exit()
        if o in ("-a", "--angles"):
            _write_out_angles = True
        if o in ("-s", "--start-mosflm"):
            _start_mosflm = True
        if o in ("-p","--pg-permutations"):
            _do_PG_permutations = True
            print a

    print "\n   XOcompare version: %s\n" % (__version__)
    print "   Extracting orientations from:\t\t%s" % inputf[0:]

    spgn = 0
    XOmat = []
    for inpf in inputf:
        gotcha = False
        for parser in (DenzoParser, XDSParser, MosflmParser):
            try:
                XOparser = parser(inpf)
                gotcha = True
            except :
                pass
            if gotcha:
                break
        if not gotcha:
            raise Exception, "Can't parse inputted orientation matrix file: %s"\
                              % inpf

        print "\n %s used to read input file: %s" % (XOparser.info, inpf)
        XOfileType = XOparser.fileType

        if not spgn:
            spgn = XOparser.spaceGroupNumber
            spg = XOparser.spaceGroupName
        else:
            spg = spg_num2symb[spgn]

        pointGroup = SPGlib[spgn][-1]
        print  "\n Space group symbol: %s,  Number: %d,  Point Group: %s" % \
                                                (spg.upper(),spgn, pointGroup)
        print ("\n       Real cell: "+3*"%10.3f"+3*"%8.2f") % tuple(XOparser.cell)
        print (" Reciprocal cell: "+3*"%10.6f"+3*"%8.2f") % tuple(XOparser.cell_r)

        Bmos = BusingLevy(XOparser.cell_r)

        # Umos =     setting matrix at current datum
        # Getting Mosflm XO
        if XOparser.fileType == "Mosflm":
            UBmos =  XOparser.U * Bmos # whitout wavelength scaling
            Umos = XOparser.U
            # XOparser.UB = UBmos * wavelength

        # Converting Denzo XO to Mosflm convention
        elif XOparser.fileType == "Denzo":
            UBmos = Qdnz2mos * XOparser.UB
            Umos = (UBmos) * Bmos.inverse()

        # Converting XDS XO to Mosflm convention
        elif XOparser.fileType == "XDS":
            UBmos = XOparser.UBxds_to_mos()/ XOparser.dict["wavelength"]
            Umos = (UBmos) * Bmos.inverse()
            A = vec3(XOparser.dict["A"])
            B = vec3(XOparser.dict["B"])
            C = vec3(XOparser.dict["C"])
            
            UBR = mat3(A, B, C).transpose()
            printmat(UBR, '\n   UBR',  "%12.6f")

        is_orthogonal(Umos)
        printmat( Umos,'\n   U',  "%12.6f")
        printmat( Bmos,'\n   B',  "%12.6f")
        printmat( UBmos,'\n  UB', "%12.6f")
        #XOmat.append(UBmos)
        
        XOmat.append(UBR)
        #for a in "A", "B", "C":
        #    print a, FMT_AXE % tuple(XOparser.dict[a])
    
    angles = 0,-55,0
    k_angles = range(-55,58,5)
    p_angles = range(-180, 180, 5)
    for _k in k_angles:
        for _p in k_angles:
            print "%4d %4d  " % (_k, _p),
            print FMT_AXE % tuple(C*GonioRotation(0, _k, _p))
    #for vect in A, B, C:
    #        print FMT_AXE % tuple(vect*GonioRotation(*angles))
            #print FMT_AXE % tuple(XOparser.dict[a])
    print len(k_angles)*len(p_angles)
