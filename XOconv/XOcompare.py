#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    
    Compare orientation matrices
    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
    
    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "23-11-2009"
__copyright__ = "Copyright (c) 2009  Pierre Legrand"
__license__ = "New BSD License"
__version__ = "0.0.1"


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

if __name__=='__main__':

    import getopt

    _debug = False
    _write_out_angles = False
    _do_PG_permutations = False
    _start_mosflm = False
    _verbose = False
    _template = "xds2mos"
    
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
            raise Exception, "Can't parse inputted orientation matrix file: %s" % inputf[0]

        print "\n %s used to read input file: %s" % (XOparser.info, inputf[0])
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

        is_orthogonal(Umos)
        printmat( Umos,'\n   U',  "%12.6f")
        printmat( Bmos,'\n   B',  "%12.6f")
        printmat( UBmos,'\n  UB', "%12.6f")
        XOmat.append(Umos)
    ############################
    Udiff = XOmat[0] * XOmat[1].inverse()
    printmat(Udiff, '\n   U*U-1',  "%12.6f")
    axis, angle = axis_and_angle(Udiff)
    print "Axis_i:  %9.6f%9.6f%9.6f" % tuple(axis),
    print "Angle_i: %10.5f degree" % (angle*R2D)

    sys.exit()
    XDS1.debut()
    XDS2.debut()
    print XDS1.UBxds_to_mos()
    print XDS2.UBxds_to_mos()
    
    sys.exit()
    matf_name = _template + ".xds2mos.umat"
    mosi_name = _template + ".xds2mos.inp"
    
    if _write_out_angles:
        print "   Writing Crystal setting angles instead of the U matrix."
        
    XDSi.debut()
    MOSi.UB = XDSi.UBxds_to_mos()
    MOSi.cell = XDSi.dict["cell"]

    
    if _do_PG_permutations:
        spgn = XDSi.dict["symmetry"]
        pointGroup = SPGlib[spgn][3]
        PGequivOperators = PGequiv[pointGroup]
        print
        print  ">>> Space Group : %s" % (SPGlib[spgn][1])
        print  ">>> Space Group number: %s" % (spgn)
        print  ">>> Point group: %s" % (pointGroup)
        print  ">>> Number of equivalent crystal ortientations: %d\n" % \
                                         (len(PGequivOperators)+1)
        allPermutedUB = getPermutUB(PGequivOperators, MOSi.UB)
        n = 0
        for _ub in allPermutedUB:
            n+=1
            print "Operator number:", n
            print _ub
        
    B = MOSi.get_B(reciprocal(MOSi.cell))
    MOSi.U = MOSi.UB * B.inverse() / XDSi.dict["wavelength"]
    
    verif = is_orthogonal(MOSi.U)
    if not verif:
        print "???  Warning: The U matrix is not orthogonal."
        print "???  Epsilon error: %.1e" % verif
    
    XDSi.dict["origin"] = XDSi.getBeamOrigin()
    XDSi.dict["omega"] = XDSi.getOmega()
    XDSi.dict["twotheta"] = XDSi.getTwoTheta()
    print "\n   Calculated Omega:    %9.2f degree" % (XDSi.dict["omega"]*r2d)
    print "   Calculated 2theta:   %9.2f degree\n" % (XDSi.dict["twotheta"]*r2d)
    
    if _write_out_angles:
        MOSi.missetingAngles = map_r2d(ThreeAxisRotation2(MOSi.U.toList(1),
                                          inversAxesOrder=1).getAngles()[0])
        MOSi.U = mat3(ex, ey, ez)
    else:    
        MOSi.missetingAngles = 0, 0, 0
        
    MOSi.write_umat(matf_name)
    
    mosDict = PARS_xds2mos(XDSi.dict)
    mosDict["matrixfile"] = matf_name
    if "mosaicity" in mosDict:
        mosDict['mosaicity_instruction'] = mosflmInpTemplate2 % mosDict
    else:
        mosDict['mosaicity_instruction'] = ""
    
    if "detector" in mosDict:
        mosDict['detector_instruction'] = mosflmInpTemplate3 % mosDict
    else:
        mosDict['detector_instruction'] = ""
            
    openWriteClose(mosi_name, mosflmInpTemplate % mosDict)
    
    if _verbose:
        Ud = MOSi.UB.decompose()[0]
        Bd = Ud.inverse() * MOSi.UB / XDSi.dict["wavelength"]
        print "   Decomposition Epsilon Ud-U = %8.1e" % diffMAT(Ud, MOSi.U)
        print "   Decomposition Epsilon Bd-B = %8.1e" % diffMAT(Bd, B)
        print 
        print str_mat(B , '   B decomposition:\n\n', "%14.7f")
    
    print str_mat(MOSi.UB, "   Mosflm UB:\n\n", "%14.7f")
    print "   New Mosflm matix file:  %s" % matf_name
    print "   New Mosflm input file:  %s\n" % mosi_name
    print "   Use -s or --start-mosflm option to start 'ipmosflm < %s'\n" % mosi_name
    
    if _start_mosflm:
        os.system("ipmosflm <  %s" % mosi_name)
