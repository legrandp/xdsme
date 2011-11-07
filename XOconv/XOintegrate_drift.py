#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Compare orientation matrices

    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).

    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "22-03-2011"
__copyright__ = "Copyright (c) 2011  Pierre Legrand"
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

def get_par(_str,match,limit=70,func=float):
    start = _str.index(match)+len(match)
    tmp = _str[start:start+limit].splitlines()[0].split()
    return map(func,tmp)

def parse_correct(infname="CORRECT.LP", infilelocation="."):
    "Extract information from XDS output CORRECT.LP and INIT.LP"
    # Extract things from CORRECT.LP
    corr = openReadClose(os.path.join(infilelocation, infname))
    ip = corr.index("PARAMETERS USING ALL IMAGES")+100
    corrp = corr[ip:ip+1400]
    corri = corr[:1500]
    corr_dict = {}
    corr_dict["rot"] = get_par(corrp,"ROTATION AXIS")
    corr_dict["beam"] = get_par(corrp,"COORDINATES (REC. ANGSTROEM)")
    corr_dict["distance"] = get_par(corrp,"DETECTOR DISTANCE (mm)")[0]
    corr_dict["origin"] = get_par(corrp,"(PIXELS) OF DIRECT BEAM")
    corr_dict["originXDS"] = get_par(corrp,"ORIGIN (PIXELS) AT ")
    corr_dict["A"] = get_par(corrp,"CELL A-AXIS")
    corr_dict["B"] = get_par(corrp,"CELL B-AXIS")
    corr_dict["C"] = get_par(corrp,"CELL C-AXIS")
    corr_dict["cell"] = get_par(corrp,"UNIT CELL PARAMETERS")
    corr_dict["mosaicity"] = get_par(corrp,"CRYSTAL MOSAICITY (DEGREES)")[0]
    iqx, iqy  = corri.index("QX=")+3, corri.index("QY=")+3
    inx, iny  = corri.index("NX=")+3, corri.index("NY=")+3
    corr_dict["pixel_size"] = float(corri[iqx:iqx+9]),float(corri[iqy:iqy+9])
    corr_dict["pixel_numb"] = int(corri[inx:inx+7]),int(corri[iny:iny+7])
    corr_dict["template"] = get_par(corri, "_DATA_FRAMES=",60,str)[0].replace("?","#")
    corr_dict["symmetry"] = int(get_par(corrp,"SPACE GROUP NUMBER")[0])
    corr_dict["detector_type"] = get_par(corri,"DETECTOR=",40,str)[0]
    corr_dict["detector_X"] = get_par(corri,"DETECTOR_X-AXIS=")
    corr_dict["detector_Y"] = get_par(corri,"DETECTOR_Y-AXIS=")
    corr_dict["phi_init"] = get_par(corri,"STARTING_ANGLE=",15)[0]
    corr_dict["num_init"] = get_par(corri,"STARTING_FRAME=",15,int)[0]
    corr_dict["delta_phi"] = get_par(corri,"OSCILLATION_RANGE=",15)[0]
    corr_dict["divergence_esd"] = get_par(corri,"BEAM_DIVERGENCE_E.S.D.=",15)[0]
    corr_dict["resolution_range"] = get_par(corri,"INCLUDE_RESOLUTION_RANGE=",20)
    corr_dict["friedel"] = get_par(corri,"FRIEDEL'S_LAW=",7,str)[0]
    corr_dict["polarization"] = get_par(corri,"FRACTION_OF_POLARIZATION=",8)[0]
    return corr_dict

def parse_integrate(infname="INTEGRATE.LP", infilelocation="."):
    "Extract information from XDS output CORRECT.LP and INIT.LP"
    # Extract things from CORRECT.LP
    integ = openReadClose(os.path.join(infilelocation, infname))
    integs = integ.split("       PROCESSING OF IMAGES ")[1:]
    integi = integ[:1500]
    all_par_dicts = []
    for integp in integs:
        par_dict = {}
        tag_num = integp[:25].split()
        #image_number = int(tag_num[0]), int(tag_num[2])
        par_dict["image_integ_start"] = int(tag_num[0])
        par_dict["rot"] = get_par(integp,"ROTATION AXIS")
        par_dict["beam"] = get_par(integp,"COORDINATES (REC. ANGSTROEM)")
        par_dict["distance"] = get_par(integp,"DETECTOR DISTANCE (mm)")[0]
        par_dict["origin"] = get_par(integp,"(PIXELS) OF DIRECT BEAM")
        par_dict["originXDS"] = get_par(integp,"ORIGIN (PIXELS) AT ")
        par_dict["A"] = get_par(integp,"CELL A-AXIS")
        par_dict["B"] = get_par(integp,"CELL B-AXIS")
        par_dict["C"] = get_par(integp,"CELL C-AXIS")
        par_dict["cell"] = get_par(integp,"UNIT CELL PARAMETERS")
        par_dict["mosaicity"] = get_par(integp,"CRYSTAL MOSAICITY (DEGREES)")[0]
        iqx, iqy  = integi.index("QX=")+3, integi.index("QY=")+3
        inx, iny  = integi.index("NX=")+3, integi.index("NY=")+3
        par_dict["pixel_size"] = float(integi[iqx:iqx+9]),float(integi[iqy:iqy+9])
        par_dict["pixel_numb"] = int(integi[inx:inx+7]),int(integi[iny:iny+7])
        par_dict["template"] = get_par(integi, "_DATA_FRAMES=",60,str)[0].replace("?","#")
        par_dict["symmetry"] = int(get_par(integp,"SPACE GROUP NUMBER")[0])
        #par_dict["detector_type"] = get_par(integi,"DETECTOR=",40,str)[0]
        #par_dict["detector_X"] = get_par(integi,"DETECTOR_X-AXIS=")
        #par_dict["detector_Y"] = get_par(integi,"DETECTOR_Y-AXIS=")
        par_dict["phi_init"] = get_par(integi,"STARTING_ANGLE=",9)[0]
        par_dict["num_init"] = get_par(integi,"STARTING_FRAME=",9,int)[0]
        par_dict["delta_phi"] = get_par(integi,"OSCILLATION_RANGE=",9)[0]
        #par_dict["divergence_esd"] = get_par(integi,"BEAM_DIVERGENCE_E.S.D.=",9)[0]
        all_par_dicts.append(par_dict)
    return all_par_dicts

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

def openReadClose(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def printmat(mat, name="", format="%12.8f"):
    if name: print "%s" % (name)
    if isinstance(mat, mat3):
        for i in 0,1,2: print 3*format % tuple(mat.getRow(i))
    else:
        for l in mat: print 3*format % tuple(l)
    
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
    
    print "\n   XOintegrate_drift version: %s\n" % (__version__)
    print "   Extracting orientations from:\t\t%s" % inputf[0:]
    
    if "CORRECT.LP" not in inputf:
        print "ERROR: Can't read CORRECT.LP."
        sys.exit()
    else:
        XOparser = XDSParser("CORRECT.LP")
        Bmos = BusingLevy(XOparser.cell_r)

        #corr_par_dict = parse_correct()
        corr_par_dict = XOparser.dict
        #UBmos = XOparser.UBxds_to_mos()/XOparser.dict["wavelength"]
        #Umos = (UBmos) * Bmos.inverse()
        #is_orthogonal(Umos)
        A = vec3(corr_par_dict["A"])
        B = vec3(corr_par_dict["B"])
        C = vec3(corr_par_dict["C"])
        UBR = mat3(A, B, C).transpose()
        printmat(UBR, '\n   Reference UBR',  "%12.6f")
        # apply symmetry operator permutation for easier comparation 
        # of INTEGRATE and CORRECT orientation matrices.
        Up = mat3(vec3(-1, 0, 0), vec3(0, 1, 0), vec3(0, 0, -1))
        UBR = Up * UBR
        printmat(UBR, '\n   Reference UBR after point-group permutation',  "%12.6f")
    if "INTEGRATE.LP" not in inputf:
        print "ERROR: Can't read INTEGRATE.LP."
        sys.exit()
    else:
        integ_par_dicts = parse_integrate()
        integ_par_dict = integ_par_dicts[0]
        A = vec3(integ_par_dict["A"])
        B = vec3(integ_par_dict["B"])
        C = vec3(integ_par_dict["C"])
        UBR = mat3(A, B, C).transpose()
        printmat(UBR, '\n   UBR',  "%12.6f")
        for integ_par_dict in integ_par_dicts[1:]:
            A = vec3(integ_par_dict["A"])
            B = vec3(integ_par_dict["B"])
            C = vec3(integ_par_dict["C"])
            UBRi = mat3(A, B, C).transpose()
            #printmat(UBRi, '\n   UBR',  "%12.6f")
            Udiff = UBR * UBRi.inverse()
            #printmat(Udiff, '\n   U*U-1',  "%12.6f")
            axis, angle = axis_and_angle(Udiff)
            #print "\n>>> DIFFERENCE_1:\n"
            print "%4d " % integ_par_dict["image_integ_start"], 
            print " Axis_i:  %9.5f%9.5f%9.5f" % tuple(axis),
            print "Angle_i: %10.5f degree" % (angle*R2D)
        