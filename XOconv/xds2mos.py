#!/usr/bin/env python
"""
    15/01/03 First version
    xds2mos.py is inspired form the wk2al.f
    program by Richard Kahn (IBS-Grenoble)
    
    Convert: Orientation matrix, Beam center, detector type
    Try to guess: the Omega and 2theta angles.
    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
    
    License: http://www.opensource.org/licenses/bsd-license.php
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "21-03-2007"
__copyright__ = "Copyright (c) 2007  Pierre Legrand"
__license__ = "New BSD License"
__version__ = "0.4.9"


import sys
import os
import math

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Converting XDS crystal orientation informations to Mosflm format.

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

   -s
   --start-mosflm
         Start "ipmosflm < dnz2mos.inp". Than clic on the "Predict" buton
         to verify predictions.
""" % _progname

detector2scanner = {
   "ADSC":              "ADSC",
   "CCDCHESS":          "MARCCD",
   "RAXIS":             "RAXIS",
   "MAR":               "SMALLMAR",
   "MAR345":            "MAR",
   "MAC":               "DIP2020",
   "SMARTCCD":          "SMART",
   "CCDD2AM":           "ESRF",
   "CCDBRANDEIS":       "B4"}
                    

# Depending on the XDS output read, different mosflm input can be generated.
mosflmInpTemplate = """
TILE      Input created from  %(title)s

DEBUG CONTROL

MATRIX        %(matrixfile)s

! The DETECTOR type should be detected automaticaly by Mosflm,
! In case of problem, try to comment/uncomment the next line.   
%(detector_instruction)s

! If the image orientation doesn't match the predictions try to
! comment/uncomment the DETECTOR line and play with the OMEGA values
! (0, 90, 180, 270)
! DETECTOR OMEGA  %(omega)12.2f

! TWOTHETA          %(twotheta)12.2f  # This is experimental...

WAVELENGTH        %(wavelength)12.5f
DISTANCE          %(distance)12.3f
BEAM              %(beam_x)12.3f %(beam_y)12.3f
PIXEL             %(pixel_x)12.5f %(pixel_y)12.5f
SYMMETRY          %(symmetry)s
TEMPLATE          %(template)s
EXTENTION         %(extention)s

%(mosaicity_instruction)s

IMAGE             %(image_numb)d PHI %(phi_i)9.3f TO %(phi_f)9.3f
GO
"""

mosflmInpTemplate2 = "MOSAICITY        %(mosaicity).3f"
mosflmInpTemplate3 = "!DETECTOR         %(detector)s"

mosaicity_conversion_factor = 4.410
# In XDS the mosaicity is discribed as:
# "the standard deviation (e.s.d) of a Gaussian modeling the rocking curve.
# If the limit of integration is set (by default) to 2% of reflection profile
# maximum, from the guassian function this correspond to a mosaic spread of
# 4.41 times the e.s.d.
    
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
            
        
    
    print "\n   xds2mos version: %s\n" % (__version__)
    print "   Extracting data from:\t\t%s" % inputf[0]

    XDSi = XDSParser()
    MOSi = MosflmParser()

    try:
        XDSi.parse(inputf[0])
    except:
        print "\nERROR! Can't parse input file: %s\n" % inputf[0]
        print _usage
        sys.exit(2)
        
    
    if "template" not in XDSi.dict:
        infilelocation, infilename = os.path.split(inputf[0])
        
        # This is necessary to get the detector and template information
        try:
            XDSi.parse(os.path.join(infilelocation, "XDS.INP"))
            print "   Detector type and image name taken from:\tXDS.INP\n"
        except:
            pass
    
    if "template" in XDSi.dict:
        _template = os.path.split(XDSi.dict['template'])[1]
        nDigit = _template.count('#')
        _template = ("%s%0"+str(nDigit)+"d") % (_template.split('#')[0],1)
                
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
