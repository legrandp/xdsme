#!/usr/bin/env python

"""
    08/05/05 First version, build from xds2mos.py
    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "14-11-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "LGPL"
__version__ = "0.2.8"

import sys
import os

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
            Going from XDS to BEST:

  A program to convert the orientation matix, extract information
  from XDS output files and write a BEST input file:

    Usage:   %s CORRECT.LP
            
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

bestTemplate = """TITLE          %(title)s
DETECTOR       %(detector)s
SITE           Not set
DIAMETER     %(diameter)10.2f
PIXEL        %(pixel)11.7f
ROTAXIS        0.00 0.00 1.00 FAST  ## What mean fast?
POLAXIS        0.00 0.00 1.00       ## How do you define this?
GAIN         %(gain)10.3f
XMOSAIC      %(xmosaic)10.3f   ## New definition!
#XCUT             0.020   ## New definition! Not realy needed: doesn't change.
PHIWIDTH     %(phiwidth)10.3f
DISTANCE     %(distance)10.2f
WAVELENGTH   %(wavelength)10.5f
POLARISATION %(polarisation)10.5f
SYMMETRY       %(symmetry)s
UB %(UB)s
CELL         %(cell)s
#RASTER           27  23  20  10   8
#SEPARATION      0.570  0.570
BEAM         %(beam_x)11.3f%(beam_y)11.3f
# end of parameter file for BEST
"""


def PARS_xds2best(xdsPar):
    "Convert XDS output parameters to Best input parameters."
    
    bestPar = {}
    bestPar["title"] = "xds2best   version: %s" % (__version__)
    if xdsPar.has_key("detector_type"):
        bestPar["detector"] = detector2scanner[xdsPar["detector_type"]]
    bestPar["diameter"] = min(xdsPar["pixel_numb"][0]*xdsPar["pixel_size"][0],
                              xdsPar["pixel_numb"][1]*xdsPar["pixel_size"][1])
    bestPar["pixel"] = xdsPar["pixel_size"][0]
    if xdsPar.has_key("mean_gain"):
        bestPar["gain"] = xdsPar["mean_gain"]
    else:
        bestPar["gain"] = 1.
     
    bestPar["xmosaic"] = xdsPar["mosaicity"]
    bestPar["phiwidth"] = xdsPar["delta_phi"]
    bestPar["distance"] = xdsPar["distance"]
    bestPar["wavelength"] = 1/vec3(xdsPar["beam"]).length()
    bestPar["polarisation"] = xdsPar["polarization"]
    bestPar["symmetry"] = spg_num2symb[xdsPar["symmetry"]]
    bestPar["UB"] = str_mat(UB, name="", format="%13.6f")[3:-1]
    bestPar["cell"] = (3*"%9.2f"+3*"%7.2f") % tuple(xdsPar["cell"])
    bestPar["beam_x"] = xdsPar["origin"][1]*xdsPar["pixel_size"][1]
    bestPar["beam_y"] = xdsPar["origin"][0]*xdsPar["pixel_size"][0]
    return bestPar

if __name__=='__main__':

    import getopt


    _start_best = False
    _do_PG_permutations = False
    _verbose = False
    _bestf_name = "xds2best.par"

    short_opt =  "hpsv"
    long_opt = ["help",  "pg-permutations", "start-best", "verbose"]

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
        if o in ("-s", "--start-best"):
            _start_best = True
        if o in ("-p","--pg-permutations"):
            _do_PG_permutations = True
        

    XDSi = XDSParser()
    

    print "\n   xds2best version: %s\n" % (__version__)
    print "   Extracting data from:\t\t%s" % inputf[0]
    
    infilelocation, infilename = os.path.split(inputf[0])

    try:
        XDSi.parse(inputf[0])
    except:
        print "\nERROR! Can't parse input file: %s\n" % inputf[0]
        print _usage
        sys.exit(2)

    try:
            XDSi.parse(os.path.join(infilelocation, "INIT.LP"))
    except:
        pass
        
    XDSi.debut()
    UB = XDSi.UBxds_to_mos()
    
    B = BusingLevy(reciprocal(XDSi.dict["cell"]))
    U = UB * B.inverse() / XDSi.dict["wavelength"]

    verif = is_orthogonal(U)
    if not verif:
        print "???  Warning: The U matrix is not orthogonal."
        print "???  Epsilon error: %.1e" % verif
    
    if "template" in XDSi.dict:
        _template = os.path.split(XDSi.dict['template'])[1].split('#')[0]
        _bestf_name = _template + "xds2best.par"
        
    openWriteClose(_bestf_name, bestTemplate % PARS_xds2best(XDSi.dict))
        
    print str_mat(UB, "\n   Mosflm type UB Matrix:\n\n","%14.7f")
    print "   New Best input parameter file:  %s" % _bestf_name
    print 
    
