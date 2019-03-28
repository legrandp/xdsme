#!/usr/bin/env python2
# -*- coding: utf-8 -*-

__version__ = "0.3.1"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "10-11-2009"
__copyright__ = "Copyright (c) 2007-2009 Pierre Legrand"
__license__ = "New BSD License http://www.opensource.org/licenses/bsd-license.php"

import sys
import os
import time
import re
import urllib

_progname = os.path.split(sys.argv[0])[1]

_usage = """
   Correction for DAC absorption on XDS INTEGRATE.HKL.

   USAGE:   %s  [OPTION]... INTEGRATE.HKL

      FILE is for one or multiple diffraction image files.

   OPTIONS:

    -h,  --help
         Print this help message.

    -l,  --length N
         Indicate the diamond path length in cm.
         For example: -l 0.15

    -e,  --energy
         OPTIONAL. Indicate the energy in keV.
         By default, the value from the header of the reflection file is read.
""" % _progname

try:
    import Numeric
except ImportError:
    try:
        import numpy as Numeric
    except ImportError:
        print "Can't import Numeric or numpy module !"
        sys.exit(2)

sin = Numeric.sin
cos = Numeric.cos
r2d = 180/Numeric.pi
d2r = Numeric.pi/180
exp = Numeric.exp
element = "C"

def get_mu_from_mucalweb(wave, element):
    energy = 12.398/wave
    _address = "http://csrri.iit.edu/cgi-bin/mucal-form?name=%s&ener=%.4f"
    f = urllib.urlopen(_address % (element,energy))
    txt = f.read()
    _i1 = txt.find("Total")+5
    mu = float(txt[_i1:txt.find("cm^2/gm",_i1)])
    return mu

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def reciprocal(cell):
    cell = Numeric.array(cell)
    sa, sb, sg = Numeric.sin(cell[3:6]/r2d)
    ca, cb, cg = Numeric.cos(cell[3:6]/r2d)
    v = cell[0]*cell[1]*cell[2]*(1-ca**2-cb**2-cg**2+2*ca*cb*cg)**0.5
    rc = Numeric.zeros(6, Numeric.Float)
    rc[0] = cell[1]*cell[2]*sa/v
    rc[1] = cell[2]*cell[0]*sb/v
    rc[2] = cell[0]*cell[1]*sg/v
    rc[3] = Numeric.arccos((cb*cg-ca)/(sb*sg)) * r2d
    rc[4] = Numeric.arccos((ca*cg-cb)/(sa*sg)) * r2d
    rc[5] = Numeric.arccos((ca*cb-cg)/(sa*sb)) * r2d
    return rc

def setmat(cell):
    cell = Numeric.array(cell,Numeric.Float)
    rcell = reciprocal(cell)
    sr = Numeric.sin(rcell[3:6]/r2d)
    cr = Numeric.cos(rcell[3:6]/r2d)
    a  = Numeric.zeros((3,3),Numeric.Float)
    a[0,0] = rcell[0]*sr[1]
    a[0,1] = rcell[1]*(cr[2]-cr[0]*cr[1])/sr[1]
    a[1,1] = ((rcell[1]*sr[0])**2-a[0,1]**2)**0.5
    a[2,0] = rcell[0]*cr[1]
    a[2,1] = rcell[1]*cr[0]
    a[2,2] = rcell[2]
    return a

def read_resolution(txt):
    idx1 = txt.index("!INCLUDE_RESOLUTION_RANGE=")
    idx2 = txt.index("\n", idx1)
    return map(float,txt[idx1+26:idx2].split())

def get_cell(txt):
    idx1 = txt.index("!UNIT_CELL_CONSTANTS=")
    idx2 = txt.index("\n", idx1)
    return map(float,txt[idx1+22:idx2].split())


def getHKLv2(inp_str):
    """ Return a list of hkls. Functional Programming style."""
    headerEnd_markup = "!END_OF_HEADER\n"
    dataEnd_markup = "!END_OF_DATA\n"
    id_endHead = inp_str.find(headerEnd_markup)+len(headerEnd_markup)
    hklOnly = lambda s: map(int, s.split()[:3])
    lines = inp_str[id_endHead:-len(dataEnd_markup)].splitlines()
    return map(hklOnly, lines)

def calc_reso(hkl,cell):
    sum = Numeric.sum
    a = setmat(cell)
    return  1/(sum(a[0,:]*hkl,1)**2 +
               sum(a[1,:]*hkl,1)**2 + 
               sum(a[2,:]*hkl,1)**2)**0.5

class XDSReflectionFile:
    """
    >>> I=XDSReflectionFile("/home/legrand/proj/benchs/xds/BENCH_xds/xds_dir.2007-10-02-10h55/INTEGRATE.HKL")
    >>> I.header["OUTPUT_FILE"]
    'INTEGRATE.HKL'
    >>> I.getData()
    >>> I.calcCorrection()
    """

    def __init__(self, fileName):
        self.fileName = fileName
        self.raw = open(fileName).read()

        headEndMarkup = "!END_OF_HEADER\n"
        dataEndMarkup = "!END_OF_DATA\n"
        indexHeadEnd = self.raw.find(headEndMarkup) + \
                                       len(headEndMarkup)
        self.rawHeader = self.raw[:indexHeadEnd-len(headEndMarkup)]
        self.rawData = self.raw[indexHeadEnd:-len(dataEndMarkup)]
        self.header = self.getHeader()
        self.data = None

    def getHeader(self):
        _re_xdsPar = r"([^ ]+[=])"
        rec_xdsPar = re.compile(_re_xdsPar)

        lpar = []
        for line in self.rawHeader.splitlines():
            l_s = rec_xdsPar.split(line[1:])
            len_s = len(l_s)
            if len_s > 1 and len_s % 2:
                for i in range(1,len_s,2):
                    lpar.append((l_s[i][:-1],l_s[i+1].strip()))
        return dict(lpar)

    def write(self):
        hklfmt = " %d %d %d %.3E %.3E %.1f %.1f %.1f %.5f %d %d %d"
        hklfmt += " %.1f %.1f %.1f %.2f %.2f %.2f %.2f %.2f \n"
        newfilename = "INTEGRATE.HKL_ABSCORRECTED"
        t0 = time.time()
        f = open(newfilename, "w")
        f.write(self.rawHeader)
        f.write("!END_OF_HEADER\n")
        rawhkl = [hklfmt % tuple(self.data[i,:]) for i in range(self.nhkl)]
        f.write("".join(rawhkl))
        f.write("!END_OF_DATA\n")
        f.close()
        print " Time taken to write the INTEGRATE.HKL file: %.2f" % (time.time() - t0)


    def getData(self):
        nitem = int(self.header["NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD"])
        self.nitem = nitem
        t0 = time.time()
        _data = map(float,self.rawData.split())
        t1 = time.time()
        print "... Converting data to float:  %.3fs" % (t1 - t0)

        if not len(_data) % nitem:
            data = [_data[i:i+nitem] for i in range(0,len(_data),nitem)]
            t2 = time.time()
            print "... Reshaping:  %.3fs" % (t2 - t1)
        else:
            print "Error in the interpretation of reflection file"
            sys.exit()
        t3 = time.time()
        self.data = Numeric.array(data)
        t4 = time.time()
        self.nhkl = self.data.shape[0]
        print "... Creating Array:  %.3fs" % (t4 - t3)
        #print self.data.shape

    def calcCorrection(self, h=0.15, mu=None):
        """Correction for diamond absorption, to apply on the raw
        integrated intensities.
        h  is is the thickness of the diamond in cm
          it is  0.13cm for the ELSA dac
          and    0.15cm for the MARACA
        Estimated density of diamond: 3.52 g/cm**3
        mu is the absorption coeficient, calculated by http://csrri.iit.edu/mucal.html
          mu/rho = 0.6     cm**2/g at 16.53 keV ==  0.75 A
          mu/rho = 0.233  cm**2/g at 33.00 keV ==  0.34 A
        """

        if self.header["OUTPUT_FILE"] != "INTEGRATE.HKL":
            print "Error. Wrong type of input file."
            print "       It must be an INTEGRATE.HKL type of file!"
        if not self.data: self.getData()

        # geting the wavelenght from INTEGRATE.HKL
        wave = float(self.header["X-RAY_WAVELENGTH"])

        if not mu:
            try:
                mu=get_mu_from_mucalweb(wave, element)
            except:
                print "Error: Can't get mu value at this energy from the mucal website.!"
                print "Try to fix proxy settings, or call calcCorrection with explicite mu value."
                print "STOP"
                sys.exit()
        print ">> At the energy of %.3f keV, using mu value of: %.4f cm**2/gm" % (12.398/wave, mu)
        rotaxis = map(float, (self.header["ROTATION_AXIS"].split()))
        #y0 = float(self.header["ORGY"])
        rho = 3.52 # in g/cm**3
        D = float(self.header["DETECTOR_DISTANCE"])

        I = self.data[:,3]
        #X = -1 * qx * (self.data[:,5] - x0)
        x0 = float(self.header["ORGX"]) 
        y0 = float(self.header["ORGY"]) 
        qx = float(self.header["QX"]) # in mm
        qy = float(self.header["QY"]) # in mm
        if abs(abs(rotaxis[0]) - 1.) < 0.002:
            print ">> Rotation axis is:  HORIZONTAL"
            Y = -1 * qx * (self.data[:,5] - x0) 
            X =  1 * qy * (self.data[:,6] - y0)
        elif abs(abs(rotaxis[1]) - 1.) < 0.002:
            print ">> Rotation axis is:  VERTICAL"
            X =  1 * qx * (self.data[:,5] - x0) 
            Y =  1 * qy * (self.data[:,6] - y0)
        else:
            print "Rotation axis in an unsupported orientation: %s" % rotaxis
            sys.exit(2)
        dphi = float(self.header["OSCILLATION_RANGE"])
        phi0 = float(self.header["STARTING_ANGLE"])
        phi = (self.data[:,7] * dphi + phi0) * d2r
        #print min(phi), max(phi), sum(phi)/len(phi)
        correc1 = exp(-h*mu*rho/cos(phi))
        correc2 = exp(-h*mu*rho/((X*sin(phi)+D*cos(phi))/(X**2+Y**2+D**2)**(0.5)))
        #correc2 = exp(-h*mu*rho/((Y*sin(phi)+D*cos(phi))/(X**2+Y**2+D**2)**(0.5)))
        correc = correc1 * correc2

        print ">> Min, Max and Mean for the following values:\n"
        listing = [(X, "X   (mm)"),
                   (Y, "Y   (mm)"),
                   (phi/d2r, "Phi (deg)"),
                   (correc, "Transmi.")]
        for val, _name in listing:
            print "%12s  %7.3f  %7.3f  %7.3f" % \
                  (_name, min(val), max(val), sum(val)/len(val))
        print ">> Listing some reflections: H, K, L, I, sigma, X, Y, Z, Transmission"
        print "   (X and Y in pixel, Z in image_number)\n"
        for i in range(20):
                j= 200*i
                print "%4d%4d%4d %10.2f %10.2f %10.0f %10.0f %10.1f" % tuple(self.data[j,0:8]),
                print "%10.4f" % correc[j]
        self.data[:,3] = I/(correc)

def get_info(hklFileName, t0=None):

    refl = opReadCl(hklFileName)
    header = headerEnd_markup = "!END_OF_HEADER\n"
    cell = get_cell(refl[:2500])
    #res = read_resolution(refl[:2500])

    if t0:
        t1 = time.time()
        print "\n... Reading:          %.3fs" % (t1 - t0)
    hkls = getHKLv2(refl)
    if t0:
        t2 = time.time()
        print "... Extracting HKLs:  %.3fs" % (t2 - t1)

    try:
        1/0 
        res = read_resolution(refl[:2500])
    except:
        import Numeric
        hkls = Numeric.array(hkls)
        if t0:
            t3 = time.time()
            print "... Creating Array:   %.3fs" % (t3 - t2)
        res = calc_reso(hkls, cell)
        if t0:
            t4 = time.time()
            print "... Calculate Reso:   %.3fs" % (t4 - t3)

    res_low, res_high =  max(res), min(res)
    return {"res_low": res_low,
            "res_high": res_high,
            "refl_numb": len(hkls),
            "cell": cell,
            "hklFileName": hklFileName}


head_fmt = """
>>> Detector distance:        %(DETECTOR_DISTANCE)s
>>> Starting angle:           %(STARTING_ANGLE)s
>>> X-ray wavelength:         %(X-RAY_WAVELENGTH)s
"""

info_fmt = """
>>> Reflection File:          %(hklFileName)s
>>> Cell parameters:          %(cell_str)s
>>> Number of reflections:    %(refl_numb)d
>>> Low  resolution limit:    %(res_low).2f
>>> High resolution limit:     %(res_high).2f
"""

def _test():
    import doctest
    doctest.testmod()

def main():
    import getopt
    short_opt =  "dtl:eh"
    long_opt = ["debug","timeit","length=","energy=","help"]
    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)

    t0 = None
    _debug = False
    _length = False
    _energy = False

    for o, a in opts:
        if o in ("-h","--help"):
            print _usage
            sys.exit(2)
        if o in ("-t", "--timeit"):
            t0 = time.time()
            get_info(inputf[0], t0)
            sys.exit()
        if o in ("-d", "--debug"):
            _anomal = True
        if o in ("-l", "--length"):
            _length = float(a)
        if o in ("-e", "--energy"):
            _energy = float(a)

    if _debug:
        _test()
    if not inputf:
        print ">> Reading input file: INTEGRATE.HKL.save"
        I = XDSReflectionFile("INTEGRATE.HKL.save")
    else:
        print ">> Reading input file: %s" % inputf[0]
        I = XDSReflectionFile(inputf[0])

    print ">> Input file type:    %s" % I.header["OUTPUT_FILE"]
    I.getData()
    print ">> Calculating corrections for diamond path of:    %.3f cm" % _length
    I.calcCorrection(_length, _energy)
    I.write()

if __name__=='__main__':
    main()
