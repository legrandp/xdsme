#!/usr/bin/env python
"""
    10/06/03 1st version legrand@embl-grenoble.fr
    
    TODO:
        - detector conversion
        - finding the good film rotation for each detector...
    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
"""

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "18-11-2005"
__copyright__ = "Copyright (c) 2005  Pierre Legrand"
__license__ = "LGPL"
__version__ = "0.4.7"

import math
import sys
import os.path

from XOconv import *

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Converting XDS crystal orientation informations to Denzo.

   A program to convert the orientation matix, extract information
   from XDS output CORRECT.LP files and write denzo pseudo dot.x files.

   USAGE:   %s  [OPTION]... CORRECT.LP
    
   OPTIONS
    
   -h
   --help
         Print this help message.
   
   -i HKLFILE
   --input-reflections HKLFILE
         Read an XDS reflection file (INTEGRATE.HKL or XDS_ASCII.HKL)
         and output some of the reflection in the .x file. Selection
         of the reflections is made throug the -n, --number option.
    
    -n
    --number
         Image number selection for the output of pseudo dot.x file.
             -n x1           (select image number x1)
             -n x1:x2        (select images number x1 to x2)
             -n x1,x2,x3     (select image number x1, x2 and x3)
         Note: This option goes with the -i, --input-reflections option.
    
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

mosaicity_conversion_factor = 4.410

# In XDS the mosaicity is discribed as:
# "the standard deviation (e.s.d) of a Gaussian modeling the rocking curve.
# If the limit of integration is set (by default) to 2% of reflection profile
# maximum, from the guassian function this correspond to a mosaic spread of
# 4.41 times the e.s.d.


FMT_dotx_mat1 = 3*"%15.8f"
FMT_dotx_mat2 = 3*"%10.6f"+"\n"
FMT_dotx_par = 4*"%12.5f"+4*"%10.4f"+6*"%2d"+"\n"
FMT_refl = ("%4d%4d%4d" +
            " 0%8.1f   100.0   1.00%6.1f %5.3f" +
            "%7.1f%7.1f 1.000    99.9\n")

FMT_dotx_end = """ 999%(boxe_size_x)4d%(boxe_size_y)4d
%(spot_profile)s
spot elliptical %(spot_radius)5.3f %(spot_radius)5.3f    0.0
background elliptical %(background_radius)5.3f %(background_radius)5.3f    0.0
goniostat single axis
goniostat orientation    0.000    0.000
motor axis  0.000000  1.000000  0.000000
profile fitting radius   22.00
resolution limits  %(resolution_limits)s
wavelength  %(wavelength).5f
monochromator   %(monochromator).3f
oscillation start  %(phi_start).2f end  %(phi_end).2f
mosaicity   %(mosaicity).3f
spindle axis%(spindle)s vertical axis%(vertical)s
unit cell    %(cell)s
crystal rotx    %(rotx).2f roty  %(roty).2f rotz  %(rotz).2f 
sector        1
raw data file '%(image_template)s'
format %(detector)s
space group %(spg)s
cassette rotx   0.0 roty    0.0 rotz    %(cassette_rotz).2f 2 theta    0.00
distance %(distance_mm)7.2f
x beam %(beam_x)8.3f  y beam %(beam_y)8.3f
film rotation   %(film_rotation).2f
crossfire y  0.000 x  0.000 xy  0.000
"""

FMT_strategy_inp = """XFILE %(xfile_name)s
COMPLETENESS 100 99 98 97 96 95
NBINS 10
>NFRAME %(nframe)d
PSFILE %(psfile)s
FORMAT %(detector)s
SPACE GROUP %(spg)s
>RESOLUTION %(resolution_limits)s
%(anomalous)s
"""

detector2denzo = {"CCDCHESS":"ccd unsupported-mar",
       "ADSC":"ccd adsc unsupported-q4",
       "RAXIS":"RAXIS", "MAR":"SMALLMAR","MAR345":"mar345 345mm 100 marpck",
       "MAC":"DIP2000", "CCDBRANDEIS":"ADSC", "CCDD2AM":"ccd unsupported-esrf",
       "UNKNOWN":"unknown"}

def detector2strategy(det_type, detector_diameter):
    if det_type == "CCDCHESS" and 166. > detector_diameter > 159.:
        return "MARCCD165"
    elif det_type == "CCDCHESS" and 227. > detector_diameter > 223.:
        return "MAR225"
    elif det_type == "CCDCHESS" and 135. > detector_diameter > 130.:
        return "MARCCD133"
    elif det_type == "ADSC" and 189. > detector_diameter > 187.:
        return "ADSCQ4R"
    elif det_type == "ADSC" and 211. > detector_diameter > 209.:
        return "ADSCQ210"
    elif det_type == "ADSC" and 316. > detector_diameter > 314.:
        return "ADSCQ315"
    elif det_type == "CCDBRANDEIS" and 201. > detector_diameter > 199.:
        return "BRAND2X2"
    elif det_type == "CCDBRANDEIS" and  96. > detector_diameter >  94.:
        return "CCDB"
    elif det_type == "MAC":
        return "DIP2"
    elif det_type == "CCDD2AM":
        return "CCDD2AM"
    else:
         return "UNKNOWN"
    #"RAXIS":"RAXIS",
    #"MAR":"SMALLMAR","MAR345":"mar345 345mm 100 marpck"

film_rotation_DB = {"MARCCD133":-90., "MARCCD165":0.,"MAR225":90.,
                    "ADSCQ4R":90., "ADSCQ210":90.,
                    "ADSCQ315":90., "DIP2":-90., "CCDD2AM":0,
                    "UNKNOWN":0.}

cassette_rotz_DB = {"MARCCD133":0., "MARCCD165":90.,"MAR225": 0.,
                    "ADSCQ4R":0., "ADSCQ210":0.,
                    "ADSCQ315":0., "DIP2":0., "CCDD2AM":0,
                    "UNKNOWN":0.}


def cercle_2D(arg):
    (x, y), (x0, y0, radius) = arg
    if ((x-x0)**2 + (y-y0)**2)**0.5 <= radius: return 1
    else: return 0

def spot_profile(boxSize, backgroundRadius, spotRadius):
    #
    x0 = boxSize[1]/2.-0.5
    y0 = boxSize[0]/2.-0.5
    size = boxSize[0]*boxSize[1]
    fX = lambda a: [a[0],]*a[1]
    lX = zip(range(boxSize[1]),boxSize[1]*[boxSize[0],])
    X = reduce(lambda a,b: a+b, map(fX, lX))
    Y = range(boxSize[0])*boxSize[1]
    matrixCorrdinate = zip(X,Y)

    dataBack = zip(matrixCorrdinate, ((x0, y0, backgroundRadius),)*size)
    dataSpot = zip(matrixCorrdinate, ((x0, y0, spotRadius),)*size)
    backgroundBox = map(cercle_2D, dataBack)
    spotBox = map(cercle_2D, dataSpot)
    #
    circles = zip(backgroundBox, spotBox)
    boxNum = map(lambda c: 1-c[0]+2*c[1],circles)
    boxFmt = boxSize[1]*(" "+"%s"*boxSize[0]+"\n")
    boxStr = boxFmt % tuple(boxNum)
    return boxStr

def write_dotx(dnz, dotx_fname, hklxy_str):
    dotx = open(dotx_fname,"w")
    dotx.write("   Generated with XDS2DNZ version %s\n" % (__version__))
    for l in range(3):
        dotx.write(FMT_dotx_mat1 % tuple(dnzDict['Amat'].getRow(l)))
        dotx.write(FMT_dotx_mat2 % tuple(dnzDict['Umat'].getRow(l)))
    #
    dotx.write(FMT_dotx_par % (dnzDict['phi_start'], dnzDict['phi_end'],
                   dnzDict['distance_pix'],
                   dnzDict['wavelength'], dnzDict['rotz'], dnzDict['roty'], dnzDict['rotx'], 
                   dnzDict['mosaicity'],
                   -1,1,-1,-1,-1,-1))
    dotx.write(hklxy_str)    
    dotx.write(FMT_dotx_end % dnz)
    dotx.close()

def d_min(detDiameter, detDistance, wavel):
    return 1/(2*math.sin(math.atan((detDiameter/2.)/detDistance)/2.)/wavel)

def read_xds_hkl(hkl_str, ftype = "INTEGRATE"):
    hkl_str = hkl_str[hkl_str.index("!END_OF_HEADER")+15:
                      hkl_str.index("!END_OF_DATA")]
    hkl_num = Numeric.array(map(float,hkl_str.split()))
    _len = len(hkl_num)
    if ftype == "INTEGRATE":hkl_num.shape = (_len/19,19)
    if ftype == "XDS_ASCII":hkl_num.shape = (_len/11,11)
    return hkl_num

def PARS_xds2dnz(xdsPar):
    "Convert XDS output parameters to Denzo parameters."

    dnz = {}
    dnz['wavelength'] = xdsPar["wavelength"]
    dnz['delta_phi'] = xdsPar["delta_phi"]
    dnz['pixel_size_x'] = xdsPar["pixel_size"][0]
    dnz['pixel_size_y'] = xdsPar["pixel_size"][1]
    dnz['distance_pix'] = xdsPar["distance"]/dnz['pixel_size_x']
    dnz['distance_mm'] = xdsPar["distance"]
    dnz['phi_init'] = xdsPar["phi_init"]
    dnz['delta_phi'] = xdsPar["delta_phi"]
    dnz['mosaicity'] = mosaicity_conversion_factor*xdsPar["mosaicity"]
    dnz['beam_x'] = xdsPar["origin"][1]*xdsPar["pixel_size"][1]
    dnz['beam_y'] = xdsPar["origin"][0]*xdsPar["pixel_size"][0]
    dnz['image_template'] = xdsPar["template"]
    dnz['spg'] = SPGlib[xdsPar["symmetry"]][0].lower()
    dnz['monochromator'] = xdsPar["polarization"]
    dnz['spot_radius'] = dnz['distance_mm']*math.tan(2.5*xdsPar["divergence_esd"]/r2d)
    dnz['diameter'] = min(xdsPar["pixel_numb"][0]*dnz['pixel_size_x'],
                       xdsPar["pixel_numb"][1]*dnz['pixel_size_y'])
    dnz['cell'] = 6*"%.3f  " % tuple(xdsPar["cell"])
    dnz['rotx'], dnz['roty'], dnz['rotz'] = rotxyz
    dnz['detector'] = detector2denzo[xdsPar["detector_type"]]
    resolution_limits = tuple(xdsPar["resolution_range"])
    if min(resolution_limits) == 0.:
        resolution_limits = (max(resolution_limits), 
               d_min(dnz['diameter'], dnz['distance_mm'], dnz['wavelength']))
    dnz['resolution_limits'] = "%.2f  %.2f" % resolution_limits
    _detector = detector2strategy(xdsPar["detector_type"],dnz['diameter'])
    dnz['film_rotation'] = film_rotation_DB[_detector]
    dnz['cassette_rotz'] = cassette_rotz_DB[_detector]
    

    # Determining denzo spot profile
    spot_radius = int(dnz['spot_radius']/(dnz['pixel_size_x'])+0.5)
    bkg_radius  = spot_radius+1
    dnz['boxe_size_x'] = (spot_radius+3)*2
    dnz['boxe_size_y'] = (spot_radius+3)*2
    dnz['background_radius'] = bkg_radius*dnz['pixel_size_x']
    dnz['spot_profile'] = spot_profile((dnz['boxe_size_x'],dnz['boxe_size_y']),
                                     bkg_radius, spot_radius)[:-1]
    return dnz

if __name__=='__main__':

    import getopt
    
    _do_PG_permutations = False
    _verbose = False
    _inputReflectionName = ""
    _template = "xds2dnz"
    _image_num_list = 1,
    _vertical = (1, 0, 0)
    _spindle =  (0, 0, 1)
    hkl_xds = ""

    short_opt =  "hpvi:n:S:V:"
    long_opt = ["help", "input-reflections=", "pg-permutations",
                "verbose", "number=", "spindle=", "vertical="]

    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)
        
    for o, a in opts:
        if o in ("-v", "--verbose"):
            _verbose = True
        if o in ("-i", "--input-reflections"):
            _inputReflectionName = a
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
        if o in ("-n", "--number"):
            if a.count(":") == 1:
                try:
                    II, IF = map(int,(a.split(":")))
                    _image_num_list = range(II, IF+1)
                except:
                    print "Error! Can't parse -n or --number option."
                    print _usage
                    sys.exit(2)
            elif a.count(","):
                try:
                    _image_num_list = map(int,(a.split(",")))
                except: pass
            else:
                try:
                    _image_num_list = [int(a),]
                except: pass
        
    print ">>> Selected Image(s) = ", _image_num_list
    
    infilelocation, infilename = os.path.split(inputf[0])

    XDSi = XDSParser()
    DNZi = DenzoParser()
    
    try:
        XDSi.parse_correct(inputf[0])
    except:
        print "\nERROR! Can't parse input file: %s\n" % inputf[0]
        print _usage
        sys.exit(2)
    
    XDSi.debut()
    Adnz = XDSi.UBxds_to_dnz()
    Udnz = DNZi.Adnz_to_Udnz(Adnz)
   
    if _inputReflectionName and os.path.isfile(_inputReflectionName):
        
        import Numeric
        
        if _inputReflectionName.count("INTEGRATE.HKL"):
            iZD, irXYD, ncol = 7, 5, 18
            hklType = "INTEGRATE"
        elif _inputReflectionName.count("XDS_ASCII.HKL"):
            iZD, irXYD, ncol = 7, 4, 10
            hklType = "XDS_ASCII"
        
        print ">>> Reading reflections in file ", _inputReflectionName
        hkl_xds = read_xds_hkl(openReadClose(_inputReflectionName), hklType)
        
        D = XDSi.dict["distance"]/XDSi.dict["pixel_size"][0]
        XY = hkl_xds[:,irXYD:irXYD+2][:,::-1]
        X, Y = XY[:,0], XY[:,1]
        X0, Y0 = XDSi.dict["origin"][1], XDSi.dict["origin"][0]
        # OBL is the cosine of incidence angle at the detector
        # In this case the Detector is supposed perpendicular to the beam.
        OBL = D/(D**2+(X-X0)**2+(Y-Y0)**2)**0.5
        OBL.shape = (OBL.shape[0],1)
        hklxy = Numeric.concatenate((hkl_xds[:,:5], OBL, XY),1)
        # it seems that denzo start counting pixels at [0,0] whereas XDS
        # start at [1,1] ???
        # hklxy[:,5:6] = hklxy[:,5:6] -1
    else:
        print ">>> Warning: Can't find reflection file: ", _inputReflectionName
    
    vertical = vec3(_vertical)
    spindle =  vec3(_spindle)
    
    rcell = reciprocal(XDSi.dict["cell"])
    cell = XDSi.dict["cell"]

    Bdnz = DNZi.get_B(rcell=rcell)
    Udnz0 = DNZi.get_U0(rcell, vertical, spindle)
    rotxyz = DNZi.UB_to_Rotxyz(Adnz, Udnz0, [0.,1.,0.])
    
    verif = is_orthogonal(Udnz)
    if not verif:
        print "???  Warning: The U matrix is not orthogonal."
        print "???  Epsilon error: %.1e" % verif
    
    pointGroup = SPGlib[XDSi.dict["symmetry"]][3]
    print  ">>> Space Group : %s" % (SPGlib[XDSi.dict["symmetry"]][1])
    print  ">>> Space Group number: %s" % (XDSi.dict["symmetry"])
    print  ">>> Point group: %s" % (pointGroup)
    print  ">>> Number of equivalent crystal ortientations: %d\n" % \
                                         (len(PGequiv[pointGroup])+1)
    _fmtrot = "  rotx%8.2f roty%8.2f rotz%8.2f"
    _fmtmat = 3*"%12.8f"+"\n"+2*(18*" "+3*"%12.8f"+"\n")
    print ">>> crystal:      " + _fmtrot % tuple(rotxyz)
    print ">>> Adnz:         " + _fmtmat % tuple(Adnz.mlist)
    
    if _do_PG_permutations:
        for equiv in PGequiv[pointGroup]:
            # Note that equiv_mat is the transpose of the normal equiv matrix!
            # because of the way mat3 matrices are contructed.
            equiv_mat = mat3(equiv[0],equiv[1],equiv[2])
            _Adnz = Adnz * equiv_mat
            _cell = reciprocal(UB_to_cellParam(_Adnz))
            
            # One way to verify that the permutation is correct is
            # that we can still extract the cell parameters from the
            # permuted UB (or Adnz) matrix
            assert abs(_cell[0] - cell[0]) < 1e-2 and \
                   abs(_cell[1] - cell[1]) < 1e-2 and \
                   abs(_cell[2] - cell[2]) < 1e-2 and \
                   abs(_cell[3] - cell[3]) < 2e-2 and \
                   abs(_cell[4] - cell[4]) < 2e-2 and \
                   abs(_cell[5] - cell[5]) < 2e-2 

            _rotxyz2 =  DNZi.UB_to_Rotxyz(_Adnz, Udnz0)
            print ">>> crystal_equiv:" + _fmtrot % tuple(_rotxyz2)
            print ">>> Adnz_equiv:   " + _fmtmat % tuple(_Adnz.mlist)
    
    dnzDict = PARS_xds2dnz(XDSi.dict)
    dnzDict['Umat'] = Udnz
    dnzDict['Amat'] = Adnz
    dnzDict['spindle'] = "%4d%4d%4d" % tuple(spindle)
    dnzDict['vertical'] = "%4d%4d%4d" % tuple(vertical)
   
    # Selecting reflections and Writing dotx files
    prefix = dnzDict['image_template'].split("/")[-1].split(".")[0].replace("#","")
    dotx_fname_fmt = prefix + "%.3d.x.xds"
    
    if hkl_xds:
        for image_num in _image_num_list:
            dotx_fname = dotx_fname_fmt % image_num
            dnzDict['phi_end'] = dnzDict['phi_init'] + \
                                 image_num * dnzDict['delta_phi']
            dnzDict['phi_start'] = dnzDict['phi_end'] - dnzDict['delta_phi']
            frame_filter = (image_num > hkl_xds[:,iZD])*\
                           (hkl_xds[:,iZD] >= (image_num-1))
            hklxy_select = Numeric.compress(frame_filter,hklxy,0)      
            nref = hklxy_select.shape[0]
            out_fmt1 = ">>> Number of reflection selected for frame %4d : %d"
            print  out_fmt1 % (image_num,nref)
            hklxy_str = nref*FMT_refl  % tuple(Numeric.ravel(hklxy_select))
            write_dotx(dnzDict, dotx_fname, hklxy_str)
            
    else:
        image_num = 1
        dotx_fname = dotx_fname_fmt % image_num
        dnzDict['phi_end'] = dnzDict['phi_init'] + \
                             image_num * dnzDict['delta_phi']
        dnzDict['phi_start'] = dnzDict['phi_end'] - dnzDict['delta_phi']
        hklxy_str = ""
        write_dotx(dnzDict, dotx_fname, hklxy_str)
        

    # Preparing STRATEGY input file.
    str_inpfname = "auto_strategy.inp"
    STR = {}
    STR['xfile_name'] = dotx_fname_fmt % _image_num_list[0]
    STR['psfile'] = str_inpfname.split(".")[0]+".ps"
    STR['spg'] = dnzDict['spg'].upper()
    STR['nframe'] = 90
    STR['resolution_limits'] = dnzDict['resolution_limits']
    STR['detector'] = detector2strategy(XDSi.dict["detector_type"],
                                                 dnzDict['diameter'])
    if XDSi.dict["friedel"] == "FALSE":
        STR['anomalous'] = "ANOMALOUS"
    else: 
        STR['anomalous'] = ""

    openWriteClose("auto_strategy.inp", FMT_strategy_inp % STR)
