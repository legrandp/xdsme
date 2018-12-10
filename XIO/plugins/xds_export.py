# -*- coding: utf-8 -*-

""" XIO plugin for the export parameters as an XDS.INP format
    See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_prepare.html
"""

__version__ = "0.4.4"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "22-11-2017"
__copyright__ = "Copyright (c) 2007-2017 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import time
import os

from pycgtypes import vec3
from pycgtypes import mat3

EX, EY, EZ = vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1)
V3FMT = "%9.6f %9.6f %9.6f"
PI = 3.1415926535897931
D2R = PI/180.

def set_detplugin_lib(dettype):
    """Find out if we can set the LIB keyword"""
    if dettype == "hdf5dec" and "XDS_LIB_HDF5DEC" in os.environ:
        return os.environ["XDS_LIB_HDF5DEC"]
    else:
        return None

def det_dist(distance, dettype):
    "Return the disance with the proper sign."
    detori = XDS_DETECTOR_DICT["orient"][dettype]
    return distance*detori[2]

def det_spindle(dettype):
    "Return the spindle axis vector."
    return V3FMT % tuple(XDS_DETECTOR_DICT["orient"][dettype][3])

def polarization(wavelength):
    "Guess the polarization fraction from the wavelength."
    if 1.5417 < wavelength < 1.5420:
        return 0.5
    else:
        return 0.99

def det_axis_x(twotheta, dettype):
    "Return the X axis vector of the Dector surface."
    detori = XDS_DETECTOR_DICT["orient"][dettype]
    return V3FMT % tuple(detori[0]*mat3().rotation(twotheta*D2R, detori[4]))

def det_axis_y(twotheta, dettype):
    "Return the Y axis vector of the Dector surface."
    detori = XDS_DETECTOR_DICT["orient"][dettype]
    return V3FMT % tuple(detori[1]*mat3().rotation(twotheta*D2R, detori[4]))

def det_beam_x(x0, y0, qx, qy, dettype):
    "Return the X position of the beam coordinate."
    _def = XDS_DETECTOR_DICT["orient"][dettype][5]
    orgx, orgy = x0/qx, y0/qy
    if _def == "XY":
        return orgx
    elif _def == "YX":
        return orgy

def det_beam_y(x0, y0, qx, qy, dettype):
    "Return the Y position of the beam coordinate."
    _def = XDS_DETECTOR_DICT["orient"][dettype][5]
    orgx, orgy = x0/qx, y0/qy
    if _def == "XY":
        return orgy
    elif _def == "YX":
        return orgx


XDS_DETECTOR_DICT = {
  "detector_name":{
    "mar":       "MAR345",
    "mar555":    "MAR345",
    "marccd":    "CCDCHESS",
    "adsc":      "ADSC",
    "raxis":     "RAXIS",
    "minicbf":   "PILATUS",
    "hdf5dec":   "EIGER",
    "mscccd":    "SATURN",
    "mscpilatus":"PILATUS",
  },
  "overload":{
    "mar":       130000,
    "mar555":    250000,
    "marccd":     65000,
    "adsc":       65000,
    "raxis":     262100, # for raxisII. raxisIV: 1000000, raxisV: 2000000.
    "minicbf":  1048500,
    "hdf5dec":  1000000,
    "mscccd":   1000000,
    "mscpilatus":1048500,
  },
  "minval":{
    "mar":        0,
    "mar555":     0,
    "marccd":     1,
    "adsc":       1,
    "raxis":      0,
    "minicbf":    0,
    "hdf5dec":    0,
    "mscccd":     1,
    "mscpilatus": 0,
  },
  "min_number_of_pixels":{
    "mar":        8,
    "mar555":     4,
    "marccd":     8,
    "adsc":       8,
    "raxis":      8,
    "minicbf":    3,
    "hdf5dec":    4,
    "mscccd":     8,
    "mscpilatus": 3,
  },
  "sensor_thickness":{
    "mar":       0,
    "mar555":    0,
    "marccd":    0,
    "adsc":      0,
    "raxis":     0,
    "minicbf":   0.32,
    "hdf5dec":   0.45,
    "mscpilatus":0.45,
  },
  "orient":{ # X_det, Y_det, distanceSign, spindle_axis, twoThetaAxis, beamdef
    "mar":       ( EX, EY,  1,  EX, -EX, "XY"),
    "mar555":    ( EX, EY,  1,  EX, -EX, "XY"),
    "marccd":    ( EX, EY,  1,  EX, -EX, "YX"),
    "adsc":      ( EX, EY,  1,  EX, -EX, "YX"),
    "raxis":     ( EX, EY,  1,  EY,  EY, "XY"),
    "minicbf":   ( EX, EY,  1,  EX,  EX, "XY"),
    "mscccd":    (-EX, EY, -1,  EY,  EY, "XY"),
    "hdf5dec":   ( EX, EY,  1,  EX,  EX, "XY"),
    "mscpilatus":( EX, EY,  1, -EY,  EX, "XY"),
  }
}

SPECIFIC_SUPPLEMENTARY_KEYWORDS = {
 "Eiger": "GAIN= 1.0\n",
 "Pilatus": "GAIN= 1.0\n",
 "PILATUS 6M": """!SPECIFIC KEYWORDS FOR PILATUS 6M
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 5500 30000
  
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
 
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13
 
! X-GEO_CORR= ../hole_mask/x_geo_corr.cbf   
! Y-GEO_CORR= ../hole_mask/y_geo_corr.cbf   

 \n""",
 "457":"""
! AS MX1 Reverse Phi
ROTATION_AXIS= -1.0 0.0 0.0
\n""",
 "928":"""
! AS MX2 Reverse Phi
ROTATION_AXIS= -1.0 0.0 0.0
\n"""
}

TEMPLATE = """! File Automaticaly generated by XIO
!       date: %s
!       Beamline SOLEIL-Proxima1
""" % (time.ctime())

TEMPLATE += """
 JOB= ALL ! XYCORR INIT COLSPOT IDXREF DEFPIX XPLAN INTEGRATE CORRECT

 NAME_TEMPLATE_OF_DATA_FRAMES= %(NAME_TEMPLATE_OF_DATA_FRAMES)s
 DATA_RANGE= %(_DRI)d %(_DRF)d
 SPOT_RANGE= %(_DRI)d %(_DRF)d
 BACKGROUND_RANGE= %(_DRI)d %(_DRB)d

 OSCILLATION_RANGE= %(OSCILLATION_RANGE).3f
 STARTING_ANGLE= %(STARTING_ANGLE).3f
 STARTING_FRAME= %(_DRI)d
 
 X-RAY_WAVELENGTH= %(X_RAY_WAVELENGTH).5f
 DETECTOR_DISTANCE= %(DETECTOR_DISTANCE).2f
 
 DETECTOR= %(DETECTOR)s
    MINIMUM_VALID_PIXEL_VALUE= %(MINIMUM_VALID_PIXEL_VALUE)d
    OVERLOAD= %(OVERLOAD)d
    DIRECTION_OF_DETECTOR_X-AXIS= %(DIRECTION_OF_DETECTOR_X-AXIS)s
    DIRECTION_OF_DETECTOR_Y-AXIS= %(DIRECTION_OF_DETECTOR_Y-AXIS)s
    NX= %(NX)d    NY= %(NY)d
    QX= %(QX).5f  QY= %(QY).5f
    ORGX= %(ORGX).2f   ORGY= %(ORGY).2f
    SENSOR_THICKNESS= %(SENSOR_THICKNESS)s

 ROTATION_AXIS= %(ROTATION_AXIS)s
 INCIDENT_BEAM_DIRECTION= 0.0 0.0 1.0
 FRACTION_OF_POLARIZATION= %(FRACTION_OF_POLARIZATION).3f
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0

 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= %(MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT)d
 STRONG_PIXEL= 8.0
 SPOT_MAXIMUM-CENTROID= 2.0

 SPACE_GROUP_NUMBER= 0
 UNIT_CELL_CONSTANTS= 0 0 0 0 0 0

 INCLUDE_RESOLUTION_RANGE= 50.0 0.0
 TRUSTED_REGION= 0.0 1.42
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000

 DELPHI= %(DELPHI).2f
 REFINE(IDXREF)= POSITION BEAM AXIS ORIENTATION CELL
 REFINE(INTEGRATE)= POSITION BEAM ORIENTATION CELL
 REFINE(CORRECT)= POSITION BEAM AXIS ORIENTATION CELL
 
 MAXIMUM_NUMBER_OF_PROCESSORS= 16
 MAXIMUM_NUMBER_OF_JOBS= 1
 RESOLUTION_SHELLS=20.0 10.0 6.0 %(_HIGH_RESOL_LIMIT).2f
 TOTAL_SPINDLE_ROTATION_RANGES=15.0 120.0 15.0
 STARTING_ANGLES_OF_SPINDLE_ROTATION=-95.0 95.0 5.0
 PROFILE_FITTING= TRUE
 STRICT_ABSORPTION_CORRECTION= TRUE
 FRIEDEL'S_LAW= TRUE
 TEST_RESOLUTION_RANGE= 20 4.0
 ! REFERENCE_DATA_SET=
 
 %(SPECIFIC_KEYWORDS)s
"""


#     Header Translator Dictionary.
#     Translate image header entries in a new dictionay
#     newdic['X_RAY_WAVELENGTH'] = float(head['Wavelength'])
#
#def TEST(x):
#  print x, type(x)
#  return round(x,2)

HTD = {
'X_RAY_WAVELENGTH':(['Wavelength'], float),
'DETECTOR_DISTANCE':(['Distance','ImageType'], det_dist),
'ROTATION_AXIS':(['ImageType'], det_spindle),
'FRACTION_OF_POLARIZATION': (['Wavelength'], polarization),
'STARTING_ANGLE':(['PhiStart'], float),
'OSCILLATION_RANGE':(['PhiWidth'], float),
'NX':(['Width'], int),
'NY':(['Height'], int),
'QX':(['PixelX'], float),
'QY':(['PixelY'], float),
'ORGX':(['BeamX','BeamY','PixelX','PixelY','ImageType'], det_beam_x),
'ORGY':(['BeamX','BeamY','PixelX','PixelY','ImageType'], det_beam_y),
'DELPHI':(['PhiWidth'], lambda x: 16*x),
'DIRECTION_OF_DETECTOR_X-AXIS':(['TwoTheta','ImageType'], det_axis_x),
'DIRECTION_OF_DETECTOR_Y-AXIS':(['TwoTheta','ImageType'], det_axis_y),
'_HIGH_RESOL_LIMIT':(['EdgeResolution'], lambda x: round(x,2)),
'SENSOR_THICKNESS':(['SensorThickness'], float)
}

#     Collect Translator Dictionary.
#     Translate collect object attributes to a new dictionay
#     newdic['SPOT_RANGE'] = list(collect.imageRanges)
#
CTD = {
'NAME_TEMPLATE_OF_DATA_FRAMES':(['xdsTemplate'], str),
'DATA_RANGE':(['imageNumbers'], lambda x: [x[0],x[-1]]),
'_DRI':(['imageNumbers'], lambda x: x[0]),
'_DRF':(['imageNumbers'], lambda x: x[-1]),
'_DRB':(['imageNumbers'], lambda x: min(x[0]+9,x[-1])),
'SPOT_RANGE':(['imageRanges'], list),
'BACKGROUND_RANGE':(['imageRanges'], lambda x: \
                                        [x[0][0], min(x[0][0]+7, x[0][1])]),
'DETECTOR':(['imageType'], lambda x: XDS_DETECTOR_DICT["detector_name"][x]),
'NAME_TEMPLATE_OF_DATA_FRAMES':(['xdsTemplate'], str),
'SENSOR_THICKNESS':(['imageType'], lambda x: \
                                   XDS_DETECTOR_DICT["sensor_thickness"][x]),
'MINIMUM_VALID_PIXEL_VALUE':(['imageType'], lambda x: \
                                             XDS_DETECTOR_DICT["minval"][x]),
'MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT':(['imageType'], lambda x: \
                               XDS_DETECTOR_DICT["min_number_of_pixels"][x]),
'OVERLOAD':(['imageType'], lambda x: XDS_DETECTOR_DICT["overload"][x]),
'_LIB':(['imageType'], set_detplugin_lib)
}
