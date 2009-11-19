
__version__ = "0.2.0"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "06-07-2004"

from xupy import *

def unBined(param):
    newparam = XParam()
    newparam.mix(param)
    newparam.QX = param.QX/2
    newparam.QY = param.QY/2
    newparam.NX = param.NX*2
    newparam.NY = param.NY*2
    return newparam


DEFAULTS = """
ROTATION_AXIS= 1.0, 0.0, 0.0
DIRECTION_OF_DETECTOR_X_AXIS = 1.0, 0.0, 0.0
DIRECTION_OF_DETECTOR_Y_AXIS = 0.0, 1.0, 0.0
FRACTION_OF_POLARIZATION = 0.98
POLARIZATION_PLANE_NORMAL = 0.0, 1.0, 0.0
INCIDENT_BEAM_DIRECTION = 0.0, 0.0, 1.0
MINIMUM_VALID_PIXEL_VALUE= 0
OVERLOAD= 65500
VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000, 30000
"""

ADSC = XParam(DEFAULTS)
ADSC.mix("""
DETECTOR =    "ADSC"
MINIMUM_VALID_PIXEL_VALUE = 1
STRONG_PIXEL= 6.0
MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 6
""")

ADSC315 = XParam(DEFAULTS)
ADSC.mix("""
DETECTOR =    "ADSC"
MINIMUM_VALID_PIXEL_VALUE = 1
STRONG_PIXEL= 8.0
MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 7
NX = NY = 3072
QX = QY = 0.10260
""")




MARCCD = XParam(DEFAULTS)
MARCCD.mix("""
DETECTOR =   "CCDCHESS"
"""
)

MARMOSAIC = XParam()
MARMOSAIC.mix(MARCCD)
MARMOSAIC.mix("""
#   Typical parameters for a MARCCD MOSAIC 225
NX = NY =      3072    # number of pixel along X and Y
QX = QY =      0.0732    # length of a pixel in X and Y
""")


ESRF_BM14_OLD = XParam()
ESRF_BM14_OLD.mix(MARCCD)
ESRF_BM14_OLD.mix("""
#   Typical parameters for a MARCCD 133
NX = NY =      2048      # number of pixel along X and Y
QX = QY =      0.06469   # length of a pixel in X and Y
""")

ESRF_BM14 = XParam()
ESRF_BM14.mix(MARMOSAIC)
ESRF_BM14.mix("""
QX = QY =      0.07319   # length of a pixel in X and Y
""")

ESRF_ID23EH1 = XParam()
ESRF_ID23EH1.mix(ADSC315)
ESRF_ID23EH1.mix("""
NX = NY =     3072
""")

ESRF_ID23_1 = ESRF_ID23EH1


BIOXHIT = XParam()
BIOXHIT.mix(MARMOSAIC)
BIOXHIT.mix("""
NX = NY =      2048      # number of pixel along X and Y
QX = QY =      0.07917   # length of a pixel in X and Y
""")

ESRF_ID14EH3_OLD = MARCCD.copy()
ESRF_ID14EH3_OLD.mix("""
NX = NY = 2048
QX = QY = 0.06469
INCIDENT_BEAM_DIRECTION = -1.0, 0.0, 0.0
ROTATION_AXIS = 0.0, 1.0, 0.0
DIRECTION_OF_DETECTOR_X_AXIS = 0.0, 1.0, 0.0
DIRECTION_OF_DETECTOR_Y_AXIS = 0.0, 0.0, -1.0
POLARIZATION_PLANE_NORMAL = 0.0, 0.0, 1.0
""")

ESRF_BM30A = XParam()
ESRF_BM30A.mix(MARCCD)
ESRF_BM30A.mix("""
#   Typical parameters for a MARCCD 165
NX = NY =      2048      # number of pixel along X and Y
QX = QY =      0.079     # length of a pixel in X and Y
""")


ESRF_ID29 = XParam()
ESRF_ID29.mix(ADSC)
ESRF_ID29.mix("""
#   Typical parameters for a ADSC Q210
NX = NY =     3072       # number of pixel along X and Y
QX = QY =     0.10386     # length of a pixel in X and Y
""")

ESRF_ID29_old = XParam()
ESRF_ID29_old.mix(ADSC)
ESRF_ID29_old.mix("""
#   Typical parameters for a ADSC Q210
NX = NY =      2048      # number of pixel along X and Y
QX = QY =      0.1024    # length of a pixel in X and Y
""")


ESRF_ID14EH1 = XParam()
ESRF_ID14EH1.mix(ADSC)
ESRF_ID14EH1.mix("""
#   Typical parameters for a ADSC Q4 or Q4R
NX = NY =      2304      # number of pixel along X and Y
QX = QY =      0.0816    # length of a pixel in X and Y
""")

ESRF_ID14EH3 = ESRF_BM30A
ESRF_ID14EH2 = ESRF_ID14EH4 = ESRF_ID14EH1
ESRF_ID14_1 = ESRF_ID14EH1
ESRF_ID14_2 = ESRF_ID14EH2
ESRF_ID14_3 = ESRF_ID14EH3
ESRF_ID14_4 = ESRF_ID14EH4

ESRF_ID29_UNBINED = unBined(ESRF_ID29) ### Marche
ESRF_ID14EH1_UNBINED = unBined(ESRF_ID14EH1)
ESRF_ID14EH2_UNBINED = unBined(ESRF_ID14EH2)
#ESRF_ID14EH3_UNBINED = unBined(ESRF_ID14EH3)
ESRF_ID14EH4_UNBINED = unBined(ESRF_ID14EH4)


MarIP345 = XParam(DEFAULTS)
MarIP345.mix("""
DETECTOR =    "MAR345"
STRONG_PIXEL= 7.0
MINIMUM_VALID_PIXEL_VALUE = 0
OVERLOAD=130000
NX = NY =      3000      # number of pixel along X and Y
QX = QY =      0.139    # length of a pixel in X and Y
""")

MarIP = XParam()
MarIP.mix(MarIP345)
MarIP.mix("""
DETECTOR =    "MAR"
STRONG_PIXEL= 7.0
ROTATION_AXIS= 0.0, 1.0, 0.0
POLARIZATION_PLANE_NORMAL= 1.0, 0.0, 0.0
""")

MarIP300 = XParam()
MarIP300.mix(MarIP)
MarIP300.mix("""
NX = NY =      3000      # number of pixel along X and Y
QX = QY =      0.139    # length of a pixel in X and Y
""")

## Aliases
adsc = ADSC
marccd = MARCCD
MarCCD = MARCCD
MAR345 = MarIP345
mar = MarIP345
