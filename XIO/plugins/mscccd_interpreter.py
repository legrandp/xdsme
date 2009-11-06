
__version__ = "0.2"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "08-04-2007"
__copyright__ = "Copyright (c) 2007 Pierre Legrand"
__license__ = "LGPL"


def getEdgeResolution(PixelX, Width, Distance, Wavelength):
            "Calculate EdgeResolution: Graeme's Method"
            from math import sin, atan
            if Distance > 0.0:
                r = 0.5 * float(PixelX) * int(Width)
                return float(Wavelength)/sin(atan(r/float(Distance)))
            else:
                return 0.

def endian(code):
    if code == 'big_endian': return '>'
    else: return '<'
    
class Interpreter:
    
    HTD = {
    # The adsc Header Translator Dictionary.
    # Potential problems:
    # - There are multiple SIZE1, SIZE2 instances.
    # = The orientation of SIZE1 and SIZE2 is unknown
    #     Not a problem as long as SIZE1 = SIZE2..

    #'ExposureTime':(['TIME'],float),
    'BeamX':(['CCD_SPATIAL_BEAM_POSITION'], lambda x: float(x.split()[0])),
    'BeamY':(['CCD_SPATIAL_BEAM_POSITION'], lambda x: float(x.split()[1])),
    'Distance':(['CCD_GONIO_VALUES'], lambda x: float(x.split()[-1])),
    'Wavelength':(['SCAN_WAVELENGTH'], float),
    'PixelX':(['CCD_SPATIAL_DISTORTION_INFO'], lambda x: float(x.split()[-2])),
    'PixelY':(['CCD_SPATIAL_DISTORTION_INFO'], lambda x: float(x.split()[-1])),
    'Width':(['SIZE1'], int),
    'Height':(['SIZE2'], int),
    'Message':(['COMMENT'], str),
    'PhiStart':(['ROTATION'], lambda x: float(x.split()[0])),
    'PhiEnd':(['ROTATION'], lambda x: float(x.split()[1])),
    'PhiWidth':(['ROTATION'], lambda x: float(x.split()[2])),
    'HeaderSize':(['HEADER_BYTES'], int),
    'EndianType':(['BYTE_ORDER'], endian),
    }
    #'EdgeResolution':(['PIXEL_SIZE','SIZE1','DISTANCE','WAVELENGTH'],
    #    getEdgeResolution),

    # Added keys from Graeme's convention.
    #'TwoTheta':(['TWO_THETA'], float),   # _FIXME_ Not really here now...
    #'SerialNumber':(['DETECTOR_SN'], str),

    SpecialRules = {
    # No special rules for now
    }
    
    Identifiers = {
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    }
    
    def getRawHeadDict(self, rawHead):
        _lis = rawHead[2:].split("}")[0].split(";\n")[:-1]
        RawHeadDict = dict([par.split("=") for par in _lis])
        RawHeadDict.update({'MESSAGE':'','TWO_THETA':'0'}) # _FIXME_ Not really here now...
        return RawHeadDict
