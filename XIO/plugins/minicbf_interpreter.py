
__version__ = "0.0.3"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "03-11-2009"
__copyright__ = "Copyright (c) 2009 Pierre Legrand"
__license__ = "LGPL"

import time

def _datetime(timestr):
    t1,t2 = timestr.split(".")
    return time.strptime(t1, "%Y/%b/%d %H:%M:%S"),float("0."+t2)
    
def _dateseconds(timestr):
    t1,msec = _datetime(timestr)
    return time.mktime(t1)+msec
    
#def _timestr():
#    pass

def getEdgeResolution(PixelX, Width, Distance, Wavelength):
            "Calculate EdgeResolution: Graeme's Method"
            from math import sin, atan
            if Distance > 0.0:
                r = 0.5 * float(PixelX) * int(Width)
                return float(Wavelength)/sin(atan(r/float(Distance)))
            else:
                return 0.

float1 = lambda x: float(x.split()[0])
float2 = lambda x: float(x.split()[0])*1e3
beamx = lambda x,y: float(x[x.find("(")+1:x.find(")")-1].split(",")[0])*float2(y)
beamy = lambda x,y: float(x[x.find("(")+1:x.find(")")-1].split(",")[1])*float2(y)


def endian(code):
    if code == 'big_endian': return '>'
    else: return '<'

_hKeys = ["Detector:", "Pixel_size", "Silicon", "Exposure_time", "Exposure_period",
"Tau", "Count_cutoff","Threshold_setting","N_excluded_pixels","Excluded_pixels:",
"Flat_field:","Trim_directory:","Wavelength","Energy_range","Detector_distance",
"Detector_Voffset","Beam_xy","Flux","Filter_transmission","Start_angle",
"Angle_increment","Detector_2theta","Polarization","Alpha","Kappa","Phi","Chi",
"Oscillation_axis","N_oscillations"]

 
class Interpreter:
    
    HTD = {
    # The adsc Header Translator Dictionary.
    # Potential problems:
    # - There are multiple SIZE1, SIZE2 instances.
    # = The orientation of SIZE1 and SIZE2 is unknown
    #     Not a problem as long as SIZE1 = SIZE2..

    'ExposureTime':(['Exposure_time'], float1),
    'BeamX':(['Beam_xy', 'Pixel_size'], beamx),
    'BeamY':(['Beam_xy', 'Pixel_size'], beamy),
    'Distance':(['Detector_distance'], float2),
    'Wavelength':(['Wavelength'], float1),
    'PixelX':(['Pixel_size'], float2),
    'PixelY':(['Pixel_size'], float2),
    'Width':(['Binary-Size-Fastest-Dimension'], int),
    'Height':(['Binary-Size-Second-Dimension'], int),
    #'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['Start_angle'], float1),
    'PhiEnd':(['Start_angle','Angle_increment'], lambda x,y: float1(x)+float1(y)),
    'PhiWidth':(['Angle_increment'], float1),
    #'EdgeResolution':(['PIXEL_SIZE','SIZE1','DISTANCE','WAVELENGTH'],
    #    getEdgeResolution),

    # Added keys from Graeme's convention.
    'TwoTheta':(['Detector_2theta'], float1),   # _FIXME_ Not really here now...
    'SerialNumber':(['Detector:'], str),
    'HeaderSize':(['HEADER_SIZE'], int),
    'DateStr':(['DATE'], str),
    'DateSeconds':(['DATE'], _dateseconds),
    #'EndianType':(['BYTE_ORDER'], endian),
    }

    SpecialRules = {
    
    # No special rules for now
    }
    
    Identifiers = {
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    #413:('ESRF','ID14EH2','ADSC Q4'),
    #420:('ESRF','ID14EH4','ADSC Q4R'),
    }
    
    def getRawHeadDict(self, rawHead):
        i1 = 28+rawHead.find("_array_data.header_contents")
        i2 = rawHead.find("_array_data.data", i1)
        i3 = rawHead.find("--CIF-BINARY-FORMAT-SECTION--",i2)+29
        i4 = i3+500
        _lis = [l[2:].strip().split(" ",1) for l in rawHead[i1:i2].splitlines() if l and l[0]=="#"]
        _lis2 =  [l[2:].strip().split(": ",1) for l in rawHead[i3:i4].splitlines() if l and l[0:2]=="X-"]
        RawHeadDict = dict([ val for val in _lis if val[0] in _hKeys])
        RawHeadDict.update(dict([ val for val in _lis2 if "Binary-" in val[0]]))
        RawHeadDict.update({'HEADER_SIZE':i3})
        RawHeadDict.update({'DATE':" ".join(_lis[1])})
        RawHeadDict.update({'MESSAGE':'','TWO_THETA':'0'}) # _FIXME_ Not really here now...
        return RawHeadDict
