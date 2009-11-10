# -*- coding: utf-8 -*-

__version__ = "0.4"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "27-10-2009"
__copyright__ = "Copyright (c) 2005-2009 Pierre Legrand"
__license__ = "LGPL"

import time
import struct
import sys

def _extrtime(time_str):
    # from str return tupple
    try:
        return time.strptime(time_str[12:27], "%m%d%H%M%Y.%S")
    except:
        return time.localtime()
    
def _dateseconds(time_str):
    # from tupple return seconds
    try:
        return time.mktime(_extrtime(time_str))
    except:
       return time.time()
        
def _datetime(time_str):
    # from str return standard str: 'Wed Oct 28 16:42:12 2009'
    return time.ctime(_dateseconds(time_str))

getDist = lambda x, y: (x or y)/1e3
divE3 = lambda x: float(x)/1e3
divE5 = lambda x: float(x)/1e5
divE6 = lambda x: float(x)/1e6

def getSerial(comment):
    serial=None
    for l in comment.splitlines():
        if l.lower().count("serial"):
	    try:
	       return l.split()[-1]
	    except:
	       pass
    
def getEdgeResolutionMARCCD(PixelX, Width, Distance, Start_Distance, Wavelength):
            "Calculate EdgeResolution: Graeme's Method"
            from math import sin, atan
            Distance = Distance or Start_Distance
            if Distance > 0.0:
                r = 0.5 * float(PixelX)/1e6 * int(Width)
                return float(Wavelength/1e5)/sin(atan(r/float(Distance)*1e3))
            else:
                return 0.

def getPhiEnd(phiEnd,phiStart,phiWidth):
    if not int(phiEnd):
        return (float(phiStart)+float(phiWidth))/1000.
    else:
        return float(phiEnd)/1000.


headerStructure = [
('tiff_stuff','1024s'),         # 
('header_type','I'),            # flag for header type  (can be used as magic number)
('header_name','16s'),          # header name (MMX)
('header_major_version','I'),   # header_major_version (n.) 
('header_minor_version','I'),   # header_minor_version (.n)
('header_byte_order','I'),      # BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel)
('data_byte_order','I'),        # BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (DEC, Intel)
('header_size','I'),            # in bytes
('frame_type','I'),             # flag for frame type
('magic_number','I'),           # to be used as a flag - usually to indicate new file
('compression_type','I'),       # type of image compression
('compression1','I'),           # compression parameter 1 
('compression2','I'),
('compression3','I'),
('compression4','I'),
('compression5','I'),
('compression6','I'),
('nheaders','I'),               # total number of headers
('pixel_count_x','I'),          # number of pixels in one line
('pixel_count_y','I'),          # number of lines in image
('bytes_per_pixel','I'),        # number of bytes per pixel
('bytes_per_line','I'),         # number of pixels between succesive rows
('bits_per_pixel','I'),         # true depth of data, in bits
('data_type','I'),              # (signed,unsigned,float...)
('saturated_pixel','I'),        # value marks pixel as saturated
('sequence','I'),               # TRUE or FAL
('nimages','I'),                # total number of images - size of each is nfast*(nslow/nimages)
('origin_location','I'),
('detector_orientation','I'),
('view_direction','I'),
('overflow_location','I'),
('overflow_8_bits','I'),
('overflow_16_bits','I'),
('junk3','500s'),
('distance','I'),
('beam_x','I'),
('beam_y','I'),
('integration_time','I'),
('exposure_time','I'),
('readout_time','I'),
('nreads','I'),
('two_theta','i'),
('start_omega','i'),		#;		/* 1000*omega angle */
('start_chi','i'),		#;			/* 1000*chi angle */
('start_kappa','i'),		#;		/* 1000*kappa angle */
('phi_start','i'),
('start_delta','i'),		#;		  /* 1000*delta angle */
('start_gamma','i'),		#;;		  /* 1000*gamma angle */
('start_xtal_to_detector','i'),	#;;   /* 1000*distance in mm (dist in um)*/
('end_twotheta','i'),		#;;		  /* 1000*two_theta angle */
('end_omega','i'),		#;;			  /* 1000*omega angle */
('end_chi','i'),		#;;		  /* 1000*chi angle */
('end_kappa','i'),		#;;			  /* 1000*kappa angle */
('phi_end','i'),
('junk7','12s'),
('axis_code','i'),
('phi_width','i'),
('detector_rotx','i'),
('detector_roty','i'),
('detector_rotz','i'),
('total_dose','i'),
('junk8','12s'),
('detector_type','i'),
('pixel_size_x','i'),
('pixel_size_y','i'),
('mean_bias','i'),
('junk9','124s'),
('lambda','I'),
('junk10','100s'),
('filetitle','128s'),
('filepath','128s'),
('filename','64s'),
('acquire_timestamp','32s'),
('header_timestamp','32s'),
('save_timestamp','32s'),
('file_comment','512s'),
('junk11','1132s')]

class Interpreter:

    HTD = {
    # The marccd Header Translator Dictionary.
    # To add:
    # + SerialNumber (for the rules) or other unique identifier

    'ExposureTime':(['exposure_time'],divE3),
    'BeamX':(['beam_x'], divE3),
    'BeamY':(['beam_y'], divE3),
    'Distance':(['distance','start_xtal_to_detector'], getDist),
    'Wavelength':(['lambda'], divE5),
    'PixelX':(['pixel_size_x'], divE6),
    'PixelY':(['pixel_size_y'], divE6),
    'Width':(['pixel_count_x'], int),
    'Height':(['pixel_count_y'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['phi_start'], divE3),
    'PhiEnd':(['phi_end','phi_start','phi_width'], getPhiEnd),
    'PhiWidth':(['phi_width'], divE3),
    'EdgeResolution':(['pixel_size_x','pixel_count_x','distance','start_xtal_to_detector','lambda'],
        getEdgeResolutionMARCCD),

    # Added keys from Graeme's convention.
    'TwoTheta':(['two_theta'], divE3),
    'SerialNumber':(['file_comment'], getSerial),  # _FIXME_ Don't know where to find SerialNumber!
                                     # Apparently their is now way to know... marcvt
                                     # can't find it, and its nowhere in the marcvt src.
    'EndianType':(['EndianType'], str),
    # 'HeaderSize':(['header_size'], int),  # Don't trust that it says 3072 instead of 4096
    'HeaderSize':(['HEADER_BYTES'], int),
    # Date and time
    #'DateTime':(['DATE'], _datetime),
    'DateStr':(['acquire_timestamp'], _datetime),
    'DateSeconds':(['acquire_timestamp'], _dateseconds),
    }

    Identifiers = {
    # ESRF info found at
    # http://www.esrf.eu/UsersAndScience/Experiments/MX/Software/PXSOFT/Denzo
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    '4':('ESRF','BM14','MarCCD 225'),
    '5':('ESRF','ID23-2','MarCCD 225'),
    '10':('EMBL_HAMBURG','???','MarCCD 225'),
    '12':('SLS','X06SA','MarCCD 225'),
    '21':('EMBL_HAMBURG','???','MarCCD 225'),
    '20':('SSRL4','BL11-1','MarCCD 325'),
    }

    SpecialRules = {
    # No special rules for now
    }

    def getRawHeadDict(self, rawHead, verbose=False):

        # Get the header endian type
        if struct.unpack('<I',rawHead[1052:1056])[0] == 1234:
            EndianType = "<"
        else:
            EndianType = ">"

        headerStructureFmt = EndianType + "".join([a[1] for a in headerStructure])
        headerStructureKeys = [a[0] for a in headerStructure]

        # uncoding from headerStructureFmt
        read_size = struct.calcsize(headerStructureFmt)
        read_unp = list(struct.unpack(headerStructureFmt, rawHead[:read_size]))  

        RawHeadDict = {}
        for l in range(len(headerStructureKeys)):
            _key = headerStructureKeys[l]
            RawHeadDict[_key] = read_unp[l]
            #if 1: #not _key.count("junk"):
            #    print "%s ->%s<-" % (_key, read_unp[l])
        RawHeadDict.update({'MESSAGE':'','HEADER_BYTES':4096,
                                         'EndianType':EndianType})

        return RawHeadDict

if __name__ == "__main__":
    from pprint import pprint
    _image = open(sys.argv[1])
    rawHead = _image.read(9000)
    h = getRawHeadDict(rawHead)

    headKeys =  [a[0] for a in headerStructure]

    l = ["beam_x","beam_y","two_theta","lambda","distance","phi_start",
            "phi_width","phi_end","pixel_size_x","pixel_size_y","pixel_count_x","pixel_count_y"]
    for k in headKeys: #l:
        print "%s:\t%s" % (k,h[k])
