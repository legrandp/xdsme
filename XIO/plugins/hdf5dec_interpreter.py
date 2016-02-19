# -*- coding: utf-8 -*-

""" XIO plugin for the Dectris HDF5 image format.
"""

__version__ = "0.0.2"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "17-02-2016"
__copyright__ = "Copyright (c) 2005-2016 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import time

def date_time(time_str):
    "from str return standard str: 'Wed Oct 28 16:42:12 2009'"
    return time.ctime(date_seconds(time_str,"%Y-%m-%dT%T"))

#def _dateseconds(timestr):
#    return time.mktime(time.strptime(timestr))

def date_seconds(time_str):
    "from tupple return seconds"
    try:
        return time.mktime(time.strptime(time_str.split(".")[0],
                                         "%Y-%m-%dT%H:%M:%S"))
    except (ValueError, TypeError), err:
        print "Warning:", err
        print "... Using time.time() instead."
        return time.time()

def get_edge_resolution(pixel_x, width, distance, wavelength):
    "Calculate EdgeResolution."
    from math import sin, atan
    distance=float(distance)
    if abs(distance) > 0.0:
        rad = 0.5 * float(pixel_x) * int(width)
        return float(wavelength)/(2*sin(0.5*atan(rad/distance)))
    else:
        return 0.

def endian(code):
    if code == 'big_endian': return '>'
    else: return '<'
      
FLOAT2 = lambda x: float(x)*1e3

# Creates a dictionary of all keys and values in the NeXus tree

class Interpreter:
    "Dummy class, container for standard Dict and Function."

    HTD = {
    # The adsc Header Translator Dictionary.
    # Potential problems:
    # - There are multiple SIZE1, SIZE2 instances.
    # = The orientation of SIZE1 and SIZE2 is unknown
    #     Not a problem as long as SIZE1 = SIZE2..

    'ExposureTime':(['TIME'],float),
    'BeamX':(['/entry/instrument/detector/beam_center_x'], float),
    'BeamY':(['/entry/instrument/detector/beam_center_y'], float),
    'Distance':(['/entry/instrument/detector/detector_distance'], FLOAT2),
    'Wavelength':(['/entry/instrument/beam/incident_wavelength'], float),
    'PixelX':(['/entry/instrument/detector/x_pixel_size'], FLOAT2),
    'PixelY':(['/entry/instrument/detector/y_pixel_size'], FLOAT2),
    'Width':(['/entry/instrument/detector/detectorSpecific/x_pixels_in_detector'], int),
    'Height':(['/entry/instrument/detector/detectorSpecific/y_pixels_in_detector'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['/entry/sample/goniometer/phi'], lambda x: x[0]),
    'PhiEnd':(['/entry/sample/goniometer/phi_end'], lambda x: x[0]),
    'PhiWidth':(['/entry/sample/goniometer/omega_range_average'], float),
    'EdgeResolution':(['/entry/instrument/detector/x_pixel_size',
       '/entry/instrument/detector/detectorSpecific/x_pixels_in_detector',
       '/entry/instrument/detector/detector_distance',
       '/entry/instrument/beam/incident_wavelength'],
           get_edge_resolution),

    # Added keys from Graeme's convention.
    'TwoTheta':(['TWOTHETA'], float),   # Example missing.
    'SerialNumber':(['/entry/instrument/detector/detector_number'], str),
    'HeaderSize':(['HEADER_BYTES'], int),
    'EndianType':(['BYTE_ORDER'], endian),
    'OscAxis':(['OSC_AXIS'], lambda x: x.lower()),
    # Date and time
    #'DateTime':(['DATE'], date_time),
    'DateStr':(['/entry/instrument/detector/detectorSpecific/data_collection_date'], str),
    'DateSeconds':(['/entry/instrument/detector/detectorSpecific/data_collection_date'], date_seconds),
    'ImageNumber':(['/entry/instrument/detector/detectorSpecific/nimages'], int),
    'SensorThickness':(['/entry/instrument/detector/sensor_thickness'], FLOAT2),
    'OverloadValue':(['/entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff'],int),
    }
    SpecialRules = {
    # No special rules for now
    }

    Identifiers = {
    # ESRF info found at
    # http://www.esrf.eu/UsersAndScience/Experiments/MX/Software/PXSOFT/Denzo
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    '401':('ALS','???','ADSC Q4'),
    }

    def __init__(self):
        self.raw_head_dict = None

    def iterate_children(self, node, nodeDict={}):
        """ iterate over the children of a neXus node """
        if node.type() == self.groupType:
            for kid in node.children():
                nodeDict = self.iterate_children(kid, nodeDict)
        else:
            nodeDict[node.path()] = node.value()
        return nodeDict

    def getRawHeadDict(self, raw_head, groupNode):
        "Intepret the ascii structure of the asdc image header."
        self.groupType = groupNode
        self.raw_head_dict = self.iterate_children(raw_head)
        #print neXus_string_tree
        # First version. Works only with local bioxsoft installation.
        self.raw_head_dict.update({'MESSAGE': '', 'TWOTHETA': '0',
                                   'HEADER_BYTES':0, 'OSC_AXIS': "phi" }) # Example missing
        return self.raw_head_dict
