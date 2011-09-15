#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Particularities of XDS - keywords
#    - Only the last keyword declatation count
#       -- except if the keyword accept multiple declatation
#          [case of:  SPOT_RANGE=
#                     EXCLUDE_RESOLUTION_RANGE=
#                     PROFILE_RANGE=]
#    - Multiple declaration on 1 line with no separator
#    - Number of values declared by 1 keyw can go from 1 to 9
#    - Some keyword names are not compatible with python name-space
#         -- the are translated: cf modified_keys
# 
import sys
import os
import re
import commands
import shutil
import fnmatch
from time import time, sleep

__version__ = "0.7.7"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "03-05-2011"
__copyright__ = "Copyright (c) 2006-2011  Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"


# Environemantal variable XDSHOME, if set, defines the place where the xds
# executables will be searched. The parallelized execs (xds_par, xscale_par)
# will be used be defaults.

def get_xdshome():
    if "XDSHOME" in os.environ.keys():
        XDSHOME = os.getenv("XDSHOME")
        if not os.path.isfile(os.path.join(XDSHOME, "xds_par")):
            #if _verbose:
            #    print "WARNING: no 'xds_par' found in $XDSHOME path (%s)." % XDSHOME
            #    print "         Using default $PATH location."
            XDSHOME = ""
    else: XDSHOME = ""
    return XDSHOME

XDSHOME = get_xdshome()
xinp = "XDS.INP"

LP_names = ["COLSPOT.LP","CORRECT.LP","DEFPIX.LP","IDXREF.LP","XSCALE.LP",
            "INIT.LP","INTEGRATE.LP","XYCORR.LP","XDS.INP","XPARM.XDS"]

multiple_keys = ("SPOT_RANGE",
                 "EXCLUDE_RESOLUTION_RANGE",
                 "PROFILE_RANGE",
                 "UNTRUSTED_RECTANGLE")

modified_keys = {
  "DIRECTION_OF_DETECTOR_Y-AXIS": "DIRECTION_OF_DETECTOR_Y_AXIS",
  "DIRECTION_OF_DETECTOR_X-AXIS": "DIRECTION_OF_DETECTOR_X_AXIS",
  "X-RAY_WAVELENGTH": "X_RAY_WAVELENGTH",
  "REFLECTING_RANGE_E.S.D.": "REFLECTING_RANGE_E_S_D_",
  "FRIEDEL'S_LAW": "FRIEDELS_LAW",
  "REFINE(IDXREF)": "REFINE_IDXREF",
  "REFINE(INTEGRATE)": "REFINE_INTEGRATE",
  "REFINE(CORRECT)": "REFINE_CORRECT",
  "BEAM_DIVERGENCE_E.S.D.": "BEAM_DIVERGENCE_E_S_D_",
  "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA": \
  "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA_BETA",
  "SPOT_MAXIMUM-CENTROID": "SPOT_MAXIMUM_CENTROID",
  "UNIT_CELL_A-AXIS": "UNIT_CELL_A_AXIS",
  "UNIT_CELL_B-AXIS": "UNIT_CELL_B_AXIS",
  "UNIT_CELL_C-AXIS": "UNIT_CELL_C_AXIS",
  "REFLECTIONS/CORRECTION_FACTOR": "REFLECTIONS_CORRECTION_FACTOR",
  "MINIMUM_I/SIGMA": "MINIMUM_I/SIGMA",
  "REFINE(IDXREF)": "REFINE_IDXREF",
  "REFINE(INTEGRATE)": "REFINE_INTEGRATE",
  "REFINE(CORRECT)": "REFINE_CORRECT",
  "X-GEO_CORR": "X_GEO_CORR",
  "Y-GEO_CORR": "Y_GEO_CORR",
  "SPOT_MAXIMUM-CENTROID": "SPOT_MAXIMUM_CENTROID" }

illegal_keys = ["SPOT_WIDTH_ALONG_X", "SPOT_WIDTH_ALONG_Y",
  "NUMBER_OF_REFLECTIONS_USED_FOR_REFINEMENT_IN_COLPROF",
  "MINIMUM_SIGNAL_TO_NOISE_FOR_LOCATING_SPOTS",
  "NUMBER_OF_FRAMES_BETWEEN_REFINEMENT_IN_COLPROF","BFRAC","WEAK",
  "MAXIMUM_RANDOM_DEVIATE_OF_INTENSITY","GAIN","RMAX"]

translate_keys = {"REFINE": "REFINE(INTEGRATE)",
                  "RESOLUTION_RANGE_FOR_ACCEPTING_REFLECTIONS": \
                  "INCLUDE_RESOLUTION_RANGE"}

modified_keys_r = {}
for k in modified_keys.keys(): modified_keys_r[modified_keys[k]] = k

xdsinp_base = """
 JOB= ALL
 DATA_RANGE= 1 46
 SPOT_RANGE= 1 9
 SPOT_MAXIMUM-CENTROID= 2.0
 BACKGROUND_RANGE= 1 10
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 5
 STRONG_PIXEL= 5.0
 OSCILLATION_RANGE= 1.0
 STARTING_ANGLE= 74.0
 STARTING_FRAME= 1
 X-RAY_WAVELENGTH= 0.978954
 NAME_TEMPLATE_OF_DATA_FRAMES= ../image_3_???.img
 DETECTOR_DISTANCE= 150.01
 DETECTOR= ADSC   MINIMUM_VALID_PIXEL_VALUE= 1   OVERLOAD= 65000
 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0
 NX= 2048    NY= 2048    QX= 0.1024    QY= 0.1024
 ORGX= 1007.8125    ORGY= 1035.15625
 ROTATION_AXIS= 1.0 0.0 0.0
 INCIDENT_BEAM_DIRECTION= 0.0 0.0 1.0
 FRACTION_OF_POLARIZATION= 0.99
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0
 ! AIR= 0.001
 SPACE_GROUP_NUMBER= 0
 UNIT_CELL_CONSTANTS= 0 0 0 0 0 0
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 7000 30000
 INCLUDE_RESOLUTION_RANGE= 45.0 0.0
 REFINE(INTEGRATE)= BEAM ORIENTATION CELL   
 DELPHI= 8.0
 MAXIMUM_NUMBER_OF_PROCESSORS= 32
 MAXIMUM_NUMBER_OF_JOBS= 1
 RESOLUTION_SHELLS=15.0 8.0 5.0 3.0
 TOTAL_SPINDLE_ROTATION_RANGES=15.0 180.0 15.0
 STARTING_ANGLES_OF_SPINDLE_ROTATION=-95.0 95.0 5.0
 TRUSTED_REGION= 0.0 1.42
 PROFILE_FITTING= TRUE
 STRICT_ABSORPTION_CORRECTION= TRUE
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA= 9
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 9
 REFINE(IDXREF)= BEAM AXIS ORIENTATION CELL
 REFINE(INTEGRATE)= BEAM ORIENTATION CELL
 REFINE(CORRECT)= DISTANCE BEAM AXIS ORIENTATION CELL
 TEST_RESOLUTION_RANGE= 20 4.5
 ! REFERENCE_DATA_SET=
"""

BravaisDico = {"a":"triclinic","m":"monoclinic","o":"orthorhombic",
"t":"tetragonal","h":"hexagonal","c":"cubic","P":"Primitive",
"R":"Rhombohedral","F":"Centered","C":"Centered","I":"Centered"}

Bravais_to_Laue = {
# Each Bravais lattice entry maps a laue_tuple.
#  Each Laue tuple contains
#    (name,first_spg_num, first_spg_name, tuple_of_spg_num)
"aP": (("-1",1,"P1",(1,)),),
"mP": (("2/m",3,"P2",(3,4)),),
"mC": (("2/m",5,"C2",(5,)),),
"mI": (("2/m",5,"C2",(5,)),),
"oP": (("mmm",16,"P222",(16,17,18,19)),),
"oC": (("mmm",21,"C222",(21,20)),),
"oF": (("mmm",22,"F222",(22,)),),
"oI": (("mmm",23,"I222",(23,24)),),
"tP": (("4/m",75,"P4",(75,76,77,78)),
       ("4/mmm",89,"P422",(89,90,91,92,93,94,95,96))),
"tI": (("4/m",79,"I4",(79,80)),
       ("4/mmm",97,"I422",(97,98))),
"hP": (("-3",143,"P3",(143,144,145)),
       ("-31m",149,"P312",(149,151,153)),
       ("-3m1",150,"P321",(150,152,154)),
       ("6/m",168,"P6",(168,169,170,171,172,173)),
       ("6/mmm",177,"P622",(177,178,179,180,181,182))),
"hR": (("-3",146,"R3",(146,)),
       ("-3m1",155,"R32",(155,))),
"cP": (("m-3",195,"P23",(195,198)),
       ("m-3m",207,"P432",(207,208,212,213))),
"cF": (("m-3",196,"F23",(196,)),
       ("m-3m",209,"F432",(209,210))),
"cI": (("m-3",197,"I23",(197,199)),
       ("m-3m",211,"I432",(211,214))),
"Unknown": ((0,0,0,(0,0)),)}

SPGlib = {
# each SGnumber entry maps (SGsymbol1,SGsymbol2,SGorder)
0:("P1","P1",1),
1:("P1","P1",1),3:("P2","P2",2),4:("P21","P2(1)",2),5:("C2","C2",4),
16:("P222","P222",4),17:("P2221","P222(1)",4),
18:("P21212","P2(1)2(1)2",4),19:("P212121","P2(1)2(1)2(1)",4),
21:("C222","C222",8),20:("C2221","C222(1)",8),22:("F222","F222",16),
23:("I222","I222",8),24:("I212121","I2(1)2(1)2(1)",8),
75:("P4","P4",4),76:("P41","P4(1)",4),77:("P42","P4(2)",4),
78:("P43","P4(3)",4),89:("P422","P422",8),90:("P4212","P42(1)2",8),
91:("P4122","P4(1)22",8),92:("P41212","P4(1)2(1)2",8),
93:("P4222","P4(2)22",8),94:("P42212","P4(2)2(1)2",8),
95:("P4322","P4(3)22",8),96:("P43212","P4(3)2(1)2",8),79:("I4","I4",8),
80:("I41","I4(1)",8),97:("I422","I422",16),98:("I4122","I4(1)22",16),
143:("P3","P3",3),144:("P31","P3(1)",3),145:("P32","P3(2)",3),
149:("P312","P312",6),150:("P321","P321",6),151:("P3112","P3(1)12",6),
152:("P3121","P3(1)21",6),153:("P3212","P3(2)12",6),
154:("P3221","P3(2)21",6),168:("P6","P6",6),169:("P61","P6(1)",6),
170:("P65","P6(5)",6),171:("P62","P6(2)",6),172:("P64","P6(4)",6),
173:("P63","P6(3)",6),177:("P622","P622",12),178:("P6122","P6(1)22",12),
179:("P6522","P6(5)22",12),180:("P6222","P6(2)22",12),
181:("P6422","P6(4)22",12),182:("P6322","P6(3)22",12),
146:("R3","R3",9),155:("R32","R32",18),195:("P23","P23",12),
198:("P213","P2(1)3",12),207:("P432","P432",24),208:("P4232","P4(2)32",24),
212:("P4332","P4(3)32",24),213:("P4132","P4(1)32",24),196:("F23","F23",48),
209:("F432","F432",96),210:("F4132","F4(1)32",96),197:("I23","I23",24),
199:("I213","I2(1)3",24),211:("I432","I432",48),214:("I4132","I4(1)32",48)}


PG2SP = {"1":(1,),"2":(3,4,5),
"222":(16,17,18,19,20,21,22,23,24),
"4":(75,76,77,78,79,80),
"422":(89,90,91,92,93,94,95,96,97,98),
"3":(143,144,145,146),
"312":(149,151,153),
"321":(150,152,154,155),
"6":(168,169,170,171,172,173),
"622":(177,178,179,180,181,182),
"23":(195,196,197,198,199),
"432":(207,208,209,210,211,212,213,214)}

def get_BravaisToSpgs():
    Bravais_to_spg = {}
    for br in Bravais_to_Laue:
        for lauespgs in Bravais_to_Laue[br]:
            if br in Bravais_to_spg:
                Bravais_to_spg[br] += lauespgs[3]
            else:
                Bravais_to_spg[br] = lauespgs[3]
    return Bravais_to_spg

# global dictionary to keep trac of the XDS.INP key order appearance
xdsKeyOrder = []

class Lattice:

    def __init__(self, cell, Bravais_type="Unknown", symmetry=None,
                      dmin=0.5, friedels_law=1):
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.

        # Minimal Bragg spacing = Maximal resolution
        # Non zerro value to avoid zerro division errors.
        self.dmin = dmin

        # Friedel'law default = true equivalent to anomal=false
        self.friedels_law = friedels_law

        # Indexing default parameters
        self.fit = None # indexing geometrical fit quality, Float
        self.character = None # Lattice character, Int (from 1 to 44,
                              # see International tables)
        self.reindexing = None # tupple of len 12

        # Scaling results
        self.rsym = None
        self.rmeas = None
        self.rmea3 = None
        self.compl = None
        self.comp3 = None
        self.isig = None

        if type(cell) == type(""):
            cell = map(float, cell.strip().split())
        if len(cell) == 6:
            self.cell = map(float,tuple(cell))
            (self.a, self.b, self.c,
             self.alpha, self.beta, self.gamma) = tuple(cell)
        elif len(cell) == 3:
            self.cell = map(float,tuple(cell))
            self.cell = cell[0], cell[1], cell[2], 90., 90., 90.
            (self.a, self.b, self.c) = tuple(cell)
        else:
            print cell, len(cell)
            print "ERROR: Number of argument incorrect for Lattice class."
            sys.exit()

        if Bravais_type in Bravais_to_Laue.keys():
            self.Bravais_type = Bravais_type

        if Bravais_type == "Unknown" and symmetry:
            self.symmetry_num = int(symmetry)
            for brav in Bravais_to_Laue:
                for Laue_spgs in Bravais_to_Laue[brav]:
                    if self.symmetry_num in Laue_spgs[3]:
                        self.Bravais_type = brav
                        break
        else: 
            if symmetry:
            # Verify that the symmetry number is 
            # compatible with the Bravais type"
               symnums = []
               for syms in Bravais_to_Laue[self.Bravais_type]:
                    symnums.extend(syms[-1])
               if symmetry not in symnums:
                   print "ERROR: Given symmetry number",
                   print " imcompatible with Bravais type."
                   sys.exit()
               else: self.symmetry_num = int(symmetry)
            else:
            # Set lowest symmetry number for the given Bravais type
                self.symmetry_num = Bravais_to_Laue[self.Bravais_type][0][1]
        #
        self.symmetry_str1 = SPGlib[self.symmetry_num][0]
        self.symmetry_str2 = SPGlib[self.symmetry_num][1]
        self.multiplicity = SPGlib[self.symmetry_num][2]

    def prt(self, fmt=6*"%7.1f"):
        return fmt % tuple(self.cell)

    def __str__(self):
        #self.cell = map(float,tuple(self.cell))
        return self.prt()

    def volum(self):
        from math import cos,pi
        d2r = pi/180
        a, b, c, al, be, ga = self.cell
        cosa, cosb, cosg = cos(al*d2r), cos(be*d2r), cos(ga*d2r)
        return a*b*c*(1-cosa**2-cosb**2-cosg**2+2*cosa*cosb*cosg)**0.5

    def idealize(self):
        _latt = self.Bravais_type[0]
        a, b, c, alpha, beta, gamma = self.cell
        if _latt == "m":
            self.alpha = 90.
            self.gamma = 90.
        if _latt == "o":
            self.alpha = 90.
            self.beta = 90.
            self.gamma = 90.
        if _latt == "t":
            a = (a+b)/2.
            self.a = a
            self.b = a
            self.alpha = 90.
            self.beta  = 90.
            self.gamma = 90.
        if _latt == "h":
            a = (a+b)/2.
            self.a = a
            self.b = a
            self.alpha = 90.
            self.beta  = 90.
            self.gamma = 120.
        if _latt == "c":
            a = (a+b+c)/3.
            self.a = a
            self.b = a
            self.c = a
            self.alpha = 90.
            self.beta  = 90.
            self.gamma = 90.

        self.cell = (self.a, self.b, self.c,
                     self.alpha, self.beta, self.gamma)


class Param:
    """A simple class to handle the parameters that fills templates"""
    def __init__(self, obj=None):
        """Constructor for the Param classes from file or string."""
        #if type(obj) == file:# needs python version >= 2.2
        from types import FileType
        if type(obj) == FileType:
            exec obj.read() in self.__dict__
            obj.close()
        if type(obj) == type(""):
            exec obj in self.__dict__

    def __getitem__(self, a):
        return self.__dict__[a]

    def __setitem__(self, a, b):
        self.__dict__[a] = b 

    def __delitem__(self, a):
        del self.__dict__[a]

    def mix(self, _par):
        """Update the current param with either 
           custom Param instance, dict or executable string."""
#        if type(_par) == dict:  # needs python version >= 2.2
        if type(_par) == type({}):
            self.__dict__.update(_par)
        elif type(_par) == type(""):
            exec _par in self.__dict__
        else:
            # if type(_par) == types.InstanceType:
            self.__dict__.update(_par.__dict__)

    def add(self, _par):
        """Complete and/or replace the current param with custom dict."""
#        if type(_par) == dict: # needs python version >= 2.2
        if type(_par) == type({}):
           self.__dict__.update(_par)

    def keys(self):
        rkeys = []
        for k  in self.__dict__.keys():
            if k not in ("__builtins__","__doc__","__module__"):
                rkeys.append(k)
        return rkeys

    def has_key(self, key):
        if key in self.keys(): return 1
        else: return 0

    def intersect(self, par2):
        return filter(self.has_key, par2.keys())


class XParam(Param):

    def add(self, par_dict):
        """Complete and/or replace the current param with custom dict."""
        for k in par_dict.keys():
            if k in multiple_keys:
                if not k in self.keys():
                    par_dict[k] = [par_dict[k]]
                else:
                    self[k].append(par_dict[k])
                    par_dict[k] = self[k]
        self.__dict__.update(par_dict)

    def copy(self):
        return xdsInp2Param(inp_str="%s" % self)

    def xds_parse(self):
        """ Parser to transforme, if possible, string variables to:
        numerical variable, or tupple of numerical or string"""
        def _parse(keyStr):
            fmt_str = "'%s',"
            fmt_num = "%s,"
            valStr = ""
            subValues = keyStr.split()
            for s in subValues:
                try:
                    float(s)
                    fmtVal = fmt_num
                except ValueError:
                    fmtVal = fmt_str
                valStr += fmtVal % (s)
            return valStr[:-1]

        for k in  self.keys():
            _to_exec = "%s= " % k
            if k in multiple_keys:
                if type(self[k]) != type(""):
                    for subv in self[k]:
                        _to_exec += "["+ _parse(subv) +"],"
                else: _to_exec += _parse(self[k])+","
            else:
                _to_exec += _parse(self[k])
            exec _to_exec in self.__dict__

    def __repr__(self):
        def _prt(self,l,_to_repr):
            line = ""
            for k in l:
                kp = k
                if k in modified_keys_r.keys(): kp = modified_keys_r[k]
                if k in self.keys():
                    if k in multiple_keys:
                        for subval in self[k]:
                            line += "%s= %s\n " % (kp, toString(subval))
                        line = line[:-1]
                    else:
                        line += "%s= %s    " % (kp, toString(self[k]))
                    _to_repr.remove(k)
            return line, _to_repr

        _repr = "!=== File generated by XUPY v%s ===!\n\n" % __version__
        _to_repr = self.keys()
        if xdsKeyOrder:
            for l1 in xdsKeyOrder:
                _line, _to_repr = _prt(self,l1,_to_repr)
                _repr += " %s\n" % _line[:-1]
        for not_print in ("__builtins__","__module__","__doc__"):
            if _to_repr.count(not_print): _to_repr.remove(not_print)
        _repr += "\n!=== Added Keywords ===!\n\n"
        # _to_repr[:], a copy of _to_repr, is used to avoid interference
        # with the _to_repr.remove(k) instruction in the _prt(l) function
        for k2 in _to_repr[:]:
            # To avoid printing of internal parameters
            if k2[0] != "_":
                _line, _to_repr = _prt(self,(k2,),_to_repr)
                _repr += " %s\n" % _line[:-1]
        return _repr    

class DataCollectInfo:
    """ A simple class to handle Data Collection Image informations."""

    def __init__(self, init):
        """ ImagesInfo can be contruct from either a string template or a
            tupple of 4 building elements (dir, prefix, nDigits, suffix).
        """
        self.image_numbers = [] # set by lookup_image_numbers()
        self.image_ranges = []  # set by lookup_image_ranges()

        if type(init) == type(""):
            self.directory, init = os.path.split(init)
            regexp = re.compile(r"(.*)_(.*)\.(.*)")
            match = regexp.match(init)
            self.prefix = match.group(1)
            numeric = match.group(2)
            self.suffix = match.group(3)
            self.nDigits = len(numeric)

        elif type(init) == type(()) and len(init) == 4:
            self.directory = init[0]
            self.prefix = init[1]
            self.nDigits = int(init[2])
            self.suffix = init[3]
        else:
            raise TypeError, "Unexpected ImagesInfo contructor argument."

        if not self.directory:
            self.directory = "."
        self.format = self.prefix+"_%0"+str(self.nDigits)+"d."+self.suffix
        expression = self.prefix +  "_([0-9]*)\." + self.suffix
        self.regexp = re.compile(expression)
        self.regexpCompress = re.compile(expression + "[\.gz|\.z|\.Z|\.bz2]*")

    def getNumber(self, imageName):
        d, imageName = os.path.split(imageName)
        match = self.regexp.match(imageName)
        return int(match.group(2))

    def getXDSTemplate(self):
        xdst = self.directory+"/"+self.prefix+"_"+\
               self.nDigits*"?"+"."+self.suffix
        if len(xdst) > 50:
            print ">>> Warning NAME_TEMPLATE_OF_DATA_FRAMES has more than 50",
            print " characters! XDS will stop."
            print ">>> Lengthy path names should be abbreviated by a symbolic",
            print " link for frames directory."
        return xdst

    def getMosflmTemplate(self):
        return  self.prefix +  "_" + self.nDigits*"#" + "." + self.suffix  

    def lookup_image_numbers(self):
        """Return a list of matching image number. Removes duplicate numbers.
        """
        images_num = []
        files = os.listdir(self.directory)
        files.sort()
        prev = 0
        for f in files:
            match = self.regexpCompress.match(f)
            if match:
                n = int(match.group(1))
                if n != prev:
                    images_num.append(int(match.group(1)))
                    prev = n
        self.image_numbers = images_num
        return images_num

    def lookup_image_ranges(self):
        """Return a range list of consecutive image number.
           For example: [[1,3],[91,93]] or [[1,90]]
        """
        seqf = []
        if self.lookup_image_numbers():
            prev, i = self.image_numbers[0]-1, self.image_numbers[0]
            for n in self.image_numbers:
                if n != prev+1:
                    seqf.append([i,prev])
                    i = n
                prev  = n
            seqf.append([i,prev])
        self.image_ranges = seqf
        return seqf

    def get_range(self, minf=None, maxf=None):
        if self.image_ranges:
            min_c, max_c = self.image_ranges[0][0], self.image_ranges[-1][-1]
            if minf: min_c = max(minf,min_c)
            if maxf: max_c = min(maxf,max_c)
            return [min_c, max_c]
        else:
            return []

    def getClosestImage(self, target):
        """Return closest existing image number to target.
        """
        target = int(target)
        if not self.image_numbers:
            self.lookup_image_numbers()
        diff = map(lambda x,y: abs(x-y), self.image_numbers, 
                               len(self.image_numbers)*(target,))
        return self.image_numbers[diff.index(min(diff))]

def toString(obj):
    """ Transforme all variables to a string format"""
    if type(obj) == type(()) or type(obj) == type([]): 
        return " ".join(map(str,obj))
    else: return str(obj)

def xdsInp2Param(inp_name="XDS.INP", inp_str=None):
    """ Translate an XDS input file to a Paramter object.
    Return the Paramter object and a list of the key order appearence"""
    global xdsKeyOrder
    if not inp_str:
        inp_str = opReadCl(inp_name)
    xdsKeyOrder = []
    # Recognition template for XDS keywords
    key_re  = r"([A-Z0-9_.\-\\'\(\)\/]+=)"
    newPar = XParam()
    allKey = []
    for line in inp_str.splitlines():
        aline = re.split("(!)",line.strip())[0]
        raw_keys = re.split(key_re,aline)[1:]
        try: nkeys = len(raw_keys)/2
        except:
            print "Can't process line:\n", aline
            print "Wrong keys are: ",raw_keys,"\nSTOP!"            
            return None
        if nkeys:
            _kOrd = []
            for i in range(0,nkeys*2,2):
                _d0 = {}
                newkey = raw_keys[i][:-1]
                if newkey in modified_keys.keys():
                    newkey = modified_keys[newkey]
                if newkey not in allKey:
                    _kOrd.append(newkey)
                    allKey.append(newkey)
                    if newkey in multiple_keys:
                          _d0[newkey] = [raw_keys[i+1].strip()]
                    else: _d0[newkey] = raw_keys[i+1].strip()
                    newPar.mix(_d0)
                elif newkey in multiple_keys:
                    values = []
                    previous = getattr(newPar,newkey)
                    if type(previous) != type(""):
                        for v in previous: values.append(v)
                    else: values.append(previous)
                    values.append(raw_keys[i+1].strip())
                    _d0[newkey] = values
                    newPar.mix(_d0)
            if _kOrd: xdsKeyOrder.append(_kOrd)
    newPar.xds_parse()
    return newPar


def latest(check_names=LP_names):
    rexp = re.compile(r"(.*)\.(\d\d\d)\Z")
    ns = [int(rexp.match(n).group(2)) for n in os.listdir(".") if \
          (rexp.match(n) and rexp.match(n).group(1) in check_names)]
    if not ns: return 0 
    else: return max(ns)

def saveLastVersion(file_names, suffix=""):
    import filecmp
    last_vnum = latest(file_names)
    if not suffix:
        for name in file_names:
            last_lnum = latest((name,))
            compare = name, name+".%03d" % last_lnum
            if os.path.isfile(compare[0]) and os.path.isfile(compare[1]):
                if not filecmp.cmp(name,name+".%03d" % last_lnum):
                     shutil.copyfile(name ,name+".%03d" % (last_vnum+1))
            elif os.path.isfile(compare[0]):
                if os.path.getsize(compare[0]) != 0:
                    shutil.copyfile(name ,name+".%03d" % (last_vnum+1))
    else:
        for name in file_names:
            if os.path.getsize(name) != 0:
                shutil.copyfile(name ,name+str(suffix))

def exec_prog(prog_name, stdinp=None, stdout= None, stderr=None):
    if not stdout : stdout = " "
    else: stdout = " > " + stdout
    if not stdinp : stdinp = " "
    else: stdinp = " < " + stdinp
    if not stderr : stderr = " "
    else: stderr = " 2>" + stderr
    os.system(prog_name+stdinp+stdout+stderr ) # use popen instead ?

def run_xds_thread(arguments):
    tpar,tinp,tout,tdir,tsave = arguments
    tmp_par = tpar.copy()
    return run_xds(tmp_par, inp_f=tinp, out_f=tout, directory=tdir, save=tsave)

def new_range(r,n):
    l = r[1]-r[0]
    new_r = [r[0]-1]
    for a in range(n): new_r.append(r[0]+int((a+1)*l/float(n)))
    return new_r 

def taskCallback(message):
    global FinishedThread
    FinishedThread += 1
    print " <==   Threads", message

def waitForAllThreadsEnd(NumberOfThreadToWait):
    global FinishedThread
    while FinishedThread != NumberOfThreadToWait:
        sleep(0.1)

def run_multi_integrate(xpar_init,inp_f=None,main_directory=None,
                            nThreads=2,init=0):    
    from ThreadPool import ThreadPool
    # Copy 
    xpar = xpar_init.copy()
    global FinishedThread
    FinishedThread = 0
    pool = ThreadPool(nThreads)
    if init:
       xpar.JOB = "INIT","DEFPIX","INTEGRATE","CORRECT"
       files_to_copy = "XPARM.XDS","X-CORRECTIONS.pck","Y-CORRECTIONS.pck"
    else:
       xpar.JOB = "DEFPIX","INTEGRATE","CORRECT"
       files_to_copy = "XPARM.XDS","BKGINIT.pck","BLANK.pck",\
                       "GAIN.pck","X-CORRECTIONS.pck","Y-CORRECTIONS.pck"
    if not main_directory: main_directory = os.getcwd()
    if os.path.isdir(main_directory):
         os.chdir(main_directory)
    else:
         print "STOP! Directory not found:",directory
         sys.exit()
    range_list = new_range(xpar.DATA_RANGE, nThreads)

    _templ = xpar.NAME_TEMPLATE_OF_DATA_FRAMES
    if type(_templ) ==  type("") and _templ[0] != "/":
        xpar.NAME_TEMPLATE_OF_DATA_FRAMES = "../" + \
                   xpar.NAME_TEMPLATE_OF_DATA_FRAMES

    if type(_templ) == type(()) and _templ[0][0] != "/":
        xpar.NAME_TEMPLATE_OF_DATA_FRAMES = "../" + \
                   xpar.NAME_TEMPLATE_OF_DATA_FRAMES[0], \
                   xpar.NAME_TEMPLATE_OF_DATA_FRAMES[1]

    hklfiles = []
    startTime = time()
    print "\n"
    for NTh in range(nThreads):
        #newdir = os.path.join(main_directory,"integrate_batch_%d" % (NTh+1))
        newdir = "integrate_batch_%d" % (NTh+1)
        if not os.path.exists(newdir): os.mkdir(newdir)
        for file in files_to_copy:
            if os.path.isfile(file):
                shutil.copyfile(file,os.path.join(newdir,file))
            else:
                print "STOP: Can't find file",file
                sys.exit()
        os.chdir(newdir)
        xpar.DATA_RANGE = range_list[NTh]+1,range_list[NTh+1]
        xpar.SPOT_RANGE = (range_list[NTh]+1,range_list[NTh]+4),
        args = (xpar,inp_f,"xds.out",os.getcwd(),0)
        pool.queueTask(run_xds_thread,args,taskCallback)
        # To avoid chdir before the xds process is started...
        print "  ==>  Threads XDS started for integration of images %4d to %4d"\
                % (xpar.DATA_RANGE[0],xpar.DATA_RANGE[1])
        sleep(2.0)
        os.chdir("..")
        hklfiles.append(os.path.join(newdir,"XDS_ASCII.HKL"))
    pool.joinAll()
    print "\n"
    while FinishedThread != nThreads:
        sleep(0.1)
    #read cell
    endTime = time()
    print "\n    Integration time: %.1f seconds" % (endTime - startTime)
    H = read_xdsascii_head(hklfiles[0])
    if H["friedels_law"] == "TRUE": H["friedels_law"] = 1
    elif H["friedels_law"] == "FALSE": H["friedels_law"] = 0
    dmin = get_maxResolution(os.path.join(newdir,"INTEGRATE.LP"))
    xlatt = Lattice(H["cell"], "Unknown", symmetry=H["sym"],
                      dmin=dmin, friedels_law=H["friedels_law"])
    run_xscale(hklfiles, "batch_merge.hkl", xlatt, save=1, out_f="xscale.out")

def run_xds(new_par, inp_f=xinp, out_f=None, directory=None, save=1):
    if directory:
        if not os.path.exists(directory):
            try: os.mkdir(directory)
            except:
                print "STOP! Can't creat xds working directory:",directory
                sys.exit()
        if os.path.isdir(directory): os.chdir(directory)
    #print "Working on images %4d to %4d in %s" % (r[0],r[1],directory)

    if inp_f:
        # Verify the presence of the specified XDS.INP file
        if not os.path.isfile(inp_f):
            print ">>> ERROR: Can't find file "+inp_f+" !"
            sys.exit()
        # If an old XDS.INP file exist, try to backup it.
        if inp_f and os.path.isfile(xinp):
            try: shutil.copyfile(xinp, xinp + "_backup")
            except: 
                print ">>> ERROR: Can't save old "+xinp+" file !"
                sys.exit()
        # Try to save the trace of the first XDS.INP template used
        if not os.path.isfile(xinp+"_init"):
            shutil.copyfile(inp_f, xinp+"_init")
        xpar = xdsInp2Param(inp_name=inp_f)
    else:
        xpar = XParam()

    xpar.mix(new_par)
    opWriteCl(xinp, "%s" % xpar)
    exec_prog(os.path.join(XDSHOME,"xds_par"), stdout=out_f, stderr="xupy.stderr")
    r = xpar.DATA_RANGE
    if save: saveLastVersion(LP_names)
    return "XDS finished for integration of images %4d to %4d" % (r[0],r[1])

def getProfilRefPar(infile="INTEGRATE.LP"): 
    p_lp = opReadCl(infile)
    if len(p_lp) <= 1500:
        print "\nERROR! Uncompleted 'INTEGRATE' step."
        sys.exit()
    #if init:
    #    st1 = p_lp.index("BASED ON SPOT PROFILE PARAMETERS")+32
    #    fit_dat_1 = p_lp[st1:st1+163].replace("DEGREES","")
    else:
        st1 = p_lp.index("* SUGGESTED VALUES FOR INPUT PARAMETERS *")
        fit_dat_1 = p_lp[st1+46:st1+165]
    return xdsInp2Param(inp_str=fit_dat_1)

def gxparm2xpar(Dir):
    shutil.copyfile(os.path.join(Dir,"GXPARM.XDS"),
                    os.path.join(Dir,"XPARM.XDS"))

def get_num(seq,num):
    L = []
    for i in num: L.append(seq[i])
    return L

facteur_names = ["Compl","ComplL","Compl3","Unique","Total","Total3",
         "Compar","Rsym","RsymL","Rsym3","Compa3","Rmeas","Rmea3","I/sig",
         "I/sigL","Misfit","Absent","AbsIav"]

def resum_scaling(lpf="CORRECT.LP", ios_threshold=2.0):
    """ Extract Scaling statistics from XSCALE.LP or CORRECT.LP."""

    lp = opReadCl(lpf)
    if len(lp) < 2000: return None
    s = XParam()

    file_type = lpf.split("/")[-1][:6]
    if file_type == "CORREC": correct = 1
    elif file_type == "XSCALE": correct = 0

    try:
        spa = lp.index("CORRECTION PARAMETERS FOR THE STANDARD ERROR OF")
        spb = lp.index(" ***********", spa+100)
        AB6 = lp[spa:spb].split()
        if correct:
            #print AB6[-6:]
            #sp1 = lp.index("      b          ISa")
            #sp2 = lp.index("  INTEGRATE.HKL   ", sp1)
            s.K1s, s.K2s = map(float, AB6[-3:-1])
        else:
            #print AB6[-50:-6]
            sp1 = lp.index("ISa0   INPUT DATA SET")
            #sp2 = lp.index("  XDS_ASCII.HKL   ", sp1)
        s.IoverSigmaAsympt =  1/((s.K1s*(s.K2s+0.0004))**0.5)
    except:
        s.IoverSigmaAsympt =  -99.9
    st2  = lp.index("  STATISTICS OF S")
    s.LowestReso = 100
    slowr = lp.index("INCLUDE_RESOLUTION_RANGE=") + 26
    s.LowestReso, s.HighestReso = lp[slowr:slowr+21].split()
    if correct:
        st3 = lp.index("NUMBER OF REJECTED MISFITS ",st2)
        st6 = lp.index("NUMBER OF UNIQUE ACCEPTED REFLECTIONS " ,st2)
        stat_g = lp[st3:st6+58].split()
        s.misfit, s.absent = get_num(stat_g,(4,10))
    else:
        st3 = lp.index("REFLECTIONS REJECTED")
        st6 = lp.index("REFLECTIONS ON OUTPUT FILE")
        s.misfit, tmp = get_num(lp[st3-24:st6].split(),(0,-1))
        s.absent = "0"
    #
    st10 = lp.index("    Corr\n\n" ,st2)+10
    st11 = lp.index("\n\n",st10)
    #st12 = lp.index("NOISE >=  3.0" ,st11)
    #st13 = lp.index("\n\n\n",st12)
    st14 = lp.index("WILSON LINE ",st11)
    #
    if correct:
        stat_tg = lp[st10:st11].splitlines()
        #stat_tg = lp[st10:st11].splitlines()[4:-2]
        #stat_tg3 = lp[st12:st13].splitlines()[4:-2]
    else:
        #st10x = lp.index("NOISE >= -2.0" ,st2)
        st10x = lp.index("= STATISTICS " ,st2)
        #st12x = lp.index("NOISE >=  4.0" ,st11)
        #stat_tg = lp[st10:st10x].splitlines()[4:-2]
        stat_tg = lp[st10:st10x].splitlines()[4:-3]
        #stat_tg3 = lp[st12:st12x].splitlines()[4:-2]

    stat_wilson = lp[st14+30:st14+75].split()
    #
    TG, TG3 = [], []
    for l in stat_tg: TG.append(l.split())
    #for l in stat_tg3: TG3.append(l.split())

    s.wilson_b, s.wilson_corr = stat_wilson[3], stat_wilson[5]
    s.reso, s.resoL =       TG[-2][0], TG[-3][0]
    s.compar, s.comparL =   TG[-1][7], TG[-2][7]
    s.total =               TG[-1][1]
    s.compl, s.complL =     TG[-1][4], TG[-2][4]
    s.rsym, s.rsymL =       TG[-1][5], TG[-2][5]
    s.rmeas, s.rmeasL =     TG[-1][9], TG[-2][9]
    s.isig, s.isigL =       TG[-1][8], TG[-2][8]
    s.anoNum, s.anoNumL  =  TG[-1][-1], TG[-2][-1]
    s.anoSig, s.anoSigL  =  TG[-1][-2], TG[-2][-2]
    s.anoCorr,s.anoCorrL =  TG[-1][-3], TG[-2][-3]
    s.unique =              TG[-1][2]
    #s.rsym3, s.rsym3L =     TG3[-1][5], TG3[-2][5]
    #s.rmeas3, s.rmeas3L =   TG3[-1][9], TG3[-2][9]
    #s.total3 =              TG3[-1][1]
    #s.compl3, s.compl3L =   TG3[-1][4], TG3[-2][4]
    #s.compar3, s.compar3L = TG3[-1][7], TG3[-2][7]
    if correct:
        stt = lp.index("   STANDARD ERROR OF REFLECTION INTENSITIES")
        stt = lp.index("--------\n", stt)
        statline = lp[stt+9:stt+86].split()
        s.LowestReso, s.HighestReso, s.iosig, s.chi2, s.rsym  = statline[0:5]
    #
    for k in s.keys():
       if type(s[k]) == str: s[k] = float(s[k].replace("%",""))
    #
    if float(s.absent):
        stnabs = lp.index("AVERAGE INTENSITY FOR", st2)
        s.AbsNum = int(lp[stnabs:stnabs+60].split()[3])
        st7 = lp.index("SYSTEMATICALLY ABSENT",st2)
        s.AbsIav = float(lp[st7+27:st7+32].strip())
    else:
        s.AbsIav = 0
        s.AbsNum = 0
    #
    reso, rsym, ios = [],[],[]
    for i in TG[:-1]:
        reso.append(float(i[0]))
        rsym.append(float(i[9][:-1]))
        ios.append(float(i[8]))
    #
    if correct:
        stcs = lp.index("SELECTED SPACE GROUP AND UNIT")
    else:
        stcs = 0
    stcell = lp.index("UNIT_CELL_CONSTANTS=", stcs)+20
    s.cell = lp[stcell:stcell+51].rstrip()
    stspg = lp.index("SPACE_GROUP_NUMBER=", stcs)+19
    s.spg_num = int(lp[stspg:stspg+5])
    s.spg_sym = SPGlib[s.spg_num][1]
    ind = 0
    for i in range(len(reso)):
        if ios[i] >= ios_threshold: ind = i
    if 0 < ind < len(reso)-1:
        M = (ios[ind] - ios_threshold)/(ios[ind] - ios[ind+1])
        s.dmin = reso[ind] + (reso[ind] - reso[ind+1]) * M
    elif ind == 0: s.dmin = reso[-1]
    else: s.dmin = reso[ind]
    #print s.dmin
    s.multiplicity = s.total/s.unique
    return s

def unpack_latticefit(_str):
    ss = _str.split()
    latt = Lattice((map(float,ss[3:9])),ss[1])
    latt.fit = float(ss[2])
    latt.character = int(ss[0])
    latt.reindexing = tuple(map(int,ss[9:]))
    return latt

def resum_idxref(idxref="IDXREF.LP"):
    list_latticesFit = []
    i_lp = opReadCl(idxref)
    st1 = i_lp.index("LATTICE-  BRAVAIS-")
    idxref_latticesFit = i_lp[st1+172:st1+5012]
    return map(unpack_latticefit, idxref_latticesFit.splitlines())

def get_xparm_cell(xp_name="XPARM.XDS"):
    return tuple(map(float,\
                opReadCl(xp_name).splitlines()[7].split()[1:]))

def select_lattices(limit = 100, idxref="IDXREF.LP"):
            selected = []
            selection = Param()
            for _latt in resum_idxref(idxref):
                 #if _latt.fit <= limit and _latt.fit != 0.0:
                 if _latt.fit <= limit:
                     selected.append(_latt)
            return selected

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def opWriteCl(filename, _str):
    f = open(filename,"w")
    r = f.write(_str)
    f.close()

_add = lambda x,y: x+y

def mean(seq):
    return reduce(_add, [s[0] for s in seq])/float(len(seq))

def wMean(seq):
    sum = reduce(_add, [s[0]*s[1] for s in seq])
    sumw = reduce(_add, [s[1] for s in seq])
    return sum/float(sumw)

def standardDeviation(seq):
    m = mean(seq)
    sum2 = reduce(_add, [ (n[0]-m)**2 for n in seq ])
    return (sum2/(len(seq)-1))**0.5

def wStandardDeviation(seq):
    m = wMean(seq)
    sum2 = reduce(_add, [ n[1]*(n[0]-m)**2 for n in seq ])
    sumw = reduce(_add, [ n[1] for n in seq ])
    return (sum2/sumw)**0.5

def read_xdsascii_head(file_name_in):
    head = {}
    head["cell"] = 0,0,0,90,90,90
    head["sym"] = 0
    head["inputfile_name"] = ""
    head["inputfile_type"] = ""
    head["merge"] = ""
    head["friedels_law"] = ""
    head["wavelength"] = 0
    head["template_name"] = ""
    if not os.path.exists(file_name_in):
        print "ERROR! Can't find file %s.\nSTOP.\n" % (file_name_in)
        sys.exit()
    raw = open(file_name_in)
    line = raw.readline()
    while line[0] == "!":
        if line.count("NAME and FORMAT") :
            line = raw.readline()
            head["inputfile_name"] = line.split()[2]
            head["inputfile_type"] = line.split()[3]
        if line.count("COMPRISES THE FOLLOWING SCALED INPUT FILES:"):
            line = raw.readline()
            head["inputfile_name"] = line.split("INPUT_FILE=")[1].strip()
            head["inputfile_type"] = "XDS_ASCII"
        #
        if line.count("UNIT_CELL_CONSTANTS="):
            head["cell"] = line[line.index("=")+1:-1].strip()
        if line.count("SPACE_GROUP_NUMBER="):
            head["sym"] = line[line.index("=")+1:-1].strip()
        if line.count("X-RAY_WAVELENGTH="):
            iw = line.index("WAVELENGTH=")+11
            head["wavelength"] = float(line[iw:iw+10].strip())
        if line.count("MERGE=") == 1:
            head["merge"] = line[line.index("MERGE=")+6:-1].split()[0].strip()
        if line.count("FRIEDEL'S_LAW="):
            head["friedels_law"] = \
                     line[line.index("FRIEDEL'S_LAW=")+14:-1].strip()
        if line.count("NAME_TEMPLATE_OF_DATA_FRAMES="):
            head["template_name"] = line[line.index("=")+1:].strip().split("??")[0].split("/")[-1]
            if head["template_name"][-1] == "_":
                head["template_name"] == head["template_name"][:-1]
        line = raw.readline()
    return head

def get_resmax_limit():
    lp = ""
    if os.path.exists("DEFPIX.LP"):
        try: lp = opReadCl("DEFPIX.LP")
        except: pass
    elif os.path.exists("INTEGRATE.LP"):
        try: lp = opReadCl("INTEGRATE.LP")
        except: pass
    else: return 1.5
    if len(lp) <= 1300: return 1.5
    st = lp.index("RESOLUTION RANGE RECORDED BY DETECTOR")
    return float(lp[st:st+80].splitlines()[0].split()[-1])


def res_bin(dmin,nbin):
   return [1./(((1./dmin)**3*bin/nbin)**(1./3)) for bin in range(1,nbin+1)]

def get_maxResolution(infile="INTEGRATE.LP"):
    p_lp = opReadCl(infile)[:1115]
    to_find = "RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)"
    resol = p_lp[p_lp.index(to_find)+len(to_find):].split()[1]
    return float(resol)

def write_xscale_resum(s,r,friedels_law,Dir=None):
    if not Dir: Dir = os.getcwd()
    resum = open(os.path.join(Dir,"xscale_resum.txt"),"w")
    print "\n\n\n\t\t\tScaling Statistics\n\n"
    #
    def print_o3(t,v1,v2):
        txt = "\t   %-22s %14s    %-s" % (t,v1,v2)
        print txt
        resum.write(txt+"\n")
    #
    print_o3("Resolution","%.2f - %.2f" % (max(r),s.reso),\
                      "(%.2f - %.2f)\n" % (s.resoL,s.reso))
    print_o3("Completeness","%.1f%%" % s.compl, "(%.1f%%)" % s.complL)
    print_o3("I/sigma(I)","%.1f " % s.isig, "(%.1f)" % s.isigL)
    print_o3("Rmeas","%.1f%%" % s.rmeas, "(%.1f%%)" % s.rmeasL)
    print_o3("Rsym","%.2f%%" % s.rsym, "(%.1f%%)" % s.rsymL)
    print_o3("Compared","%d   " % s.compar, "(%d)" % s.comparL)
    print_o3("Measured","%d   " % s.total,"")
    print_o3("Unique","%d   " % s.unique,"")
    print_o3("Multiplicity","%.1f " % (s.total/s.unique),"")
    print_o3("Rejected misfits","%d   " % s.misfit,"")
    if s.absent:
        print_o3("   with <Iabs>/<I>","%.1f%%" % s.AbsIav,"")
    if friedels_law == "FALSE":
        print_o3("Anomalous contrib.","%.1f " % s.anom, "")
    print_o3("Wilson scaling (B/Corr)","%.1f  " % s.wilson_b,\
                                 "%.2f" % s.wilson_corr)
    print_o3("Estimated Res_max","%.2f" % s.dmin,"")
    print 
    resum.close()

def run_xscale(files, hklout, lattice, Dir=None,
               out_f=None, nbin=12, merge="TRUE", save=0):
    xscale_inp = """   ! RESOLUTION_SHELLS= %s
    SPACE_GROUP_NUMBER= %d
    UNIT_CELL_CONSTANTS= %s 
    OUTPUT_FILE=%s\n    FRIEDEL'S_LAW= %s MERGE= %s
    STRICT_ABSORPTION_CORRECTION= FALSE
    """
    if lattice.friedels_law: friedels_law = "TRUE"
    else: friedels_law = "FALSE"
    #spg = lattice.symmetry_num
    #cell = lattice.cell
    #print """cell=%s  spg=%d   Dir=%s
    #dmin=%s   nbin=%d   friedels_law=%s   merge=%s
    #""" % (cell,spg,Dir,dmin,nbin,friedels_law,merge)
    inp_line = "INPUT_FILE=%s  XDS_ASCII 100 %s\n"
    resbin = ""
    if not Dir: Dir = os.getcwd()
    for reso in res_bin(lattice.dmin,nbin-2): resbin = resbin + "%5.2f" % reso
    resbin = "15. 10. " + resbin
    #
    script = xscale_inp % (resbin,lattice.symmetry_num, lattice,
                           hklout, friedels_law, merge)
    for file in files:
         script = script + inp_line % (file, lattice.dmin)
    open(os.path.join(Dir,"XSCALE.INP"),"wb").write(script)
    #
    exec_prog("cd %s; %s" % (Dir,os.path.join(XDSHOME,"xscale_par")),
                              stdout=out_f, stderr="xupy.stderr")
    if save: saveLastVersion("XSCALE.LP")
    s = resum_scaling(lpf=os.path.join(Dir,"XSCALE.LP"))
    if not s:
        print "\nERROR while running XSCALE"
        sys.exit()

    r = 100, lattice.dmin
    write_xscale_resum(s,r,friedels_law)

def guess_imageType(image_name):
    """Return a the image detector type and compress type.
       Do not handel swaping..."""
    image = open(image_name,"r")
    head = image.read(9000)
    image.close()
    imageType = []
    imageCompress = "NO"
    if  head[:15] == "{\nHEADER_BYTES=" and head.count(";\nDETECTOR_SN=") and \
        head.count(";\nPIXEL_SIZE="): imageType.append("ADSC")
    if  head.count("\nSCANNER") and head.count("mar research") and \
        head.count("\nPROGRAM"): imageType.append("MarIP345")
    if  head[:3] == "II*" and \
        head[1028:1031] == "MMX" : imageType.append("MarCCD")
    if head.count("CCP4 packed image"): imageCompress = "PCK"
    if len(imageType) == 1: return imageType[0], imageCompress
    elif len(imageType) == 0: return "Unknown", imageCompress
    elif len(imageType) > 1:
        print "ERROR: Can't choose the detector type between:", imageType
        sys.exit()

def get_number_of_processors():
    platf = None
    try:
        if "linux" in sys.platform:
            platf = int(commands.getoutput("egrep -c '^processor' /proc/cpuinfo"))
        else:
            #"darwin" in sys.platform:
            # or [Free|Net|Open]BSD and MacOS X
            platf = int(commands.getoutput("sysctl -n hw.ncpu"))
    except:
        platf = 4
    return platf
