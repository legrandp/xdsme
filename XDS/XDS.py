#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 Changes in version 0.4.4 6-11-2009.
 - Improved index searching strategy.
 - Fix XIO compatibility with pilatus det.
 - Improve detector specific definitions in XIO exports.
 Changes in version 0.4.3 28-10-2009.
 - Fixe compatibility problems with Python2.2 and Python2.3
 - Add proper error message when image file can't be read.
 Changes in version 0.4.0 27-07-2009
 - Strategy: after indexing ask to choose which lattice and rerun idxref+xplan
 - Reference option for adding reference dataset (XPLAN parsing need modif).
 - during INTEGRATE, display overloads and mean strong refl/images.
 - Add compatibility for subprocess.Popen (python > 2.4.0) and Popen2 for
   previous versions

 TODO-0: If just the space group is selected and not the cell:
         try to find the proper cell if it is not ambigous
         (like P21212, P2122,,,),
 TODO-1: Add in the CORRECT summary the Rmrgd-F, image range, and total
         overloaded refl.
 TODO-2: Start multiple COLSPOT with different thresholds+ multiple IDXREF...
"""

__version__ = "0.4.5"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "23-11-2009"
__copyright__ = "Copyright (c) 2006-2009 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import os
import sys
import re

if sys.version_info <= (2, 4, 0):
    from popen2 import Popen3
else:
    from subprocess import Popen, PIPE

from pycgtypes import vec3
from xupy import XParam, xdsInp2Param, opWriteCl, \
                 saveLastVersion, LP_names, xdsinp_base, \
                 SPGlib, Lattice, resum_scaling, \
                 get_BravaisToSpgs
import XIO

_progname = os.path.split(sys.argv[0])[1]
_usage = """
   Running XDS automaticaly...

   USAGE:   %s  [OPTION]... FILES

      FILES is for one or multiple diffraction image files.

   OPTIONS:

    -h,  --help
         Print this help message.

    -1,-2,-3,-4,-5
         Go directly to a particular step:
         -1: XYCOOR + INIT
         -2: COLSPOT
         -3: IDXREF
         -4: DEFPIX + INTEGRATE
         -5: CORRECT

    -i,  --xds-input
         Give direct XDS Keyword input.
         For example: -i "DETECTOR_DISTANCE= 167.0 JOB= IDXREF AIR= 0.002"

    -a,  --anomal
         Distinguishes Friedel paires for scaling, strategy and completeness
         statistics. Default is no anoulous contribution.

    -A,  --Anomal
         Like -a, but also set "STRICT_ABSORPTION_CORRECTION" to True.
         It usualy gives better scaling statistics with redunduncy > 2.

    -r,  --high-resolution
         Set a high resolution cutoff. Default is 0 (no cutoff).

    -R,  --low-resolution
         Set a low resolution cutoff. Default is 9999.

    -b,  --beam-center-optimize-i
         Starting from the initial given values, search and optimize the beam
         center coordinates (given by -x, -y or extracted form the header).
         Best solution is chosen after i-score ranking.

    -B,  --beam-center-optimize-z
         Like -b/--beam-center-optimize-i, but best solution is chosen with after
         a z-score ranking.

    -d,  --distance
         Set the detector to crystal distance.

    -c,  --cell
         Set the expected cell.
         For example: -c "79 79 38 90 90 90"

    -f,  --reference FILE
         Defines a reference data set used during the XPLAN and CORRECT steps.
         For example: -f ../ref/XDS_ASCII.HKL

    -O,  --oscillation
         Set frame oscillation range in degree.
         For example: -c 0.5 

    -p,  --project
         Set the project name. The default is the prefix taken from
         image names. The working directory will be: xds_process_"project"  

    -s,  --spg
         Set the expected space group using either the space group number
         or simple string.
         For example: -s 18 or -s P21212

    -S, --strategy
         Force to go for calculating strategy (XPLAN) and then stops.

    -x,  --beam-x
         Set a new value for ORGX: X-coordinates (in pixels) of the
         detector origin.

    -y,  --beam-y
         Set a new value for ORGY: Y-coordinates (in pixels) of the
         detector origin.

    -v,  --verbose
         Turn on verbose output.

    -w, --wavelength
         Set the x-ray wavelength.

    --slow,
         Set parameters to process either more accurately.

    --weak,
         Set parameters to index in case of weak spots.

""" % _progname

_fmt_hello = """
    Simplified XDS Processing\n
    Diffraction Setup Parameters:\n
  Detector distance:             %(DETECTOR_DISTANCE)8.2f mm
  X-ray wavelength:            %(X_RAY_WAVELENGTH)10.4f A
  Oscillation range:           %(OSCILLATION_RANGE)10.4f degree\n
  Beam coordinate X:             %(ORGX)8.1f pixel
                  Y:             %(ORGY)8.1f pixel
  Resolution range:           %(INCLUDE_RESOLUTION_RANGE)11s
  Image range:                %(DATA_RANGE)11s
"""
#         RMSd spot position    %%() pixels   %%() degree

_fmt_final_stat = """
      Refined Parameters and Scaling Statistics
      =========================================\n
      Image Range   %(image_start)5d  to  %(image_last)5d

      Space group   number    %(spg_num)d
                    symbol    %(spg_sym)s

      Cell parameters     %(cell)s

      Resolution           %(LowestReso)8.2f -%(reso)6.2f\
    (%(resoL).2f - %(reso).2f)

      Completeness                    %(compl)5.1f%%   (%(complL).1f%%)
      I/sigma(I)                     %(isig)6.1f    (%(isigL).1f)
      Rmeas                          %(rmeas)6.1f%%   (%(rmeasL).1f%%)
      Rsym                          %(rsym)7.2f%%   (%(rsymL).1f%%)
      Multiplicity               %(multiplicity)10.1f
      Compared                   %(compar)10d    (%(comparL)d)
      Measured                   %(total)10d
      Unique                     %(unique)10d
      Rejected misfits           %(misfit)10d
      Wilson scaling (B/Corr)    %(wilson_b)10.1f    (%(wilson_corr).2f)
"""

_fmt_AbsIav = "      Systematic absente reflection measured %(AbsNum)6d \
 with <Iabs>/<I> =  %(AbsIav).1f%%\n"

_fmt_anomal = """
      Anomalous pairs measured   %(anoNum)10d
      SigAno                     %(anoSig)10.3f    (%(anoSigL).3f)
      Anomalous Correlation      %(anoCorr)10.1f%%   (%(anoCorrL).1f%%)
"""

test_str = """
    s.anoNum, s.anoNumL  =  TG[-1][-1], TG[-2][-1]
    s.anoSig, s.anoSigL  =  TG[-1][-2], TG[-2][-2]
    s.anoCorr,s.anoCorrL =  TG[-1][-3], TG[-2][-3]
    if friedels_law == "FALSE":
        print_o3("Anomalous contrib.","%.1f " % s.anom, "")
    print_o3("Estimated Res_max","%.2f" % s.dmin,"")
"""

STEPMARK = re.compile(r"^( [*]{5} (\w{4,}) [*]{5}  )")
integrate_step = re.compile(r" PROCESSING OF IMAGES ")
integrate_mosaicity = re.compile(r"CRYSTAL MOSAICITY \(DEGREES\)")
integrate_strong = re.compile(r"REFLECTIONS ACCEPTED FOR REFINEMENT")
_rF, _rI = r"[\ ]+([0-9\.]+)", r"[\ ]+([\d]+) "
scalegr =   r" "+_rI+r"  (\d)"+_rF+r"  ....... "+4*_rI+2*_rF
sc = re.compile(scalegr)

XDSHOME = os.getenv('XDS')

def unpack_latticefit2(lattice_string):
    "From lattice_string to Lattice object."
    lats = lattice_string[2:].split()
    latt = Lattice((map(float, lats[3:9])), lats[1])
    latt.fit = float(lats[2])
    latt.character = int(lats[0])
    #latt.reindexing = tuple(map(int,ss[9:]))
    return latt

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)

def makeXDSImageLinks(imagename_list, linkDir="img_links",
                       prefix="image", start_num=1):
    """All image names in the imagename_list are supposed to be part
    of one continous sequence of collected images.
    Todo:
        - How to safely modulate PhiStart outside the [-180,180] range ?
    """
    link_list = []
    if linkDir not in os.listdir("."):
        try:
            _mkdir(linkDir)
        except Exception, e:
            print "Error\n", e
            sys.exit(0)
    #
    linkDir = os.path.abspath(linkDir)
    _collect = {}
    _osc_range = []
    for _image in imagename_list:
        _I = XIO.Image(_image)
        if VERBOSE:
            print _image
        # How to safely modulate PhiStart outside the [-180,180] range?
        if VERBOSE:
            print "\tPhiStart %8.2f" % _I.header['PhiStart']
        if VERBOSE:
            print "\tPhiWidth %8.2f" % _I.header['PhiWidth']
        _collect[_I.header['PhiStart']] = _image
        _osc_range.append(_I.header['PhiWidth'])

    if max(_osc_range) != min(_osc_range):
        print "Error. Image list contains different oscillation range!"
        sys.exit(0)
    #
    _osc_starts = _collect.keys()
    _osc_starts.sort()
    _osc1 = _osc_starts[0]
    _oscR = _osc_range[0]
    for _osc in _osc_starts:
        _num =  start_num+ (_osc-_osc1)/_oscR
        link_name = os.path.join(linkDir, prefix+"_%04.0f.img" % _num)
        #print _osc,_osc1,_oscR,1+(_osc-_osc1)/_oscR, _num,link_name
        if os.path.lexists(link_name) and os.path.islink(link_name):
            if VERBOSE:
                print "==> Removing existing link: %s" % link_name
            os.remove(link_name)
        os.symlink(os.path.abspath(_collect[_osc]), link_name)
        link_list.append(link_name)
    return link_list

class XDSLogParserException(Exception):
    """This level of exception raises a recoverable error which can be fixed."""

class XDSExecError(Exception):
    ""

class XDSLogParser:
    """ Parser for the xds *.LP files.
    """
    def __init__(self, filename, run_dir="", verbose=False, raiseErrors=True):
        self.results = {}
        self.info = "XDS Parser"
        self.fileType = "XDS"
        self.verbose = verbose
        #
        if not run_dir:
            run_dir = "./"
        self.run_dir = run_dir
        #
        full_filename = os.path.join(self.run_dir, filename)
        #
        try:
            fp = open(full_filename, "r")
            self.lp = fp.read()
            fp.close()
        except:
            raise IOError, "Can't read file: %s" % full_filename
        # Catch Errors:
        _err = self.lp.find(" !!! ERROR !!!" )
        _err_type = None
        _err_level = None
        if _err != -1:
            _err_msg = self.lp[_err:]
            if _err_msg.count("INSUFFICIENT PERCENTAGE (< 70%) OF INDEXED REF"):
                _err_type = "IDXREF. Percentage of indexed"
                _err_type += " reflections bellow 70%.\n"
                _err_level = "WARNING"
            if _err_msg.count("INSUFFICIENT NUMBER OF ACCEPTED SPOTS."):
                _err_type = "IDXREF. Insufficient number of accepted spots."
                _err_level = "FATAL"
        if _err_level in ("FATAL", "ERROR") and raiseErrors:
            raise XDSExecError, _err_type
        #
        if self.verbose and _err != -1:
            print "\n %s in %s" % (_err_level, _err_type)
        #
        if full_filename.count("INIT.LP"): self.parse_init()
        elif full_filename.count("COLSPOT.LP"): self.parse_colspot()
        elif full_filename.count("IDXREF.LP"): self.parse_idxref()
        elif full_filename.count("XPLAN.LP"): self.parse_xplan()
        elif full_filename.count("INTEGRATE.LP"): self.parse_integrate()
        elif full_filename.count("CORRECT.LP"): self.parse_correct()
        else:
            raise IOError, "Don't know how to parse file: %s" % full_filename

    def get_par(self, match, limit=75, func=None,
                      multiLine=False, start=0, before=False):
        "Extract parameters from XDS .LP lines."
        try:
            if before:
                limit = start
                start = self.lp.index(match)-start
            else:
                start = self.lp.index(match, start) + len(match)
        except Exception, e:
            raise e
            #return 0
        if multiLine:
            _raw = self.lp[start:start+limit].split()
        else:
            _raw = self.lp[start:start+limit].splitlines()[0].split()
        if not func:
            for vType in (int, float, str):
                try:
                    vType(_raw[0])
                    func = vType
                except ValueError:
                    pass
                if func: break
        if not func:
            raise ValueError, "get_par function can't process value '%s'" % _raw
        pars = map(func, _raw)
        if len(pars) == 1: return pars[0]
        else: return pars

    def _get_lattices_table(self):
        list_latticesFit = []
        st1 = self.lp.index("LATTICE-  BRAVAIS-   QUALITY")
        _table = self.lp[st1:st1+6000].splitlines()[3:47]
        #for sol in _table: print sol
        return map(unpack_latticefit2, _table)

    def _get_index_origins_table(self):
        st0 = self.lp.index(" DL\n  ORIGIN\n")+14
        st1 = self.lp.index(" SELECTED:     INDEX_ORIGIN=")-2
        return map(lambda s: \
                  map(float, s.split()), self.lp[st0:st1].splitlines())

    def parse_init(self):
        "Parse INIT.LP"
        R_d, R_f = self.results, self.get_par
        #
        R_d["background_range"] = R_f("BACKGROUND_RANGE=")
        R_d["mean_gain"] = R_f("MEAN GAIN VALUE")
        R_d["min_gain"] = R_f("MINIMUM GAIN VALUE IN TABLE")
        R_d["max_gain"] = R_f("MAXIMUM GAIN VALUE IN TABLE")
        R_d["mean_background"] = R_f("BACKGROUND COUNTS IN A DATA IMAGE PIXEL")
        #
        prp =  "  Looking at images %(background_range)s\n"
        prp += "  Mean Gain:        %(mean_gain).1f\n"
        prp += "  Min table gain:   %(min_gain).2f\n"
        prp += "  Max table gain:   %(max_gain).2f\n"
        prp += "  Mean Backgroud:   %(mean_background).1f\n"
        if self.verbose:
            print prp % R_d
        return R_d, prp

    def parse_colspot(self):
        "Parse COLSPOT.LP"
        R_d, R_f = self.results, self.get_par
        #
        R_d["strong_pixels"] = R_f("EXTRACTED FROM IMAGES")
        R_d["weak_spots_ignored"] = R_f("WEAK SPOTS OMITTED")
        R_d["out_of_center_spots"] = R_f("SPOT MAXIMUM OUT OF CENTER")
        R_d["spot_number"] = self.get_spot_number()
        R_d["time"] = R_f("elapsed wall-clock time", 11)

        prp =  "  Number of spots found:    %(spot_number)10d\n"
        prp += "  Out of center rejected:   %(out_of_center_spots)10d\n"
        prp += "  Weak spots rejected:      %(weak_spots_ignored)10d\n"
        prp += "  Number of spots accepted: %(spot_number)10d\n"
        if self.verbose:
            print prp % R_d
        return R_d, prp

    def parse_idxref(self):
        "Parse IDXREF.LP"
        R_d, R_f = self.results, self.get_par
        #
        rexp1 = r".* (\d+) OUT OF\ +(\d+) SPOTS INDEXED\..*"
        rexp2 = r".* QX=\ +([\d|\.]+)\ +QY=\ +([\d|\.]+)"
        #if "! ERROR !" in self.lp:
        #    raise XDSLogParserException, "Error while parsing XDS logfile"
        nIS, nTS = map(int, re.match(rexp1, self.lp, re.DOTALL).groups())
        qx, qy = map(float, re.match(rexp2, self.lp, re.DOTALL).groups())
        meanPixel = (qx+qy)/2
        R_d["indexed_spots"] = nIS
        R_d["total_spots"] = nTS
        R_d["indexed_percentage"] = 100.*nIS/nTS
        #
        st0 = self.lp.index("START OF INTEGRATION *****")
        st1 = "STANDARD DEVIATION OF SPOT    POSITION (PIXELS)"
        st2 = "STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)"
        st3 = "UNIT CELL PARAMETERS"
        st4 = "SPACE GROUP NUMBER"
        st5 = "COORDINATES (PIXELS) OF DIRECT BEAM"
        #
        R_d["oscillation_range"] = R_f("OSCILLATION_RANGE=") 
        R_d["xy_spot_position_ESD"] = R_f(st1, start=st0)
        R_d["z_spot_position_ESD"] = R_f(st2, start=st0)
        R_d["index_origin_table"] = self._get_index_origins_table()
        R_d["lattices_table"] = self._get_lattices_table()
        R_d["refined_cell"] = R_f(st3, start=st0)
        R_d["refined_cell_str"] = 6*"%.2f " % \
                                           tuple(R_d["refined_cell"])
        R_d["space_group_number"] = R_f(st4, start=st0)
        R_d["direct_beam_pixels"] = R_f(st5, start=st0)
        R_d["direct_beam_mm"] = R_d["direct_beam_pixels"][0]*qx, \
                                           R_d["direct_beam_pixels"][1]*qy
        R_d["bmx"], R_d["bmy"] = R_d["direct_beam_mm"]
        R_d["bpx"], R_d["bpy"] = R_d["direct_beam_pixels"]

        origin_t = R_d["index_origin_table"]
        origin_n = len(origin_t)
        quality_t = [x[3] for x in origin_t if x[3] < 2.]
        #R_d["index_score"] = reduce(lambda a,b: a+b, quality_t)/len(quality_t)
        max_ot = min(origin_n, 5)
        R_d["shift_pixel"] = origin_t[0][4]*meanPixel
        R_d["shift_mm"] = origin_t[0][4]
        
        prp = """  Unit cell parameters:   %(refined_cell_str)s
  Space group number:     %(space_group_number)s
  Indexed spots:          %(indexed_percentage).1f%% (%(indexed_spots)d/%(total_spots)d)
  Spot prediction ESD:       %(xy_spot_position_ESD).2f   pixels and  %(z_spot_position_ESD).2f degrees
  Refined beam position (in mm):      (%(bmx)9.3f, %(bmy)9.3f)
  Refined beam position (in pixels):  (%(bpx)9.2f, %(bpy)9.2f)
  Shift in beam position: %(shift_mm)9.2f mm  (%(shift_pixel).1f pixels)
"""
        prp2 = "  Size of the origine index table: %(origin_n)7d\n" % vars()
        ppa, ppb = "\n\tQuality:       ", "\n\tShift (mm):    "
        ppc, ppd = "\n\tShift (pixels):", "\n\tBeam X (mm):   "
        ppe, ppf = "\n\tBeam Y (mm):   ", "\n\tIndex Origin:  "
        for i in range(max_ot):
            ppa += "%9.2f," % (origin_t[i][3])
            ppb += "%9.2f," % (origin_t[i][4]*meanPixel)
            ppc += "%9.1f," % (origin_t[i][4])
            ppd += "%9.1f," % (origin_t[i][5]*qx)
            ppe += "%9.1f," % (origin_t[i][6]*qy)
            ppf += "%3d%3d%3d," % tuple(origin_t[i][0:3])
        prp2 += "  Origin ranking for the best %d solutions: " % max_ot
        prp2 += ppa[:-1] + ppb[:-1] + ppc[:-1]
        prp2 += ppd[:-1] + ppe[:-1] + ppf[:-1] + "\n"
        #prp += " Index origin score: %.2f\n" % (R_d["index_score"])
        if self.verbose == 1:
            print (prp + prp2) % R_d
        elif self.verbose == 2:
            print prp % R_d
        return R_d, prp

    def parse_integrate(self):
        "Parse INTEGRATE.LP"
        R_d, R_f = self.results, self.get_par
        R_d["reflections"] = R_f("REFLECTIONS SAVED ON FILE",
                                  start=9, func=int, before=True)
        R_d["divergence"] = R_f("BEAM_DIVERGENCE_E.S.D.= ")
        R_d["mosaicity"] = R_f("REFLECTING_RANGE_E.S.D.= ")
        prp =  "\n  Number of reflection integrated:      %(reflections)d\n"
        prp += "  Estimated divergence:                 %(divergence).3f\n"
        prp += "  Estimated mosaicity:                  %(mosaicity).3f\n"
        if self.verbose:
            print prp % R_d
        return R_d, prp

    def parse_xplan(self):
        "Parse XPLAN.LP"
        R_d, R_f = self.results, self.get_par
        R_d["spacegroup"] = R_f("SPACE_GROUP_NUMBER=")
        R_d["unitcell"] = 6*" %8.2f" % tuple(R_f("UNIT_CELL_CONSTANTS="))
        R_d["friedels_law"] = R_f("FRIEDEL'S_LAW=")[0]
        st0 = self.lp.index(72*"*")
        st1 = self.lp.index(72*"*", st0+72)
        st2 = self.lp.index(72*"*", st1+72)
        #
        prp =  "  Friedel's law:   %(friedels_law)s\n"
        prp += "  Spacegroup:      %(spacegroup)d\n"
        prp += "  Unitcell:        %(unitcell)s\n"
        if self.verbose:
            print prp % R_d
            print
            print self.lp[st0:st2]
        return R_d, prp

    def parse_correct(self):
        "Parse CORRECT.LP"
        R_d, R_f = self.results, self.get_par

        R_d["RMSd_spotPosition"] = R_f("SPOT    POSITION (PIXELS)")
        R_d["RMSd_spindlePosition"] = R_f("SPINDLE POSITION (DEGREES)")
        R_d["Mosaicity"] = R_f("CRYSTAL MOSAICITY (DEGREES)")
        r = R_f(" "+"-"*74+"\n")
        R_d["I_sigma"], R_d["Rsym"] = r[2], r[4]
        R_d["Compared"], R_d["Total"] = r[6], r[7]
        ### Select Diffraction range.
        sp1 = self.lp.index("RESOLUTION RANGE  I/Sigma")
        sp2 = self.lp.index(10*"-", sp1)
        _table = self.lp[sp1:sp2].splitlines()[3:-1]
        _table = [ map(float, l.split()[1:3]) for l in _table ]
        R_d["HighResCutoff"] = self.get_proper_resolition_range(_table)
        prp = ""
        if R_d["Mosaicity"]:
            prp += "  RMSd spot position in pixel:      %(RMSd_spotPosition)9.2f\n"
            prp += "  RMSd spot position in degree:     %(RMSd_spindlePosition)9.2f\n"
            prp += "  Refined Mosaicity:                %(Mosaicity)9.2f\n\n"
        prp += "  Rsym:                             %(Rsym)9.1f\n"
        prp += "  I/sigma:                          %(I_sigma)9.1f\n"
        if R_d["HighResCutoff"]:
            prp += "  Suggested high resolution cutoff: %(HighResCutoff)9.2f\n"
        prp += "  Compared reflections:                 %(Compared)d\n"
        prp += "  Total number of measures:             %(Total)d\n"
        if self.verbose:
            print prp % R_d
        return R_d, prp

    def get_proper_resolition_range(self, res_table):
        "High res is selected when at least 3 values of I/sigma are below 1."
        _highT, _high = [], None
        for res, IoS in res_table:
            if IoS < 1.:
                _highT.append(res)
                if not _high and len(_highT) == 3:
                    _high = _highT[0]
            else:
                _highT = []
            #print "%8.3f  %8.3f  %s" % (res, IoS, IoS >= 1.)
        if not _high and len(_highT) >= 1:
            _high = _highT[0]
        #print "Suggested high resolution cut-off: %.2f" % _high
        return _high

    def get_spot_number(self):
        "Read the number of spot directly from SPOT.XDS"
        _execstr = "wc -l %s/SPOT.XDS" % self.run_dir
        if sys.version_info <= (2, 4, 0):
            s = os.popen(_execstr)
            wc = s.read(); s.close()
        else:
            wc = Popen([_execstr], stdout=PIPE, shell=True).communicate()[0]
        return int(wc.split()[0])

MIN_SPOT_NUMBER = 200
LATTICE_GEOMETRIC_FIT_CUTOFF = 50
FRAMES_PER_COLSPOT_SEQUENCE = 16 # number of frames per sequence in COLSPOT.
STEPS = "INIT", "COLSPOT", "IDXREF", "INTEGRATE", "CORRECT"
SPOTFILENAME = "SPOT.XDS"

class XDS:
    "Main class for runing xds step by step."

    def __init__(self, obj=None, linkToImages=True):
        """Constructor for the Param classes from file or string."""
        #
        self.linkToImages = linkToImages
        self.__cancelled = 0
        self.__lastOutp = 0
        self.mode = []
        if XDSHOME:
            self.__execfile = os.path.join(XDSHOME,"xds_par")
        else:
            self.__execfile = "xds_par"
        self.running = 0
        self.outp = []
        self.run_dir = "."
        self.status = None
        self.inpParam = XParam()
        self.collectDir = "./"
        self.link_name_to_image = "img"
        #
        if type(obj) == file:
            exec obj.read() in self.inpParam.__dict__
            obj.close()
        if type(obj) == str:
            exec obj in self.inpParam.__dict__

    def set_collect_dir(self, dirname):
        "Set the collect directory"
        if os.path.isdir(dirname):
            self.collectDir = dirname
        else:
            raise XIO.XIOError, "Can't find %s directory" % dirname

    def _creat_process(self, _execstr):
        "Return a process with pipe redirected IO."
        if sys.version_info <= (2, 4, 0):
            self.wait_value = -1
            return Popen3(_execstr)
        else:
            self.wait_value = None
            return Popen(_execstr, stdin=PIPE, stdout=PIPE, 
                              stderr=PIPE, bufsize=1, close_fds=True,
                              universal_newlines=True)

    def cancel (self):
        ""
        self.__cancelled = 1

    def getOutp(self):
        if not self.__cancelled:
            nLine = len(self.outp)
            diff = nLine - self.__lastOutp
            self.__lastOutp = nLine
            if diff:
                return "".join(self.outp[-diff:])[:-1]
            else: return ""

    #def copyInitFilesTo(self, newDir):
    #    initFileToCopy = ""

    def saveLastLP(self, file_names, suffix=""):
        "Use xupy function to copy all LP files that have changed."
        saveLastVersion(file_names, suffix="")

    def run(self, run_dir=None, rsave=None, verbose=True):
        "Control the runing of the xds process and parse the output."
        self.__cancelled = 0
        self.running = 1
        self.step = 0
        self.step_name = ""
        self.outp = []
        self.init_dir = os.getcwd()
        if run_dir:
            self.run_dir = run_dir
        if not self.run_dir:
            self.run_dir = "."#
        result = 0
        #
        if self.run_dir:
            if not os.path.exists(self.run_dir):
                try:
                    os.mkdir(self.run_dir)
                except:
                    print
                    raise XIO.XIOError, \
                       "STOP! Can't creat xds working directory:",self.run_dir
            #
            if os.path.isdir(self.run_dir):
                os.chdir(self.run_dir)
                if self.linkToImages:  
                    if not os.path.exists(self.link_name_to_image):
                        os.system("ln -sf %s %s" % (self.collectDir, \
                                                    self.link_name_to_image))
                        #os.system("ln -sf .. %s" % (self.link_name_to_image))
                    #else:
                    #    raise XIO.XIOError, \
                    #     "STOP! Can't creat link %s in working directory: %s" \
                    #     % (self.link_name_to_image, self.run_dir)
        #
        #
        opWriteCl("XDS.INP", "%s" % self.inpParam)
        #
        xdsProcess = self._creat_process(self.__execfile)
        _init_parse = True
        overloaded_spots = 0
        while self.running:
            self.status = xdsProcess.poll() 
            if self.status != self.wait_value:
                self.running = 0
                break
            if self.__cancelled:
                os.kill(xdsProcess.pid, 9)
                break
            if self.wait_value == -1:
                nl = xdsProcess.fromchild.readline()
            else:
                nl = xdsProcess.stdout.readline()
                #nl = xdsProcess.communicate()
            # inline parsing of stdout
            if self.step_name == "INTEGRATE":
                if _init_parse:
                    print "    Processing    Mean #Strong  ",
                    print "Estimated   Overloaded"
                    print "    Image Range   refl./image   ",
                    print "Mosaicity   reflections\n"
                    _talbInt = []
                    _init_parse = False
                if integrate_step.search(nl):
                    print nl[44:50]+" - "+nl[56:-1],
                    nimages = int(nl[56:-1]) - int(nl[44:50]) + 1
                elif integrate_strong.search(nl):
                    print "%11.0f" % (float(nl.split()[0])/nimages),
                elif integrate_mosaicity.search(nl):
                    print " %11.3f" % float(nl.split()[3]),
                    print " %11d" %  overloaded_spots
                    overloaded_spots = 0
                hit = sc.search(nl)
                if hit:
                    _talbInt = hit.groups()
                    overloaded_spots += int(hit.groups()[3])
            sm = STEPMARK.match(nl)
            if sm:
                self.step += 1
                self.step_name = sm.group(2)
                #if VERBOSE:
                if verbose:
                    print "\n --->  Running job: %20s\n" % self.step_name
            if nl: self.outp.append(nl)
        self.step += 1
        self.step_name = "FINISHED"
        if self.__cancelled: result = -1
        if rsave:
            saveLastVersion(LP_names)
        #if VERBOSE:
        #    print "End of XDS run"
        os.chdir(self.init_dir)
        return 1

    def spots_resolution_cutoff(self, res_cutoff, verbose=False):
        "Read the SPOT.XDS file and filter spots using a resolution cutoff."
        from math import atan2, sin
        import shutil
        #
        spotsFileName = os.path.join(self.run_dir, "SPOT.XDS")
        # Save the SPOT file and open a new one
        shutil.copy(spotsFileName, spotsFileName+".bck")
        spots = open(spotsFileName+".bck").readlines()
        newspots = open(os.path.join(self.run_dir, SPOTFILENAME),"w")
        # Get parameters for the resol calculation
        xo, yo = self.inpParam["ORGX"], self.inpParam["ORGY"]
        rx, ry = self.inpParam["QX"], self.inpParam["QY"]
        D = self.inpParam["DETECTOR_DISTANCE"]
        # the resolution calculation function
        resolCal = lambda s, D, xo, yo, rx, ry: \
                   0.5/sin(atan2(((rx*(float(s[: 10])  -xo))**2 +
                                  (ry*(float(s[10:20])-yo))**2)**0.5,D)/2.)
        filtredSpots = [s for s in spots \
	                  if resolCal(s,D,xo,yo,rx,ry) >= res_cutoff]
        #
        newspots.writelines(filtredSpots)
        ni, nf = len(spots), len(filtredSpots)
        if verbose:
            print ">> Selected spots with %.2f resolution cutoff:" % \
	                                                     (res_cutoff),
            print "%d / %d (%.1f%%)" % (nf, ni, nf*100./ni)
        newspots.close()

    def run_init(self):
        "Runs the 2 first steps: XYCORR and INIT"
        self.inpParam["TRUSTED_REGION"] = [0, 1.20]
        self.inpParam["JOB"] = "XYCORR", "INIT"
        i1, i2 = self.inpParam["DATA_RANGE"]
        if "slow" in self.mode:
            self.inpParam["BACKGROUND_RANGE"] =  i1, min(i2, i1+11)
        else:
            self.inpParam["BACKGROUND_RANGE"] =  i1, min(i2, i1+3)
        self.run(rsave=True)
        res = XDSLogParser("INIT.LP", run_dir=self.run_dir, verbose=1)
        return res.results

    def run_colspot(self):
        "Runs the COLSPOT step."
        self.inpParam["JOB"] = "COLSPOT",
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 1
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 8
        _trial = 0

        if "slow" in self.mode:
            FRAMES_PER_COLSPOT_SEQUENCE = 32
        elif "fast" in self.mode:
            FRAMES_PER_COLSPOT_SEQUENCE = 4
        else:
            FRAMES_PER_COLSPOT_SEQUENCE = 16
        if "weak" in self.mode:
            self.inpParam["STRONG_PIXEL"] = 4.5
            FRAMES_PER_COLSPOT_SEQUENCE = 32
        # Selecting spot range(s),
        # self.inpParam["SPOT_RANGE"] is set to Collect.imageRanges by the
        # xds export function XIO
        # 
        _dc = XIO.Collect("foo_001.bar")
        _dc.imageNumbers = _dc._ranges_to_sequence(self.inpParam["SPOT_RANGE"])
        #
        min_fn, max_fn = self.inpParam["DATA_RANGE"] 
        dPhi = self.inpParam["OSCILLATION_RANGE"]
        _fpcs = FRAMES_PER_COLSPOT_SEQUENCE
        _2fpcs = 1 + 2 * FRAMES_PER_COLSPOT_SEQUENCE

        if (max_fn - min_fn + 1) >= _2fpcs:
            # use two range ex: i-i+2, f-2,f
            # with f at maximum 90 degre distance
            max_frame = min(max_fn, min_fn + int(89./dPhi + _fpcs))
            spot_ranges = ((min_fn, min_fn + _fpcs - 1),
                          (max_frame - _fpcs + 1, max_frame))
        else:
            spot_ranges = (min_fn, min(min_fn + _2fpcs - 1, max_fn)),
        # Restrict to matching collected images...
        #print "T1", spot_ranges
        self.inpParam["SPOT_RANGE"] = _dc.lookup_imageRanges(False, \
                                              mask_range=spot_ranges)
        #
        self.run(rsave=True)
        _rs = "  Image range(s) for spot collection: "
        for _r in self.inpParam["SPOT_RANGE"]:
            _rs += ("  [%d - %d]," % tuple(_r))
        print _rs[:-1] + "\n"

        res = XDSLogParser("COLSPOT.LP", run_dir=self.run_dir, verbose=1)
        while res.results["spot_number"] < MIN_SPOT_NUMBER and _trial < 4:
            _trial += 1
            self.inpParam["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] -= 1
            self.inpParam["STRONG_PIXEL"] -= 1.5
            #self.inpParam["SPOT_MAXIMUM_CENTROID"] += 1
            print "Insuficiant number of spot (minimum set to %d)." % \
                                                         MIN_SPOT_NUMBER
            print "Recollecting spots. Trial number %d" % _trial
            self.run(rsave=True)
            res = XDSLogParser("COLSPOT.LP", run_dir=self.run_dir, verbose=1)
        return res.results

    def run_idxref(self, beam_center_search=False, ranking_mode="ZSCORE"):
        "Runs the IDXREF step. Can try to search for better beam_center."
        #
        self.inpParam["JOB"] = "IDXREF",
        self.inpParam["TRUSTED_REGION"] = [0, 0.98]
        self.run(rsave=True)
        try:
            res = XDSLogParser("IDXREF.LP", run_dir=self.run_dir, verbose=1)
        except XDSExecError, err:
            print "ERROR in", err

        RD = res.results
        qx, qy = self.inpParam["QX"], self.inpParam["QY"]
        dist = self.inpParam["DETECTOR_DISTANCE"]
        det_X = vec3(self.inpParam["DIRECTION_OF_DETECTOR_X-AXIS"])
        det_Y = vec3(self.inpParam["DIRECTION_OF_DETECTOR_Y-AXIS"])
        det_Z = det_X.cross(det_Y)
        det_params = dist, det_X, det_Y, det_Z, qx, qy

        #RD["indexed_percentage"] < 70. or \
        #if beam_center_search or RD["xy_spot_position_ESD"] > 2. or \
        #  RD["z_spot_position_ESD"] > 2*self.inpParam["OSCILLATION_RANGE"]:
        if beam_center_search:
            TestResults = [RD]
            print " Number of possible beam coordinates: %d" % \
                              len(RD["index_origin_table"])
            maxTestOrigin = min(60, len(RD["index_origin_table"]))
            origins = RD["index_origin_table"][:maxTestOrigin]
            for origin in origins:
                # We first need to calculate the beam_origin from the
                # beam_coordinate and beam_vector given in the table
                beam = vec3(origin[7:10])
                beam_origin = get_beam_origin(origin[5:7], beam, det_params)
                self.inpParam["ORGX"] = beam_origin[0]
                self.inpParam["ORGY"] = beam_origin[1]
                self.inpParam["INCIDENT_BEAM_DIRECTION"] = tuple(beam)
                #print "DEBUG:  %7.1f %7.1f  - %7.1f %7.1f" % \
                #  (coorx, coory, self.inpParam["ORGX"], self.inpParam["ORGY"])
                print "   Testing beam coordinate: (%.2fmm, %.2fmm) = " % \
                                           (origin[5]*qx, origin[6]*qy), 
                print "  %.1f, %.1f" % (origin[5], origin[6])
                self.run(rsave=True, verbose=False)
                try:
                    TestResults.append(XDSLogParser("IDXREF.LP", 
                                           run_dir=self.run_dir, 
                                           verbose=0, raiseErrors=True).results)
                except XDSExecError, err:
                    print "\t\tError in", err
            print "\n"
            # Need to lookup in the results for the beam-center giving 
            best_index_rank = rank_indexation(TestResults, ranking_mode)
            #for o in origins:
            #    print origins.index(o), o[:-3]
            best_origin = origins[best_index_rank[ranking_mode]-1]
            if VERBOSE:
                print best_index_rank
                fmt = "%4i%4i%4i%7.2f%7.2f%8.1f%8.1f%9.5f%9.5f%9.5f"
                print "best_index_rank", best_index_rank[ranking_mode]
                print "best_origin", fmt % tuple(best_origin[:10])
            best_beam = vec3(best_origin[7:10])
            best_beam_coor = best_origin[5:7]
            best_beam_orig = get_beam_origin(best_beam_coor,
                                             best_beam, det_params)
            self.inpParam["ORGX"], self.inpParam["ORGY"] = best_beam_orig
            if VERBOSE:
                print best_beam_orig, best_beam_coor
            self.inpParam["INCIDENT_BEAM_DIRECTION"] = tuple(best_beam)
            self.run(rsave=True)
            res = XDSLogParser("IDXREF.LP", run_dir=self.run_dir)
        return res.results

    def check_fileout(self, fileout):
        "Checking normal terminaison."
        if not os.path.exists(os.path.join(self.run_dir, fileout)):
            err = "Abnormal terminaison. Can't locate file: '%s'" % fileout
            print err
            raise Exception(err)

    def run_xplan(self, image_ranges, ridx=None):
        "Running the strategy."
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 8
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 1

        select_strategy(ridx, self.inpParam)
        print "\n Starting strategy calculation."
        self.inpParam["JOB"] = "IDXREF",
        self.run(rsave=True)
        res = XDSLogParser("IDXREF.LP", run_dir=self.run_dir, verbose=2)
        # Select just the internal circle of the detector.
        self.inpParam["TRUSTED_REGION"] = 0.0, 1.0
        self.inpParam["JOB"] = "DEFPIX", "XPLAN"
        self.run(rsave=True)
        res = XDSLogParser("XPLAN.LP", run_dir=self.run_dir, verbose=1)      
        sys.exit(0)

    def run_integrate(self, image_ranges, ridx=None):
        "Running INTEGRATE."
        self.inpParam["TRUSTED_REGION"] = [0, 1.0]
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 8
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 1
        if self.mode == "slow":
            self.inpParam["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA"] = 13
            self.inpParam["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA"] = 13

        "Runs the 2 first steps: DEFPIX and INTEGRATE"
        if len(image_ranges) == 1:
            self.inpParam["JOB"] = "DEFPIX", "INTEGRATE"
            self.run(rsave=True)
            res = XDSLogParser("INTEGRATE.LP", run_dir=self.run_dir, verbose=1)
            self.check_fileout("INTEGRATE.HKL")
        else:
            #print "\n Error in the INTEGRATE step:"
            print "\n Image range:", image_ranges
            print " Multi-sweep integration not yet implemanted. Sorry.\n"
            sys.exit(0)
        return res.results

    def run_correct(self):
        "Runs the last step: CORRECT"
        self.inpParam["JOB"] = "CORRECT",
        self.run(rsave=True)
        res = XDSLogParser("CORRECT.LP", run_dir=self.run_dir, verbose=1)
        L, H = self.inpParam["INCLUDE_RESOLUTION_RANGE"]
        newH = res.results["HighResCutoff"]
        if newH > H:
            print "   ->  Rerunning CORRECT with a new high resolution cutoff."
            self.inpParam["INCLUDE_RESOLUTION_RANGE"] = L, newH
            self.run(rsave=True)
        s = resum_scaling(lpf=os.path.join(self.run_dir,"CORRECT.LP"))
        s["image_start"], s["image_last"] = self.inpParam["DATA_RANGE"]
        if not s:
            print "\nERROR while running CORRECT"
            sys.exit()
        print _fmt_final_stat % vars(s)
        if s.absent:
            print _fmt_AbsIav % vars(s)
        if self.inpParam["FRIEDEL'S_LAW"] == "FALSE":
            print _fmt_anomal % vars(s)
        return res.results

    def run_scaleLaueGroup(self):
        """Runs the CORRECT step with reindexation for all the selected Laue
           group
        1 - Get the selected bravais Lattices from IDXREF
        2 - Filtrate the equivalents (same geometry and reindexation)
        3 - For each one of the selected lattices:
                    in a seperated dir,
                    for all the laue symmetry compatible with
                    the bravais lattice geometry run the CORRECT scaling 
        4 - Rank all the scaling from the parsing of all the CORRECT.LP
        """
        return 1 #res.resutls

def rank_indexation(indexations, ranking_mode="ISCORE"):
    "Rank indexations obtained using different beam-center coordinates."
    
    best_beam_center = None
    rank_items = ["indexed_percentage", "xy_spot_position_ESD",
                  "z_spot_position_ESD", "quality_contrast","i_score"]
    rank_table = {}
    for items in rank_items:
        rank_table[items] = []
    #
    prp = " Indexed spots:          %(indexed_percentage).1f%%"
    prp += "    (%(indexed_spots)d/%(total_spots)d)\n"
    prp += " Spot prediction ESD:       %(xy_spot_position_ESD).2f "
    prp += "pixels and  %(z_spot_position_ESD).2f degrees"        
    n = 0
    i_score = []
    for indexation in indexations:
        n += 1
        print " Test indexation number: %d" % n
        print prp % indexation
        #
        origin_t = indexation["index_origin_table"]
        quality_contrast = origin_t[1][3] - origin_t[0][3]
        indexation["quality_contrast"] = quality_contrast
        indexation["i_score"] = indexation["indexed_percentage"]/(
                                   2*indexation["xy_spot_position_ESD"] + 
                                   indexation["z_spot_position_ESD"]/ \
                                      indexation["oscillation_range"])
        i_score.append(indexation["i_score"])
        #
        for items in rank_items:
            rank_table[items].append(indexation[items])
        #
        print " Contrast in the quality of indexation: ", quality_contrast
        p4, p6 = "\n\tQuality:       ", "\n\tShift (pixels):"
        p7, p8 = "\n\tBeam X (pixel):", "\n\tBeam Y (pixel):"
        p9 = "\n\tIndex Origin:  "
        for i in range(min(len(origin_t), 5)):
            p4 += "%9.2f," % (origin_t[i][3])
            p6 += "%9.1f," % (origin_t[i][4])
            p7 += "%9.1f," % (origin_t[i][5])
            p8 += "%9.1f," % (origin_t[i][6])
            p9 += "%3d%3d%3d," % tuple(origin_t[i][0:3])
        #
        print p4[:-1] + p6[:-1] + p7[:-1] + p8[:-1] +  p9[:-1] + "\n"
    #
    z_table = {}
    print "%22s: " % "Test number", " %3d"*n % tuple(range(1, n+1))
    for item in rank_table:
        isorted = rank_table[item][:]
        if item in ["indexed_percentage", "quality_contrast", "i_score"]:
            reverse = True
        else:
            reverse = False
        isorted.sort(reverse=reverse)
        #
        rank = [isorted.index(i) + 1 for i in rank_table[item]]
        print "%22s: " % item,
        print " %3d"*len(rank) % tuple(rank)
        z_table[item] = rank
    #
    z_score = []
    for idq in range(len(z_table["quality_contrast"])):
        z_score.append(z_table["quality_contrast"][idq] +
                       z_table["xy_spot_position_ESD"][idq] +
                       z_table["z_spot_position_ESD"][idq])
    print "%22s: " % "z_score",
    print " %3d"*len(z_score) % tuple(z_score)
    z_best_index = z_score.index(min(z_score))
    i_best_index = i_score.index(max(i_score))
    best_beam_center = {}
    best_beam_center["ISCORE"] = \
        indexations[i_best_index]["index_origin_table"][0][5:7]
    best_beam_center["ZSCORE"] = \
        indexations[z_best_index]["index_origin_table"][0][5:7]
    if ranking_mode == "ISCORE":
        zflag, iflag = "   ", "***"
    else:
        iflag, zflag = "   ", "***"
    _best =  best_beam_center[ranking_mode]
    fmt1 =  "%s Best  %s_score rank: %3d  for Solution #%-3d"
    fmt2 = " beamx=%7.1f beamy=%7.1f"
    print 
    print fmt1 % (iflag, "I", 1, i_best_index+1),
    print fmt2 % tuple(best_beam_center["ISCORE"])
    print fmt1 % (zflag, "Z", min(z_score), z_best_index+1),
    print fmt2 % tuple(best_beam_center["ZSCORE"])
    return {"ISCORE":i_best_index, "ZSCORE": z_best_index}

def get_beam_origin(beam_coor, beam_vec, det_parameters):
    "Calculate beam_origin from beam_coordinate."
    dist, det_x, det_y, det_z, qx, qy = det_parameters
    beamOx, beamOy, beamOz = beam_coor[0]*qx, beam_coor[1]*qy, beam_vec*det_z
    return (beamOx - beam_vec*det_x*dist/beamOz)/qx, \
           (beamOy - beam_vec*det_y*dist/beamOz)/qy


#def resolution2trustedRegion(high_res, dist, beam_center, pixel_size, npixel):
    # Usefull for the IDXREF stage. One can use the TRUSTED_REGION keyword to
    # cut unwanted spots at low or high resolution.
    # different mode can be used. Internal, external or midle.
    # Internal: set the smallest RMAX radius,
    # External: set the biggest RMAX radius and Midle...

#def write_autoPar(adpPar):
#    ""
#    link_name_to_image = "img"
#    newdir = adpPar["prefix"] + "adp_process"
#    #
#    if not os.path.exists(newdir):
#        try: os.mkdir(newdir)
#        except:
#            raise XIO.XIOError, \
#                 "STOP! Can't creat adp working directory:", newdir
#    if os.path.isdir(newdir):
#        img_dir = os.path.abspath(adpPar["img_dir"])
#        os.chdir(newdir)
#        if not os.path.exists(link_name_to_image) or \
#                           os.path.islink(link_name_to_image):
#            os.system("ln -sf %s %s" % (img_dir, link_name_to_image))
#            adpPar["img_dir"] = link_name_to_image
#        #
#        keys = adpPar.keys()
#        keys.sort()
#        paramStr = "".join(["%s = %s\n" % (k, adpPar[k]) for k in keys])
#        opWriteCl("auto.par", paramStr)
#    os.chdir("..")

def parse_spacegroup(spginp):
    "Try to interpret spg input string from command line."
    try:
        spg_int = int(spginp)
    except ValueError:
        spg_int = 0
        spginp_up = spginp.upper()
        for spgn in SPGlib:
            if spginp_up in SPGlib[spgn]:
                spg_int = spgn
                break
    if spg_int:
        spg_info = SPGlib[spg_int]
        spg_str = "  Imposed Space group:  %s,  number %d" % \
               (spg_info[1], spg_int)
    else:
        raise Exception, "\nERROR: Unrecognised space group: %s\n" % spginp
    return spg_int, spg_info, spg_str

def select_strategy(idxref_results, xds_par):
    "Interactive session to select strategy parameters."
    sel_spgn = xds_par["SPACE_GROUP_NUMBER"]
    sel_ano =  xds_par["FRIEDEL'S_LAW"]
    #print xds_par["UNIT_CELL_CONSTANTS"] 
    valid_inp = False
    bravais_to_spgs = get_BravaisToSpgs()
    # Select LATTICE
    while not valid_inp:
        defSel = 1
        if sel_spgn != 0:
            # choose the lattice solution according to the selected spg.
            i = 0
            for LAT in idxref_results["lattices_table"]:
                if LAT.fit <= LATTICE_GEOMETRIC_FIT_CUTOFF:
                    i += 1
                    if sel_spgn in bravais_to_spgs[LAT.Bravais_type]:
                        defSel = i
        selection = raw_input("\n Select a solution number [%d]: " % defSel)
        # If the selection is not compatible with the spg, set not valid
        _sel = selection.split()
        selnum = 1
        try:
            if len(_sel) == 1:
                selnum = int(_sel[0])
                valid_inp = True
            elif len(_sel) == 0:
                selnum = defSel
                valid_inp = True
            else:
                raise Exception, "Invalid selection input."
        except Exception, err:
            print "\n ERROR. ", err
    selLat = idxref_results["lattices_table"][selnum-1]
    if sel_spgn == 0:
        sel_spgn = selLat.symmetry_num
    valid_inp = False
    # Select SPACEGROUP
    print " Possible spacegroup for this lattice are:\n"
    for spgsymb in bravais_to_spgs[selLat.Bravais_type]:
        print "  %15s, number: %3d" % (SPGlib[spgsymb][1], spgsymb)
    while not valid_inp:
        selection = raw_input("\n Select the spacegroup [%s, %d]: "
                             % (SPGlib[sel_spgn][1], sel_spgn))
        _sel = selection.split()
        try:
            if len(_sel) == 1:
                sel_spgn, _spg_info, _spg_str = parse_spacegroup(_sel[0])
                # selSpgS = _spg_info[1]
                valid_inp = True
            elif len(_sel) == 0:
                valid_inp = True
            else:
                raise Exception, "Invalid selection input."
            if sel_spgn not in bravais_to_spgs[selLat.Bravais_type]:
                valid_inp = False
                msg = "Inconsistant combinaison of Bravais lattice"
                msg += " and spacegroup.\n For this Bravais Lattice"
                msg += " (%s), spacegroup should be one of these:\n\n" % \
                        (selLat.Bravais_type)
                for spgsymb in bravais_to_spgs[selLat.Bravais_type]:
                    msg += "  %15s, number: %3d\n" % \
                                 (SPGlib[spgsymb][1], spgsymb)
                raise Exception, msg
        except Exception, err:
            print "\n ERROR. ", err
    valid_inp = False
    # Select ANOMALOUS
    while not valid_inp:
        if sel_ano == "TRUE":
            txt3 = "N/y"
        else:
            txt3 = "Y/n"
        selection = raw_input(" Anomalous [%s]: " % txt3)
        try:
            _ans =  selection.strip()
            if _ans == "":
                valid_inp = True    
            elif _ans[0] in "Yy":
                xds_par["FRIEDEL'S_LAW"] = "FALSE"
                valid_inp = True
            elif _ans[0] in "Nn":
                xds_par["FRIEDEL'S_LAW"] = "TRUE"
                valid_inp = True
            else:
                raise Exception, "Invalid answer [Y/N]."
        except Exception, err:
            print "\n ERROR. ", err
    print "\n Selected  cell paramters:  ", selLat
    if sel_spgn > 2:
        selLat.idealize()
        print " Idealized cell parameters: ", selLat.prt()
        xds_par["UNIT_CELL_CONSTANTS"] = selLat.prt()
    xds_par["SPACE_GROUP_NUMBER"] = sel_spgn
    return xds_par


if __name__ == "__main__":

    import getopt

    short_opt =  "123456aAbBc:d:f:i:O:p:s:Sr:R:x:y:vw:SF"
    long_opt = ["anomal",
                "Anomal",
                "beam-x=",
                "beam-y=",
                "spg=",
                "strategy",
                "high-resolution="
                "low-resolution="
                "cell=",
                "distance",
                "reference=",
                "oscillation",
                "project",
                "beam-center-optimize-i",
                "beam-center-optimize-z",
                "xds-input=",
                "verbose",
                "wavelength=",
                "slow", "weak"]

    if len(sys.argv) == 1:
        print _usage
        sys.exit(2)
    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print _usage
        sys.exit(2)

    WARNING = ""
    VERBOSE = False
    DEBUG = False
    WEAK = False
    ANOMAL = False
    _strict_corr = False
    _beam_x = 0
    _beam_y = 0
    _spg = 0
    _strategy = False
    _res_high = 0
    _distance = 0
    _oscillation = 0
    _project = ""
    _wavelength = 0
    _res_low = 45
    _reference = False
    _beam_center_optimize = False
    _beam_center_ranking = "ZSCORE"
    _cell = ""
    _xds_input = ""
    _beam_in_mm = False
    SLOW = False
    _fast = False
    _step = 1

    for o, a in opts:
        if o == "-v":
            VERBOSE = True
        if o in ("-a", "--anomal"):
            ANOMAL = True
        if o in ("-A", "--Anomal"):
            ANOMAL = True
            _strict_corr = True
        if o[1] in "123456":
            _step = int(o[1])
        if o in ("-s", "--spg"):
            _spg, _spg_info, _spg_str = parse_spacegroup(a)
        if o in ("-i", "--xds-input"):
            _xds_input = a
        if o in ("-c", "--cell"):
            _cell = a
        if o in ("-d", "--distance"):
            _distance = float(a)
        if o in ("-f", "--reference"):
            if os.path.isfile(a):
                _reference = str(a)
            else:
                print "\n  ERROR: %s is not a regular file." % a
                print "  STOP!\n"
                sys.exit()
        if o in ("-O", "--oscillation"):
            _oscillation = float(a)
        if o in ("-p", "--project"):
            _project = str(a)
        if o in ("-S", "--strategy"):
            _strategy = True
        if o in ("-w", "--wavelength"):
            _wavelength = float(a)
        if o in ("-r", "--high-resolution"):
            _res_high = float(a)
        if o in ("-R", "--low-resolution"):
            _res_low = float(a)
        if o in ("-x", "--beam_x"):
            if "mm" in a:
                _beam_in_mm = True
                a = a.replace("mm","")
            _beam_x = float(a)
        if o in ("-y", "--beam_y"):
            if "mm" in a:
                _beam_in_mm = True
                a = a.replace("mm","")
            _beam_y = float(a)
        if o in ("-b", "--beam-center-optimize-i"):
            _beam_center_optimize = True
            _beam_center_ranking = "ISCORE"
        if o in ("-B", "--beam-center-optimize-z"):
            _beam_center_optimize = True
            _beam_center_ranking = "ZSCORE"
        if o in ("--slow"):
            SLOW = True
        if o in ("--weak"):
            WEAK = True
        if o in ("-h", "--help"):
            print _usage
            sys.exit()

    if not os.path.isfile(inputf[0]):
        print "\nERROR. Can't open file: %s\n" % inputf[0] 
        sys.exit(2)
    _coll = XIO.Collect(inputf[0])
    if not _project:
        newDir = "xds_process_" + _coll.prefix
    else:
        newDir = "xds_process_" + _project
    #
    _linkimages = False
    if not _coll.isContinuous(inputf):
        _linkimages = True
        link_dir = "img_links"
        inputf = makeXDSImageLinks(inputf,
                                    os.path.join(newDir,link_dir),
                                    "image")
        #collect.setDirectory(link_dir)
        #collect.prefix = prefix

    try:
        collect = XIO.Collect(inputf)
        collect.interpretImage()
        collect.image.info()
        collect.lookup_imageRanges(forceCheck=False)

    except XIO.XIOError, _mess:
        print _mess
        print "\nError: Can't access to file(s) %s.\nStop." % inputf
        sys.exit(2)

    imgDir = collect.directory
    newPar = collect.export("xds")

   # In case no beam origin is defined, take the detector center.
    if newPar["ORGX"] == 0:
        newPar["ORGX"] = newPar["NX"]/2.
    if newPar["ORGY"] == 0:
        newPar["ORGY"] = newPar["NY"]/2.

    #write_autoPar(collect.export("adp"))

    # Update some default values.
    # Defined by XDS.export xds:
    #    newPar["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = 9
    newPar["STRONG_PIXEL"] = 7

    newrun = XDS()

    if _beam_in_mm:
        _beam_x = _beam_x / newPar["QX"]
        _beam_y = _beam_y / newPar["QY"]
    if ANOMAL:
        newPar["FRIEDEL'S_LAW"] = "FALSE"
    else:
        newPar["FRIEDEL'S_LAW"] = "TRUE"
    if _strict_corr:
        newPar["STRICT_ABSORPTION_CORRECTION"] = "TRUE"
    if _beam_x:
        newPar["ORGX"] = _beam_x
    if _beam_y:
        newPar["ORGY"] = _beam_y
    if _spg and _cell:
        newPar["SPACE_GROUP_NUMBER"] = _spg
        newPar["UNIT_CELL_CONSTANTS"] = _cell
    elif _spg and not _cell:
        WARNING = "  WARNING: Spacegroup is defined but not cell."
        WARNING += " Waiting for indexation for setting cell."
    elif _cell and not _spg:
        WARNING = "  WARNING: Cell is defined but not spacegroup,"
        WARNING += " setting spacegroup to P1."
        newPar["SPACE_GROUP_NUMBER"] = 1
        newPar["UNIT_CELL_CONSTANTS"] = _cell
    if _distance:
        newPar["DETECTOR_DISTANCE"] = _distance
    if _reference:
        newPar["REFERENCE_DATA_SET"] = "../"+_reference
    if _oscillation:
        newPar["OSCILLATION_RANGE"] = _oscillation
    if _wavelength:
        newPar["X_RAY_WAVELENGTH"] = _wavelength
    if _xds_input:
        newPar.update(xdsInp2Param(inp_str=_xds_input))
    if _res_high or _res_low != 35:
        newPar["INCLUDE_RESOLUTION_RANGE"] = _res_low, _res_high

    if _linkimages:
        collect.setDirectory(link_dir)
    else:
        collect.setDirectory(newrun.link_name_to_image)

    newPar["NAME_TEMPLATE_OF_DATA_FRAMES"] = collect.xdsTemplate

    newPar["DELPHI"] = 8 * newPar["OSCILLATION_RANGE"]
    newrun.inpParam.mix(xdsInp2Param(inp_str=xdsinp_base))
    newrun.inpParam.mix(newPar)
    newrun.set_collect_dir(os.path.abspath(imgDir))
    newrun.run_dir = newDir

    if SLOW:
        newrun.inpParam["DELPHI"] = 12 * newPar["OSCILLATION_RANGE"]
        newrun.mode.append("slow")
    if WEAK:
        newrun.mode.append("weak")

    print _fmt_hello % vars(newrun.inpParam)
    if WARNING:
        print WARNING
    if _spg:
        print _spg_str

    #newrun.run()
    R1 = R2 = R3 = R4 = R5 = None
    if _step > 1:
        print "\n Starting at step: %d (%s)\n" % (_step, STEPS[_step-1])
    if _step <= 1:
        R1 = newrun.run_init()
    if _step <= 2:
        R2 = newrun.run_colspot()
    if _step <= 3:
        if newrun.mode == "slow" and _res_high:
            print "   Applying a SPOT RESOLUTION CUTOFF: %.2f A" % _res_high
            #newrun.spots_resolution_cutoff(_res_high)
        R3 = newrun.run_idxref(_beam_center_optimize, _beam_center_ranking)
    if R3:
        i = 0
        _selected_cell = []
        print "    TABLE OF POSSIBLE LATTICES:\n"
        print " num  Symm  quality  mult     a      b      c",
        print "   alpha   beta  gamma"
        print " "+"-"*67
        fmt_lat = "%3d) %5s %7.2f  %4d  %s"
        for LAT in R3["lattices_table"]:
            if LAT.fit <= LATTICE_GEOMETRIC_FIT_CUTOFF:
                i += 1
                print fmt_lat % (i, LAT.symmetry_str1,
                            LAT.fit, LAT.multiplicity, LAT)
        # If not multiple possible solutions (like P2, or P1...) try to define
        # unitcell from spacegroup.
        #if _spg and not _cell:
        if (len(collect.imageRanges) > 1) or _strategy:
            newrun.run_xplan(collect.imageRanges, ridx=R3)
    if _step <= 4:
        R4 = newrun.run_integrate(collect.imageRanges)
    if _step <= 5:
        R5 = newrun.run_correct()

    #import time 
    #t0 = time.time()
    #newrun.spots_resolution_cutoff(3.,"S1.XDS")
    #t1 = time.time()
    #print "New: time = %10.3f" % (t1-t0)


