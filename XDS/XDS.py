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
 TODO-1: Add in the CORRECT summary the Rmrgd-F, and total
         overloaded refl.
 TODO-2: Start multiple COLSPOT with different thresholds+ multiple IDXREF...
 TODO-3: Generating plots !
"""

__version__ = "0.4.6"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "02-01-2010"
__copyright__ = "Copyright (c) 2006-2010 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import os
import sys
import re
import math

if sys.version_info <= (2, 4, 0):
    from popen2 import Popen3
else:
    from subprocess import Popen, PIPE

from XOconv.pycgtypes import mat3
from XOconv.pycgtypes import vec3
from XOconv.XOconv import reciprocal, UB_to_cellParam, BusingLevy, \
                          volum, r2d, cosd, sind, ex, ey, ez

from pointless import pointless, is_pointless_installed
from xupy import XParam, xdsInp2Param, opWriteCl, \
                 saveLastVersion, LP_names, xdsinp_base, \
                 SPGlib, Lattice, resum_scaling, \
                 get_BravaisToSpgs, get_number_of_processors
import XIO

PROGNAME = os.path.split(sys.argv[0])[1]
USAGE = """
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

""" % PROGNAME

FMT_HELLO = """
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

FMT_FINAL_STAT = """
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

FMT_ABSENCES = "      Systematic absente reflection measured %(AbsNum)6d \
 with <Iabs>/<I> =  %(AbsIav).1f%%\n"

FMT_ANOMAL = """
      Anomalous pairs measured   %(anoNum)10d
      SigAno                     %(anoSig)10.3f    (%(anoSigL).3f)
      Anomalous Correlation      %(anoCorr)10.1f%%   (%(anoCorrL).1f%%)
"""

STEPMARK = re.compile(r"^( [*]{5} (\w{4,}) [*]{5}  )")
INTEGRATE_STEP_RE = re.compile(r" PROCESSING OF IMAGES ")
INTEGRATE_MOSAICITY_RE = re.compile(r"CRYSTAL MOSAICITY \(DEGREES\)")
INTEGRATE_STRONG_RE = re.compile(r"REFLECTIONS ACCEPTED FOR REFINEMENT")
rrf, rri = r"[\ ]+([0-9\.]+)", r"[\ ]+([\d]+) "
SCALE_RE = re.compile(r" "+rri+r"  (\d)"+rrf+r"  ....... "+4*rri+2*rrf)

XDS_HOME = os.getenv('XDS')

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

def make_xds_image_links(imagename_list, dir_name="img_links",
                       prefix="image", start_num=1):
    """All image names in the imagename_list are supposed to be part
    of one continous sequence of collected images.
    Todo:
        - How to safely modulate PhiStart outside the [-180,180] range ?
    """
    link_list = []
    if dir_name not in os.listdir("."):
        try:
            _mkdir(dir_name)
        except Exception, err:
            print "Error\n", err
            sys.exit(0)
    #
    dir_name = os.path.abspath(dir_name)
    collect_im = {}
    osc_ranges = []
    for _image in imagename_list:
        image = XIO.Image(_image)
        if VERBOSE:
            print _image
        # How to safely modulate PhiStart outside the [-180,180] range?
        if VERBOSE:
            print "\tPhiStart %8.2f" % image.header['PhiStart']
        if VERBOSE:
            print "\tPhiWidth %8.2f" % image.header['PhiWidth']
        collect_im[image.header['PhiStart']] = _image
        osc_ranges.append(image.header['PhiWidth'])

    if max(osc_ranges) != min(osc_ranges):
        print "Error. Image list contains different oscillation range!"
        sys.exit(0)
    #
    osc_starts = collect_im.keys()
    osc_starts.sort()
    for _osc in osc_starts:
        _num =  start_num+ (_osc-osc_starts[0])/osc_ranges[0]
        link_name = os.path.join(dir_name, prefix+"_%04.0f.img" % _num)
        if os.path.lexists(link_name) and os.path.islink(link_name):
            if VERBOSE:
                print "==> Removing existing link: %s" % link_name
            os.remove(link_name)
        os.symlink(os.path.abspath(collect_im[_osc]), link_name)
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
            if _err_msg.count(" CANNOT READ IMAGE "):
                _err_type = "Some images connot be read"
                _err_level = "WARNING"
            # IDXREF ERROR Messages:
            elif _err_msg.count("INSUFFICIENT PERCENTAGE (<"):
                _err_type = "IDXREF. Percentage of indexed"
                _err_type += " reflections bellow 70%.\n"
                _err_level = "WARNING"
            elif _err_msg.count("INSUFFICIENT NUMBER OF ACCEPTED SPOTS."):
                _err_type = "IDXREF. INSUFFICIENT NUMBER OF ACCEPTED SPOTS."
                _err_level = "FATAL"
            elif _err_msg.count("CANNOT INDEX REFLECTIONS"):
                _err_type = "IDXREF. CANNOT INDEX REFLECTIONS."
                _err_level = "FATAL"
            else:
                print "\n %s \n" % (self.lp[_err:-1])
                sys.exit()
        if _err_level in ("FATAL", "ERROR") and raiseErrors:
            raise XDSExecError, (_err_level, _err_type)

        if self.verbose and _err != -1:
            print "\n %s in %s" % (_err_level, _err_type)

        if full_filename.count("INIT.LP"):
            self.parse_init()
        elif full_filename.count("COLSPOT.LP"):
            self.parse_colspot()
        elif full_filename.count("IDXREF.LP"):
            self.parse_idxref()
        elif full_filename.count("XPLAN.LP"):
            self.parse_xplan()
        elif full_filename.count("INTEGRATE.LP"):
            self.parse_integrate()
        elif full_filename.count("CORRECT.LP"):
            self.parse_correct()
        else:
            raise IOError, "Don't know how to parse file: %s" % full_filename

    def get_par(self, match, limit=75, func=None,
                      multi_line=False, start=0, before=False):
        "Extract parameters from XDS .LP lines."
        try:
            if before:
                limit = start
                start = self.lp.index(match)-start
            else:
                start = self.lp.index(match, start) + len(match)
        except Exception, err:
            raise err
        if multi_line:
            _raw = self.lp[start:start+limit].split()
        else:
            _raw = self.lp[start:start+limit].splitlines()[0].split()
        if not func:
            for var_type in (int, float, str):
                try:
                    var_type(_raw[0])
                    func = var_type
                except ValueError:
                    pass
                if func:
                    break
        if not func:
            raise ValueError, "get_par function can't process value '%s'" % _raw
        pars = map(func, _raw)
        if len(pars) == 1:
            return pars[0]
        else: return pars

    def _get_lattices_table(self):
        "Extract lattice table"
        st1 = self.lp.index("LATTICE-  BRAVAIS-   QUALITY")
        _table = self.lp[st1:st1+6000].splitlines()[3:47]
        return map(unpack_latticefit2, _table)

    def _get_index_origins_table(self):
        "Extract origin table"
        st0 = self.lp.index(" DL\n  ORIGIN\n")+14
        st1 = self.lp.index(" SELECTED:     INDEX_ORIGIN=")-2
        return map(lambda s: \
                  map(float, s.split()), self.lp[st0:st1].splitlines())

    def parse_init(self):
        "Parse INIT.LP"
        rdi, gpa = self.results, self.get_par
        #
        rdi["background_range"] = gpa("BACKGROUND_RANGE=")
        rdi["mean_gain"] = gpa("MEAN GAIN VALUE")
        rdi["min_gain"] = gpa("MINIMUM GAIN VALUE IN TABLE")
        rdi["max_gain"] = gpa("MAXIMUM GAIN VALUE IN TABLE")
        rdi["mean_background"] = gpa("BACKGROUND COUNTS IN A DATA IMAGE PIXEL")
        #
        prp =  "  Looking at images %(background_range)s\n"
        prp += "  Mean Gain:        %(mean_gain).1f\n"
        prp += "  Min table gain:   %(min_gain).2f\n"
        prp += "  Max table gain:   %(max_gain).2f\n"
        prp += "  Mean Backgroud:   %(mean_background).1f\n"
        if self.verbose:
            print prp % rdi
        return rdi, prp

    def parse_colspot(self):
        "Parse COLSPOT.LP"
        rdi, gpa = self.results, self.get_par
        #
        rdi["strong_pixels"] = gpa("EXTRACTED FROM IMAGES")
        rdi["weak_spots_ignored"] = gpa("WEAK SPOTS OMITTED")
        rdi["out_of_center_spots"] = gpa("SPOT MAXIMUM OUT OF CENTER")
        rdi["spot_number"] = self.get_spot_number()
        rdi["time"] = gpa("elapsed wall-clock time", 11)

        prp =  "  Number of spots found:    %(spot_number)10d\n"
        prp += "  Out of center rejected:   %(out_of_center_spots)10d\n"
        prp += "  Weak spots rejected:      %(weak_spots_ignored)10d\n"
        prp += "  Number of spots accepted: %(spot_number)10d\n"
        if self.verbose:
            print prp % rdi
        return rdi, prp

    def parse_idxref(self):
        "Parse IDXREF.LP"
        rdi, gpa = self.results, self.get_par
        #
        rexp1 = r".* (\d+) OUT OF\ +(\d+) SPOTS INDEXED\..*"
        rexp2 = r".* QX=\ +([\d|\.]+)\ +QY=\ +([\d|\.]+)"
        #if "! ERROR !" in self.lp:
        #    raise XDSLogParserException, "Error while parsing XDS logfile"
        nis, nts = map(int, re.match(rexp1, self.lp, re.DOTALL).groups())
        qx, qy = map(float, re.match(rexp2, self.lp, re.DOTALL).groups())
        meanPixel = (qx+qy)/2
        rdi["indexed_spots"] = nis
        rdi["total_spots"] = nts
        rdi["indexed_percentage"] = 100.*nis/nts
        #
        st0 = self.lp.index("START OF INTEGRATION *****")
        st1 = "STANDARD DEVIATION OF SPOT    POSITION (PIXELS)"
        st2 = "STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)"
        st3 = "UNIT CELL PARAMETERS"
        st4 = "SPACE GROUP NUMBER"
        st5 = "COORDINATES (PIXELS) OF DIRECT BEAM"
        #
        rdi["oscillation_range"] = gpa("OSCILLATION_RANGE=")
        rdi["xy_spot_position_ESD"] = gpa(st1, start=st0)
        rdi["z_spot_position_ESD"] = gpa(st2, start=st0)
        rdi["index_origin_table"] = self._get_index_origins_table()
        rdi["lattices_table"] = self._get_lattices_table()
        rdi["refined_cell"] = gpa(st3, start=st0)
        rdi["refined_cell_str"] = 6*"%.2f " % \
                                           tuple(rdi["refined_cell"])
        rdi["space_group_number"] = gpa(st4, start=st0)
        rdi["direct_beam_pixels"] = gpa(st5, start=st0)
        rdi["direct_beam_mm"] = rdi["direct_beam_pixels"][0]*qx, \
                                           rdi["direct_beam_pixels"][1]*qy
        rdi["bmx"], rdi["bmy"] = rdi["direct_beam_mm"]
        rdi["bpx"], rdi["bpy"] = rdi["direct_beam_pixels"]

        origin_t = rdi["index_origin_table"]
        origin_n = len(origin_t)
        quality_t = [x[3] for x in origin_t if x[3] < 2.]
        #rdi["index_score"] = reduce(lambda a,b: a+b, quality_t)/len(quality_t)
        max_ot = min(origin_n, 5)
        rdi["shift_pixel"] = origin_t[0][4]
        rdi["shift_mm"] = origin_t[0][4]*meanPixel
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
        #prp += " Index origin score: %.2f\n" % (rdi["index_score"])
        if self.verbose == 1:
            print (prp + prp2) % rdi
        elif self.verbose == 2:
            print prp % rdi
        return rdi, prp

    def parse_integrate(self):
        "Parse INTEGRATE.LP"
        rdi, gpa = self.results, self.get_par
        rdi["reflections"] = gpa("REFLECTIONS SAVED ON FILE",
                                  start=9, func=int, before=True)
        rdi["divergence"] = gpa("BEAM_DIVERGENCE_E.S.D.= ")
        rdi["mosaicity"] = gpa("REFLECTING_RANGE_E.S.D.= ")
        prp =  "\n  Number of reflection integrated:      %(reflections)d\n"
        prp += "  Estimated divergence:                 %(divergence).3f\n"
        prp += "  Estimated mosaicity:                  %(mosaicity).3f\n"
        if self.verbose:
            print prp % rdi
        return rdi, prp

    def parse_xplan(self):
        "Parse XPLAN.LP"
        rdi, gpa = self.results, self.get_par
        rdi["spacegroup"] = gpa("SPACE_GROUP_NUMBER=")
        rdi["unitcell"] = 6*" %8.2f" % tuple(gpa("UNIT_CELL_CONSTANTS="))
        rdi["friedels_law"] = gpa("FRIEDEL'S_LAW=")[0]
        st0 = self.lp.index(72*"*")
        st1 = self.lp.index(72*"*", st0+72)
        st2 = self.lp.index(72*"*", st1+72)
        #
        prp =  "  Friedel's law:   %(friedels_law)s\n"
        prp += "  Spacegroup:      %(spacegroup)d\n"
        prp += "  Unitcell:        %(unitcell)s\n"
        if self.verbose:
            print prp % rdi
            print
            print self.lp[st0:st2]
        return rdi, prp

    def parse_correct(self):
        "Parse CORRECT.LP"
        rdi, gpa = self.results, self.get_par

        sp1 = self.lp.index("b              INPUT DATA SET")
        sp2 = self.lp.index("  INTEGRATE.HKL   ", sp1)
        K1s, K2s = map(float, self.lp[sp1+30: sp2].split())
        print "  Variance estimate scaling (K1, K2): %6.3f, %.3e" % \
                                                   (4*K1s, (K2s/4+0.0001))
        rdi["IoverSigmaAsympt"] =  1/((K1s*(K2s+0.0004))**0.5)
        print "  Upper theoritical limit of I/sigma: %8.3f" % \
                                                   rdi["IoverSigmaAsympt"]
        rdi["RMSd_spotPosition"] = gpa("SPOT    POSITION (PIXELS)")
        rdi["RMSd_spindlePosition"] = gpa("SPINDLE POSITION (DEGREES)")
        rdi["Mosaicity"] = gpa("CRYSTAL MOSAICITY (DEGREES)")
        r = gpa(" "+"-"*74+"\n")
        rdi["I_sigma"], rdi["Rsym"] = r[2], r[4]
        rdi["Compared"], rdi["Total"] = r[6], r[7]
        ### Select Diffraction range.
        sp1 = self.lp.index("RESOLUTION RANGE  I/Sigma")
        sp2 = self.lp.index(10*"-", sp1)
        _table = self.lp[sp1:sp2].splitlines()[3:-1]
        _table = [ map(float, l.split()[1:3]) for l in _table ]
        rdi["HighResCutoff"] = self.get_proper_resolition_range(_table)
        prp = ""
        if rdi["Mosaicity"]:
            prp += "  RMSd spot position:     %(RMSd_spotPosition)9.2f pix,"
            prp += " %(RMSd_spindlePosition)6.2f deg.\n"
            prp += "  Refined Mosaicity:                %(Mosaicity)9.2f\n\n"
        prp += "  Rsym:                             %(Rsym)9.1f\n"
        prp += "  I/sigma:                          %(I_sigma)9.1f\n"
        if rdi["HighResCutoff"]:
            prp += "  Suggested high resolution cutoff: %(HighResCutoff)9.2f\n"
        prp += "  Compared reflections:                 %(Compared)d\n"
        prp += "  Total number of measures:             %(Total)d\n"
        if self.verbose:
            print prp % rdi
        return rdi, prp

    def get_proper_resolition_range(self, res_table):
        "High res is selected when at least 3 values of I/sigma are below 1."
        high_n, high_hit = [], None
        for res, IoS in res_table:
            if IoS < 1.:
                high_n.append(res)
                if not high_hit and len(high_n) == 3:
                    high_hit = high_n[0]
            else:
                high_n = []
            #print "%8.3f  %8.3f  %s" % (res, IoS, IoS >= 1.)
        if not high_hit and len(high_n) >= 1:
            high_hit = high_n[0]
        #print "Suggested high resolution cut-off: %.2f" % high_hit
        return high_hit

    def get_spot_number(self):
        "Read the number of spot directly from SPOT.XDS"
        _execstr = "wc -l %s/SPOT.XDS" % self.run_dir
        if sys.version_info <= (2, 4, 0):
            spot_file = os.popen(_execstr)
            wc_out = spot_file.read()
            spot_file.close()
        else:
            wc_out = Popen([_execstr], stdout=PIPE, shell=True).communicate()[0]
        return int(wc_out.split()[0])

MIN_SPOT_NUMBER = 200
LATTICE_GEOMETRIC_FIT_CUTOFF = 50
FRAMES_PER_COLSPOT_SEQUENCE = 16 # number of frames per sequence in COLSPOT.
JOB_STEPS = "INIT", "COLSPOT", "IDXREF", "INTEGRATE", "CORRECT"
SPOTFILENAME = "SPOT.XDS"

class XDS:
    "Main class for runing xds step by step."

    def __init__(self, obj=None, link_to_images=True):
        """Constructor for the Param classes from file or string."""
        #
        self.link_to_images = link_to_images
        self.__cancelled = 0
        self.__lastOutp = 0
        self.mode = []
        if XDS_HOME:
            self.__execfile = os.path.join(XDS_HOME,"xds_par")
        else:
            self.__execfile = "xds_par"
        self.running = 0
        self.outp = []
        self.run_dir = "."
        self.status = None
        self.inpParam = XParam()
        self.collect_dir = "./"
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
            self.collect_dir = dirname
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

    def cancel(self):
        "Cancel the job."
        self.__cancelled = 1

    def get_outp(self):
        "Collect the latest output."
        if not self.__cancelled:
            nLine = len(self.outp)
            diff = nLine - self.__lastOutp
            self.__lastOutp = nLine
            if diff:
                return "".join(self.outp[-diff:])[:-1]
            else:
                return ""

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
            self.run_dir = "."
        result = 0
        if self.run_dir:
            if not os.path.exists(self.run_dir):
                try:
                    os.mkdir(self.run_dir)
                except OSError, err:
                    raise XIO.XIOError, \
                     ("\nSTOP! Can't create xds working directory: %s\n" % \
                                                              self.run_dir)
            if os.path.isdir(self.run_dir):
                os.chdir(self.run_dir)
                if self.link_to_images:  
                    if not os.path.exists(self.link_name_to_image):
                        os.system("ln -sf %s %s" % (self.collect_dir, \
                                                    self.link_name_to_image))
                        #os.system("ln -sf .. %s" % (self.link_name_to_image))
                    #else:
                    #    raise XIO.XIOError, \
                    #     "STOP! Can't creat link %s in working directory: %s" \
                    #     % (self.link_name_to_image, self.run_dir)
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
                lines = xdsProcess.fromchild.readline()
            else:
                lines = xdsProcess.stdout.readline()
                #lines = xdsProcess.communicate()
            # ilines parsing of stdout
            if self.step_name == "INTEGRATE":
                if _init_parse:
                    print "    Processing    Mean #Strong  ",
                    print "Estimated   Overloaded"
                    print "    Image Range   refl./image   ",
                    print "Mosaicity   reflections\n"
                    table_int = []
                    _init_parse = False
                if INTEGRATE_STEP_RE.search(lines):
                    print lines[44:50]+" - "+lines[56:-1],
                    nimages = int(lines[56:-1]) - int(lines[44:50]) + 1
                elif INTEGRATE_STRONG_RE.search(lines):
                    print "%11.0f" % (float(lines.split()[0])/nimages),
                elif INTEGRATE_MOSAICITY_RE.search(lines):
                    print " %11.3f" % float(lines.split()[3]),
                    print " %11d" %  overloaded_spots
                    overloaded_spots = 0
                hit = SCALE_RE.search(lines)
                if hit:
                    table_int = hit.groups()
                    overloaded_spots += int(hit.groups()[3])
            sm = STEPMARK.match(lines)
            if sm:
                self.step += 1
                self.step_name = sm.group(2)
                #if VERBOSE:
                if verbose:
                    print "\n --->  Running job: %20s\n" % self.step_name
            if lines:
                self.outp.append(lines)
        self.step += 1
        self.step_name = "FINISHED"
        if self.__cancelled:
            result = -1
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
                   0.5/sin(atan2(((rx*(float(s[:10])  -xo))**2 +
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
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        #self.inpParam["TRUSTED_REGION"] = [0, 1.20]
        self.inpParam["JOB"] = "XYCORR", "INIT"
        i1, i2 = self.inpParam["DATA_RANGE"]
        #if "slow" in self.mode:
        if SLOW:
            self.inpParam["BACKGROUND_RANGE"] =  i1, min(i2, i1+11)
        else:
            self.inpParam["BACKGROUND_RANGE"] =  i1, min(i2, i1+3)
        self.run(rsave=True)
        res = XDSLogParser("INIT.LP", run_dir=self.run_dir, verbose=1)
        return res.results

    def run_colspot(self):
        "Runs the COLSPOT step."
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        self.inpParam["JOB"] = "COLSPOT",
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 1
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = NUMBER_OF_PROCESSORS
        _trial = 0

        frames_per_colspot_sequence = FRAMES_PER_COLSPOT_SEQUENCE
        if "slow" in self.mode:
            frames_per_colspot_sequence = 32
        elif "fast" in self.mode:
            frames_per_colspot_sequence = 4
        else:
            frames_per_colspot_sequence = 16
        if "weak" in self.mode:
            self.inpParam["STRONG_PIXEL"] = 4.5
            frames_per_colspot_sequence = 32
        # Selecting spot range(s),
        # self.inpParam["SPOT_RANGE"] is set to Collect.imageRanges by the
        # xds export function XIO
        cfo = XIO.Collect("foo_001.bar")
        cfo.imageNumbers = cfo._ranges_to_sequence(self.inpParam["SPOT_RANGE"])
        #
        min_fn, max_fn = self.inpParam["DATA_RANGE"] 
        dPhi = self.inpParam["OSCILLATION_RANGE"]
        _fpcs = frames_per_colspot_sequence
        _2fpcs = 1 + 2 * frames_per_colspot_sequence

        if (max_fn - min_fn + 1) >= _2fpcs:
            # use two range ex: i-i+2, f-2,f
            # with f at maximum 90 degre distance
            max_frame = min(max_fn, min_fn + int(89./dPhi + _fpcs))
            spot_ranges = ((min_fn, min_fn + _fpcs - 1),
                          (max_frame - _fpcs + 1, max_frame))
        else:
            spot_ranges = (min_fn, min(min_fn + _2fpcs - 1, max_fn)),
        # Restrict to matching collected images...
        self.inpParam["SPOT_RANGE"] = cfo.lookup_imageRanges(False, \
                                              mask_range=spot_ranges)
        self.run(rsave=True)
        _rs = "  Image range(s) for spot collection: "
        for sub_range in self.inpParam["SPOT_RANGE"]:
            _rs += ("  [%d - %d]," % tuple(sub_range))
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
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        self.inpParam["JOB"] = "IDXREF",
        # this prevent bad spot to be included.
        saved_trusted_region = self.inpParam["TRUSTED_REGION"]
        if saved_trusted_region[1] > 0.98:
            self.inpParam["TRUSTED_REGION"] = [0, 0.98]
        self.run(rsave=True)
        try:
            res = XDSLogParser("IDXREF.LP", run_dir=self.run_dir, verbose=1)
        except XDSExecError, err:
            print " !!! ERROR in", err[1], "\n"
            if err[0] == "FATAL":
                sys.exit()

        RD = res.results
        qx, qy = self.inpParam["QX"], self.inpParam["QY"]
        dist = self.inpParam["DETECTOR_DISTANCE"]
        det_x = vec3(self.inpParam["DIRECTION_OF_DETECTOR_X-AXIS"])
        det_y = vec3(self.inpParam["DIRECTION_OF_DETECTOR_Y-AXIS"])
        det_z = det_x.cross(det_y)
        det_params = dist, det_x, det_y, det_z, qx, qy

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
        # Set back the Trusted_region to larger values.
        self.inpParam["TRUSTED_REGION"] = saved_trusted_region
        return res.results

    def check_fileout(self, fileout):
        "Checking normal terminaison."
        if not os.path.exists(os.path.join(self.run_dir, fileout)):
            err = "Abnormal terminaison. Can't locate file: '%s'" % fileout
            print err
            raise Exception(err)

    def run_xplan(self, ridx=None):
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        "Running the strategy."
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = NUMBER_OF_PROCESSORS
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 1

        select_strategy(ridx, self.inpParam)
        print "\n Starting strategy calculation."
        self.inpParam["JOB"] = "IDXREF",
        self.run(rsave=True)
        res = XDSLogParser("IDXREF.LP", run_dir=self.run_dir, verbose=2)
        # Select just the internal circle of the detector.
        self.inpParam["JOB"] = "DEFPIX", "XPLAN"
        self.run(rsave=True)
        res =  XDSLogParser("XPLAN.LP", run_dir=self.run_dir, verbose=1)      
        return res.results

    def run_integrate(self, image_ranges):
        "Running INTEGRATE."
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 8
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 1
        if "slow" in self.mode:
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

    def run_pre_correct(self):
        """Runs a first pass of CORRECT to evaluate high_res and
           point group.
        """
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        # run pointless on INTEGRATE.HKL        
        if not is_pointless_installed():
            print "!!  Warning. Pointless program doesn't seems to be installed."
            print "  -> Skipping pointless analysis."
            likely_spg = [["P1", 0],]
        else:
            print "     Pointless analysis on the INTEGRATE.HKL file"
            print "     "+44*"="
            try:
                likely_spg, new_cell = pointless(dir_name=self.run_dir,
                                                 hklinp="INTEGRATE.HKL")
            except:
                print "  -> ERROR. While running Pointless... skiping this step."
                likely_spg = [["P1", 0],]

        self.inpParam["JOB"] = "CORRECT",
        if not SPG:
            # run first CORRECT in P1 with the cell used for integration.
            # read the cell parameters from the XPARM.XDS file
            self.inpParam["SPACE_GROUP_NUMBER"] = 1
            xparm_file = os.path.join(self.run_dir, "XPARM.XDS")
            #xparm_file = "XPARM.XDS"
            self.inpParam["UNIT_CELL_CONSTANTS"] = \
               map(float, (open(xparm_file,'r').readlines()[7]).split()[1:])
        # run CORRECT
        self.run(rsave=True)
        res = XDSLogParser("CORRECT.LP", run_dir=self.run_dir, verbose=1)
        L, H = self.inpParam["INCLUDE_RESOLUTION_RANGE"]
        newH = res.results["HighResCutoff"]
        if newH > H and not RES_HIGH:
            H = newH
        if SPG:
            spg_choosen = SPG
        else:
            spg_choosen = likely_spg[0][1]
            lattice = Lattice(new_cell, symmetry=spg_choosen)
            lattice.idealize()
            #reidx_mat = likely_spg[0][-1]
            #new_cell = new_reidx_cell(self.inpParam["UNIT_CELL_CONSTANTS"],
            #                          reidx_mat)
            self.inpParam["UNIT_CELL_CONSTANTS"] = lattice.cell
        return (L, H), spg_choosen

    def run_correct(self, res_cut=(1000, 0), spg_num=0):
        "Runs the last step: CORRECT"
        if res_cut[1]:
            print "   ->  New high resolution limit: %.2f Å" % res_cut[1]
            self.inpParam["INCLUDE_RESOLUTION_RANGE"] = res_cut
        if spg_num:
            print "   ->  Usging spacegroup: %s  #%d" % \
                                   (SPGlib[spg_num][1], spg_num)
        lattice = Lattice(self.inpParam["UNIT_CELL_CONSTANTS"],
                          symmetry=spg_num)
        lattice.idealize()
        self.inpParam["UNIT_CELL_CONSTANTS"] = lattice.cell
        self.inpParam["JOB"] = "CORRECT",
        self.inpParam["SPACE_GROUP_NUMBER"] = spg_num
        self.run(rsave=True)
        res = XDSLogParser("CORRECT.LP", run_dir=self.run_dir, verbose=1)
        s = resum_scaling(lpf=os.path.join(self.run_dir,"CORRECT.LP"))
        s["image_start"], s["image_last"] = self.inpParam["DATA_RANGE"]
        if not s:
            print "\nERROR while running CORRECT"
            sys.exit()
        print FMT_FINAL_STAT % vars(s)
        if s.absent:
            print FMT_ABSENCES % vars(s)
        if self.inpParam["FRIEDEL'S_LAW"] == "FALSE":
            print FMT_ANOMAL % vars(s)

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

    prp = " Indexed spots:          %(indexed_percentage).1f%%"
    prp += "    (%(indexed_spots)d/%(total_spots)d)\n"
    prp += " Spot prediction ESD:       %(xy_spot_position_ESD).2f "
    prp += "pixels and  %(z_spot_position_ESD).2f degrees"        
    nind = 0
    i_score = []
    for indexation in indexations:
        nind += 1
        print " Test indexation number: %d" % nind
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
        pp4, pp6 = "\n\tQuality:       ", "\n\tShift (pixels):"
        pp7, pp8 = "\n\tBeam X (pixel):", "\n\tBeam Y (pixel):"
        pp9 = "\n\tIndex Origin:  "
        for i in range(min(len(origin_t), 5)):
            pp4 += "%9.2f," % (origin_t[i][3])
            pp6 += "%9.1f," % (origin_t[i][4])
            pp7 += "%9.1f," % (origin_t[i][5])
            pp8 += "%9.1f," % (origin_t[i][6])
            pp9 += "%3d%3d%3d," % tuple(origin_t[i][0:3])
        #
        print pp4[:-1] + pp6[:-1] + pp7[:-1] + pp8[:-1] +  pp9[:-1] + "\n"
    #
    z_table = {}
    print "%22s: " % "Test number", " %3d"*nind % tuple(range(1, nind+1))
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
        
def new_reidx_cell(init_cell, reidx_mat):
    "Applies the reindexing card to initial cell parameters and return a new cell"
    UB = BusingLevy(reciprocal(init_cell))
    REIDX = mat3(reidx_mat)
    return reciprocal(UB_to_cellParam(UB*REIDX))
    
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
        print USAGE
        sys.exit(2)
    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print USAGE
        sys.exit(2)

    NUMBER_OF_PROCESSORS = get_number_of_processors()
    WARNING = ""
    VERBOSE = False
    DEBUG = False
    WEAK = False
    ANOMAL = False
    STRICT_CORR = False
    BEAM_X = 0
    BEAM_Y = 0
    SPG = 0
    STRATEGY = False
    RES_HIGH = 0
    _distance = 0
    _oscillation = 0
    _project = ""
    _wavelength = 0
    RES_LOW = 50
    _reference = False
    _beam_center_optimize = False
    _beam_center_ranking = "ZSCORE"
    _cell = ""
    XDS_INPUT = ""
    _beam_in_mm = False
    SLOW = False
    _fast = False
    STEP = 1

    for o, a in opts:
        if o == "-v":
            VERBOSE = True
        if o in ("-a", "--anomal"):
            ANOMAL = True
        if o in ("-A", "--Anomal"):
            ANOMAL = True
            STRICT_CORR = True
        if o[1] in "123456":
            STEP = int(o[1])
        if o in ("-s", "--spg"):
            SPG, _spg_info, _spg_str = parse_spacegroup(a)
        if o in ("-i", "--xds-input"):
            XDS_INPUT = a
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
            STRATEGY = True
        if o in ("-w", "--wavelength"):
            _wavelength = float(a)
        if o in ("-r", "--high-resolution"):
            RES_HIGH = float(a)
        if o in ("-R", "--low-resolution"):
            RES_LOW = float(a)
        if o in ("-x", "--beam_x"):
            if "mm" in a:
                _beam_in_mm = True
                a = a.replace("mm","")
            BEAM_X = float(a)
        if o in ("-y", "--beam_y"):
            if "mm" in a:
                _beam_in_mm = True
                a = a.replace("mm","")
            BEAM_Y = float(a)
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
            print USAGE
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
        print "Discontinous naming scheme, creating ling."
        _linkimages = True
        link_dir_name = "img_links"
        inputf = make_xds_image_links(inputf,
                                    os.path.join(newDir,link_dir_name),
                                    "image")
        #collect.setDirectory(link_dir_name)
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

    # Update some default values defined by XIO.export_xds:
    # In case no beam origin is defined, take the detector center.
    if newPar["ORGX"] == 0:
        newPar["ORGX"] = newPar["NX"]/2.
    if newPar["ORGY"] == 0:
        newPar["ORGY"] = newPar["NY"]/2.
    # This is to correct the starting angle in case first image is not 1.
    newPar["STARTING_ANGLE"] = newPar["STARTING_ANGLE"] - \
              newPar["OSCILLATION_RANGE"]*(newPar["DATA_RANGE"][0] - 1)
    newPar["STRONG_PIXEL"] = 7

    newrun = XDS()

    if _beam_in_mm:
        BEAM_X = BEAM_X / newPar["QX"]
        BEAM_Y = BEAM_Y / newPar["QY"]
    if ANOMAL:
        newPar["FRIEDEL'S_LAW"] = "FALSE"
    else:
        newPar["FRIEDEL'S_LAW"] = "TRUE"
    if STRICT_CORR:
        newPar["STRICT_ABSORPTION_CORRECTION"] = "TRUE"
    if BEAM_X:
        newPar["ORGX"] = BEAM_X
    if BEAM_Y:
        newPar["ORGY"] = BEAM_Y
    if SPG and _cell:
        newPar["SPACE_GROUP_NUMBER"] = SPG
        newPar["UNIT_CELL_CONSTANTS"] = _cell
    elif SPG and not _cell:
        WARNING = "  WARNING: Spacegroup is defined but not cell."
        WARNING += " Waiting for indexation for setting cell."
    elif _cell and not SPG:
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
    #if XDS_INPUT:
    #    newPar.update(xdsInp2Param(inp_str=XDS_INPUT))
    if RES_HIGH or (RES_LOW != 50):
        newPar["INCLUDE_RESOLUTION_RANGE"] = RES_LOW, RES_HIGH

    if _linkimages:
        collect.setDirectory(link_dir_name)
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

    print FMT_HELLO % vars(newrun.inpParam)
    if WARNING:
        print WARNING
    if SPG:
        print _spg_str

    #newrun.run()
    R1 = R2 = R3 = R4 = R5 = None
    if STEP > 1:
        print "\n Starting at step: %d (%s)\n" % (STEP, JOB_STEPS[STEP-1])
    if STEP <= 1:
        R1 = newrun.run_init()
    if STEP <= 2:
        R2 = newrun.run_colspot()
    if STEP <= 3:
        if SLOW and RES_HIGH:
            print "   Applying a SPOT RESOLUTION CUTOFF: %.2f A" % RES_HIGH
            #newrun.spots_resolution_cutoff(RES_HIGH)
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
        if (len(collect.imageRanges) > 1) or STRATEGY:
            newrun.run_xplan(ridx=R3)
    if STEP <= 4:
        R4 = newrun.run_integrate(collect.imageRanges)
    if STEP <= 5:
        (h, l), spgn  = newrun.run_pre_correct()
        newrun.run_correct((h, l), spgn)



