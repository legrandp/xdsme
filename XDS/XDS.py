#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 TODO-0: If just the space group is selected and not the cell:
         try to find the proper cell if it is not ambigous
         (like P21212, P2122,,,),
 TODO-1: Add in the CORRECT summary the Rmrgd-F, and total
         overloaded refl.
 TODO-2: Start multiple COLSPOT with different thresholds+ multiple IDXREF.
 TODO-3: Generating plots !
"""

__version__ = "0.5.5.2"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "23-06-2017"
__copyright__ = "Copyright (c) 2006-2017 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import os
import sys
import re

if sys.version_info <= (2, 4, 0):
    from popen2 import Popen3
else:
    from subprocess import Popen, PIPE

from XOconv.pycgtypes import mat3
from XOconv.pycgtypes import vec3
from XOconv.XOconv import reciprocal, UB_to_cellParam, BusingLevy, XDSParser

from pointless import pointless, is_pointless_installed, run_aimless, \
                 run_xdsconv
from xupy import XParam, xdsInp2Param, opWriteCl, \
                 saveLastVersion, LP_names, xdsinp_base, \
                 SPGlib, Lattice, resum_scaling, \
                 get_BravaisToSpgs, get_number_of_processors, \
                 EXCLUDE_ICE_RING, gxparm2xpar, getProfilRefPar
import XIO
from CChalf_xdsme_dev import CalculateAimlessHighRes
                             

PROGNAME = os.path.split(sys.argv[0])[1]
USAGE = """
   Running XDS automatically...

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

    -a,  --anomal
         Distinguishes Friedel pairs for scaling, strategy and completeness
         statistics. Default is no anomalous contribution.

    -A,  --Anomal
         Like -a, but also set "STRICT_ABSORPTION_CORRECTION" to True.
         It usually gives better scaling statistics with redundancy > 2.

    -b,  --beam-center-optimize-i
         Starting from the initial given values, search and optimize the beam
         center coordinates (given by -x, -y or extracted form the header).
         Best solution is chosen after i-score ranking.

    -B,  --beam-center-optimize-z
         Like -b/--beam-center-optimize-i, but best solution is chosen with
         after a z-score ranking.

    -c,  --cell
         Set the expected cell.
         For example: -c "79 79 38 90 90 90"

    -d,  --distance
         Set the detector to crystal distance.

    -f,  --reference FILE
         Defines a reference data set used during the XPLAN & CORRECT steps.
         For example: -f ../ref/XDS_ASCII.HKL

    -F, --first-frame
         Specify the first frame to be used in the DATA_RANGE (see also -L)

    -i,  --xds-input
         Give direct XDS Keyword input.
         For example: -i "DETECTOR_DISTANCE= 167.0 JOB= IDXREF AIR= 0.002"

    -I,  --ice
         Exclude resolution ranges where ice-rings occurs (3.897, 3.669,
         3.441, 2.671, 2.249, 2.249, 1.948, 1.918, 1.883, 1.721 A).

    -L, --last-frame
         Specify the last frame to be used in the DATA_RANGE (see also -F).
         This can be useful in case of radiation damage.

    -O,  --oscillation
         Set frame oscillation range in degree.
         For example: -c 0.5

    -n, --nthreads
         Set the maximum number of threads to use. Default is to use the
         maximum available.
         For example: -n 4

    -M,  --orientation-matrix
         Input crystal orientation matrix.
         For example: -M XPARM.XDS

    -p,  --project NAME
         Set the project name. The default is the prefix taken from
         image names. The working directory will be: xds_process_"project"

    -r,  --high-resolution
         Set a high resolution cutoff. Default is 0 (no cutoff).

    -R,  --low-resolution
         Set a low resolution cutoff. Default is 50 A.

    -s,  --spg
         Set the expected space group using either the space group number
         or simple string.
         For example: -s 18 or -s P21212

    -S, --strategy
         Force to go for calculating strategy (XPLAN) and then stops.

    -x,  --beam-x
         Set a new value for ORGX: X-coordinates (in pixels) of the
         detector origin. It may be given in mm if the value is directly
         ended by "mm", (e.g. -x 109.1mm).

    -y,  --beam-y
         Set a new value for ORGY: Y-coordinates (in pixels) of the
         detector origin. It may be given in mm if the value is directly
         ended by "mm", (e.g. -y 106.4mm).

    -W,  --beam-center-swap
         From the header recorded X and Y beam-center coordinate values,
         try the 8 possible permutations and select the best one based on
         z-score ranking. This very useful if indexing fails, the convention
         for recording these values may not be identical from one synchrotron
         to another.

    -v,  --verbose
         Turn on verbose output.

    -w, --wavelength
         Set the x-ray wavelength.

    --slow,
         Set parameters to process either more accurately.

    --weak,
         Set parameters to index in case of weak spots.

    --brute,
         Try hard to index. To be used in resistant cases.

    --invert,
         Invert the rotation axis. For example, used at the Australian Synchrotron (-1, 0, 0).

    --optimize, --O[1-3] 
         After a first integration, gather differente parameters to optimize
         profiles and post-refined diffraction prediction and re-run from
         step 4 (DEFPIX, INTEGRATE, CORRECT). Default for --O or --optimize
         --O1 Use postrefined geometry.
         --O2 Update profiles as suggested in INTEGRATE.LP.
         --O3 Combine O1 and O2 optimizations.
         
    -E, --exec PATH
         Path for the directory containing the executables, if different from
         or not in the default path.

""" % PROGNAME

FMT_HELLO = """
    Diffraction Setup Parameters:\n
  Detector distance:             %(DETECTOR_DISTANCE)8.2f mm
  X-ray wavelength:            %(X_RAY_WAVELENGTH)10.4f A
  Oscillation range:           %(OSCILLATION_RANGE)10.4f degree\n
  Beam coordinate X:             %(ORGX)8.1f pixel
                  Y:             %(ORGY)8.1f pixel
  Image range:                %(DATA_RANGE)11s
"""

FMT_FINAL_STAT = """
      Refined Parameters and Scaling Statistics
      =========================================\n
      Name template    %(name)s
      Data range   %(image_start)5d  to  %(image_last)5d

      Space group   number    %(spg_num)d
                    symbol    %(spg_sym)s

      Cell parameters     %(cell)s

      Resolution           %(LowestReso)8.2f -%(reso)6.2f\
    (%(resoL).2f - %(reso).2f)

      Completeness                    %(compl)5.1f%%   (%(complL).1f%%)
      I/sigma(I)                     %(isig)6.2f    (%(isigL).2f)
      Rmeas                          %(rmeas)6.1f%%   (%(rmeasL).1f%%)
      Rsym                          %(rsym)7.2f%%   (%(rsymL).1f%%)
      Multiplicity               %(multiplicity)10.1f
      Compared                   %(compar)10d    (%(comparL)d)
      Measured                   %(total)10d
      Unique                     %(unique)10d
      Rejected misfits           %(misfit)10d
      Wilson scaling (B/Corr)    %(wilson_b)10.1f    (%(wilson_corr).2f)
"""

FMT_ABSENCES = "      Systematic absent reflections measured %(AbsNum)6d \
 with <Iabs>/<I> =  %(AbsIav).1f%%\n"

FMT_ANOMAL = """
      Anomalous pairs measured   %(anoNum)10d
      SigAno                     %(anoSig)10.3f    (%(anoSigL).3f)
      Anomalous Correlation      %(anoCorr)10.1f%%   (%(anoCorrL).1f%%)
"""

STEPMARK = re.compile(r"^( [*]{5} (\w{4,}) [*]{5} )")
INTEGRATE_STEP_RE = re.compile(r" PROCESSING OF IMAGES ")
INTEGRATE_MOSAICITY_RE = re.compile(r"CRYSTAL MOSAICITY \(DEGREES\)")
INTEGRATE_STRONG_RE = re.compile(r"REFLECTIONS ACCEPTED FOR REFINEMENT")
RRF, RRI = r"[\ ]+([0-9\.]+)", r"[\ ]+([\d]+) "
SCALE_RE = re.compile(r" "+RRI+r"  (\d)"+RRF+r"  ....... "+4*RRI+2*RRF)

def _get_omatrix(_file):
    P = XDSParser(_file)
    return P.dict["symmetry"], P.dict["cell"], \
                              [P.dict["A"], P.dict["B"], P.dict["C"]]

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
                       prefix="image", start_num=1, rotationAxis=False):
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
        image = XIO.Image(_image, rotationAxis=rotationAxis)
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
    """This level of exception raises a recoverable error which can be fixed.
    """

class XDSExecError(Exception):
    ""

class XDSLogParser:
    """ Parser for the xds *.LP files.
    """
    def __init__(self, filename="", run_dir="",
                 verbose=False, raiseErrors=True):
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
        if filename:
            try:
                fp = open(full_filename, "r")
                self.lp = fp.read()
                fp.close()
            except:
                raise IOError, "Can't read file: %s" % full_filename
        else:
            self.lp = ""
        # Catch Errors:
        _err = self.lp.find(" !!! ERROR " )
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
                _err_type += " reflections bellow limit.\n"
                _err_level = "WARNING"
            elif _err_msg.count("ERROR IN REFINE !!! RETURN"):
                _err_type = "IDXREF. Can't refine cell paramters."
                _err_level = "FATAL"
            elif _err_msg.count("USELESS DATA SET"):
                _err_type = "INTEGRATE:  USELESS DATA SET."
                _err_type += " Not enough images or bad diffraction ?"
                _err_level = "FATAL"
            elif _err_msg.count("SOLUTION IS INACCURATE"):
                _err_type = "IDXREF. Solution is inaccurate.\n"
                _err_level = "WARNING"
            elif _err_msg.count("INSUFFICIENT NUMBER OF ACCEPTED SPOTS."):
                _err_type = "IDXREF. INSUFFICIENT NUMBER OF ACCEPTED SPOTS."
                _err_level = "FATAL"
            elif _err_msg.count("CANNOT INDEX REFLECTIONS"):
                _err_type = "IDXREF. CANNOT INDEX REFLECTIONS."
                _err_level = "FATAL"
            elif _err_msg.count("CANNOT CONTINUE WITH A TWO-DIMENSIONAL"):
                _err_type = "IDXREF. CANNOT INDEX REFLECTIONS."
                _err_level = "FATAL"
            else:
                print "\n %s \n" % (self.lp[_err:-1])
                sys.exit()
        if _err_level in ("FATAL", "ERROR") and raiseErrors:
            raise XDSExecError, (_err_level, _err_type)

        if self.verbose and _err != -1:
            print "\n  !!! %s in %s" % (_err_level, _err_type)

        if full_filename.count("INIT.LP"):
            self.parse_init()
        elif full_filename.count("COLSPOT.LP"):
            self.parse_colspot()
        elif full_filename.count("IDXREF.LP"):
            self.parse_idxref()
        elif full_filename.count("XPLAN.LP"):
            self.parse_xplan()
        elif full_filename.count("DEFPIX.LP"):
            self.parse_defpix()
        elif full_filename.count("INTEGRATE.LP"):
            self.parse_integrate()
        elif full_filename.count("CORRECT.LP"):
            self.parse_correct()
        else:
            if filename:
                raise IOError, "Don't know how to parse file: %s" % \
                               full_filename

    def get_par(self, match, limit=75, func=None, multi_line=False,
                      start=0, before=False, match_end=None):
        "Extract parameters from XDS .LP lines."
        try:
            if before:
                limit = start
                start = self.lp.index(match)-start
            else:
                start = self.lp.index(match, start) + len(match)
        except Exception, err:
            raise err
        if match_end:
            end = self.lp.index(match_end, start + 1)
        else:
            end = start+limit
        if multi_line:
            _raw = self.lp[start:end].split()
        else:
            _raw = self.lp[start:end].splitlines()[0].split()
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
            raise ValueError, "get_par function can't process value '%s'" \
                               % _raw
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
        try:
            rdi["mean_gain"] = gpa("MEAN GAIN VALUE")
            rdi["min_gain"] = gpa("MINIMUM GAIN VALUE IN TABLE")
            rdi["max_gain"] = gpa("MAXIMUM GAIN VALUE IN TABLE")
        except ValueError: 
            rdi["mean_gain"] = rdi["min_gain"] = rdi["max_gain"] = 1.
        rdi["mean_background"] = gpa("KGROUND COUNTS IN A DATA IMAGE PIXEL")
        #
        prp =  "  Looking at images %(background_range)s\n"
        prp += "  Mean Gain:        %(mean_gain).1f\n"
        prp += "  Min table gain:   %(min_gain).2f\n"
        prp += "  Max table gain:   %(max_gain).2f\n"
        prp += "  Mean Background:  %(mean_background).2f\n"
        if self.verbose:
            print prp % rdi
        return rdi, prp

    def parse_colspot(self):
        "Parse COLSPOT.LP"
        rdi, gpa = self.results, self.get_par
        #
        rdi["strong_pixels"] = gpa("EXTRACTED FROM IMAGES")
        rdi["weak_spots_ignored"] = gpa("WEAK SPOTS OMITTED")
        rdi["out_of_center_spots"] = gpa("CLOSE TO UNTRUSTED REGION")
        rdi["spot_number"] = self.get_spot_number()
        rdi["time"] = gpa("elapsed wall-clock time", 11)

        prp = "  Close to untrusted region: %(out_of_center_spots)10d\n"
        prp += "  Weak spots rejected:       %(weak_spots_ignored)10d\n"
        prp += "  Number of spots accepted:  %(spot_number)10d\n"
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
        st6 = "SUBTREE    POPULATION\n"
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

        subtrees = gpa(st6, multi_line=True, func=int, match_end="\n\n ")
        rdi["substrees"] = [subtrees[i] for i in range(1, len(subtrees), 2)]
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
  Shift in beam position: %(shift_mm)9.2f mm  (%(shift_pixel).1f pixels)\n"""
        prp2 = "  Size of the origin index table: %(origin_n)7d\n" % vars()
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

    def parse_defpix(self):
        "Parse DEFPIX.LP"
        rdi, gpa = self.results, self.get_par
        rdi["value_range"] = gpa("TRUSTED_DETECTOR_PIXELS= ")
        prp = "  Value range for trusted detector pixels: %(value_range)s"
        if self.verbose:
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

        try:
            sp1 = self.lp.index("   INPUT DATA SET")
            sp2 = self.lp.index("  INTEGRATE.HKL   ", sp1)
            K1s, K2s = map(float, self.lp[sp1+18: sp2].split())[:2]
            rdi["IoverSigmaAsympt"] =  1/((K1s*(K2s+0.0004))**0.5)
        except:
            try:
                    sp1 = self.lp.index("a        b          ISa") + 23
                    rdi["IoverSigmaAsympt"] = float(self.lp[sp1:sp1+31].split()[2])
            except:
                    rdi["IoverSigmaAsympt"] =  0.0
        #print "  Variance estimate scaling (K1, K2): %8.3f, %12.3e" % \
        #                                       (4*K1s, (K2s/4+0.0001))
        rdi["RMSd_spotPosition"] = gpa("SPOT    POSITION (PIXELS)")
        rdi["RMSd_spindlePosition"] = gpa("SPINDLE POSITION (DEGREES)")
        rdi["Space_group_num"] = gpa("SPACE GROUP NUMBER ")
        rdi["Cell_ref"] = gpa("UNIT CELL PARAMETERS")
        rdi["Mosaicity"] = gpa("CRYSTAL MOSAICITY (DEGREES)")
        r = gpa(" "+"-"*74+"\n")
        rdi["I_sigma"], rdi["Rsym"] = r[2], r[4]
        rdi["Compared"], rdi["Total"] = r[6], r[7]
        ### Select Diffraction range.
        sp1 = self.lp.index("RESOLUTION RANGE  I/Sigma")
        sp2 = self.lp.index(10*"-", sp1)
        _table = self.lp[sp1:sp2].splitlines()[3:-1]
        _table = [ map(float, l[:26].split()[1:3]) for l in _table ]

        rdi["HighResCutoff"] = self.get_proper_resolition_range(_table)
        prp = ""
        if rdi["Mosaicity"]:
            prp += "  RMSd spot position:     %(RMSd_spotPosition)19.2f pix,"
            prp += "%(RMSd_spindlePosition)6.2f deg.\n"
            prp += "  Refined Mosaicity:       %(Mosaicity)29.2f deg.\n\n"
        prp += "  Rsym:                             %(Rsym)9.1f\n"
        prp += "  I/sigma:                          %(I_sigma)9.1f\n"
        if rdi["HighResCutoff"]:
            prp += "  Suggested high resolution cutoff: %(HighResCutoff)9.2f"
        prp += "\n  Compared reflections:                 %(Compared)d\n"
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

    def run_exec_str(self, execstr):
        if sys.version_info <= (2, 4, 0):
            spot_file = os.popen(execstr)
            outp = spot_file.read()
            spot_file.close()
        else:
            outp = Popen([execstr], stdout=PIPE, shell=True).communicate()[0]
        return outp

    def get_xds_version(self):
        "Get the version of XDS"
        _execstr = "cd /tmp; %s | grep VERSION" % \
                                         os.path.join(XDS_PATH,"xds_par")
        wc_out = self.run_exec_str(_execstr)
        return wc_out.strip()[24:-12].replace(")","")

    def get_spot_number(self):
        "Read the number of spot directly from SPOT.XDS"
        _execstr = "wc -l %s/SPOT.XDS" % self.run_dir
        wc_out = self.run_exec_str(_execstr)
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
        if XDS_PATH:
            self.__execfile = os.path.join(XDS_PATH, "xds_par")
        else:
            self.__execfile = "xds_par"
        self.running = 0
        self.outp = []
        self.run_dir = "."
        self.status = None
        self.inpParam = XParam()
        self.collect_dir = "./"
        self.link_name_to_image = "img"
        self.running_processes = []
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

    def run(self, run_dir=None, rsave=None, verbose=True, async=False):
        "Control the runing of the xds process and parse the output."
        self.__cancelled = 0
        self.running = 1
        self.step = 0
        self.step_name = ""
        self.outp = []
        self.init_dir = os.getcwd()
        self.async = async
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
                    print "ERROR: %s" % err
                    raise XIO.XIOError, \
                     ("\nSTOP! Can't create xds working directory: %s\n" % \
                                                              self.run_dir)
            if os.path.isdir(self.run_dir):
                os.chdir(self.run_dir)
                if self.link_to_images:
                    if not os.path.exists(self.link_name_to_image):
                        os.system("ln -sf '%s' %s" % (self.collect_dir, \
                                                    self.link_name_to_image))
                        #os.system("ln -sf .. %s" % (self.link_name_to_image))
                    #else:
                    #    raise XIO.XIOError, \
                    #     "STOP! Can't creat link %s in working directory: %s" \
                    #     % (self.link_name_to_image, self.run_dir)
        opWriteCl("XDS.INP", "%s" % self.inpParam)
        #
        # self.running_processes
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

    #def run_idxref_optimize(self, number_of_test=4, verbose=False):
    #    "Run COLSPOT + DXREF with different spot search paramters"
    #    min_pixels = [4, 7, 10, 15]
    #    strong_pixel = [11, 9, 7, 5]

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
        # default is min of 3 degrees or 8 images.
        dPhi = self.inpParam["OSCILLATION_RANGE"]
        if BRUTE:
            bkgr =  i1, i1+40
        elif SLOW or WEAK:
            bkgr =  i1, min(i2, min(i1+15, i1+int(7./dPhi)))
        else:
            bkgr =  i1, min(i2, min(i1+7, i1+int(3./dPhi)))
        self.inpParam["BACKGROUND_RANGE"] = bkgr
        self.run(rsave=True)
        res = XDSLogParser("INIT.LP", run_dir=self.run_dir, verbose=1)
        if res.results["mean_background"] < 1.:
            print "  WARNING: INIT has found a very LOW mean background.\n" + \
                  "  -> Setting FIXED_SCALE_FACTOR for INTEGRATE step."
            self.inpParam["DATA_RANGE_FIXED_SCALE_FACTOR"] = i1, i2, 1.
        return res.results

    def run_colspot(self):
        "Runs the COLSPOT step."
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        self.inpParam["JOB"] = "COLSPOT",
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = 1
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = NUMBER_OF_PROCESSORS
        _trial = 0

        # DEFAULT=3.2 deg., SLOW=6.4 deg., FAST=1.6 deg.
        dPhi = self.inpParam["OSCILLATION_RANGE"]
        frames_per_colspot_sequence = FRAMES_PER_COLSPOT_SEQUENCE
        if "slow" in self.mode:
            frames_per_colspot_sequence = int(round(6.4/dPhi, 0))
        elif "fast" in self.mode:
            frames_per_colspot_sequence = int(round(1.6/dPhi, 0))
        elif BRUTE:
            frames_per_colspot_sequence = int(round(60./dPhi, 0))
            self.inpParam["VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS"] = \
                 5000, 30000
            self.inpParam["STRONG_PIXEL"] = 4.5
        else:
            frames_per_colspot_sequence = int(round(3.2/dPhi, 0))
        if "weak" in self.mode:
            self.inpParam["STRONG_PIXEL"] = 4.5
            self.inpParam["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] -= 1
            frames_per_colspot_sequence = int(round(12.8/dPhi, 0))
        # Selecting spot range(s),
        # self.inpParam["SPOT_RANGE"] is set to Collect.imageRanges by the
        # xds export function XIO
        cfo = XIO.Collect("foo_001.bar", rotationAxis=INVERT)
        cfo.imageNumbers = cfo._ranges_to_sequence(self.inpParam["SPOT_RANGE"])
        #
        min_fn, max_fn = self.inpParam["DATA_RANGE"]
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
        if BRUTE:
            self.inpParam["SPOT_RANGE"] = (min_fn, int(89./dPhi + _fpcs)),
        self.run(rsave=True)
        _rs = "  Image range(s) for spot collection: "
        for sub_range in self.inpParam["SPOT_RANGE"]:
            _rs += ("  [%d - %d]," % tuple(sub_range))
        print _rs[:-1] + "\n"

        res = XDSLogParser("COLSPOT.LP", run_dir=self.run_dir, verbose=1)
        while res.results["spot_number"] < MIN_SPOT_NUMBER and _trial < 4:
            _trial += 1
            min_pixels = int(self.inpParam["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"])
            self.inpParam["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = max(min_pixels-1, 1)
            self.inpParam["STRONG_PIXEL"] -= 1.
            #self.inpParam["SPOT_MAXIMUM_CENTROID"] += 1
            print "Insuficiant number of spot (minimum set to %d)." % \
                                                         MIN_SPOT_NUMBER
            print "Recollecting spots. Trial number %d" % _trial
            self.run(rsave=True)
            res = XDSLogParser("COLSPOT.LP", run_dir=self.run_dir, verbose=1)
        return res.results

    def run_idxref(self, beam_center_search=False, ranking_mode="ZSCORE",
                         beam_center_swap=False):
        "Runs the IDXREF step. Can try to search for better beam_center."
        res = None
        test_results = []
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
            if err[0] == "FATAL" and not (beam_center_swap or beam_center_search):
                sys.exit()
        except Exception, err:
            print err
            sys.exit()

        qx, qy = self.inpParam["QX"], self.inpParam["QY"]
        dist = self.inpParam["DETECTOR_DISTANCE"]
        det_x = vec3(self.inpParam["DIRECTION_OF_DETECTOR_X-AXIS"])
        det_y = vec3(self.inpParam["DIRECTION_OF_DETECTOR_Y-AXIS"])
        det_z = det_x.cross(det_y)
        det_params = dist, det_x, det_y, det_z, qx, qy

        #RD["indexed_percentage"] < 70. or \
        #if beam_center_search or RD["xy_spot_position_ESD"] > 2. or \
        #  RD["z_spot_position_ESD"] > 2*self.inpParam["OSCILLATION_RANGE"]:
        if res:
            test_results.append(res.results)
        if beam_center_swap:
            x, y = self.inpParam["ORGX"], self.inpParam["ORGY"]
            mx, my = self.inpParam["NX"] - x, self.inpParam["NY"] - y
            origins = [[y, x], [mx, my], [my, mx],
                       [ x, my], [y, mx], [mx, y], [my, x]]
            for origin in origins:
                self.inpParam["ORGX"] = origin[0]
                self.inpParam["ORGY"] = origin[1]
                print "   Testing beam coordinate: (%.2fmm, %.2fmm) = " % \
                                           (origin[0]*qx, origin[1]*qy),
                print "  %.1f, %.1f" % (origin[0], origin[1])
                self.run(rsave=True, verbose=False)
                try:
                    test_results.append(XDSLogParser("IDXREF.LP",
                                           run_dir=self.run_dir,
                                           verbose=0, raiseErrors=True).results)
                except XDSExecError, err:
                    print "\t\tError in", err
        if beam_center_search:
            RD = res.results
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
                    test_results.append(XDSLogParser("IDXREF.LP",
                                           run_dir=self.run_dir,
                                           verbose=0, raiseErrors=True).results)
                except XDSExecError, err:
                    print "\t\tError in", err
        if beam_center_search or beam_center_swap:
            print "\n"
            # Need to lookup in the results for the beam-center giving
            best_index_rank = rank_indexation(test_results, ranking_mode)
            #for o in origins:
            #    print origins.index(o), o[:-3]
            best_origin = origins[best_index_rank[ranking_mode]-1]
            if VERBOSE:
                print best_index_rank
                #fmt = "%4i%4i%4i%7.2f%7.2f%8.1f%8.1f%9.5f%9.5f%9.5f"
                print "best_index_rank", best_index_rank[ranking_mode]
                #print "best_origin", fmt % tuple(best_origin)
            if beam_center_search:
                best_beam = vec3(best_origin[7:10])
                best_beam_coor = best_origin[5:7]
                best_beam_orig = get_beam_origin(best_beam_coor,
                                             best_beam, det_params)
                self.inpParam["ORGX"], self.inpParam["ORGY"] = best_beam_orig
                self.inpParam["INCIDENT_BEAM_DIRECTION"] = tuple(best_beam)
            else:
                self.inpParam["ORGX"], self.inpParam["ORGY"] = best_origin
            # Running again with updated best parameters
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
        if BRUTE:
           self.inpParam["DELPHI"] = 20.
        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        self.inpParam["MAXIMUM_NUMBER_OF_PROCESSORS"] = NUMBER_OF_PROCESSORS
        self.inpParam["MAXIMUM_NUMBER_OF_JOBS"] = 1
        if ("slow" in self.mode) or BRUTE:
            self.inpParam["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA_BETA"] = 13
            self.inpParam["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA"] = 13

        "Runs the 2 first steps: DEFPIX and INTEGRATE"
        self.inpParam["JOB"] = "DEFPIX",
        self.run(rsave=True)
        res = XDSLogParser("DEFPIX.LP", run_dir=self.run_dir, verbose=1)

        if len(image_ranges) >= 1:
            self.inpParam["JOB"] = "INTEGRATE",
            self.run(rsave=True)
            res = XDSLogParser("INTEGRATE.LP", run_dir=self.run_dir, verbose=1)
            self.check_fileout("INTEGRATE.HKL")
        #else:
        #    #print "\n Error in the INTEGRATE step:"
        #    print "\n Image range:", image_ranges
        #    print " Multi-sweep integration not yet implemanted. Sorry.\n"
        #    sys.exit(0)
        return res.results

    def run_pre_correct(self, cutres=True):
        """Runs a first pass of CORRECT to evaluate high_res and
           point group.
        """
        def _get_cell(_file):
	    _txt_file = open(_file,'r').readlines()
	    if "XPARM.XDS" in _txt_file[0]:
	       return map(float, (_txt_file[3]).split()[1:])
            else:
	       return map(float, (_txt_file[7]).split()[1:])

        if XDS_INPUT:
            self.inpParam.mix(xdsInp2Param(inp_str=XDS_INPUT))
        # run pointless on INTEGRATE.HKL
        if not is_pointless_installed():
            print "!!  Warning. Pointless program not installed."
            print "  -> Skipping pointless analysis."
            likely_spg = [["P1", 0],]
            new_cell = False
        else:
            print "     Pointless analysis on the INTEGRATE.HKL file"
            print "     "+44*"="
            try:
                likely_spg, new_cell = pointless(dir_name=self.run_dir,
                                                 hklinp="INTEGRATE.HKL")
            except:
                raise
                print "  -> ERROR. While running Pointless. Skipped"
                likely_spg = [["P1", 0],]
                new_cell = False
        self.inpParam["JOB"] = "CORRECT",
        if not SPG:
            # run first CORRECT in P1 with the cell used for integration.
            # read the cell parameters from the XPARM.XDS file
            self.inpParam["SPACE_GROUP_NUMBER"] = 1
            try:
                xparm_file = os.path.join(self.run_dir, "XPARM.XDS")
                self.inpParam["UNIT_CELL_CONSTANTS"] = _get_cell(xparm_file)
            except:
                os.chdir("..")
                self.inpParam["UNIT_CELL_CONSTANTS"] = _get_cell(xparm_file)
        # run CORRECT
        self.run(rsave=True)
        res = XDSLogParser("CORRECT.LP", run_dir=self.run_dir, verbose=1)
        print "  Upper theoretical limit of I/sigma: %8.3f" % \
                                             res.results["IoverSigmaAsympt"]

        L, H = self.inpParam["INCLUDE_RESOLUTION_RANGE"]
        newH = res.results["HighResCutoff"]
        if cutres is True and newH > H and not RES_HIGH:
            H = newH
        if SPG:
            spg_choosen = SPG
        else:
            spg_choosen = likely_spg[0][1]
            # Re-order pointless cell-axes in case of orthorombic SPG.
            spgSplit = likely_spg[0][0].split()
            # if cell is coming from pointless, it need reordering
            # in orthorombic cases
            if new_cell:
                a, b, c, A, B, G = new_cell
                if spg_choosen == 18:
                    if spgSplit[1] == "2":
                        new_cell = [b, c, a, A, B, G]
                    elif spgSplit[2] == "2":
                        new_cell = [a, c, b, A, B, G]
                elif spg_choosen == 17:
                    if spgSplit[1] == "21":
                        new_cell = [b, c, a, A, B, G]
                    elif spgSplit[2] == "21":
                        new_cell = [a, c, b, A, B, G]
            else:
                new_cell = self.inpParam["UNIT_CELL_CONSTANTS"]
            lattice = Lattice(new_cell, symmetry=spg_choosen)
            lattice.idealize()
            self.inpParam["UNIT_CELL_CONSTANTS"] = lattice.cell
            #reidx_mat = likely_spg[0][-1]
            #new_cell = new_reidx_cell(self.inpParam["UNIT_CELL_CONSTANTS"],
        return (L, H), spg_choosen

    def run_correct(self, res_cut=(1000, 0), spg_num=0):
        "Runs the last step: CORRECT"
        if res_cut[1]:
            print "   ->  New high resolution limit: %.2f Å" % res_cut[1]
            self.inpParam["INCLUDE_RESOLUTION_RANGE"] = res_cut
        if spg_num:
            print "   ->  Using spacegroup: %s  #%d" % \
                                   (SPGlib[spg_num][1], spg_num)
        lattice = Lattice(self.inpParam["UNIT_CELL_CONSTANTS"],
                          symmetry=spg_num)
        lattice.idealize()
        self.inpParam["UNIT_CELL_CONSTANTS"] = lattice.cell
        self.inpParam["JOB"] = "CORRECT",
        self.inpParam["SPACE_GROUP_NUMBER"] = spg_num
        self.run(rsave=True)
        res = XDSLogParser("CORRECT.LP", run_dir=self.run_dir, verbose=1)
        print "  Upper theoritical limit of I/sigma: %8.3f" % \
                                             res.results["IoverSigmaAsympt"]        
        s = resum_scaling(lpf=os.path.join(self.run_dir,"CORRECT.LP"))
        if not s:
            print "\nERROR while running CORRECT"
            sys.exit()
        s["image_start"], s["image_last"] = self.inpParam["DATA_RANGE"]
        s["name"] = os.path.basename(self.inpParam["NAME_TEMPLATE_OF_DATA_FRAMES"])
        print s.last_table
        print FMT_FINAL_STAT % vars(s)
        if s.absent:
            print FMT_ABSENCES % vars(s)
        if self.inpParam["FRIEDEL'S_LAW"] == "FALSE":
            print FMT_ANOMAL % vars(s)
        if s.compl > 85.:
            if RUN_AIMLESS: run_aimless(self.run_dir)
            if RUN_XDSCONV: run_xdsconv(self.run_dir)

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
    """Applies the reindexing card to initial cell parameters and
    return a new cell"""
    UB = BusingLevy(reciprocal(init_cell))
    REIDX = mat3(reidx_mat)
    return reciprocal(UB_to_cellParam(REIDX*UB))

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
    spg_found = False
    try:
        spg_int = int(spginp)
        spg_found = True
    except ValueError:
        #spg_int = 0
        spginp_up = spginp.upper()
        for spgn in SPGlib:
            if spginp_up in SPGlib[spgn]:
                spg_int = spgn
                spg_found = True
                break
    if spg_found:
        if spg_int == 0:
             spg_int = 1
        spg_info = SPGlib[spg_int]
        spg_str = "  Imposed Space group:  %s,  number %d" % \
               (spg_info[1], spg_int)
    else:
        raise Exception, "\nERROR: Unrecognised space group: %s\n" % spginp
    return spg_int, spg_info, spg_str

def select_strategy(idxref_results, xds_par):
    "Interactive session to select strategy parameters."
    sel_spgn = SPG #xds_par["SPACE_GROUP_NUMBER"]
    sel_ano =  xds_par["FRIEDEL'S_LAW"]
    #print xds_par["UNIT_CELL_CONSTANTS"]
    valid_inp = False
    bravais_to_spgs = get_BravaisToSpgs()
    # Select LATTICE
    while not valid_inp:
        def_sel = 1
        if sel_spgn != 0:
            # choose the lattice solution according to the selected spg.
            i = 0
            for LAT in idxref_results["lattices_table"]:
                if LAT.fit <= LATTICE_GEOMETRIC_FIT_CUTOFF:
                    i += 1
                    if sel_spgn in bravais_to_spgs[LAT.Bravais_type]:
                        def_sel = i
        selection = raw_input("\n Select a solution number [%d]: " % def_sel)
        # If the selection is not compatible with the spg, set not valid
        _sel = selection.split()
        selnum = 1
        try:
            if len(_sel) == 1:
                selnum = int(_sel[0])
                valid_inp = True
            elif len(_sel) == 0:
                selnum = def_sel
                valid_inp = True
            else:
                raise Exception, "Invalid selection input."
        except Exception, err:
            print "\n ERROR. ", err
    sel_lat = idxref_results["lattices_table"][selnum-1]
    if sel_spgn == 0:
        sel_spgn = sel_lat.symmetry_num
    valid_inp = False
    # Select SPACEGROUP
    print " Possible spacegroup for this lattice are:\n"
    for spgsymb in bravais_to_spgs[sel_lat.Bravais_type]:
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
            if sel_spgn not in bravais_to_spgs[sel_lat.Bravais_type]:
                valid_inp = False
                msg = "Inconsistant combinaison of Bravais lattice"
                msg += " and spacegroup.\n For this Bravais Lattice"
                msg += " (%s), spacegroup should be one of these:\n\n" % \
                        (sel_lat.Bravais_type)
                for spgsymb in bravais_to_spgs[sel_lat.Bravais_type]:
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
    print "\n Selected  cell paramters:  ", sel_lat
    if sel_spgn > 2:
        sel_lat.idealize()
        print " Idealized cell parameters: ", sel_lat.prt()
        xds_par["UNIT_CELL_CONSTANTS"] = sel_lat.prt()
    xds_par["SPACE_GROUP_NUMBER"] = sel_spgn
    return xds_par


if __name__ == "__main__":

    import getopt

    short_opt =  "123456aAbBc:d:E:f:F:i:IL:O:M:n:p:s:Sr:R:x:y:vw:WSF"
    long_opt = ["anomal",
                "Anomal",
                "beam-x=",
                "beam-y=",
                "ice",
                "invert",
                "spg=",
                "strategy",
                "high-resolution=",
                "low-resolution=",
                "last-frame",
                "first-frame",
                "cell=",
                "distance",
                "reference=",
                "oscillation",
                "orientation-matrix=",
                "nthreads=",
                "project=",
                "exec=",
                "beam-center-optimize-i",
                "beam-center-optimize-z",
                "beam-center-swap",
                "xds-input=",
                "verbose",
                "optimize",
                "O1","O2","O3","O",
                "wavelength=",
                "slow", "weak", "brute",
                "CChalf="]

    if len(sys.argv) == 1:
        print USAGE
        sys.exit(2)
    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        # print help information and exit:
        print USAGE
        sys.exit(2)

    DIRNAME_PREFIX = "xdsme_"
    NUMBER_OF_PROCESSORS = min(32, get_number_of_processors())
    # Use a maximum of 16 proc. by job. Change it if you whant another limit.
    WARNING = ""
    VERBOSE = False
    DEBUG = False
    WEAK = False
    ANOMAL = False
    ICE = False
    STRICT_CORR = False
    BEAM_X = 0
    BEAM_Y = 0
    SPG = 0
    STRATEGY = False
    RES_HIGH = 0
    DISTANCE = 0
    OSCILLATION = 0
    ORIENTATION_MATRIX = False
    PROJECT = ""
    WAVELENGTH = 0
    RES_LOW = 50
    FIRST_FRAME = 0
    LAST_FRAME = 0
    REFERENCE = False
    _beam_center_optimize = False
    _beam_center_ranking = "ZSCORE"
    _beam_center_swap = False
    CELL = ""
    XDS_INPUT = ""
    _beam_in_mm = False
    SLOW = False
    FAST = False
    BRUTE = False
    STEP = 1
    OPTIMIZE = 0
    INVERT = False
    XDS_PATH = ""
    RUN_XDSCONV = True
    RUN_AIMLESS = True
    CChalf = None

    for o, a in opts:
        if o == "-v":
            VERBOSE = True
        if o in ("-E", "--exec"):
            XDS_PATH = a
            xdsf = os.path.join(XDS_PATH, "xds_par")
            if not (os.path.isfile(xdsf) and os.access(xdsf, os.X_OK)):
                print "WARNING: no 'xds_par' exec found in path %s." % XDS_PATH
                print "         Using default $PATH location."
                XDS_PATH = ""
            else:
                os.environ['PATH'] = "%s%s%s" % (XDS_PATH, os.pathsep,
                                                           os.environ["PATH"])
        if o in ("-a", "--anomal"):
            ANOMAL = True
            STRICT_CORR = False
        if o in ("-A", "--Anomal"):
            ANOMAL = True
            STRICT_CORR = True
        if o in ("-I", "--ice"):
            ICE = True
        if o[1] in "123456":
            STEP = int(o[1])
        if o in ("-s", "--spg"):
            SPG, _spg_info, _spg_str = parse_spacegroup(a)
        if o in ("-i", "--xds-input"):
            XDS_INPUT = a
        if o in ("-c", "--cell"):
            CELL = a
        if o in ("-d", "--distance"):
            DISTANCE = float(a)
        if o in ("-f", "--reference"):
            if os.path.isfile(a):
                REFERENCE = str(a)
            else:
                print "\n  ERROR: Can't open reference file %s." % a
                print "  STOP!\n"
                sys.exit()
        if o in ("-F", "--first-frame"):
            FIRST_FRAME = int(a)
        if o in ("-L", "--last-frame"):
            LAST_FRAME = int(a)
        if o in ("--CChalf"):
            CChalf = float(a)
        if o in ("-O", "--oscillation"):
            OSCILLATION = float(a)
        if o in ("-M", "--orientation-matrix"):
            if os.path.isfile(a):
                ORIENTATION_MATRIX = str(a)
            else:
                print "\n  ERROR: Can't open orientation matrix file %s." % a
                print "  STOP!\n"
                sys.exit()
        if o in ("-n","--nthreads"):
            NUMBER_OF_PROCESSORS = int(a)
        if o in ("-p", "--project"):
            PROJECT = str(a)
        if o in ("-S", "--strategy"):
            STRATEGY = True
        if o in ("-w", "--wavelength"):
            WAVELENGTH = float(a)
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
        if o in ("-W", "--beam-center-swap"):
            _beam_center_swap = True
        if o in ("--optimize", "--O"):
            OPTIMIZE = 3
            STEP = 4
        if "--O" in o and len(o) == 4:
            STEP = 4
            try:
                OPTIMIZE = int(o[-1])
            except:
                pass
            if OPTIMIZE > 3:
                OPTIMIZE = 3
        if o in ("--slow"):
            SLOW = True
        if o in ("--brute"):
            BRUTE = True
        if o in ("--weak"):
            WEAK = True
        if o == "--invert":
            INVERT = True
        if o in ("-h", "--help"):
            print USAGE
            sys.exit()

    if not inputf:
        print "\nFATAL ERROR. No image file specified.\n"
        sys.exit(2)
    elif not os.path.isfile(inputf[0]):
        print "\nFATAL ERROR. Image file %s not found.\n" % inputf[0]
        sys.exit(2)
    else:
        # TODO cycle over input_file with try/except to avoid XIOError
        _coll = XIO.Collect(inputf[0])
    if not PROJECT:
        newDir = DIRNAME_PREFIX + _coll.prefix
    else:
        newDir = DIRNAME_PREFIX + PROJECT
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
        collect = XIO.Collect(inputf, rotationAxis=INVERT)
        collect.interpretImage()
        collect.image.info()
        collect.lookup_imageRanges(forceCheck=False)

    except XIO.XIOError, _mess:
        print _mess
        print "\nError: Can't access to file(s) %s.\nStop." % inputf
        sys.exit(2)

    imgDir = collect.directory
    newPar = collect.export("xds")
    #import pprint
    #pprint.pprint(newPar)
    # Update some default values defined by XIO.export_xds:
    # In case no beam origin is defined, take the detector center.
    if newPar["ORGX"] == 0:
        newPar["ORGX"] = newPar["NX"]/2.
    if newPar["ORGY"] == 0:
        newPar["ORGY"] = newPar["NY"]/2.
    # This is to correct the starting angle in case first image is not 1.
    newPar["STARTING_ANGLE"] = newPar["STARTING_ANGLE"] - \
              newPar["OSCILLATION_RANGE"]*(newPar["DATA_RANGE"][0] - 1)
    newPar["STRONG_PIXEL"] = 6
    newPar["RESOLUTION_SHELLS"] = 15.0, 7.0, newPar["_HIGH_RESOL_LIMIT"]
    newPar["TEST_RESOLUTION_RANGE"] = 20, newPar["_HIGH_RESOL_LIMIT"]+1.5

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
    else:
        newPar["STRICT_ABSORPTION_CORRECTION"] = "FALSE" 
    if BEAM_X:
        newPar["ORGX"] = BEAM_X
    if BEAM_Y:
        newPar["ORGY"] = BEAM_Y
    if FIRST_FRAME:
        newPar["DATA_RANGE"][0] = FIRST_FRAME
    if LAST_FRAME:
        newPar["DATA_RANGE"][1] = LAST_FRAME
    if ICE:
        newPar.update(EXCLUDE_ICE_RING)
    if SPG and CELL:
        newPar["SPACE_GROUP_NUMBER"] = SPG
        newPar["UNIT_CELL_CONSTANTS"] = CELL
    elif SPG and not CELL:
        WARNING = "  WARNING: Spacegroup is defined but not cell."
        WARNING += " Waiting for indexation for setting cell."
    elif CELL and not SPG:
        WARNING = "  WARNING: Cell is defined but not spacegroup,"
        WARNING += " setting spacegroup to P1."
        newPar["SPACE_GROUP_NUMBER"] = 1
        newPar["UNIT_CELL_CONSTANTS"] = CELL
    if DISTANCE:
        newPar["DETECTOR_DISTANCE"] = DISTANCE
    if REFERENCE:
        if REFERENCE[0] == "/" or REFERENCE[0] == "~":
            newPar["REFERENCE_DATA_SET"] = REFERENCE
        else:
            newPar["REFERENCE_DATA_SET"] = "../"+REFERENCE
    if OSCILLATION:
        newPar["OSCILLATION_RANGE"] = OSCILLATION
    if ORIENTATION_MATRIX:
        try:
            _spg, cell, omat = _get_omatrix(ORIENTATION_MATRIX)
            SPG, _spg_info, _spg_str = parse_spacegroup(_spg)
            newPar["SPACE_GROUP_NUMBER"] = SPG
            newPar["UNIT_CELL_CONSTANTS"] = cell
            newPar["UNIT_CELL_A_AXIS"] = omat[0]
            newPar["UNIT_CELL_B_AXIS"] = omat[1]
            newPar["UNIT_CELL_C_AXIS"] = omat[2]
        except:
            print "\nERROR Can't import orientation matrix from: %s" % \
                                  ORIENTATION_MATRIX
            sys.exit()
    if WAVELENGTH:
        newPar["X_RAY_WAVELENGTH"] = WAVELENGTH
    if OPTIMIZE in [1, 3]:
        gxparm2xpar(newDir)
    if OPTIMIZE >= 2:
        profiles = getProfilRefPar(os.path.join(newDir, "INTEGRATE.LP"))
        newPar.update(profiles)
        print "UPDATED PROFILES: %s" % profiles
    #if XDS_INPUT:
    #    newPar.update(xdsInp2Param(inp_str=XDS_INPUT))
    if "_HIGH_RESOL_LIMIT" in newPar:
        newPar["INCLUDE_RESOLUTION_RANGE"] = RES_LOW, \
                                             newPar["_HIGH_RESOL_LIMIT"]
    if RES_HIGH:
        newPar["INCLUDE_RESOLUTION_RANGE"] = RES_LOW, RES_HIGH

    if _linkimages:
        collect.setDirectory(link_dir_name)
    else:
        collect.setDirectory(newrun.link_name_to_image)

    newPar["NAME_TEMPLATE_OF_DATA_FRAMES"] = collect.xdsTemplate

    if "SPECIFIC_KEYWORDS" in newPar.keys():
        specific_keys = newPar["SPECIFIC_KEYWORDS"]
        del newPar["SPECIFIC_KEYWORDS"]
    else:
        specific_keys = ""
    newrun.inpParam.mix(xdsInp2Param(inp_str=xdsinp_base+specific_keys))
    newrun.inpParam.mix(newPar)
    newrun.set_collect_dir(os.path.abspath(imgDir))
    newrun.run_dir = newDir
    #print newPar
    # Setting DELPHI as a fct of OSCILLATION_RANGE, MODE and NPROC
    _MIN_DELPHI = 5. # in degree
    _DELPHI = NUMBER_OF_PROCESSORS * newrun.inpParam["OSCILLATION_RANGE"]
    while _DELPHI < _MIN_DELPHI:
        _DELPHI *= 2
    newrun.inpParam["DELPHI"] = _DELPHI

    if SLOW:
        newrun.inpParam["DELPHI"] *= 2
        newrun.mode.append("slow")
    if WEAK:
        newrun.mode.append("weak")

    if XDS_PATH: print ">> XDS_PATH set to: %s" % XDS_PATH
    print "\n    Simplified XDS Processing"
    print "\n      xds   version: %18s" % XDSLogParser().get_xds_version()
    print "      xdsme version: %18s" % __version__
    print FMT_HELLO % vars(newrun.inpParam)
    print "  Selected resolution range:       %.2f - %.2f A" % \
                                           newPar["INCLUDE_RESOLUTION_RANGE"]
    print "  Number of processors available:    %3d\n" % NUMBER_OF_PROCESSORS

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
        if RES_HIGH:
            print "   Applying a SPOT RESOLUTION CUTOFF: %.2f A" % RES_HIGH
            # July 2013: spot resolution cutoff is now included in xds
            #newrun.spots_resolution_cutoff(RES_HIGH, verbose=True)
        R3 = newrun.run_idxref(_beam_center_optimize,
                               _beam_center_ranking,
                               _beam_center_swap)
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
        # If not multiple possible solutions (like P2, or P1...)try to define
        # unitcell from spacegroup.
        #if _spg and not _cell:
        if (len(collect.imageRanges) > 1) or STRATEGY:
            newrun.run_xplan(ridx=R3)
    if STEP <= 4:
        R4 = newrun.run_integrate(collect.imageRanges)
    if STEP <= 5:
        if CChalf is not None:
            (l, h), spgn  = newrun.run_pre_correct(cutres=False)
            RUN_XDSCONV=False
            newrun.run_correct((l, h), spgn)
            Newh=CalculateAimlessHighRes(filename="%s_aimless.log"%_coll.prefix, \
                                         run_dir=newrun.run_dir, verbose=1, CChalf=CChalf)
            RUN_XDSCONV=True
            if Newh is not None:
                newrun.run_correct((l, Newh), spgn)
            else:
                print "==>  Processing with full resolution range  <=="
                newrun.run_correct((l, h), spgn)
        else:
            (l, h), spgn  = newrun.run_pre_correct(cutres=True)
            newrun.run_correct((l, h), spgn)        