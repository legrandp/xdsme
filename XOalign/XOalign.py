#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2012, Pierre Legrand
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# Neither the name of the <ORGANIZATION> nor the names of its contributors may
# be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

# 25th April 2005
# Original fortran code written by:
#   Phil Evans
#   MRC LMB, Cambridge
#   October 1995

# DOTO:
#  - Add DOCTEST in functions and all...
#  - Verify the datum function with fortran gonset (gonion inversAxesOrder???)
#  - more tests...

__version__ = "0.2.9"
__author__ = "Pierre Legrand (pierre legrand \at synchrotron-soleil fr)"
__date__ = "16-11-2012"
__copyright__ = "Copyright (c) 2004-2012 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import os
import sys
import getopt
import math

from pycgtypes import vec3, mat3
from ThreeAxisRotation2 import ThreeAxisRotation2

_progname = os.path.split(sys.argv[0])[1]
_usage = """

   XOalign - Calculate possible 3-axis goniometer settings
             to realigne crystal axes.\n

   USAGE:   %s  [OPTION]... FILE

      FILE is the file containing the crystal orientation information.
            Supported types are:

                 - Denzo dot.x reflection file.
                 - Mosflm matrix files (.umat or .mat).
                 - XDS output files (CORRECT.LP, IDXREF.LP, XPARM.XDS)

             NOTE:
                 The Mosflm matrix file does not contain the space groupe
                 information. In that case, when using the -p otion, for
                 adding point group permutations, one need to specify the
                 crystal symmetry using the -s option.

   OPTIONS:

    -d
    --debug
         Turn on debug mode.

    -D angle1,angle2,angle3
    --datum angle1,angle2,angle3
         Specify the goniometer "datum" (setting) corresponding to the
         inputted crystal orientation. The three angles are given in degree.
         For exemple: -D -15.4,121.55,0
         Default datum is 0,0,0.

    -h
    --help
         Print this help message.

    -p
    --pg-permutations
         Use the allowed point group permutations.
         Given the crystal symmetry information (taken from the Denzo or XDS
         crystal orientation file, or given by the -s option), the equivalent
         crystal orientations are used to calculate all possible goniometers
         settings. This is the default.

    -n
    --no-pg-permutations
         Do not use the allowed point group permutations.

    -m mode_name
    --mode mode_name
         Specify the crystal frame to laboratory frame alignement mode.
         There are at present two possible mode:

             Mode "MAIN": Crystal vector 1 is aligned along the first
                          goniometer axis (Omega). Crystal vector 2
                          is placed in the plane containing the beam and
                          and the first goniometer axis.
             Mode "CUSP": Crystal vector 1 is aligned perpendicular to
                          both the beam vector and the first
                          goniometer axis (Omega). Crystal vector 2
                          is placed in the plane containing crystal vector 1
                          and the first goniometer axis.

         Default mode is "MAIN"
         Exemple: -m maine or -m cusp

    -s name_or_number
    --space-group name_or_number
         Specify the crystal space group. The space group name or number
         can be used. For exemple: -s p21212 or -s 18 or -s I213.

    -V crystal_vector1
    --aligned-crystal-vector-1 crystal_vector1
         Specify the first crystal vector to be aligned. The "gonset"
         notation is used to define these vectors:

             These may be given as a principle axis:
                 "a", "b", "c", "a*", "b*", or "c*"
             or as a reciprocal space vector enclosed in brackets:
                 "(h k l)",
             or a real space vector in square brackets:
                 "[a b c]"

         For exemple: -V a or -V "a*" or -V "[1 0 0]" or -V "(0 1 1)"
         Default first angle is: "a*"

    -W crystal_vector2
    --aligned-crystal-vector-2 crystal_vector2
         Specify the second crystal vector to be aligned.
         For exemple: -W b or -W "b*" or -W "[0 1 0]" or -W "(1 1 0)"
         Default second angle is: "b*"

    -O  xo,yo,zo
    --rotation-axis-1  xo,yo,zo
         Specify the first rotation axis of the goniometer, corresponding
         to the "Omega" axis on a "Kappa" goniometer.
         It is important that there is no space between comas.
         For exemple: -O 0.,0.,1 (this is the default).

    -K  xk,yk,zk
    --rotation-axis-2  xk,yk,zk
         Specify the second rotation axis of the goniometer, corresponding
         to the "Kappa" axis on a "Kappa" goniometer.
         It is important that there is no space between comas.
         For exemple: -K -0.76604444311897801,0.,0.64278760968653936
          (this is the default, which correspond to an 50deg alpha angle).

    -P  xp,yp,zp
    --rotation-axis-3  xp,yp,zp
         Specify the third rotation axis of the goniometer, corresponding
         to the "Phi" axis on a "Kappa" goniometer.
         It is important that there is no space between comas.
         For exemple: -P 0.,0.,1 (this is the default).

   AUTHORS
       Python version:
           Pierre Legrand  Proxima1 SOLEIL
           E-mail: pierre legrand \at synchrotron-soleil fr
       Original fortran version:
           Phil Evans MRC LMB, Cambridge

""" % _progname

from XOconv import mat3T, printmat, is_orthogonal, spg_num2symb, BusingLevy, \
                   SPGlib, map_r2d, PGequiv, openWriteClose, openReadClose,  \
                   rootSquareSum, random_3axes, kappaVector, SPGlib2
from XOconv import MosflmParser, DenzoParser, XDSParser

VERBOSE = True
r2d = 180/math.pi
radian2degree = lambda a: a*r2d
degree2radian = lambda a: a/r2d

ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)
X, Y, Z = ex, ey, ez
Qdnz2mos = mat3T(ez, ex, ey)

class CrystalVector(vec3):
    """ Define a crystal vector to represent reciprocal or direct space vectors

    NOTE that it can accept fractional coordinates like
         CrystalVector("(1.2 1.22 4.9)")

    NOTE that as it inherit from the Vector class, CrystalVectors support the
    usual arithmetic operations ('v1', 'v2': vectors, 's': scalar):

    -  'v1+v2'           (addition)
    -  'v1-v2'           (subtraction)
    -  'v1*v2'           (scalar product)
    -  's*v1', 'v1*s'    (multiplication with a scalar)
    -  'v1/s'            (division by a scalar)

    BUT: after these arithmetic operations, a standard Vector is returned.

    >>> print CrystalVector("a*")
    ( 1  0  0 )
    >>> print CrystalVector("c")
    [ 0  0  1 ]
    >>> xv = CrystalVector("(1 1 4)")
    >>> xv
    <CrystalVector ( 1  1  4 )>
    >>> xv.setPrintPrecision(4)
    >>> print xv
    ( 1.0000  1.0000  4.0000 )
    >>> assert xv.space == 'rec'
    >>> assert xv.inverse_space == 'dir'
    >>> assert xv*xv == 18.0
    >>> assert xv+xv == vec3(2.0,2.0,8.0)
    >>> assert xv*4 == vec3(4.0,4.0,16.0)
    """

    def __init__(self, initStr, printPrecision=0):
        "CrystalVector Constructor."

        if not type(initStr) == type(""):
            raise ValueError, 'A string value is needed.'

        axesVectors = {'a':ex, 'b':ey, 'c':ez}

        initStr = initStr.strip()
        vect = initStr[1:-1].split()
        nval = len(vect)
        axeN = initStr[0].lower()
        axeOK = axeN == 'a' or axeN == 'b' or axeN == 'c'

        if initStr[0] == "[" and initStr[-1] == "]" and nval == 3:
            self.space = 'dir'
            self.x, self.y, self.z = map(float, vect)
        elif initStr[0] == "(" and initStr[-1] == ")" and nval == 3:
            self.space = 'rec'
            self.x, self.y, self.z = map(float, vect)
        elif initStr.count("*") == 1 and len(initStr) == 2 and axeOK:
            self.space = 'rec'
            self.x, self.y, self.z = axesVectors[axeN]
        elif initStr.count("*") == 0 and len(initStr) == 1 and axeOK:
            self.space = 'dir'
            self.x, self.y, self.z = axesVectors[axeN]
        else:
            raise ValueError, 'Unrecognised string constructor: %s' % initStr

        if  self.space == 'dir':
            self.inverse_space = 'rec'
            self.space_name = 'real space'
        elif self.space == 'rec':
            self.inverse_space = 'dir'
            self.space_name = 'reciprocal space'

        self.init_name = initStr
        self.setPrintPrecision(printPrecision)

    def __str__(self):
        return self.prt_fmt % (self.x, self.y, self.z)

    def __repr__(self):
        return ('<' + self.__class__.__name__ + ' ' + str(self) + '>')

    def setPrintPrecision(self, printPrecision):
        if type(printPrecision) == type(1):
            self.printPrecision = printPrecision
            fmt1 = 3*(' %.'+str(printPrecision)+'f ')
            if self.space == 'dir':
                self.prt_fmt = '[%s]' % fmt1
            elif self.space == 'rec':
                self.prt_fmt = '(%s)' % fmt1
        else:
            raise ValueError, 'A int value is needed.'

def gnsmod(mode, beamVector, Goniometer):
    """Set laboratory frame mode.
    mode           = 'MAIN' or'CUSP'
    beamVector     = Vector of the beam
    Goniometer     = ThreeAxisRotation2 object


    return 2 laboratory frame Vectors l1, l2
    If mode = 'MAIN', then l1 = e1, l2 = beamVector
    If mode = 'CUSP'  then l1 = e1 x beamVector, l2 = e1
    """

    if _mode == 'MAIN':
        return Goniometer.e1, beamVector
    elif _mode == 'CUSP':
        return Goniometer.e1.cross(beamVector), Goniometer.e1
    else:
        raise Exception, 'Unrecognized mode'

def gnsnow(orthogMatrices, U0mos, Goniometer, beamVector):
    """ Calculate various angles from current setting.

    orthogMatrices = {'rec':B, 'dir':Bm1t}
    B      Cell orthogonalization matrix B (multiplies h)
    Bm1t   To orthogonalize real space vector (B**-1)T
    U0mos  True orthogonalization matrix, at zero goniostat angles
    D      Datum matrix D: current orientation matrix U = D U0mos

    Returns:
    Angles (in degrees) at current datum between:
         /  crystal axes in reciprocal (a*,b*,c*), real space (a,b,c)
    and /__ spindle (e3), beamVector at current datum, and Omega axis (e1)
    """
    e1, e2, e3 = Goniometer.rotationAxes
    angles_XtalAxis = {'rec':{'e3':[], 'beamVector':[], 'e1':[]},
                       'dir':{'e3':[], 'beamVector':[], 'e1':[]}}


    # angleD: Calculate the angle between v1 and v2 in degree
    # Different results than the Vector associated method: v1.angle(v2)
    # (because of the use of abs(v1*v2)?)
    # Return the result in degree.
    angleD = lambda v1,v2: math.acos(abs(v1*v2))*r2d

    # Loop abc axes (real & reciprocal)
    for space in orthogMatrices.keys():
        for axis in ex, ey, ez:

            # Orthogonalize (w1 = B*axis) and then
            # rotate by true setting matrix U0 + Normalize
            # Normalize needed: to avoid math.acos exception:
            #    ValueError: math domain error

            v = U0mos*(orthogMatrices[space]*axis).normalize()

            # Angle between vector & spindle (e3) = acos(v . e3)
            angles_XtalAxis[space]['e3'].append(angleD(v,e3))

            # v rotated by curent datum

            Dv = Goniometer.tensor * v

            # Angle between Dv and beam beamVector
            angles_XtalAxis[space]['beamVector'].append(angleD(Dv,beamVector))

            # Angle between Dv and e1

            angles_XtalAxis[space]['e1'].append(angleD(Dv,e1))

    return angles_XtalAxis

def now(cell, Umos, U0mos, phi123, orthogMatrices, Goniometer, beamVector,
                                                   datum, desiredSetting):
    """now() prints out various information on the input orientation and
       desired settings.
    """
    AXES_NAME = {'dir':['a', 'b', 'c'], 'rec':['a*', 'b*', 'c*']}
    print (" Cell dimensions: " + 3*"%8.3f" + 3*"%8.2f") % tuple(cell)
    print (" Setting matrix at current Datum [U]:\n" +
              3*("  ( "+3*"%9.5f"+" )\n")) % tuple(Umos.mlist)
    print " Setting angles Phix, Phiy, Phiz : %8.3f%8.3f%8.3f" % \
                          tuple(phi123)

    res = gnsnow(orthogMatrices, U0mos, Goniometer, beamVector)
    fmt = "\n "+3*"%8.3f"+"   Relative zero position (Datum)"+\
          "  Omega  Kappa  Phi"
    print fmt % tuple(datum)
    space = 'rec'
    for key in camAxesKeys:
        axesAngles = res[space][key]
        nearest = AXES_NAME[space][axesAngles.index(min(axesAngles))]
        print (23*" "+"%s   Reciprocal axis nearest to %s") % \
                                            (nearest,camAxes[key])

    fmt1 = '  Angles between crystal axes and'
    fmt2 = '    "       "       "      "  and'
    fmt = {'e3':fmt1, 'beamVector':fmt2, 'e1':fmt2}
    for space in "rec", "dir":
        print
        print (3*"      %-2s") % tuple(AXES_NAME[space])
        for axis in camAxesKeys:
            print " "+3*"%8.3f" % tuple(res[space][axis]),
            print fmt[axis],camAxes[axis]

    print "\n\n Desired setting:"
    for vi in 0, 1:
        print " %s Crystal vector %d, %s (v%d)" % \
               (desiredSetting[vi], vi+1, desiredSetting[vi].space_name, vi+1)
    #<<<<<<<<<<<<<<<<< NOW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def solve(desiredSetting, VL, orthogMatrices, U0, Goniometer,
                              PG_permutions=[[X,Y,Z]]):
    """ Solve for datum position
    Uses:
        VL              Laboratory vectors to match
        orthogMatrices  Dict containing {'rec':Bm1t, 'dir':B}
        desiredSetting  List of the two vectors defining the desired setting
        U0              Zero-angle setting matrix
        PG_permutions   Point group permutations

    Return:
        datumSolutions (list of datum solutions in degrees)

    Local variables
        U0m1    [U0]**-1
        C       base vectors in crystal frame
        Cm1     [C]**-1
        L       base vectors in laboratory frame
        SL      L matrix with permuted signs
        CU      [C]**-1 [U0]**-1
    """

    # Orthogonalize crystal vectors h -> yv
    # YV crystal frame vectors defining desired setting (orthogonalized)

    YV = []
    for vi in 0, 1:
        # if reciproc: B*vi.vector elif direct: B1mt*vi.vector
        B = orthogMatrices[desiredSetting[vi].space]
        YV.append(B*desiredSetting[vi])

    #  referential_permutations sign permutations for four permutations of
    #        parallel/antiparallel (rotation axis & beam)
    #    y1 // e1, y2 // beamVector;  y1 anti// e1, y2 // beamVector
    #    y1 // e1, y2 anti// beamVector;  y1 anti// e1, y2 anti// beamVector

    referential_permutations = [ ex,  ey,  ez], \
                               [-ex, -ey,  ez], \
                               [ ex, -ey, -ez], \
                               [-ex,  ey, -ez]

    # Construct orthogonal frame in crystal space defined by vectors
    # y1 = yv(,1) and y2 = yv(,2). The base vectors form the COLUMNS
    # of matrix C
    #  - 1st vector along YV1
    #  - 3rd vector perpendicular to YV1 & YV2
    #  - 2nd vector completes the RH set

    C1 = YV[0].normalize()
    C3 = YV[0].cross(YV[1]).normalize()
    C2 = C3.cross(C1).normalize()
    C = mat3(C1, C2, C3)

    # Similarly, construct laboratory frame matrix defined by VL1, VL2
    #  1st vector along VL1
    #  3rd vector perpendicular to VL1 & VL2
    #  2nd vector completes the RH set

    L1 = VL[0].normalize()
    L3 = VL[0].cross(VL[1]).normalize()
    L2 = L3.cross(L1).normalize()
    L = mat3(L1, L2, L3)

    # The matrix which superimposes the two frames is L C**-1. This
    # is then the orientation matrix DU, including the required datum
    # position D.  The datum matrix is then  D = L C**-1 U**-1
    # We need to consider four possible alignments & values of D,
    # corresponding to y1 parallel & antiparallel to vl1, and y2
    # parallel & antiparallel to vl2. To get these, we change
    # the signs of l1 & l2, or l2 & l3, or both, as set out in array LS.

    datumSolutions = []
    independentSolution = []

    # Loop four solutions for D, each with possibly two sets of goniostat
    # angles

    Cm1 = C.inverse()
    B = orthogMatrices['rec']
    UB0 = U0 * B

    # Loop for introducing Point Group permutations:
    for PG_operator in PG_permutions:
        # We use the operator in reciprocal space so it need a transposition.
        PG_operator = mat3T(PG_operator)
        # printmat(B.dot(PG_operator),'B.dot(PG_operator)')

        # The PG operator is applied on the UB0 matrix
        U0perm = UB0 * PG_operator * B.inverse()
        if _debug:
            print
            printmat(PG_operator, 'PG_operator', "%12.6f")
            printmat(U0, 'U0', "%12.6f")
            printmat(U0perm, 'U0perm', "%12.6f")

        if not is_orthogonal(U0perm):
            print "!!!!!!!!!!!!!!! ERROR: Improper U0perm matrice!!"

        #  CU = C**-1 . U**-1
        CU = Cm1 * U0perm.inverse()
        for referential_permut in referential_permutations:

            SL = L * mat3T(referential_permut)
            _DM = SL * CU

            # angles: there are zero or two possible sets of goniostat angles
            # which correspond to D, although one or both may be inaccessible

            try:
                _2s = ThreeAxisRotation2(_DM.mlist,
                                        Goniometer.rotationAxes,
                                        inversAxesOrder=0).getAngles()
                twoSolutions = (_2s[0][2],_2s[0][1],_2s[0][0]), \
                               (_2s[1][2],_2s[1][1],_2s[1][0])

                for oneSolution in twoSolutions:

                    solutionKey = "%6.2f%6.2f%6.2f" % tuple(oneSolution)
                    if solutionKey not in independentSolution:
                        independentSolution.append(solutionKey)
                        datumSolutions.append(map(radian2degree, oneSolution))
                        #print ">>>>>>>>>",solutionKey,"number ",
                        #print len(independentSolution)
                    else:
                        if _debug:
                            print ">>>>>>>>> Redundant solution"
            except ValueError:
                pass

    return datumSolutions

def print_solutions(datumSolutions, v1v2, axes_names=("Omega","Kappa","Phi")):
    print "\n Solutions for Datum positions:  %s,  %s\n" % v1v2
    print (15*" "+3*"%10s") % axes_names

    if datumSolutions:
        nsol = 1
        for sol in datumSolutions:
            print "%10d    " % nsol,
            print 3*"%10.3f" % (sol[2], sol[1], sol[0])
            nsol += 1
    else:
        print ' *** No solutions found ***'

def _test0():
    import doctest
    doctest.testmod()


def main(GoniometerAxes, inputFile, mode, v1, v2, datum, beam, spgn=None):

    global XOfileType

    # The goniostat consists of axes e1 carrying e2 carrying e3
    # (eg Omega, Kappa, Phi). GoniometerAxes = e1, e2, e3
    # gnsdtm calculate D from datum and GoniometerAxes.
    # D == Goniometer.tensor == datum matrix
    Goniometer = ThreeAxisRotation2((0.,0.,0.), GoniometerAxes,
                                                inversAxesOrder=0)

    # DATUM definition in degree.
    setDatum = tuple(datum) # in degree
    Goniometer.setAngles(map(degree2radian,setDatum))

    # MODE = 'CUSP' or 'MAIN'
    VL = gnsmod(mode, beam, Goniometer)

    # V1/V2
    V1 = CrystalVector(v1, printPrecision=4)
    V2 = CrystalVector(v2, printPrecision=4)
    if "%s" % V1  != "%s" % V2:
        desiredSetting = V1, V2
    else:
        print "ERROR: The two crystal aligned vector can't be identical!"
        sys.exit(2)

    if type(inputf[0]) == str:
        gotcha = False
        for parser in (DenzoParser, XDSParser, MosflmParser):
            try:
                XOparser = parser(inputf[0])
                gotcha = True
            except :
                pass
            if gotcha:
                break
        if not gotcha:
            raise Exception, "Can't parse input orientation matrix file: %s" %\
                inputf[0]
    elif hasattr(inputf, "fileType"):
        XOparser = inputf


    if VERBOSE:
        print "\n %s used to read input file: %s" % (XOparser.info, inputf[0])
    XOfileType = XOparser.fileType

    if not spgn:
        spgn = XOparser.spaceGroupNumber
        spg = XOparser.spaceGroupName
    else:
        spg = spg_num2symb[spgn]

    pointGroup = SPGlib[spgn][-1]
    if VERBOSE:
        print  "\n Space group symbol: %s,  Number: %d,  Point Group: %s" % \
                                            (spg.upper(),spgn, pointGroup)
        print ("\n       Real cell: "+3*"%10.3f"+3*"%8.2f") % \
                                                       tuple(XOparser.cell)
        print (" Reciprocal cell: "+3*"%10.6f"+3*"%8.2f") % \
                                                       tuple(XOparser.cell_r)

    Bmos = BusingLevy(XOparser.cell_r)

    # Umos =     setting matrix at current datum
    # Getting Mosflm XO
    if XOparser.fileType == "Mosflm":
        UBmos =  XOparser.U * Bmos # whitout wavelength scaling
        Umos = XOparser.U
        # XOparser.UB = UBmos * wavelength

    # Converting Denzo XO to Mosflm convention
    elif XOparser.fileType == "Denzo":
        UBmos = Qdnz2mos * XOparser.UB
        Umos = (UBmos) * Bmos.inverse()

    # Converting XDS XO to Mosflm convention
    elif XOparser.fileType == "XDS":
        UBmos = XOparser.UBxds_to_mos()/ XOparser.dict["wavelength"]
        Umos = (UBmos) * Bmos.inverse()

    is_orthogonal(Umos)
    if VERBOSE:
        printmat( Umos,'\n   U',  "%12.6f")
        printmat( Bmos,'\n   B',  "%12.6f")
        printmat( UBmos,'\n  UB', "%12.6f")

    Bm1t = Bmos.inverse().transpose()
    orthogMatrices = {'rec':Bmos, 'dir':Bm1t}

    # There are two ways to write the crystal orientation in mosflm:
    # using the U matrix or the phi1, phi2, phi3 misseting angles (in degree).
    # == gnspxg: decompose Umos into 3 "missetting" angles phi123
    phi123r = ThreeAxisRotation2(Umos.mlist, inversAxesOrder=1).getAngles()
    # keep only the first solution and convert it in degre.
    phi123 = map_r2d(phi123r[0])

    # This new Umos def - equivalent to the previous - is needed to
    # stay exactly compatible with gonset.
    # Umos2 = ThreeAxisRotation(phi123r[0], inversAxesOrder=1).tensor
    # printmat(Umos2, 'Umos2')
    # print "\nDiff Umos/Umos2 =======>> %8.1e" %  diffMAT(Umos2,Umos)

    # U0mos =    zero-angle setting matrix
    #           [U] = [D] [U0] so [U0] = [D]**-1 [U]
    U0mos = Goniometer.tensor.inverse() * Umos # == gnsszr
    if VERBOSE:
        print ("\n phixyz: "+3*"%8.2f"+"\n") % tuple(phi123)
        now(XOparser.cell, Umos, U0mos, phi123, orthogMatrices, Goniometer,
                             beam, setDatum, desiredSetting)

    PG_permutions = [[X,Y,Z]]
    if _do_PG_permutations:
        PG_permutions.extend(PGequiv[pointGroup])

    datumSolutions = solve(desiredSetting, VL, orthogMatrices,
                              U0mos, Goniometer, PG_permutions)

    # Ajouter la verification directement au moment du calcul
    # Permutation of the UB matrix... (U matrix) only?
    # How the permutation is handled for hexagonal.
    #
    n = 0
    #if 0:
    if _debug:
        for newdatum in datumSolutions:
            n += 1
            Goniometer.setAngles(map(degree2radian, newdatum))

            newUmos = Goniometer.tensor * U0mos
            printmat(newUmos, "\nNew Umos matrix. Solution numb. %3d" % n,\
                       "%12.6f")
            phi123 = ThreeAxisRotation2(newUmos.mlist,
                                        inversAxesOrder=1).getAngles()
            phi123 = map_r2d(phi123[0])

            print ("\n phixyz: "+3*"%8.2f"+"\n") % tuple(phi123)

    return datumSolutions

def run_gonset(gonioAxes, inputFile, mode, v1, v2, datum, beamVect, fileType):

    if fileType not in ("Mosflm", "Denzo"):
         print "\nSorry, gonset is unable to accept "
         print "orientation file of the %s type!\n" % fileType
         sys.exit(2)

    # write "GNSDEF" file, containing the goniometer rotation axes vectors
    axesName = "Omega", "Kappa", "Phi"
    gnsdef = "# Source:\n %10.5f%10.5f%10.5f\n" % tuple(beamVect)
    gnsdef += "#\n# Gonimeter Axes\n"
    for i in range(3):
        gnsdef += "%-10s" % axesName[i]
        gnsdef += "%10.5f%10.5f%10.5f\n" % tuple(gonioAxes[i])
    openWriteClose("GNSDEF",gnsdef)

    # Write gonset input file:
    gonsetinp = "now\nmode %s\n"% (mode)
    if fileType == "Denzo":
        gonsetinp += "denzo %s\n" % (inputFile)
    elif fileType == "Mosflm":
        gonsetinp += "refix %s\n" % (inputFile)
    gonsetinp += "datum %10.4f%10.4f%10.4f\n" % tuple(datum)
    gonsetinp += "v1 %s\nv2 %s\nnow\nsolve\nend\n" % (v1, v2)
    openWriteClose("gonset.inp",gonsetinp)

    # Run gonset
    #os.system("gonset_pl < gonset.inp > gonset_pl.log")
    os.system("gonset < gonset.inp > gonset.log")

    # Extract gonset results
    gonsetlog = openReadClose("gonset.log")
    try:
        idx1 = gonsetlog.index(" Omega     Kappa       Phi")+27
    except:
        print " *** No solutions for gonset ***"
        sys.exit()

    idx2 = gonsetlog.index("Data line--- end") - 2
    line2sol = lambda s: map(float, s.split()[1:])
    solutions = map(line2sol, gonsetlog[idx1:idx2].splitlines())

    return solutions

def compareSolutions(solutions1, solutions2, _epsilon=0.1):
    "Check if both solution matchs, with a tolerated difference of _epsilon."
    l = 0
    allMatchs = True
    for s1 in solutions1:
        match = False
        minRMSdiff = 1000
        for s2 in solutions2:
            l += 1
            vecDiff =  vec3(s1[2],s1[1],s1[0]) - vec3(s2)
            RMSdiff = rootSquareSum(vecDiff)/3.
            #print vecDiff, RMSdiff

        if RMSdiff < minRMSdiff:
            minRMSdiff = RMSdiff

        if RMSdiff < _epsilon:
                match = True
                break
    if match:
            solutions2.remove(s2)
            print "Good match for solution: %9.3f%9.3f%9.3f" % tuple(s1),
            print " Minimum RMSdiff = %.3f" % minRMSdiff
    else:
            allMatchs = False
            print "Warning: no match for solution: %9.3f%9.3f%9.3f" %tuple(s1),
            print " Minimum RMSdiff = %.3f" % minRMSdiff
    print l
    return allMatchs

def _random_crystal_axes():
    xaxes = ["a", "b", "c", "a*", "b*", "c*", "[0 1 1]", "(2 0 1)"]
    _ta1 = random.choice(xaxes)
    xaxes.remove(_ta1)
    if len(_ta1) == 1: xaxes.remove(_ta1+"*")
    if len(_ta1) == 2: xaxes.remove(_ta1[0])

    _ta2 = random.choice(xaxes)
    print " Testing align crystal axes : %s, %s" % (_ta1, _ta2)
    return _ta1, _ta2

def _random_gonioAxes():
    _randomGoniometerAxes = random_3axes()
    print " Testing Goniometer Axes"
    return _randomGoniometerAxes

def _random_mode():
    _randomMode = random.choice(("MAIN", "CUSP"))
    print " Testing Mode: %s" % _randomMode
    return _randomMode

def _random_datum():
    _randomDatum = random.uniform(-180.,180.), \
             random.uniform(-180.,180.), \
             random.uniform(-180.,180.)
    print " Testing Datum : %7.1f%7.1f%7.1f" % _randomDatum
    return _randomDatum


if __name__ == '__main__':

    camAxes = {'e3':'Phi axis','beamVector':'BEAM at Datum','e1':'Omega axis'}
    camAxesKeys = 'e3', 'beamVector', 'e1'

    # Definition equivalent to GNSDEF
    # StandardAxes = ex, ey, ez
    NoniusKappaAxes = [ez, kappaVector(50./r2d), ez]
    CrystalLogicKappaAxes = [ez, kappaVector(49.64/r2d), ez]
    EMBLminiKappaAxes = [ez, kappaVector(-24/r2d), ez]
    EulerAxes = [ez, ex, ez]

    from XOalign_sitedef import GONIOMETER_AXES, GONIOMETER_AXES_NAMES, \
                                    GONIOMETER_DATUM, GONIOMETER_NAME

    # Definitions:
    #    GONIOMETER_NAME = "SOLEIL PROXIMA-1 CrystalLogic"
    #    GONIOMETER_AXES = CrystalLogicKappaAxes
    #    GONIOMETER_AXES_NAMES = ("Omega","Kappa","Phi")
    #    GONIOMETER_DATUM = (0,0,0)  # in degree
    
    print "\n XOalign: Calculate possible 3-axis goniometer settings to",
    print "realigne crystal axes."
    print " Version: %s\n" % (__version__)
    print "   Goniometer definition used (from XOalign_sitedef.py): "
    print "     Name: %s\nAxes:" % GONIOMETER_NAME
    for name, axis in zip(GONIOMETER_AXES_NAMES, GONIOMETER_AXES):
        print "%14s:   %8.5f%8.5f%8.5f" % tuple([name]+list(axis))
    print "     Datum:  %8.2f%8.2f%8.2f  (in degree)" % GONIOMETER_DATUM
    # Default parameters
    _debug = False
    _test = False
    _verbose = False
    VERBOSE = True

    _do_PG_permutations = True
    _space_group_numb = 0

    # Definition equivalent to GNSDEF
    _beam_vector = ex

    # Default orientation
    _v1, _v2 = "", ""
    _mode = 'MAIN'

    XOfileType = "Denzo"

    short_opt =  "dD:hnm:O:K:P:ps:tvV:W:"
    long_opt = ["debug", "test", "datum=", "help", "no-pg-permutations",
                "pg-permutations", "mode=", "space-group="
                "aligned-crystal-vector-1=", "aligned-crystal-vector-2=",
                "rotation-axis-1=","rotation-axis-2=","rotation-axis-3="]

    try:
        opts, inputf = getopt.getopt(sys.argv[1:], short_opt, long_opt)
    except getopt.GetoptError:
        print _usage
        sys.exit(2)

    for o, a in opts:
        if o == "-v":
            _verbose = True
        elif o in ("-h", "--help"):
            print _usage
            sys.exit()
        elif o in ("-d", "--debug"):
            _debug = True
        elif o in ("-t", "--test"):
            import random
            #_debug = True
            _test = True
        elif o in ("-p","--pg-permutations"):
            _do_PG_permutations = True
        elif o in ("-n","--no-pg-permutations"):
            _do_PG_permutations = False
        elif o in ("-s","--space-group"):
            try:
                _space_group_numb = int(a)
            except:
                pass
            if not _space_group_numb:
                try:
                    _space_group_numb = SPGlib2[a.lower()]
                except:
                    pass
            if not _space_group_numb:
                print "\n  Error. Can't unknown space group.\n"
                print _usage
                sys.exit(2)
        elif o in ("-V","--aligned-crystal-vector-1"):
            _v1 = a
        elif o in ("-W","--aligned-crystal-vector-2"):
            _v2 = a
        elif o in ("-O","--rotation-axis-1"):
            GONIOMETER_AXES[0] = vec3(map(float, a.split(","))).normalize()
        elif o in ("-K","--rotation-axis-2"):
            GONIOMETER_AXES[1] = vec3(map(float, a.split(","))).normalize()
        elif o in ("-P","--rotation-axis-3"):
            GONIOMETER_AXES[2] = vec3(map(float, a.split(","))).normalize()
        elif o in ("-D","--datum"):
            GONIOMETER_DATUM = map(float, a.split(","))
        elif o in ("-m","--mode"):
            if a.upper() in ["MAIN", "CUSP"]:
                _mode = a.upper()
            else:
                print "\n  Error. Unknown mode option.\n"
                print _usage
                sys.exit(2)

    gotcha = False
    for parser in (DenzoParser, XDSParser, MosflmParser):
        try:
            XOparser = parser(inputf[0])
            gotcha = True
        except :
            pass
        if gotcha:
            break
    if not gotcha:
        raise Exception, "Can't parse input orientation matrix file: %s" %\
            inputf[0]
    if not _space_group_numb:
        _space_group_numb = XOparser.spaceGroupNumber

    if  _space_group_numb in [143, 144, 145, 149, 150, 151, 152, 153, 154, 168,
          169, 170, 171, 172, 173, 177, 178, 179, 180, 181, 182, 146, 155]:
        print "WARNING: The symmetry permutation option is not yet working for"
        print "         the Trigonal and Hexagonal space groups."
        _do_PG_permutations = False

    if _debug and sys.version_info[:3] > (2,2,0):
        _test0()

    if _test:
        GONIOMETER_DATUM = _random_datum()
        GONIOMETER_AXES = _random_gonioAxes()
        _v1, _v2 =  _random_crystal_axes()
        _mode =_random_mode()
    all_solutions = {}
    if _v1 and _v2:
        allSets = ((_v1, _v2),)
    elif _v1 == "a*":
        allSets = (("a*", "b*"), ("a*", "c*"))
    elif _v1 == "b*":
        allSets = (("b*", "a*"), ("b*", "c*"))
    elif _v1 == "c*":
        allSets = (("c*", "a*"), ("c*", "b*"))
    else:
        allSets = (("a*", "b*"), ("a*", "c*"),
                   ("b*", "a*"), ("b*", "c*"),
                   ("c*", "a*"), ("c*", "b*"))

    for v1v2 in allSets:
        _v1, _v2 = v1v2
        all_solutions[v1v2] = main(GONIOMETER_AXES, inputf[0], _mode, _v1,
                              _v2, GONIOMETER_DATUM, _beam_vector, _space_group_numb)
        print_solutions(all_solutions[v1v2], v1v2, GONIOMETER_AXES_NAMES)
        if not _verbose:
            VERBOSE = False
        # Compare results of Gonset and XOalign
        if _debug or _test:
            print
            GonsetSolutions = run_gonset(GONIOMETER_AXES, inputf[0], _mode,
                                         _v1, _v2, GONIOMETER_DATUM, _beam_vector,
                                         XOfileType)
            allMatchs = compareSolutions(GonsetSolutions,
                                         all_solutions[v1v2], 0.02)

            if not allMatchs:
                print "Error !!! Results do not comape well with gonset!!!!!"
    independant_solutions = {}
    for sols in all_solutions:
        for sol in all_solutions[sols]:
            key = "%9.3f%9.3f" % (sol[1],sol[0])
            if key not in independant_solutions:
                independant_solutions[key] = [sols]
            else:
                if sols not in independant_solutions[key]:
                    independant_solutions[key].append(sols)
    n = 1
    print "\n\n Independent Solutions for %s goniometer:\n" % GONIOMETER_NAME
    print "       %9s%9s      Settings" % GONIOMETER_AXES_NAMES[1:3]
    for isol in independant_solutions:
            print "%4d   %s   %s" % (n, isol, independant_solutions[isol])
            n += 1
