#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2005, Pierre Legrand
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA 02111-1307 USA

"""
    19/10/05 First version

    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).

    TODO:
      - 
"""

__author__ = "Pierre Legrand (pierre legrand \at synchrotron-soleil fr)"
__date__ = "21-11-2009"
__copyright__ = "Copyright (c) 2007-2009  Pierre Legrand"
__license__ = "New BSD"
__version__ = "0.3.4"


import sys
import os.path
import math

from pycgtypes import vec3
from pycgtypes import mat3
from ThreeAxisRotation2 import *


r2d = 180/math.pi
cosd = lambda a: math.cos(a/r2d)
sind = lambda a: math.sin(a/r2d)
map_r2d = lambda l: map(lambda x: x*r2d, l)
map_d2r = lambda l: map(lambda x: x/r2d, l)
str2floats = lambda s: map(float, s.split())
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)

Qdnz2mos = mat3( ez, ex, ey).transpose()
Qdnz2xds = mat3(-ey, ex,-ez).transpose()
Qmos2xds = mat3( ez,-ey, ex).transpose()
Qmos2dnz = mat3( ey, ez, ex).transpose()

DNZAxes = ey, -ex, -ez

if sys.version_info[:3] < (2,2,0):
    True = 1
    False = 0

_debug = False
#_debug = True

def mat3T(*args):
    if len(args) == 3:
        return mat3(args[0], args[1], args[2]).transpose()
    elif len(args) == 1:
        return mat3(args[0][0], args[0][1], args[0][2]).transpose()
    else:
        raise Exception, "Number of argument should be 1 or 3!"

class ParserError(Exception):
    """This level of exception raises a recoverable error which can be fixed.
    """

class DenzoParser:

    def __init__(self, filename=None):
        self.DNZAxes = ey, -ex, -ez
        self.verticalAxis = vec3(1, 0, 0)
        self.spindleAxis = vec3(0, 0, 1)
        self.motorAxis = [0.,1.,0.]
        self.info = "Denzo Parser"
        self.fileType = "Denzo"
        if filename:
            self.parse(filename)
            self.spaceGroupName = self.spg.upper()
            self.spaceGroupNumber = SPGlib2[self.spg.lower()]

    def parse(self, filename):
        "Denzo x-file parser"
        try:
            xfile = open(filename,"r").read().splitlines()
            xhead = xfile[:7]
            xtail = xfile[-30:]

        except:
            raise ParserError, "Error, Can't read file: %s" % filename

        if xhead[0][:6] == "HEADER":
            xhead =  xhead[1:]
        self.title = xhead[0]
        mats = map(str2floats, xhead[1:4])
        self.UB = mat3(mats[0][:3], mats[1][:3], mats[2][:3]).transpose()
        self.U =  mat3(mats[0][3:], mats[1][3:], mats[2][3:]).transpose()

        if len(xhead[4].split()) == 4: line1, line2 = xhead[4], xhead[5][:40]
        else: line1, line2 = xhead[4][:48], xhead[4][48:88]
        self.phi0, self.phi1, self.xtod, self.wavel = str2floats(line1)
        self.rotz, self.roty, self.rotx, self.mosaic = str2floats(line2)
        self.crystal_setting = self.rotz, self.roty, self.rotx

        # Extract reciprocal unit cell vectors
        self.Ar = vec3(self.UB.getColumn(0))
        self.Br = vec3(self.UB.getColumn(1))
        self.Cr = vec3(self.UB.getColumn(2))

        # Extract reciprocal cell parameters
        self.cell_r = UB_to_cellParam(self.UB)
        self.cell = reciprocal(self.cell_r)

        # Calculate direct unit cell vectors
        self.volum_r = self.Ar.cross(self.Br)*self.Cr
        self.volum = 1/self.volum_r

        self.A = self.Br.cross(self.Cr)*self.volum
        self.B = self.Cr.cross(self.Ar)*self.volum
        self.C = self.Ar.cross(self.Br)*self.volum

        for line in xtail:
            lineSplit = line.split()
            if line.upper().count("SPACE GROUP"):
                self.spg = lineSplit[2]
            elif line.upper().count("SPINDLE AXIS"):
                self.spindleAxis = map(int,lineSplit[2:5])
                self.verticalAxis = map(int,lineSplit[7:])
            elif line.upper().count("MOTOR AXIS"):
                self.motorAxis = map(float,lineSplit[2:5])
            elif line.upper().count("DISTANCE"):
                self.distance = float(lineSplit[1])
            elif line.upper().count("X BEAM"):
                self.beam_x = float(lineSplit[2])
                self.beam_y = float(lineSplit[5])
            elif line.upper().count("SECTOR"):
                self.sector = int(lineSplit[1])
            elif line.upper().count("RAW DATA FILE"):
                self.template = str(lineSplit[-1]).replace("'","")
            elif line.upper().count("UNIT CELL"):
                self.cell2 = map(float,lineSplit[2:])
                # Verify that the cell extracted from UB correspond
                # to the cell read from the xfile tail
                assert abs(self.cell2[0] - self.cell[0]) < 1e-2 and \
                       abs(self.cell2[1] - self.cell[1]) < 1e-2 and \
                       abs(self.cell2[2] - self.cell[2]) < 1e-2 and \
                       abs(self.cell2[3] - self.cell[3]) < 2e-2 and \
                       abs(self.cell2[4] - self.cell[4]) < 2e-2 and \
                       abs(self.cell2[5] - self.cell[5]) < 2e-2 

        # Verify that the calculation method for UB_to_Rotxyz works correctly
        _rotx, _roty, _rotz = self.UB_to_Rotxyz()
        assert abs(_rotx - self.rotx) < 2e-2 and \
               abs(_roty - self.roty) < 2e-2 and \
               abs(_rotz - self.rotz) < 2e-2

        # Verify that the calculation method for Adnz_to_Udnz works correctly
        _U = self.Adnz_to_Udnz()
        print diffMAT(_U, self.U)
        assert diffMAT(_U, self.U) < 5e-6

    def get_U0(self, rcell=None, vertical=None, spindle=None, clean=False):
        "Calculate denzo U0 from spindle, verctical"
        if not rcell: rcell = self.cell_r
        if not vertical: vertical = self.verticalAxis
        if not spindle: spindle = self.spindleAxis

        Bmat = self.get_B(rcell)
        vertical = vec3(vertical)
        spindle = vec3(spindle)

        U0y = (Bmat * spindle).normalize()
        U0xi = Bmat * vertical
        U0x = (U0xi - (U0xi * U0y) * U0y).normalize()

        U0 = mat3(U0x, U0y, U0x.cross(U0y)).transpose()

        # cleaning... Just cosmetic, not realy needed.
        if clean: U0 = cleanU0(U0)
        return U0

    def get_B(self, rcell=None):
        """ Denzo Orthogonalisation matrix
          b* is aligned with spindle axis (2-nd coordinate)
          a* is in the plane perpendicular to the beam (1-st,2-nd coords)
        """
        if not rcell:
            rcell = self.cell_r
        sr = map(sind, rcell[3:6])
        cr = map(cosd, rcell[3:6])
        B  = mat3()
        B13 = (cr[1]-cr[2]*cr[0])/sr[2]
        B[0,0] = rcell[0] * sr[2]
        B[1,0] = rcell[0] * cr[2]
        B[1,1] = rcell[1]
        B[0,2] = rcell[2] * B13
        B[1,2] = rcell[2] * cr[0]
        B[2,2] = rcell[2] * (sr[0]**2 - B13**2)**0.5
        
        return B

    def UB_to_Rotxyz(self, UB=None, U0=None, vertical=None,
                                    spindle=None, motorAxis=None):
        
        if not UB: UB = self.UB
        if not vertical: vertical = self.verticalAxis
        if not spindle: spindle = self.spindleAxis
        if not motorAxis: motorAxis = self.motorAxis
        
        rcell = UB_to_cellParam(UB)
        if not U0: U0 = self.get_U0(rcell, vertical, spindle)
        B = self.get_B(rcell)
                
        U0B = U0 * B
        # U = UB.B**-1
        U = UB * U0B.inverse()
        verif = is_orthogonal(U)
        if not verif:
            print "???  Warning: The U matrix is not orthogonal."
            print "???  Epsilon error: %.1e" % verif

        if  motorAxis == [0., -1., 0.]:  self.DNZAxes =  -ey,  ex, -ez
        elif motorAxis == [0., 1., 0.]:  self.DNZAxes =   ey, -ex, -ez
            
        angles_pair = ThreeAxisRotation2(U.mlist, self.DNZAxes).getAngles()
        return map_r2d(angles_pair[0])
         
    def Adnz_to_Udnz(self, Adnz=None):
        """!!!! Not well that Adnz and Udnz are permutated!!!
        This function only tries to follow dot.x header internal rules...
        """
        if not Adnz:
            Adnz = self.UB
            cell_r = self.cell_r
        else:
            cell_r = UB_to_cellParam(Adnz)
        BL = BusingLevy(cell_r)
        return Qdnz2mos * (Adnz * BL.inverse())
         
    def Adnz_to_rotxyz(self, Adnz, Bdnz, vertical, spindle):
        U2a = (Bdnz * spindle).normalize()
        U1a = Bdnz * vertical
        U1b = (U1a - (U1a * U2a) * U2a).normalize()
        Udnz0 = mat3(U1b, U2a, U1b.cross(U2a)).transpose()
        Adnz0 = Udnz0 * Bdnz
        RMAT = Adnz * Adnz0.inverse()
        dnzrmat = ThreeAxisRotation2(RMAT.mlist, self.DNZAxes)
        rotxyz = dnzrmat.getAngles()
        return rotxyz, Udnz0
        
#    def RotXYZ2UB(self):
#         
#        return self.get_U0(self.cell_r)

class XDSParser:

    def __init__(self, filename=None):
        
        # XDS Dictionary
        self.dict = {}
        self.info = "XDS Parser"
        self.fileType = "XDS"
        
        if filename:
            self.parse(filename)
            if self.dict.has_key("symmetry"):
                self.spaceGroupNumber = self.dict["symmetry"]
                self.spaceGroupName = spg_num2symb[self.dict["symmetry"]]
            if self.dict.has_key("cell"):
                self.cell = self.dict["cell"]
                self.cell_r = reciprocal(self.dict["cell"])
    
    def parse(self, filename):
        if filename.count("XDS.INP"): self.parse_xdsinp(filename)
        elif filename.count("INIT.LP"): self.parse_init(filename)
        elif filename.count("CORRECT.LP"): self.parse_correct(filename)
        elif filename.count("IDXREF.LP"): self.parse_idxref(filename)
        elif filename.count("XPARM.XDS"): self.parse_xparm(filename)
        else:
            raise ParserError, "Error, Can't parse file: %s" % filename
        
    def get_par(self, _str,match,limit=70,func=float):
        start = _str.index(match)+len(match)
        tmp = _str[start:start+limit].splitlines()[0].split()
        return map(func,tmp)

    def parse_xdsinp(self, infname="XDS.INP"):
        xinp = openReadClose(infname)
        self.dict["detector_type"] = self.get_par(xinp,"DETECTOR=",40,str)[0]
        self.dict["template"] = self.get_par(xinp, "_OF_DATA_FRAMES=",60,str)[0].replace("?","#")

    def parse_init(self, infname="INIT.LP", infilelocation="."):
        # Extract things from INIT.LP
        try:
            init = openReadClose(os.path.join(infilelocation, infname))
            self.dict["mean_gain"] = self.get_par(init,"MEAN GAIN VALUE",40)[0]
        except:
            pass

    def parse_correct(self, infname="CORRECT.LP", infilelocation="."):
        "Extract information from XDS output CORRECT.LP and INIT.LP"
        # Extract things from CORRECT.LP
        corr = openReadClose(os.path.join(infilelocation, infname))
        ip = corr.index("PARAMETERS USING ALL IMAGES")+100
        corrp = corr[ip:ip+1400]
        corri = corr[:1500]
        self.dict["rot"] = self.get_par(corrp,"ROTATION AXIS")
        self.dict["beam"] = self.get_par(corrp,"COORDINATES (REC. ANGSTROEM)")
        self.dict["distance"] = self.get_par(corrp,"DETECTOR DISTANCE (mm)")[0]
        self.dict["origin"] = self.get_par(corrp,"(PIXELS) OF DIRECT BEAM")
        self.dict["originXDS"] = self.get_par(corrp,"ORIGIN (PIXELS) AT ")
        self.dict["A"] = self.get_par(corrp,"CELL A-AXIS")
        self.dict["B"] = self.get_par(corrp,"CELL B-AXIS")
        self.dict["C"] = self.get_par(corrp,"CELL C-AXIS")
        self.dict["cell"] = self.get_par(corrp,"UNIT CELL PARAMETERS")
        self.dict["mosaicity"] = self.get_par(corrp,"CRYSTAL MOSAICITY (DEGREES)")[0]
        iqx, iqy  = corri.index("QX=")+3, corri.index("QY=")+3
        inx, iny  = corri.index("NX=")+3, corri.index("NY=")+3
        self.dict["pixel_size"] = float(corri[iqx:iqx+9]),float(corri[iqy:iqy+9])
        self.dict["pixel_numb"] = int(corri[inx:inx+7]),int(corri[iny:iny+7])
        self.dict["template"] = self.get_par(corri, "_DATA_FRAMES=",60,str)[0].replace("?","#")
        self.dict["symmetry"] = int(self.get_par(corrp,"SPACE GROUP NUMBER")[0])
        self.dict["detector_type"] = self.get_par(corri,"DETECTOR=",40,str)[0]
        self.dict["detector_X"] = self.get_par(corri,"DETECTOR_X-AXIS=")
        self.dict["detector_Y"] = self.get_par(corri,"DETECTOR_Y-AXIS=")
        self.dict["phi_init"] = self.get_par(corri,"STARTING_ANGLE=",15)[0]
        self.dict["num_init"] = self.get_par(corri,"STARTING_FRAME=",15,int)[0]
        self.dict["delta_phi"] = self.get_par(corri,"OSCILLATION_RANGE=",15)[0]
        self.dict["divergence_esd"] = self.get_par(corri,"BEAM_DIVERGENCE_E.S.D.=",15)[0]
        self.dict["resolution_range"] = self.get_par(corri,"INCLUDE_RESOLUTION_RANGE=",20)
        self.dict["friedel"] = self.get_par(corri,"FRIEDEL'S_LAW=",7,str)[0]
        self.dict["polarization"] = self.get_par(corri,"FRACTION_OF_POLARIZATION=",8)[0]

    def parse_idxref(self, infname="IDXREF.LP"):

        idxr = openReadClose(infname)
        idxrp = idxr[idxr.index("FRACTION PARAMETERS USED AT START OF INTEGRAT"):]
        idxri = idxr[:1200]
        self.dict["rot"] = self.get_par(idxrp,"ROTATION AXIS")
        self.dict["beam"] = self.get_par(idxrp,"COORDINATES (REC. ANGSTROEM)")
        self.dict["distance"] = self.get_par(idxrp,"DETECTOR DISTANCE (mm)")[0]
        self.dict["origin"] = self.get_par(idxrp,"(PIXELS) OF DIRECT BEAM")
        self.dict["originXDS"] = self.get_par(idxrp,"ORIGIN (PIXELS) AT ")
        self.dict["A"] = self.get_par(idxrp,"CELL A-AXIS")
        self.dict["B"] = self.get_par(idxrp,"CELL B-AXIS")
        self.dict["C"] = self.get_par(idxrp,"CELL C-AXIS")
        self.dict["cell"] = self.get_par(idxrp,"UNIT CELL PARAMETERS")
        iqx, iqy  = idxri.index("QX=")+3, idxri.index("QY=")+3
        self.dict["pixel_size"] = float(idxri[iqx:iqx+9]),float(idxri[iqy:iqy+9])
        inx, iny  = idxri.index("NX=")+3, idxri.index("NY=")+3
        self.dict["pixel_numb"] = int(idxri[inx:inx+8]),int(idxri[iny:iny+8])
        self.dict["template"] = self.get_par(idxri, "_DATA_FRAMES=",60,str)[0].replace("?","#")
        self.dict["symmetry"] = int(self.get_par(idxrp,"SPACE GROUP NUMBER")[0])
        self.dict["detector_X"] = self.get_par(idxri,"DETECTOR_X-AXIS=")
        self.dict["detector_Y"] = self.get_par(idxri,"DETECTOR_Y-AXIS=")
        self.dict["phi_init"] = self.get_par(idxri,"STARTING_ANGLE=",15)[0]
        self.dict["num_init"] = self.get_par(idxri,"STARTING_FRAME=",15,int)[0]
        self.dict["delta_phi"] = self.get_par(idxri,"OSCILLATION_RANGE=",15)[0]

    def parse_xparm(self, infname="XPARM.XDS"):

        xparm = openReadClose(infname).splitlines()
        if len(xparm[0].split()) == 3:
            self.dict["rot"] = map(float,xparm[0].split())
        if len(xparm[0].split()) == 6:
            self.dict["rot"] = map(float,xparm[0].split()[3:])
        self.dict["beam"] = map(float,xparm[1].split()[1:])
        self.dict["distance"] = float(xparm[3].split()[0])
        self.dict["originXDS"] = map(float,xparm[3].split()[1:])
        self.dict["A"] = map(float,xparm[8].split())
        self.dict["B"] = map(float,xparm[9].split())
        self.dict["C"] = map(float,xparm[10].split())
        self.dict["cell"] = map(float,xparm[7].split()[1:])
        self.dict["pixel_size"] = map(float,xparm[2].split()[2:])
        self.dict["pixel_numb"] = map(float,xparm[2].split()[:2])
        self.dict["symmetry"] = int(xparm[7].split()[0])
        self.dict["num_init"], self.dict["phi_init"], self.dict["delta_phi"] = \
                   map(float,xparm[0].split()[:3])
        self.dict["detector_X"] = map(float,xparm[4].split())
        self.dict["detector_Y"] = map(float,xparm[5].split())

    def debut(self):
        "Do simple cristallographic calculations from XDS initial parameters"
        
        A = vec3(self.dict["A"])
        B = vec3(self.dict["B"])
        C = vec3(self.dict["C"])
        
        volum = A.cross(B)*C
        Ar = B.cross(C)/volum
        Br = C.cross(A)/volum
        Cr = A.cross(B)/volum
        UBxds = mat3(Ar,Br,Cr)

        BEAM = vec3(self.dict["beam"])
        wavelength = 1/BEAM.length()
        
        self.dict["cell_volum"] = volum
        self.dict["wavelength"] = wavelength
        self.dict["Ar"] = Ar
        self.dict["Br"] = Br
        self.dict["Cr"] = Cr
        self.dict["UB"] = UBxds
        

    def getOmega(self):
        """Calculate an Omega value (in radian) wich defines how the fast (X) and
        slow (Y) axis of detector files are  orientated toward the camera frame.
        The calculation of this omega value is supposed to reflect the Mosflm
        definition...

        But it seems that I get different values from the mosflm defaults... This
        may be due to: A) My missanderstanding of the mosflm documentation, B)
        Some tricks in the image reading routines.

        Nonetheless, this calculated value works for translating
        correctly the beam coordinates from XDS to mosflm [at least in the tested
        cases of MARCCD, MAR345 and ADSC detector images].

        Reference:
        http://www.ccp4.ac.uk/dist/x-windows/Mosflm/doc/mosflm_user_guide.html#a3
        """
        # Xd = CAMERA_y = beam
        # Yd = CAMERA_z = rot
        Xd =  vec3(self.dict["beam"]).normalize()
        Yd =  vec3(self.dict["rot"]).normalize()
        CAMERA_x = Xd.cross(Yd)
        CAMERA = mat3(CAMERA_x, Xd, Yd).transpose()
            
        # This is the definition of the fast:X and slow:Y axis for the detector files.
        XDSdetector_X = vec3(self.dict["detector_X"])
        XDSdetector_Y = vec3(self.dict["detector_Y"])
                
        # Now this axes are translated in the mosflm Camera frame
        Xs = XDSdetector_X*CAMERA
        Ys = XDSdetector_Y*CAMERA
        
        # Both angles should be identical.
        omegaX = Xd.angle(Xs)
        omegaY = Yd.angle(Ys)
        
        if _debug:
            print "DEBUG: X xds: fast =",XDSdetector_X
            print "DEBUG: Y xds: slow =",XDSdetector_Y
            print "DEBUG: Xs: fast =",XDSdetector_X,"->", Xs
            print "DEBUG: Ys: slow =",XDSdetector_Y,"->", Ys
            print "DEBUG: Xd: ", Xd
            print "DEBUG: Yd: ", Yd
            print "DEBUG: OmegaX:   %8.2f" % (omegaX*r2d)
            print "DEBUG: OmegaY:   %8.2f" % (omegaY*r2d)
            
        return omegaX

    def getTwoTheta(self):
        """Tries to calculate the 2theta angle (in radian) of the detector.
        I am not completely sure of this calculation. How 2theta is precisely
        geometricaly defined in mosflm ?
        I need to look in the mosflm code where it is taken into account.
        """
        BEAM = vec3(self.dict["beam"])
        ROT  = vec3(self.dict["rot"]).normalize()
        camY = ROT.cross(BEAM)

        XDSdetector_X = vec3(self.dict["detector_X"]).normalize()
        XDSdetector_Y = vec3(self.dict["detector_Y"]).normalize()
        #XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)

        #print beam.angle(XDSdetector_Z)*r2d
        if abs(ROT * XDSdetector_X) - 1 <= 0.05:
            detecorVector = -XDSdetector_Y
            #print 1
        elif abs(ROT * XDSdetector_Y) - 1 <= 0.05:
            detecorVector = XDSdetector_X
            #print 2
        else:
            raise Exception, "Can't calculate TwoTheta angle"    
        return camY.angle(detecorVector)    

    def getBeamOrigin(self):
        """Calculate the direct beam coordinates on the detector from beamOrigin."""

        distance = self.dict["distance"]
        beam = vec3(self.dict["beam"])

        XDSdetector_X = vec3(self.dict["detector_X"])
        XDSdetector_Y = vec3(self.dict["detector_Y"])
        XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)

        # Calculate the direct beam coordinates on the detector
        beamOx = self.dict["originXDS"][0]*self.dict["pixel_size"][0]
        beamOy = self.dict["originXDS"][1]*self.dict["pixel_size"][1]
        beamOz = beam*XDSdetector_Z

        beamX = beamOx + beam*XDSdetector_X*distance/beamOz
        beamY = beamOy + beam*XDSdetector_Y*distance/beamOz
        beamXp = beamX/self.dict["pixel_size"][0]
        beamYp = beamY/self.dict["pixel_size"][1]
        
        if _debug:
            if "origin" in self.dict.keys():
                print "\nDEBUG: BEAM center read from XDS in pixel:",
                print  "%9.2f %9.2f" % tuple(self.dict["origin"])
                # When given by XDS, verifies that my calculation is correct
                assert (self.dict["origin"][0] - beamXp) < 0.05
                assert (self.dict["origin"][1] - beamYp) < 0.05
            print "DEBUG: BEAM center calculated in pixel:\t%9.2f %9.2f" % (beamXp,beamYp)
            print "DEBUG: BEAM center calculated in mm:\t\t%9.2f %9.2f\n" % (beamX,beamY)
        return beamX, beamY

    def get_beam_coordinate(self):
        """Calculate the direct beam coordinates on the detector
           from beam origin."""

        distance = self.dict["distance"]
        beam = vec3(self.dict["beam"])

        XDSdetector_X = vec3(self.dict["detector_X"])
        XDSdetector_Y = vec3(self.dict["detector_Y"])
        XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)

        # Calculate the direct beam coordinates on the detector
        beamOx = self.dict["originXDS"][0]*self.dict["pixel_size"][0]
        beamOy = self.dict["originXDS"][1]*self.dict["pixel_size"][1]
        beamOz = beam*XDSdetector_Z

        beamX = beamOx + beam*XDSdetector_X*distance/beamOz
        beamY = beamOy + beam*XDSdetector_Y*distance/beamOz
        beamXp = beamX/self.dict["pixel_size"][0]
        beamYp = beamY/self.dict["pixel_size"][1]

        if _debug:
            if "origin" in self.dict.keys():
                print "\nDEBUG: BEAM center read from XDS in pixel:",
                print  "%9.2f %9.2f" % tuple(self.dict["origin"])
                # When given by XDS, verifies that my calculation is correct
                assert (self.dict["origin"][0] - beamXp) < 0.05
                assert (self.dict["origin"][1] - beamYp) < 0.05
            print "DEBUG: BEAM center calculated in pixel:\t%9.2f %9.2f" % (beamXp,beamYp)
            print "DEBUG: BEAM center calculated in mm:\t\t%9.2f %9.2f\n" % (beamX,beamY)
        return beamX, beamY

    def get_beam_origin(self):
        """Calculate the direct beam origin from the
           detector beam coordinate."""

        distance = self.dict["distance"]
        beam = vec3(self.dict["beam"])

        XDSdetector_X = vec3(self.dict["detector_X"])
        XDSdetector_Y = vec3(self.dict["detector_Y"])
        XDSdetector_Z = XDSdetector_X.cross(XDSdetector_Y)

        # Calculate the direct beam coordinates on the detector
        beamCx = self.dict["origin"][0]*self.dict["pixel_size"][0]
        beamCy = self.dict["origin"][1]*self.dict["pixel_size"][1]
        beamCz = beam*XDSdetector_Z

        beamX = beamCx - beam*XDSdetector_X*distance/beamCz
        beamY = beamCy - beam*XDSdetector_Y*distance/beamCz
        beamXp = beamX/self.dict["pixel_size"][0]
        beamYp = beamY/self.dict["pixel_size"][1]

        if _debug:
            if "origin" in self.dict.keys():
                print "\nDEBUG: BEAM center read from XDS in pixel:",
                print  "%9.2f %9.2f" % tuple(self.dict["origin"])
                # When given by XDS, verifies that my calculation is correct
                assert (self.dict["origin"][0] - beamXp) < 0.05
                assert (self.dict["origin"][1] - beamYp) < 0.05
            print "DEBUG: BEAM center calculated in pixel:\t%9.2f %9.2f" % (beamXp,beamYp)
            print "DEBUG: BEAM center calculated in mm:\t\t%9.2f %9.2f\n" % (beamX,beamY)
        return beamXp, beamYp


    def UBxds_to_mos(self):
        """ Convert the XDS direct space Orientation Matrix to a mosflm OM
        
        Mosflm CAMERA coordinate frame has orthonormal axes with:
        
          z // rotation axis
          y perpendicular to z and to the beam
          x perpendicular to y and z (along the beam)
        
        For more details see the mosflm documentation:
        http://www.ccp4.ac.uk/dist/x-windows/Mosflm/doc/mosflm_user_guide.html#a3
        """
        
        if "UB" not in self.dict.keys():
            self.debut()
            
        BEAM = vec3(self.dict["beam"])
        ROT = vec3(self.dict["rot"])
        UBxds = self.dict["UB"]
        
        CAMERA_z = ROT.normalize()
        CAMERA_y = CAMERA_z.cross(BEAM).normalize()
        CAMERA_x = CAMERA_y.cross(CAMERA_z)
        CAMERA = mat3(CAMERA_x,CAMERA_y,CAMERA_z).transpose()
        
        return  CAMERA * UBxds * self.dict["wavelength"]
         

    def UBxds_to_dnz(self):
        """ Convert the XDS direct space Orientation Matrix to a mosflm OM
        
        Denzo CAMERA coordinate frame has orthonormal axes with:

          y // to the rotation (spindel) axis
          z // to the beam
          x perpendicular to z and to the beam
        
        For more details see the denzo documentation:
        http://www.ccp4.ac.uk/dist/x-windows/Mosflm/doc/mosflm_user_guide.html#a3
        """
            
        if "UB" not in self.dict.keys():
            self.debut()
            
        BEAM = vec3(self.dict["beam"])
        ROT = vec3(self.dict["rot"])
        UBxds = self.dict["UB"]
                
        CAMERA_y = ROT.normalize()
        CAMERA_x = CAMERA_y.cross(BEAM).normalize()
        CAMERA_z = CAMERA_x.cross(CAMERA_y)
        CAMERA = mat3(CAMERA_x,CAMERA_y,CAMERA_z).transpose()
                    
        return  CAMERA * UBxds

class MosflmParser:

    def __init__(self, filename=None):
        self.dict = {}
        self.info = "Mosflm Parser"
        self.fileType = "Mosflm"
        if filename:
            self.parse_umat(filename)
            self.cell_r = reciprocal(self.cell)
            self.B = self.get_B()
            self.spaceGroupNumber = 1
            self.spaceGroupName = 'P1'

    def parse_umat(self, infname):
        try:
            mosFile = map(float,open(infname).read().split())
        except:
            raise Exception, "Error! Can't parse Mosflm matrices file."
        
        self.UB = mat3(mosFile[:9])
        self.Ur = mat3(mosFile[12:21])
        # Ur = reference orientation matrix (corresponding to setting
        # angles "missetingAngles" given below)
        
        self.cell = mosFile[21:27]
        self.missetingAngles = mosFile[27:30]
        self.dummy = mosFile[9:12]
        
        MatXYZ = ThreeAxisRotation2(map_d2r(self.missetingAngles)).tensor
        self.U = MatXYZ * self.Ur

    def write_umat(self, filename):
        "Writes a mosflm formated crystal oritentation matrice file"
        
        matf = open(filename,"w")
        fmt_cell = 6 * "%12.4f" + "\n"
        fmt_angles = 3 * "%12.3f" + "\n"
        matf.write(str_mat(self.UB, name="", format="%12.8f"))
        matf.write(fmt_angles % (0,0,0))
        matf.write(str_mat(self.U, name="", format="%12.8f"))
        matf.write(fmt_cell % tuple(self.cell))
        matf.write(fmt_angles % tuple(self.missetingAngles))
        matf.close()
    
    def UB_to_wavelength(self, UB=None, cell=None):
        """Extracting the wavelength from UBmos and the cell parameters
           UBmos = U * B * wavelength"""
        
        if not UB: UB = self.UB 
        if not cell: cell = self.cell
        
        wcell = reciprocal(UB_to_cellParam(UB))
        w3 = map(lambda t: t[0]/t[1], zip(cell[:3], wcell[:3]))
        if (abs(w3[0] - w3[1]) < 1e-4) and (abs(w3[0] - w3[2]) < 1e-4):
            return (w3[0]+w3[1]+w3[2])/3.
        else:
            raise Exception, 'Error in extracting wavelength from UB and cell!'

    def get_B(self, rcell=None):
        " Mosflm Orthogonalisation matrix (Busing & Levy)"
        if not rcell:  rcell = self.cell_r
        return BusingLevy(rcell)

def openReadClose(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def openWriteClose(filename,text):
    f = open(filename,'w')
    f.write(text)
    f.close()

def diffMAT(m1,m2):
    diffsq = 0
    for a in (m1-m2):
        diffsq += a[0]**2+a[1]**2+a[2]**2
    return diffsq**0.5

def rootSquareSum(mdiff):
    if isinstance(mdiff, mat3):
        return reduce(lambda x, y: x+y, map(lambda x: x**2, mdiff.mlist))**0.5
    else:
        return reduce(lambda x, y: x+y, map(lambda x: x**2, mdiff))**0.5

def printmat(mat, name="", format="%12.8f"):
    if name: print "%s" % (name)
    if isinstance(mat, mat3):
        for i in 0,1,2: print 3*format % tuple(mat.getRow(i))
    else:
        for l in mat: print 3*format % tuple(l)

def str_mat(mat, name="", format="%12.6f"):
    "Return a string representation of a 3x3 matrix"
    Pformat = name + 3*(3*format+"\n")
    return Pformat % tuple(mat.toList(1))

def volum(cell):
    """
    Calculate the cell volum from either:
     - the 6 standard cell parameters (a, b, c, alpha, beta, gamma)
     - or the 3 vectors A, B, C
    """
    if (len(cell) == 6) and (type(cell[0]) == float):
        # expect a, b, c, alpha, beta, gamma (angles in degree).
        ca, cb, cg = map(cosd, cell[3:6])
        return cell[0]*cell[1]*cell[2]*(1-ca**2-cb**2-cg**2+2*ca*cb*cg)**0.5
    elif (len(cell) == 3) and isinstance(cell[0], vec3):
        # expect vectors of the 3 cell parameters
        A, B, C = cell
        return A*B.cross(C)
    else:
        print "error in volum()"
        return "Can't parse input arguments."  

def reciprocal(cell):
    "Calculate the 6 reciprocal cell parameters: a*, b*, c*, alpha*, beta*..."
    sa, sb, sg = map(sind, cell[3:6])
    ca, cb, cg = map(cosd, cell[3:6])
    v = volum(cell)
    rc = (cell[1]*cell[2]*sa/v, 
          cell[2]*cell[0]*sb/v,
          cell[0]*cell[1]*sg/v,
          math.acos((cb*cg-ca)/(sb*sg)) * r2d,
          math.acos((ca*cg-cb)/(sa*sg)) * r2d,
          math.acos((ca*cb-cg)/(sa*sb)) * r2d)
    return rc

def UB_to_cellParam(UB):
    """Return an array containing the cell parameters with angles en degree
    >>> ub = mat3(0.0045624910668708527, 0.0013380296069423175, -0.0019732516590096985,
                  0.0014703215926493108, 0.0037937417049515054, 0.0057564982133741704,
                 7.3231240428790203e-05, -0.002607820316488004, 0.007361827462991322)
    >>> print UB_to_cellParam(ub)
    """
    Ar = vec3(UB.getColumn(0))
    Br = vec3(UB.getColumn(1))
    Cr = vec3(UB.getColumn(2))
    return (Ar.length(), Br.length(), Cr.length(),
            Br.angle(Cr)*r2d, Cr.angle(Ar)*r2d, Ar.angle(Br)*r2d)

def is_orthogonal(R, epsilon=5e-6, debug=True):
    """Check if a tensor is orthogonal:R.transpose() = R.inverse()
       Usefull to check proper caclculation of U orientation matrices."""
    Mdiff = diffMAT(R.transpose(), R.inverse())
    print "==debug:  Mdiff= %.2e    abs(1-det)= %.2e" % \
                                      (Mdiff, abs(1-R.determinant()))
    #print "\n Rotation error =====>> %7.1e\n" % Mdiff
    if Mdiff < epsilon:
        return True
    else:
        return Mdiff


def BusingLevy(rcell):
    cosr = map(cosd, rcell[3:6])
    sinr = map(sind, rcell[3:6])
    Vr = volum(rcell)
    BX = ex*rcell[0]
    BY = rcell[1]*(ex*cosr[2] + ey*sinr[2])
    c = rcell[0]*rcell[1]*sinr[2]/Vr
    cosAlpha = (cosr[1]*cosr[2] - cosr[0])/(sinr[1]*sinr[2])
    BZ = vec3([rcell[2]*cosr[1],
                -1*rcell[2]*sinr[1]*cosAlpha,
                1/c])
    return mat3(BX,BY,BZ)

def cleanU0(Umat,resid=3e-6):
    "cleaning... Just cosmetic, not realy needed."
    sin120 = math.sin(2*pi/3)
    for i in range(3):
        for j in range(3):
            v = Umat[i,j]
            vsign = cmp(v, 1e-30)
            if abs(v) < resid: Umat[i,j] = 0.
            elif abs(v - 0.5) < resid: Umat[i,j] = 0.5*vsign
            elif abs(v - sin120) < resid: Umat[i,j] = sin120*vsign
    return Umat

def getPermutUB(PGoperators, UB, _epsilonCell=1e-2):
    "Apply Point group permutations to UB, return a list of permutated UB"

    cell = reciprocal(UB_to_cellParam(UB))
    permutedList = []

    for equiv in PGoperators:
        # Note that equiv_mat is the transpose of the normal equiv matrix!
        # because of the way mat3 matrices are contructed.
        equiv_mat = mat3(equiv[0],equiv[1],equiv[2])
        new_UB = UB * equiv_mat

        # One simple way to verify that the permutation is correct is
        # to extract the cell parameters from the permuted UB matrix
        if _debug:
            new_cell = reciprocal(UB_to_cellParam(new_UB))
            assert abs(new_cell[0] - cell[0]) < _epsilonCell and \
                   abs(new_cell[1] - cell[1]) < _epsilonCell and \
                   abs(new_cell[2] - cell[2]) < _epsilonCell and \
                   abs(new_cell[3] - cell[3]) < _epsilonCell and \
                   abs(new_cell[4] - cell[4]) < _epsilonCell and \
                   abs(new_cell[5] - cell[5]) < _epsilonCell 

        permutedList.append(new_UB)
    return permutedList

def getPermutU(PGoperators, U, _epsilonCell=1e-4, debug=True):
    "Apply Point group permutations to U, return a list of permutated U"

    permutedList = []
    for equiv in PGoperators:
        # Note that equiv_mat is the transpose of the normal equiv matrix!
        # because of the way mat3 matrices are contructed.
        equiv_mat = mat3(equiv[0],equiv[1],equiv[2])
        new_U = equiv_mat.transpose() * U
        if is_orthogonal(new_U):
            permutedList.append(new_U)
        else:
            print "Internal Error: permuted U matrix is not orthogonal !"
            print new_U
    return permutedList

def _test():

    ub = mat3(0.00195366,  0.01690921, -0.00112061,
              0.00676745, -0.00440955, -0.00560790,
             -0.00585598,  0.00054534, -0.01059838)
    ub = ub/0.93100
    u = mat3(0.2132794,   0.9671680,  -0.1381950,
             0.7387952,  -0.2522162,  -0.6249549,
            -0.6392915,   0.0311922,  -0.7683315)

    cell = (103.7050, 53.2511, 78.8808, 90.0000, 101.4632, 90.0000)

    print diffMAT(ub.decompose()[0], u)
    print ub.decompose()[0] - u
    print cell
    print reciprocal(UB_to_cellParam(ub))

spg_num2symb = {1:'P1',3:'P2',4:'P21',5:'C2',16:'P222',17:'P2221',18:'P21212',
19:'P212121',20:'C2221',21:'C222',22:'F222',23:'I222',24:'I212121',75:'P4',
76:'P41',77:'P42',78:'P43',79:'I4',80:'I41',89:'P422',90:'P4212',91:'P4122',
92:'P41212',93:'P4222',94:'P42212',95:'P4322',96:'P43212',97:'I422',98:'I4122',
143:'P3',144:'P31',145:'P32',146:'R3',149:'P312',150:'P321',151:'P3112',
152:'P3121',153:'P3212',154:'P3221',155:'R32',168:'P6',169:'P61',170:'P65',
171:'P62',172:'P64',173:'P63',177:'P622',178:'P6122',179:'P6522',180:'P6222',
181:'P6422',182:'P6322',195:'P23',196:'F23',197:'I23',198:'P213',199:'I213',
207:'P432',208:'P4232',209:'F432',210:'F4132',211:'I432',212:'P4332',
213:'P4132',214:'I4132'}

SPGlib2 = {'p1':1,'p2':3,'p21':4,'c2':5,'p222':16,'p2221':17,'p21212':18,
'p212121':19,'c2221':20,'c222':21,'f222':22,'i222':23,'i212121':24,'p4':75,
'p41':76,'p42':77,'p43':78,'i4':79,'i41':80,'p422':89,'p4212':90,'p4122':91,
'p41212':92,'p4222':93,'p42212':94,'p4322':95,'p43212':96,'i422':97,'i4122':98,
'p3':143,'p31':144,'p32':145,'r3':146,'p312':149,'p321':150,'p3112':151,
'p3121':152,'p3212':153,'p3221':154,'r32':155,'p6':168,'p61':169,'p65':170,
'p62':171,'p64':172,'p63':173,'p622':177,'p6122':178,'p6522':179,'p6222':180,
'p6422':181,'p6322':182,'p23':195,'f23':196,'i23':197,'p213':198,'i213':199,
'p432':207,'p4232':208,'f432':209,'f4132':210,'i432':211,'p4332':212,
'p4132':213,'i4132':214}

SPGlib = {
# each SGnumber entry maps (SGsymbol1,SGsymbol2,SGorder,PGsymbol)
 1: ('P1', 'P1', 1, '1'),
 3: ('P2', 'P2', 2, '2'),
 4: ('P21', 'P2(1)', 2, '2'),
 5: ('C2', 'C2', 4, '2'),
 16: ('P222', 'P222', 4, '222'),
 17: ('P2221', 'P222(1)', 4, '222'),
 18: ('P21212', 'P2(1)2(1)2', 4, '222'),
 19: ('P212121', 'P2(1)2(1)2(1)', 4, '222'),
 20: ('C2221', 'C222(1)', 8, '222'),
 21: ('C222', 'C222', 8, '222'),
 22: ('F222', 'F222', 16, '222'),
 23: ('I222', 'I222', 8, '222'),
 24: ('I212121', 'I2(1)2(1)2(1)', 8, '222'),
 75: ('P4', 'P4', 4, '4'),
 76: ('P41', 'P4(1)', 4, '4'),
 77: ('P42', 'P4(2)', 4, '4'),
 78: ('P43', 'P4(3)', 4, '4'),
 79: ('I4', 'I4', 8, '4'),
 80: ('I41', 'I4(1)', 8, '4'),
 89: ('P422', 'P422', 8, '422'),
 90: ('P4212', 'P42(1)2', 8, '422'),
 91: ('P4122', 'P4(1)22', 8, '422'),
 92: ('P41212', 'P4(1)2(1)2', 8, '422'),
 93: ('P4222', 'P4(2)22', 8, '422'),
 94: ('P42212', 'P4(2)2(1)2', 8, '422'),
 95: ('P4322', 'P4(3)22', 8, '422'),
 96: ('P43212', 'P4(3)2(1)2', 8, '422'),
 97: ('I422', 'I422', 16, '422'),
 98: ('I4122', 'I4(1)22', 16, '422'),
 143: ('P3', 'P3', 3, '3'),
 144: ('P31', 'P3(1)', 3, '3'),
 145: ('P32', 'P3(2)', 3, '3'),
 146: ('R3', 'R3', 9, '3'),
 149: ('P312', 'P312', 6, '312'),
 150: ('P321', 'P321', 6, '321'),
 151: ('P3112', 'P3(1)12', 6, '312'),
 152: ('P3121', 'P3(1)21', 6, '321'),
 153: ('P3212', 'P3(2)12', 6, '312'),
 154: ('P3221', 'P3(2)21', 6, '321'),
 155: ('R32', 'R32', 18, '321'),
 168: ('P6', 'P6', 6, '6'),
 169: ('P61', 'P6(1)', 6, '6'),
 170: ('P65', 'P6(5)', 6, '6'),
 171: ('P62', 'P6(2)', 6, '6'),
 172: ('P64', 'P6(4)', 6, '6'),
 173: ('P63', 'P6(3)', 6, '6'),
 177: ('P622', 'P622', 12, '622'),
 178: ('P6122', 'P6(1)22', 12, '622'),
 179: ('P6522', 'P6(5)22', 12, '622'),
 180: ('P6222', 'P6(2)22', 12, '622'),
 181: ('P6422', 'P6(4)22', 12, '622'),
 182: ('P6322', 'P6(3)22', 12, '622'),
 195: ('P23', 'P23', 12, '23'),
 196: ('F23', 'F23', 48, '23'),
 197: ('I23', 'I23', 24, '23'),
 198: ('P213', 'P2(1)3', 12, '23'),
 199: ('I213', 'I2(1)3', 24, '23'),
 207: ('P432', 'P432', 24, '432'),
 208: ('P4232', 'P4(2)32', 24, '432'),
 209: ('F432', 'F432', 96, '432'),
 210: ('F4132', 'F4(1)32', 96, '432'),
 211: ('I432', 'I432', 48, '432'),
 212: ('P4332', 'P4(3)32', 24, '432'),
 213: ('P4132', 'P4(1)32', 24, '432'),
 214: ('I4132', 'I4(1)32', 48, '432')}

X, Y, Z = ex, ey, ez
PGequiv={'1':[],'2':[[-X,Y,-Z],],'222':[[-X,-Y,Z],
[X,-Y,-Z],[-X,Y,-Z]],'4':[[-Y,X,Z],[-X,-Y,Z],[Y,-X,Z]],
'422':[[-Y,X,Z],[-X,-Y,Z],[Y,-X,Z],[X,-Y,-Z],[-X,Y,-Z],
[Y,X,-Z],[-Y,-X,-Z]],'3':[[-Y,X-Y,Z],[-X+Y,-X,Z]],
'312':[[-Y,X-Y,Z],[-X+Y,-X,Z],[-Y,-X,-Z],[-X+Y,Y,-Z],[X,X-Y,-Z]],
'321':[[-Y,X-Y,Z],[-X+Y,-X,Z],[X-Y,-Y,-Z],[-X,-X+Y,-Z],[Y,X,-Z]],
'6':[[X-Y,X,Z],[-Y,X-Y,Z],[-X,-Y,Z],[-X+Y,-X,Z],[Y,-X+Y,Z]],
'622':[[X-Y,X,Z],[-Y,X-Y,Z],[-X,-Y,Z],[-X+Y,-X,Z],[Y,-X+Y,Z],
[X-Y,-Y,-Z],[-X,-X+Y,-Z],[Y,X,-Z],[-Y,-X,-Z],[-X+Y,Y,-Z],[X,X-Y,-Z]],
'23':[[Z,X,Y],[Y,Z,X],[-Y,-Z,X],[Z,-X,-Y],[-Y,Z,-X],[-Z,-X,Y],
[-Z,X,-Y],[Y,-Z,-X],[-X,-Y,Z],[X,-Y,-Z],[-X,Y,-Z]],
'432':[[-Y,X,Z],[-X,-Y,Z],[Y,-X,Z],[X,-Z,Y],[X,-Y,-Z],[X,Z,-Y],
[Z,Y,-X],[-X,Y,-Z],[-Z,Y,X],[Z,X,Y],[Y,Z,X],[-Y,-Z,X],[Z,-X,-Y],
[-Y,Z,-X],[-Z,-X,Y],[-Z,X,-Y],[Y,-Z,-X],[Y,X,-Z],[-Y,-X,-Z],
[-X,Z,Y],[-X,-Z,-Y],[Z,-Y,X],[-Z,-Y,-X]]}

def getPGequiv(symbol, add_identity=False, invert=False):
    OPs = []
    if add_identity:
        OPs = [[X, Y, Z],]
    OPs += PGequiv[symbol]
    if invert:
       invOPs = []
       for OP in OPs:
           invOPs.append([-1*OP[0],-1*OP[1],-1*OP[2]])
       OPs = invOPs
    return OPs
        

if __name__ == '__main__':
    _test()
