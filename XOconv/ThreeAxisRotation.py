# This module defines a class that represent coordinate rotation by
# a three axis rotation system. Handy for goniosta systems, Euler angles or Kappa
# angles manipulation.
#
# Need the ScientificPython module by Konrad Hinsen
# http://starship.python.net/crew/hinsen/scientific.html
# NOTE: The getAngles method has been integrated to ScientificPython since the
#       2.4.9 version (in the Rotation class with method name ...).

# TODO:
#  - Add the posibility to input/output angles in radian (defaults) or degree.
#  - Add doctest.

__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "16-05-2005"
__copyright__ = "Copyright (c) 2004-2005  Pierre Legrand"
__license__ = "GPL"

print "New"

import Numeric
from Scientific.Geometry.VectorModule import *
from Scientific.Geometry.TensorModule import *
from Scientific.Geometry.Transformation import *
import math

pi = math.pi
r2d = 180/pi

# Utility functions
def mod_angle(angle, mod):
    return (angle + mod/2.) % mod - mod/2

def kappaVector(alpha):
    return Vector([-1*math.sin(alpha), 0, math.cos(alpha)])

def angleFromSineAndCosine2(y, x):
    try: return math.atan2(y, x)
    except: raise TypeError, 'ThreeAxisRotation:FAILURE in atan2'

# Some Axes Definitions
StandardAxes = ex, ey, ez
NoniusKappaAxes = ez, kappaVector(50/r2d), ez
EulerAxes = ez, ex, ez
DenzoAxes = ez, ey, ex

class ThreeAxisRotation(Rotation):
    
    """Rotational transformation using three successive axes.
    This is a subclass of Rotation (which is a subclass of Transformation).

    Constructor:

      Rotation(init, rotationAxes=(ex, ey, ez), inversAxesOrder=0), where:
     
      - init can be either:    
         - a tensor object containing the rotation matrix,
         - a Rotation object,
         - a set of three angles (a1, a2, a3) in radian
         
      - rotationAxes is a tupple of 3 vectors (e1, e2, e3) defining the three
        rotation axes. Default axes are (ex, ey, ez).
      
      The resulting rotation is obtained by applying successivly:
        - a1 angle around e1 axis, and then
        - a2 angle around e2 axis, and finaly
        - a3 angle around e3.

      - If  inversAxesOrder is set to 1, the the rotations are done in the
        reverse order (a3 angle around e3, a2 around e2 and a1 around e1).
     
     Rq: ThreeAxisRotation((a1,a2,a3),rotationAxes=(e1,e2,e3),inversAxesOrder=0)
         is equivalent to:
         ThreeAxisRotation((a3,a2,a1),rotationAxes=(e3,e2,e1),inversAxesOrder=1)
      
    """
    
    def __init__(self, init, rotationAxes=(ex, ey, ez), inversAxesOrder=0):
        self.inversAxesOrder = inversAxesOrder
        self.e1 = Vector(rotationAxes[0]).normal()
        self.e2 = Vector(rotationAxes[1]).normal()
        self.e3 = Vector(rotationAxes[2]).normal()
        if inversAxesOrder:
            self.e1 = Vector(rotationAxes[2]).normal()
            self.e3 = Vector(rotationAxes[0]).normal()
        self.rotationAxes = self.e1, self.e2, self.e3
	if hasattr(init,'is_tensor') == 1 :
	    self.tensor = init
	elif hasattr(init,'is_rotation') == 1:
	    self.tensor = init.tensor
	elif len(init) == 3:
            print "len3"
            a1, a2, a3 = init
            print a1, a2, a3
            if inversAxesOrder: a3, a2, a1 = init
            R1 = Rotation(self.e1, a1)
            R2 = Rotation(self.e2, a2)
            R3 = Rotation(self.e3, a3)
	    self.tensor = (R1*R2*R3).tensor
	elif len(init) == 9 and type(init) == list:
            "Used for compatibility with other tensor types, like cgtypes.mat3"
            "In that case, use mat3.mlist."
            self.tensor = Tensor([[init[0], init[1], init[2]],
                                  [init[3], init[4], init[5]],
                                  [init[6], init[7], init[8]]])
	else:
	    raise TypeError, 'no valid arguments'
        
    is_3axes_rotation = 1
    
    def addAngles(self, angles):
        """Add angles to the current rotation.
         Does it works???."""
        a1, a2, a3 = angles
        if self.inversAxesOrder: a3, a2, a1 = angles
        R1 = Rotation(self.e1, a1)
        R2 = Rotation(self.e2, a2)
        R3 = Rotation(self.e3, a3)
        #self.tensor = (R1*R2*R3*self).tensor
        self.tensor = (self*R1*R2*R3).tensor
        
    def setAngles(self, angles):
        "Set angles defining a new rotation"
        a1, a2, a3 = angles
        if self.inversAxesOrder: a3, a2, a1 = angles
        R1 = Rotation(self.e1, a1)
        R2 = Rotation(self.e2, a2)
        R3 = Rotation(self.e3, a3)
        self.tensor = (R1*R2*R3).tensor

    def getAngles(self, tolerance=1e-7):
        """
        Get angles from the matrix given a set of 3 unit vectors describing
        a 3-axis goniometer. Basicly this a reimplementation of the David
        Thomas's algorithm [1] described by Gerard Bricogne in [2]:
        
        [1] "Modern Equations of Diffractometry. Goniometry." D.J. Thomas
        Acta Cryst. (1990) A46 Page 321-343.
        
        [2] "The ECC Cooperative Programming Workshop on Position-Sensitive
        Detector Software." G. Bricogne,
        Computational aspect of Protein Crystal Data Analysis,
        Proceedings of the Daresbury Study Weekend (23-24/01/1987)
        Page 122-126
        """
        
        # We are searching for the three angles a1, a2, a3
        e1, e2, e3 = self.rotationAxes
        # If 2 consecutive axes are parallel: decomposition is not meaningful
        if (e1.cross(e2)).length() < tolerance or \
           (e2.cross(e3)).length() < tolerance :
            raise ValueError, 'Consecutive parallel axes. To many solutions'
        ROT = Rotation(self.tensor)
        w = ROT(e3)
        
        # Solve the equation : _a.cosx + _b.sinx = _c
        _a = e1*e3 - (e1*e2)*(e2*e3)
        _b = e1*(e2.cross(e3))
        _c = e1*w - (e1*e2)*(e2*e3)
        _norm = (_a**2 + _b**2)**0.5
        
        # Checking for possible errors in initial Rot matrix
        if _norm == 0: raise ValueError, 'ThreeAxisRotation:FAILURE 1, norm = 0'
        if abs(_c/_norm) > 1+tolerance:
            raise ValueError, 'ThreeAxisRotation:FAILURE 2' + \
                 'malformed rotation Tensor (non orthogonal?) %.8f' % (_c/_norm)
        _c_norm = _c/_norm
        if _c_norm > 1:
            if _c_norm < 1 + tolerance: _c_norm = 1
            else: raise ValueError, 'Step1: No solution'
        _th = angleFromSineAndCosine2(_b/_norm, _a/_norm)
        _xmth = math.acos((_c_norm))
                
        # a2a and a2b are the two possible solutions to the equation.
        a2a = mod_angle((_th + _xmth), 2*pi)
        a2b = mod_angle((_th - _xmth), 2*pi)
        
        solutions = []
        # for each solution, find the two other angles (a1, a3).
        for a2 in (a2a, a2b):
            R2 = Rotation(e2, a2)
            v =  R2(e3)
            v1 = v - (v*e1)*e1
            w1 = w - (w*e1)*e1
            norm = ((v1*v1)*(w1*w1))**0.5
            if norm == 0: 
                # in that case rotation 1 and 3 are about the same axis
                # so any solution for rotation 1 is OK
                a1 = 0.
            else:
                cosa1 = (v1*w1)/norm
                sina1 = v1*(w1.cross(e1))/norm
                a1 = mod_angle(angleFromSineAndCosine2(sina1, cosa1), 2*pi)
                
            R3 = Rotation(e2, -1*a2)*Rotation(e1, -1*a1)*ROT
            # u = normalized test vector perpendicular to e3
            # if e2 and e3 are // we have an exception before.
            # if we take u = e1^e3 then it will not work for Euler and Kappa axes.
            u = (e2.cross(e3)).normal()
            cosa3 = u*R3(u)
            sina3 = u*(R3(u).cross(e3))
            a3 =  mod_angle(angleFromSineAndCosine2(sina3, cosa3), 2*pi)
            
            if self.inversAxesOrder: solutions.append(Numeric.array([a3, a2, a1]))
            else: solutions.append(Numeric.array([a1, a2, a3]))
            
        # Gives the closest solution to 0,0,0 first
        if Numeric.add.reduce(solutions[0]**2) > Numeric.add.reduce(solutions[1]**2):
            solutions = [solutions[1], solutions[0]]
        return solutions

    def asQuaternion(self):
        from Scientific.Geometry.Quaternion import Quaternion
        axis, angle = self.axisAndAngle()
        sin_angle_2 = math.sin(0.5*angle)
        cos_angle_2 = math.cos(0.5*angle)
        return Quaternion(cos_angle_2, sin_angle_2*axis[0],
                          sin_angle_2*axis[1], sin_angle_2*axis[2])

# Tests functions
def random_3axes():
    r = random.uniform
    return (Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normal(),
            Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normal(),
            Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normal())

def test_internal_coherance(axes_i):
    sum = Numeric.add.reduce
    r = random.uniform
    a1, a2, a3 = r(-180,180),r(-180,180),r(-180,180)
    angles_i = Numeric.array([a1, a2, a3])
    angles_r = angles_i[::-1]
    axes_r = axes_i[2], axes_i[1], axes_i[0]
    
    rotdnz_new = ThreeAxisRotation(angles_i/r2d, rotationAxes=axes_i, inversAxesOrder=1)
    rotdnz_inv = ThreeAxisRotation(angles_r/r2d, rotationAxes=axes_r, inversAxesOrder=0)
    diffmat1 = sum(sum(Numeric.absolute(rotdnz_new.tensor - rotdnz_inv.tensor)))
    
    sol1_new, sol2_new = rotdnz_new.getAngles()
    sol1_inv, sol2_inv = rotdnz_inv.getAngles()
    diffinit = min(sum(Numeric.absolute(sol1_new*r2d-angles_i)),
                   sum(Numeric.absolute(sol2_new*r2d-angles_i)))
    diffinv  = min(sum(Numeric.absolute(sol1_inv[::-1]*r2d-angles_i)),
                   sum(Numeric.absolute(sol2_inv[::-1]*r2d-angles_i)))
    
    #print rotdnz_new.asQuaternion()
    if diffmat1 > 1e-12 or diffinit > 2e-11 or diffinv > 2e-11:
        print "I"+3*"%10.2f" % tuple(Numeric.array(angles_i))
        print "N"+3*"%10.2f" % tuple(Numeric.array(sol1_new)*r2d)
        print "N"+3*"%10.2f" % tuple(Numeric.array(sol2_new)*r2d)
        print 'Matrix difference \t%.1e' % diffmat1
        print 'Init Angle difference \t%.1e' % diffinit
        print 'Inv  Angle difference \t%.1e' % diffinv

    
if __name__ == '__main__':
    import random
    
    for a in range(300):
        axes = random_3axes()
        test_internal_coherance(axes)
    for axes in (StandardAxes, DenzoAxes, NoniusKappaAxes, EulerAxes, random_3axes()):
        print axes
        for a in range(300):
            test_internal_coherance(axes)
