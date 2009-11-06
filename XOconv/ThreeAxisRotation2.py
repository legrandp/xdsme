# This module defines a class that represent coordinate rotation by
# a three axis rotation system. Handy for goniosta systems, Euler angles or Kappa
# angles manipulation.
# 
# TODO:
#  - Add the posibility to input/output angles in radian (defaults) or degree.
#  - Add doctest.

__version__ = "0.2.5"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "27-09-2006"
__copyright__ = "Copyright (c) 2004-2006  Pierre Legrand"
__license__ = "GPL"

from pycgtypes import  vec3
from pycgtypes import  mat3
import math

pi = math.pi
r2d = 180/pi
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)

Vector = vec3

# Utility functions
square = lambda a: a**2
addReduceSq = lambda l: reduce(lambda a, b: a + b, map(square, l))

def Rotation2(axis, angle):
    return mat3().rotation(angle, axis)

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

class ThreeAxisRotation2:
    
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
         (a1, a2, a3 in radian).

      - If  inversAxesOrder is set to 1, the the rotations are done in the
        reverse order (a3 angle around e3, a2 around e2 and a1 around e1).
     
     Rq: ThreeAxisRotation((a1,a2,a3),rotationAxes=(e1,e2,e3),inversAxesOrder=0)
         is equivalent to:
         ThreeAxisRotation((a3,a2,a1),rotationAxes=(e3,e2,e1),inversAxesOrder=1)
      
    """
    
    def __init__(self, init, rotationAxes=(ex, ey, ez), inversAxesOrder=0):
        self.inversAxesOrder = inversAxesOrder
        self.e1 = Vector(rotationAxes[0]).normalize()
        self.e2 = Vector(rotationAxes[1]).normalize()
        self.e3 = Vector(rotationAxes[2]).normalize()
        if inversAxesOrder:
            self.e1 = Vector(rotationAxes[2]).normalize()
            self.e3 = Vector(rotationAxes[0]).normalize()
        self.rotationAxes = self.e1, self.e2, self.e3
	if hasattr(init,'is_tensor') == 1:
	    self.tensor = init
	elif hasattr(init,'is_rotation') == 1:
	    self.tensor = init.tensor
	elif len(init) == 3:
            a1, a2, a3 = init
            if inversAxesOrder: a3, a2, a1 = init
            R1 = Rotation2(self.e1, a1)
            R2 = Rotation2(self.e2, a2)
            R3 = Rotation2(self.e3, a3)
	    self.tensor = R1*R2*R3
	elif len(init) == 9 and \
            (type(init) == type([]) or type(init) == type(())):
            "Used for compatibility with other tensor types, like cgtypes.mat3"
            "In that case, use mat3.mlist."
            self.tensor = mat3(list(init))
	else:
	    raise TypeError, 'no valid arguments'
        
    is_3axes_rotation = 1
    
    def addAngles(self, angles):
        """Add angles to the current rotation.
         Does it works???."""
        a1, a2, a3 = angles
        if self.inversAxesOrder: a3, a2, a1 = angles
        R1 = Rotation2(self.e1, a1)
        R2 = Rotation2(self.e2, a2)
        R3 = Rotation2(self.e3, a3)
        #self.tensor = (R1*R2*R3*self).tensor
        self.tensor = self.tensor*R1*R2*R3
        
    def setAngles(self, angles):
        "Set angles defining a new rotation"
        a3, a2, a1 = angles
        if self.inversAxesOrder: a1, a2, a3 = angles
        R1 = Rotation2(self.e1, a1)
        R2 = Rotation2(self.e2, a2)
        R3 = Rotation2(self.e3, a3)
        self.tensor = R1*R2*R3

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
        w = self.tensor * e3
        
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
            v =  Rotation2(e2, a2) * e3
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
                
            R3 = Rotation2(e2, -1*a2)*Rotation2(e1, -1*a1)*self.tensor
            # u = normalized test vector perpendicular to e3
            # if e2 and e3 are // we have an exception before.
            # if we take u = e1^e3 then it will not work for Euler and Kappa axes.
            u = (e2.cross(e3)).normalize()
            cosa3 = u * R3 * u
            sina3 = u * ((R3*u).cross(e3))
            a3 =  mod_angle(angleFromSineAndCosine2(sina3, cosa3), 2*pi)
            
            if self.inversAxesOrder: solutions.append([a3, a2, a1])
            else: solutions.append([a1, a2, a3])
            
        # Gives the closest solution to 0,0,0 first
        if addReduceSq(solutions[0]) > addReduceSq(solutions[1]):
            solutions = [solutions[1], solutions[0]]
        return solutions

#    def asQuaternion(self):
#        from Scientific.Geometry.Quaternion import Quaternion
#        axis, angle = self.axisAndAngle()
#        sin_angle_2 = math.sin(0.5*angle)
#        cos_angle_2 = math.cos(0.5*angle)
#        return Quaternion(cos_angle_2, sin_angle_2*axis[0],
#                          sin_angle_2*axis[1], sin_angle_2*axis[2])

# Tests functions
def random_3axes():
    import random
    r = random.uniform
    return (Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normalize(),
            Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normalize(),
            Vector(r(-1.,1.),r(-1.,1.),r(-1.,1.)).normalize())


def test_1(axes_i):
    import random
    
    r = random.uniform
    a1, a2, a3 = r(-pi,pi),r(-pi,pi),r(-pi,pi)
    angles_i = [a1, a2, a3]
    angles_r = [a3, a2, a1]
    axes_r = axes_i[2], axes_i[1], axes_i[0]
    rotdnz_new = ThreeAxisRotation2(angles_i, rotationAxes=axes_i, inversAxesOrder=0)
    rotdnz_inv = ThreeAxisRotation2(angles_r, rotationAxes=axes_r, inversAxesOrder=1)
    diffmat1 = addReduceSq((rotdnz_new.tensor - rotdnz_inv.tensor).mlist)

    sol1_new, sol2_new = rotdnz_new.getAngles()
    sol1_inv, sol2_inv = rotdnz_inv.getAngles()
    sol1_inv2, sol2_inv2 = sol1_inv[:], sol2_inv[:]
    sol1_inv2.reverse()
    sol2_inv2.reverse()
    
    diffinit = min(addReduceSq(map(diff, zip(sol1_new,angles_i))),
                   addReduceSq(map(diff, zip(sol2_new,angles_i))))
    diffinv  = min(addReduceSq(map(diff, zip(sol1_inv2,angles_i))),
                   addReduceSq(map(diff, zip(sol2_inv2,angles_i))))
    
    #if 1:
    if diffmat1 > 2e-24 or diffinit > 2e-24 or diffinv > 2e-24:
        print "0"+3*"%10.2f" % tuple(angles_i)
        print "N"+3*"%10.2f" % tuple(sol1_new)
        print "N"+3*"%10.2f" % tuple(sol2_new)
        print "I"+3*"%10.2f" % tuple(sol1_inv2)
        print "I"+3*"%10.2f" % tuple(sol2_inv2)
        print 'Matrix difference \t%.1e' % diffmat1
        print 'Init Angle difference \t%.1e' % diffinit
        print 'Inv  Angle difference \t%.1e' % diffinv

def test_2(axes_i):
    from ThreeAxisRotation import ThreeAxisRotation
    import Numeric, random
    
    r = random.uniform
    a1, a2, a3 = r(-pi,pi),r(-pi,pi),r(-pi,pi)
    angles_i = [a1, a2, a3]
    rot1 = ThreeAxisRotation (angles_i, rotationAxes=axes_i, inversAxesOrder=1)
    rot2 = ThreeAxisRotation2(angles_i, rotationAxes=axes_i, inversAxesOrder=1)
    rot1t = mat3(list(Numeric.ravel(rot1.tensor.array)))
    diffmat1 = math.sqrt(addReduceSq((rot1t - rot2.tensor).mlist))

    sol1A, sol1B = rot1.getAngles()
    sol2A, sol2B = rot2.getAngles()
    
    diff1I = math.sqrt(min(addReduceSq(map(diff, zip(sol1A,angles_i))), 
                           addReduceSq(map(diff, zip(sol1B,angles_i)))))
    diff2I = math.sqrt(min(addReduceSq(map(diff, zip(sol2A,angles_i))), 
                           addReduceSq(map(diff, zip(sol2B,angles_i)))))
    diffAA = math.sqrt(addReduceSq(map(diff, zip(sol2A,sol1A))))
    diffBB = math.sqrt(addReduceSq(map(diff, zip(sol2B,sol1B))))
    
    #if 1:
    ae = 1e-12
    if diffmat1 > 1e-12 or diffAA > ae or diffBB > ae or diff1I > ae or diff2I > ae:
        print "0 "+3*"%10.2f" % tuple(angles_i)
        print "1A"+3*"%10.2f" % tuple(sol1A)
        print "2A"+3*"%10.2f" % tuple(sol2A)
        print "1B"+3*"%10.2f" % tuple(sol1B)
        print "2B"+3*"%10.2f" % tuple(sol2B)
        print 'Matrix difference \t%.1e' % diffmat1
        print '1I Angle difference \t%.1e' % diff1I
        print '2I Angle difference \t%.1e' % diff1I
        print 'AA Angle difference \t%.1e' % diffAA
        print 'BB Angle difference \t%.1e' % diffBB

if __name__ == '__main__':

    import random
    diff = lambda p: p[0]-p[1]
    
    ntest = 1000
    for a in range(ntest):
        axes = random_3axes()
        test_1(axes)

    for axes in (StandardAxes, DenzoAxes, NoniusKappaAxes, EulerAxes, random_3axes()):
        print axes
        for a in range(ntest):
            test_1(axes)
