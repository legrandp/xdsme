# -*- coding: utf-8 -*-

"""
This is a Numeric/numpy free port of the method:
Scientific.Geometry.Transformation.Rotation.axisAndAngle(self)
From Konrad Hinsen ScientificPython
http://dirac.cnrs-orleans.fr/plone/software/scientificpython
"""

__author__ = "Pierre Legrand (pierre legrand \at synchrotron-soleil fr)"
__date__ = "23-11-2009"
__copyright__ = "Copyright (c) 2009  Pierre Legrand"
__version__ = "0.1.0"

from pycgtypes import  vec3
from pycgtypes import  mat3
import math

R2D = 180/math.pi

def asymmetrical_part(mat_3):
    "Return the asymmetrical part."
    if len(mat_3.mlist) == 9 and len(mat_3) == 3:
        return 0.5*(mat_3 - mat_3.transpose())
    else:
        raise ValueError('Not yet implemented')   

def dyadic_product(vector1, vector2):
    "Dyadic product of two vectors."
    matr1, matr2 = mat3(), mat3()
    matr1.setColumn(0, vector1)
    matr2.setRow(0, vector2)
    return matr1*matr2

def trace(ma3):
    "Return the trace of the matrix."
    return ma3[0, 0]+ma3[1, 1]+ma3[2, 2]

def angle_from_sine_and_cosine(ylen, xlen):
    "Indirection to atan2 with check."
    try:
        return math.atan2(ylen, xlen)
    except:
        raise TypeError, 'AxisAndAngle:FAILURE in atan2'

def axis_and_angle(mat_3):
    """From a rotation matrix return a corresponding rotation as an
       axis (a normalized vector) and angle (in radians).
       The angle is in the interval (-pi, pi]
    """
    asym = -asymmetrical_part(mat_3)
    axis = vec3(asym[1, 2], asym[2, 0], asym[0, 1])
    sine = axis.length()
    if abs(sine) > 1.e-10:
        axis = axis/sine
        projector = dyadic_product(axis, axis)
        cosine = trace((mat_3-projector))/(3.-axis*axis)
        angle = angle_from_sine_and_cosine(sine, cosine)
    else:
        tsr = 0.5*(mat_3+mat3(1))
        diag = tsr[0, 0], tsr[1, 1], tsr[3, 3] 
        i = tsr.index(max(diag))
        axis = vec3(tsr.getRow(i)/(tsr[i, i])**0.5)
        angle = 0.
        if trace(tsr) < 2.:
            angle = math.pi
    return axis, angle

# Test code
if __name__ == '__main__':

    from Scientific.Geometry import Vector ##.Transformation import *
    from Scientific.Geometry.Transformation import Rotation
    from random import random
    #
    Q = mat3(0.36, 0.48, -0.8, -0.8, 0.6, 0, 0.48, 0.64, 0.60)
    axis_q, angle_q = axis_and_angle(Q)
    print "Axis_q:  %9.6f%9.6f%9.6f" % tuple(axis_q),
    print "Angle_q: %10.5f" % (angle_q*R2D)
    #
    for iii in range(1e6):
        axis_i = list(vec3([random(), random(), random()]).normalize())
        angle_i = 3*random()
        rme = mat3().rotation(angle_i, vec3(axis_i))
        axis_1, angle_1 = axis_and_angle(rme)

        v = Vector(axis_i)
        r = Rotation(v, angle_i)
        axis_2, angle_2 = r.axisAndAngle()
        axis_d = (axis_1 - vec3(tuple(axis_2))).length()
        angle_d = abs(angle_1 - angle_2)
        if (angle_d  > 1e-13) or (axis_d > 1e-13):
            print "Angle_d:  %.3e" % (angle_d*R2D),
            print "  Axis_length_diff:  %.3e" % axis_d 
            print "Axis_i:  %9.6f%9.6f%9.6f" % tuple(axis_i),
            print "Angle_i: %10.5f" % (angle_i*R2D)
            print "Axis_1:  %9.6f%9.6f%9.6f" % tuple(axis_1),
            print "Angle_1: %10.5f" % (angle_1*R2D)
            print "Axis_2:  %9.6f%9.6f%9.6f" % tuple(axis_2),
            print "Angle_2: %10.5f" % (angle_2*R2D)

