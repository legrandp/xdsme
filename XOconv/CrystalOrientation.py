#!/usr/bin/env python2
"""
    19/10/05 First version
    
    Uses the mat3 and vec3 classes from Python Computer Graphics Kit v1.2.0
    module by Matthias Baas (see http://cgkit.sourceforge.net).
    
    Uses the ScientificPython module by Konrad Hinsen
    http://starship.python.net/crew/hinsen/scientific.html
    
    TODO:
      - Convert ThreeAxisOritentation to use mat3 & vec3...
      - 
"""

import sys
import os
import math

#try:
#    from cgtypes import vec3
#    from cgtypes import mat3
#except ImportError:
#    try:
#        from pycgtypes import vec3
#        from pycgtypes import mat3
#    except ImportError:
#        print "Can't find cgtypes or pycgtype modules... Please check your"
#        print "$PYTHONPATH variable if it is installed somewhere or instal" 
#        print "one of these module for vector and matrix operations."
#        sys.exit()

from pycgtypes import vec3, mat3

class mat3(mat3):
    
    def dot(self, other):
        return self * other
    
    def __str__(self):
        fmt="%9.5f"
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mlist
        return ('['+fmt%m11+', '+fmt%m12+', '+fmt%m13+']\n'+
                '['+fmt%m21+', '+fmt%m22+', '+fmt%m23+']\n'+
                '['+fmt%m31+', '+fmt%m32+', '+fmt%m33+']')
class vec3(vec3):
    
    def normal(self):
        return self.normalize()
        
# Vector = vec3
# Tensor = mat3    
r2d = 180/math.pi
cosd = lambda a: math.cos(a/r2d)
sind = lambda a: math.sin(a/r2d)
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)


def diffMAT(m1,m2):
    diffsq = 0
    for a in (m1-m2):
        diffsq += a[0]**2+a[1]**2+a[2]**2
    return diffsq**0.5

def volum(cell):
    """
    
    """
    if len(cell) == 6 and type(cell[0]) == float:
        # expect a, b, c, alpha, beta, gamma (angles in degree).
        ca, cb, cg = map(cosd, cell[3:6])
        return cell[0]*cell[1]*cell[2]*(1-ca**2-cb**2-cg**2+2*ca*cb*cg)**0.5
         
    elif len(cell) == 3 and isinstance(cell[0], vec3):
        # expect vectors of the 3 cell parameters
        a, b, c = cell
        return a*b.cross(c)
    else:
        print "error in volum()"
        return ""  

def reciprocal(cell):
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
    ... 
    """
    Ar = vec3(UB.getColumn(0))
    Br = vec3(UB.getColumn(1))
    Cr = vec3(UB.getColumn(2))
    return (Ar.length(), Br.length(), Cr.length(),
            Br.angle(Cr)*r2d, Cr.angle(Ar)*r2d, Ar.angle(Br)*r2d)


def BusingLevy(rcell):
    cosr = map(cosd, rcell[3:6])
    sinr = map(sind, rcell[3:6])
    Vr = volum(rcell)
    X = ex*rcell[0]
    Y = rcell[1]*(ex*cosr[2] + ey*sinr[2])
    c = rcell[0]*rcell[1]*sinr[2]/Vr
    cosAlpha = (cosr[1]*cosr[2] - cosr[0])/(sinr[1]*sinr[2])
    Z = vec3([rcell[2]*cosr[1],
                -1*rcell[2]*sinr[1]*cosAlpha,
                1/c])
    return mat3(X,Y,Z)

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

if __name__ == '__main__':
    _test()
