#!/usr/bin/env python

import sys
from XOconv import *

pi = math.pi
r2d = 180/pi
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)

Vector = vec3
def kappaVector(alpha=49.64/r2d):
    return Vector([-1*math.cos(alpha), 0, math.sin(alpha)])

kappa = Rotation2(kappaVector(), 30./r2d)

if __name__ == "__main__":

    omega_angle = raw_input(" OMEGA angle at start of the collect in deg [0.0]? ")
    kappa_angle = raw_input(" KAPPA angle at start of the collect in deg [0.0]? ")

    if not omega_angle:
        omega_angle = 0.
    else:
        omega_angle = float(omega_angle)
    if not kappa_angle:
        kappa_angle = 0.
    else:
        kappa_angle = float(kappa_angle)

    print "\n ***  PHI Rotation Vector for XDS ***"
    print "      With OMEGA angle= %8.2f deg" %  omega_angle
    print "      With KAPPA angle= %8.2f deg" %  kappa_angle
    omega = Rotation2(ex, omega_angle/r2d)
    kappa = Rotation2(kappaVector(), kappa_angle/r2d)
    phi = ex*kappa*omega
    CrystalLogicKappaAxes = [ex, kappaVector(49.64/r2d), ex]
    Goniometer = ThreeAxisRotation2([omega_angle/r2d,kappa_angle/r2d,0.],
                                CrystalLogicKappaAxes,
				inversAxesOrder=0)
    #phi = ex*Goniometer.tensor
    print "      For image collection using PHI axis with a KAPPA angle != 0,"
    print "      in XDS use:"
    print "\n      ROTATION_AXIS= %.5f %.5f %.5f\n" % (phi[0], phi[1], phi[2])
