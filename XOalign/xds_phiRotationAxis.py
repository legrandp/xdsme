#!/usr/bin/env python2

# Script used for the SOLEIL_PROXIMA-1 Kappa, in the case where
# data are collected with the phi axis when Kappa and Omega are non-null.
# The script outputs the ROTATION_AXIS= vector to use with XDS.
# Initial version Nov 2009

__version__ = "2015_07_25"
__author__ = "Pierre Legrand (pierre legrand \at synchrotron-soleil fr)"
__copyright__ = "Copyright (c) 2009-2015 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
from XOconv import *

usage = "%s [omage,kappa]|[omega kappa]\n "

pi = math.pi
r2d = 180/pi
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)

Vector = vec3
def kappaVector(alpha=49.64/r2d):
    return Vector([-1*math.cos(alpha), 0, math.sin(alpha)])

kappa = Rotation2(kappaVector(), 30./r2d)

if __name__ == "__main__":

    progname = sys.argv[0].split("/")[-1]
    if ("-h" in sys.argv) \
        or  ("--help" in sys.argv) or ("-help" in sys.argv):
        print usage % (progname, "|".join(options))
        sys.exit(1)

    if len(sys.argv) == 2:
        #try:
            omega_angle,kappa_angle = map(float,sys.argv[1].split(","))
        #except:
        #    print "Error: Can't parse input angles. Format needed: -20.3,-30."
        #    sys.exit(1)
    elif len(sys.argv) == 3:
        try:
            omega_angle = float(sys.argv[1])
            kappa_angle = float(sys.argv[2])
        except:
            print "Error: Can't parse input angles. Format needed: -20.3 -30."
            sys.exit(1)
    else:
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
    print "  For OMEGA angle= %8.2f deg" %  omega_angle
    print "      KAPPA angle= %8.2f deg" %  kappa_angle
    omega = Rotation2(ex, -omega_angle/r2d)
    kappa = Rotation2(kappaVector(), kappa_angle/r2d)
    phi = ex*kappa*omega
    CrystalLogicKappaAxes = [ex, kappaVector(49.64/r2d), ex]
    Goniometer = ThreeAxisRotation2([omega_angle/r2d,kappa_angle/r2d,0.],
                                CrystalLogicKappaAxes,
                                inversAxesOrder=0)
    #phi = ex*Goniometer.tensor
    print " For image collection using PHI axis with a KAPPA angle != 0 in XDS.INP use:"
    print "\n ROTATION_AXIS= %.5f %.5f %.5f" % (phi[0], phi[1], phi[2])
