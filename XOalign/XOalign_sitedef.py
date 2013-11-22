from math import cos, sin, pi
from pycgtypes import vec3, mat3

r2d = 180/pi
ex, ey, ez = vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)
def kappaVector(alpha):
    return vec3([-sin(alpha), 0, cos(alpha)])

def kappaVectorVertical(alpha):
    return vec3([0, -cos(alpha), sin(alpha)])

# Frame definition used: Cambridge as used in Mosflm
# X = Beam_vector direction of the X-ray photons
# Z = Spindle axis, such that looking down this axis
#     towards the sample, positive phi is anti-clockwise
# Y = Defined to give a right handed coordinate system


# -- SOLEIL's PX1 CrystalLogic Goniometer definitions --
GONIOMETER_NAME = "SOLEIL PROXIMA-1 CrystalLogic"
GONIOMETER_AXES_NAMES = ("Omega","Kappa","Phi")
GONIOMETER_AXES = [ez, kappaVector(49.64/r2d), ez]
GONIOMETER_DATUM = (0,0,0)  # in degree

# -- DLS's MiniKappa Goniometer definitions --
#GONIOMETER_NAME = "DLS's MiniKappa"
#GONIOMETER_AXES_NAMES = ("Omega","Kappa","Phi")
#GONIOMETER_AXES = [ez, kappaVector(-24/r2d), ez]
#GONIOMETER_DATUM = (0,0,45)  # in degree
#GONIOMETER_AXES = [[0.00211, 0.00143, 1.], [0.28907, 0.28990, 0.91236], [0.00691, -0.00364, 0.99997]]
#GONIOMETER_DATUM = (0,0,0)  # in degree

# -- DLS's MiniKappa Goniometer definitions --
#GONIOMETER_NAME = "Petra's P14 MiniKappa"
#GONIOMETER_AXES_NAMES = ("Omega","Kappa","Phi")
#GONIOMETER_AXES = [ez, kappaVector(-24/r2d), ez]
#GONIOMETER_AXES = [-ey, kappaVectorVertical(24/r2d), -ey]
#GONIOMETER_DATUM = (0,0,27.2)  # in degree

# -- SLS's Prigo Goniometer definitions --
#GONIOMETER_NAME = "SLS's Prigo"
#GONIOMETER_AXES_NAMES = ("Omega","Chi","Phi")
#GONIOMETER_AXES = [ez, ex, ez]
#GONIOMETER_DATUM = (0,0,90)  # in degree

