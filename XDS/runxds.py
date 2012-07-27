#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Little utility that help automate tests and optimization with XDS.
"""

__version__ = "0.4.5"
__date__ = "12-06-2012"
__author__ = "Pierre Legrand (pierre.legrand _at_ synchrotron-soleil.fr)"
__copyright__ = "Copyright (c) 2004-2012 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import shutil

from xupy import saveLastVersion, xdsInp2Param, \
                 getProfilRefPar, run_xds, LP_names


if __name__ == '__main__':

    xp = {}
    #xp = xdsInp2Param()
    #xp = XParam()
    #xp.SPOT_RANGE = [2, 12],[44, 54]
    #print xp.SPOT_RANGE
    #xp.SPOT_RANGE = "2 12", "44 54"
    #xp.DATA_RANGE = 8, 8
    #xp.INCLUDE_RESOLUTION_RANGE= 70.0, 2.0
    #xp.NBX = 3
    #xp.NBY = 3
    if "-a" in sys.argv:
        sys.argv.remove('-a')
        xp.update(getProfilRefPar())
        xp["JOB"] = "DEFPIX", "INTEGRATE", "CORRECT"
        xp["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA_BETA"] = 15
        xp["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA"] = 15
        xp["REFINE_INTEGRATE"] = "ORIENTATION", "BEAM", "CELL" #"DISTANCE",
        shutil.copyfile("GXPARM.XDS","XPARM.XDS")
    if "-p" in sys.argv:
        sys.argv.remove('-p')
        xp.update(getProfilRefPar())
        xp["JOB"] = "DEFPIX", "INTEGRATE", "CORRECT"
        xp["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA_BETA"] = 15
        xp["NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA"] = 15
        xp["REFINE_INTEGRATE"] = "ORIENTATION", "BEAM", "CELL" #"DISTANCE",
    if "-i" in sys.argv:
        optid = sys.argv.index("-i")
        _xds_input = sys.argv[optid+1]
        xp.update(xdsInp2Param(inp_str=_xds_input))
        sys.argv.remove('-i')
        sys.argv.remove(_xds_input)
    if "-norun" in sys.argv:
        saveLastVersion(LP_names)
        sys.exit()
    else:
        ARGV = sys.argv[1:]
        while ARGV:
            ARG = ARGV.pop()
            xp.update(xdsInp2Param(inp_str=ARG))
        run_xds(xp, inp_f="XDS.INP")
        saveLastVersion(LP_names)
