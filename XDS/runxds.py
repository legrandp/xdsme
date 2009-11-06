#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys, os
from xupy import *

__version__ = "0.4.3"
__date__ = "21-10-2004"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__copyright__ = "Copyright (c) 2004 Pierre Legrand"
__license__ = "LGPL"
            
if __name__ == '__main__':

    xp = XParam()
    #xp = xdsInp2Param()
    #xp.SPOT_RANGE = [2, 12],[44, 54]
    #print xp.SPOT_RANGE
    #xp.SPOT_RANGE = "2 12", "44 54"
    #xp.DATA_RANGE = 8, 8
    #xp.INCLUDE_RESOLUTION_RANGE= 70.0, 2.0
    #xp.NBX = 3
    #xp.NBY = 3
    if "-a" in sys.argv:
        sys.argv.remove('-a')
        xp.mix(getProfilRefPar())
        xp.JOB = "DEFPIX", "INTEGRATE", "CORRECT"
        xp.REFINE_INTEGRATE= "ORIENTATION", "BEAM",  "CELL" #"DISTANCE",
        shutil.copyfile("GXPARM.XDS","XPARM.XDS")
    if "-i" in sys.argv:
        sys.argv.remove('-i')
        _xds_input = sys.argv[-1]
        xp.update(xdsInp2Param(inp_str=_xds_input))
    if "-norun" in sys.argv:
        saveLastVersion(LP_names)
        sys.exit()
#    if "-adp" in sys.argv:
#        
    else:
        par_str = ""
        _argv = sys.argv[1:]
        while _argv:
            par_str += "%s\n" % _argv.pop()
        xp.mix(par_str)
        run_xds(xp, inp_f="XDS.INP")
        saveLastVersion(LP_names)
