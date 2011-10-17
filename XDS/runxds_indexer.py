#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Little utility that help automate tests and optimization with XDS.
"""

__version__ = "0.1.0"
__date__ = "11-10-2011"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__copyright__ = "Copyright (c) 2011 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import shutil
from XDS import XDSLogParser
from xupy import saveLastVersion, xdsInp2Param, \
                 getProfilRefPar, run_xds, LP_names

gnuplot_template = """set dgrid3d %d,%d
set pm3d
splot 'indexer.log' u 1:2:3 with lines
"""

OUT_FILE = open("index.log","a")
def log(txt):
    OUT_FILE.write(txt)
    sys.stdout.write(txt)

MIN_SPOT_SIZE_LIST_1 = range(2, 36,2)
STRONG_PIXEL_LIST_1 = [2,3,4,5,6,7,8,9,10,11,12] #range(2,12,1)

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
    #if "-a" in sys.argv:
    #    sys.argv.remove('-a')
    #    xp.update(getProfilRefPar())
    #    xp["JOB"] = "DEFPIX", "INTEGRATE", "CORRECT"
    #    xp["REFINE_INTEGRATE"] = "ORIENTATION", "BEAM", "CELL" #"DISTANCE",
    #    shutil.copyfile("GXPARM.XDS","XPARM.XDS")
    if "-i" in sys.argv:
        optid = sys.argv.index("-i")
        _xds_input = sys.argv[optid+1]
        xp.update(xdsInp2Param(inp_str=""))
        sys.argv.remove('-i')
        sys.argv.remove(_xds_input)
    #if "-norun" in sys.argv:
    #    saveLastVersion(LP_names)
    #    sys.exit()
    ARGV = sys.argv[1:]
    while ARGV:
        ARG = ARGV.pop()
        xp.update(xdsInp2Param(inp_str=ARG))

    open("gnuplot.inp").write("gnuplot_template" % (len(STRONG_PIXEL_LIST_1), len(MIN_SPOT_SIZE_LIST_1)))
    xp["JOB"] = "COLSPOT", "IDXREF"
    for spot_size in MIN_SPOT_SIZE_LIST_1:
        for strong_pixel in STRONG_PIXEL_LIST_1:
            xp["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = spot_size
            xp["STRONG_PIXEL"] = strong_pixel
            run_xds(xp, inp_f="XDS.INP", out_f="xds_indexer.log")
            saveLastVersion(LP_names)
            res = XDSLogParser("IDXREF.LP", verbose=False).results
            log( "%4d %5.1f" % (spot_size, strong_pixel))
            log( "%(indexed_percentage)6.1f %(indexed_spots)9d %(total_spots)9d" % res)
            log( "%(xy_spot_position_ESD)6.2f %(z_spot_position_ESD)6.2f" % res)
            log( "  %(refined_cell_str)s\n" % res)

