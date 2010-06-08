#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Little utility that help automate tests and optimization with XDS.
"""

__version__ = "0.0.1"
__date__ = "06-06-2010"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__copyright__ = "Copyright (c) 2010 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import os

from xupy import  xdsInp2Param

if __name__ == '__main__':
    # run_xds(new_par, inp_f=xinp, out_f=None, directory=None, save=1):
    inp_f1, inp_f2 = sys.argv[1:3]
    if inp_f1:
        # Verify the presence of the specified XDS.INP file
        if not os.path.isfile(inp_f1):
            print ">>> ERROR: Can't find file "+inp_f1+" !"
            sys.exit()
        xpar1 = xdsInp2Param(inp_name=inp_f1)
    if inp_f2:
        # Verify the presence of the specified XDS.INP file
        if not os.path.isfile(inp_f2):
            print ">>> ERROR: Can't find file "+inp_f2+" !"
            sys.exit()
        xpar2 = xdsInp2Param(inp_name=inp_f2)
    for k in xpar1.keys():
        if k in xpar2.keys():
            if xpar1[k] != xpar2[k]:
                print "%s:   %s   %s" % (k, xpar1[k], xpar2[k])
    #print xpar1.keys()
