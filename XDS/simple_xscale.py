#!/usr/bin/env python2
# -*- coding: utf-8 -*-

__version__ = "0.0.2"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "01-10-2018"
__copyright__ = "Copyright (c) 2018 Pierre Legrand"
__license__ = "LGPL"

import os
from subprocess import Popen, PIPE
import distutils.spawn

XSCALE_PROGNAME = "xscale_par"
if not distutils.spawn.find_executable(XSCALE_PROGNAME):
    print("Executable xscale_par not found in PATH.")

XSCALE_INP_TEMPL = """OUTPUT_FILE= XSCALE.HKL
STRICT_ABSORPTION_CORRECTION= TRUE
  MERGE=FALSE
  INPUT_FILE= %s
    CORRECTIONS= NONE
"""

def run_exec_str(execstr):
    outp = Popen([execstr], stdout=PIPE, shell=True).communicate()[0]
    return outp

def simple_xscale(inpfn="XDS_ASCII.HKL", run_dir=None):
    "Create a simple XSCALE.INP file and run xscale."
    init_dir = os.getcwd()
    if run_dir:
        os.chdir(run_dir)
    if not distutils.spawn.find_executable(XSCALE_PROGNAME):
        print("Executable xscale_par not found in PATH.")
    xinp = open("XSCALE.INP","w")
    xinp.write(XSCALE_INP_TEMPL % inpfn)
    xinp.close()
    lp = run_exec_str(XSCALE_PROGNAME)
    os.chdir(init_dir)
    return lp

if __name__ == '__main__':
    import sys
    import os.path
    lpstr = ""
    if len(sys.argv) >= 2:
        inpf = sys.argv[1]
        if os.path.isfile(inpf):
            lpstr = simple_xscale(inpf)
        else:
            print "Error. Can't access %s file." % (inpf)
    else:
        if os.path.isfile("XDS_ASCII.HKL"):
            lpstr = simple_xscale()
        else:
            print "Error. Can't access XDS_ASCII.HKL file."
    if lpstr:
        lp = open("xscale.log", "w")
        lp.write(lpstr)

