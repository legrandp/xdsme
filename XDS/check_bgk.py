#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Little utility to analyse background stats of a datacollection with XDS.
"""

__version__ = "0.0.1"
__date__ = "01-02-2018"
__author__ = "Pierre Legrand (pierre.legrand _at_ synchrotron-soleil.fr)"
__copyright__ = "Copyright (c) 2018 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import os
import XIO

from xupy import saveLastVersion, xdsInp2Param, \
                 getProfilRefPar, run_xds, LP_names, opWriteCl


DIRNAME_PREFIX = "backgroud_xds_"

def mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

def make_xds_input(image_name):
    try:
        image_name_real = os.path.realpath(image_name)
        datacoll = XIO.Collect(image_name_real)
        datacoll.interpretImage()
        datacoll.lookup_imageRanges()
    except XIO.XIOError:
        print "\nError while trying to acceess %s.\nSorry." % filename
        sys.exit()
    newPar = datacoll.export("xds")
    newPar.pop("SPECIFIC_KEYWORDS")
    newPar["JOB"] = "XYCORR", "INIT"
    newPar["BACKGROUND_RANGE"] =  newPar["DATA_RANGE"]
    TMPL = datacoll.get_export_template("xds")
    newDir = DIRNAME_PREFIX + datacoll.prefix
    return TMPL, newPar, newDir

if __name__ == '__main__':
        tmpl, par, new_dir = make_xds_input(sys.argv[1])
        mkdir(new_dir)
        os.chdir(new_dir)
        opWriteCl("XDS.INP", tmpl[:-25] % par)
        run_xds(par, inp_f="XDS.INP")
        saveLastVersion(LP_names)
