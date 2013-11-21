#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Initial version 02/04/03 legrand@embl-grenoble.fr
"""

__version__ = "0.3.1"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "26-03-2013"
__copyright__ = "Copyright (c) 2006-2013  Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import os
import fnmatch
import xupy

facteur_repr = [['LowR','LowestReso',"%5.1f"],
                ['HighR','HighestReso',"%5.2f"],
                ['Cmpl','compl',"%5.1f"],
                ['CmplL','complL',"%5.1f"],
                ['Unique','unique',"%7d"],
                ['Total','total',"%7d"],
                ['Compar','compar',"%7d"],
                ['Rsym','rsym',"%6.2f"],
                ['RsymL','rsymL',"%6.1f"],
                ['Rmeas','rmeas',"%6.1f"],
                ['Isig','isig',"%5.1f"],
                ['IsigL', 'isigL',"%5.1f"],
                ['Misfit','misfit',"%6d"],
                ['CC1/2L','cchalfL',"%6.1f"],
                ['AnoC','anoCorr',"%4.f"],
                ['AnoS','anoSig',"%5.2f"],
                ['Abs','absent',"%4d"],
                ['IAbs','AbsIav',"%5.1f"],
                ['ISa','IoverSigmaAsympt','%6.2f']]

#                ['Total3','total3',"%6d"],
#                ['Cmpl3','compl3',"%5.1f"],
#                ['Compa3','compar3',"%6d"],
#                ['Rsym3','rsym3',"%6.1f"],
#                ['Rmea3','rmeas3',"%6.1f"],

def list_lp(pattern='*'):
    result = []
    names = os.listdir(".")
    for name in names:
            if fnmatch.fnmatch(name, pattern) and os.path.isfile(name):
                    result.append(name)
    result.sort()
    return result

if __name__ == '__main__':

    digit_template = ".???"
    if "-old" in sys.argv:
        digit_template = ".*"
        sys.argv.remove("-old")

    facteur_names = "    "
    for a in facteur_repr:
        fmt = a[2].split(".")[0].replace("d","")+"s"
        facteur_names += ((fmt % a[0]) + " ")
    print facteur_names
    if len(sys.argv) >= 2:
        file_list = sys.argv[1:]
    else:
        file_list = list_lp("CORRECT.LP"+digit_template)
        file_list+= list_lp("XSCALE.LP"+digit_template)
    #print "   "+len(facteur_names)*"%7s" % tuple(facteur_names)
    for lp in file_list:
        stat = xupy.resum_scaling(lpf=lp)
        print "%1s%3s" % (lp[0],lp.split(".")[-1]),
        if not stat: print "Error in CORRECT."; continue
        for k in facteur_repr:
            print k[2] % stat[k[1]],
        print
