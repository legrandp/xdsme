#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1.2"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "06-06-2010"
__copyright__ = "Copyright (c) 2010 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import os
from xml.dom import minidom

def is_pointless_installed():
    tmpfile = "/tmp/xdsme_pointless.xml"
    cmline = "rm -f %s; pointless XMLOUT %s > /dev/null 2>&1"
    os.system(cmline % (tmpfile, tmpfile))
    if os.path.isfile(tmpfile):
        return True
    else:
        return False
    os.system("rm -f %s" % tmpfile)

def process_pointless_xml():
    xml_inp = "XDS_pointless.xml"
    logf_name = "XDS_pointless.log"
    pname = 'a', 'b', 'c', 'alpha', 'beta', 'gamma'
    likely_spacegroups = []
    prob_max = 0
    zone_list = []
    get_elem = lambda n, m, f: f(
                n.getElementsByTagName(m)[0].childNodes[0].data.strip())
    # Reading final choosen cell paramters from the logfile
    logf = open(logf_name)
    logf_raw = logf.read()
    logf_raw.find(" * Dataset ID, project")
    logf.close()
    id1 = logf_raw.find(" * Dataset ID, project")
    id2 = logf_raw.find(" * Number of Columns", id1)
    new_cell_par = map(float, logf_raw[id1:id2].splitlines()[5].split())
    # Reading informations from the XML file
    try:
        dom = minidom.parse(xml_inp)
        cell = dom.getElementsByTagName('cell')[0]
        spg_list = dom.getElementsByTagName('SpacegroupList')[0]
        lattice_sym = dom.getElementsByTagName('LatticeSymmetry')[0]
        new_cell = lattice_sym.getElementsByTagName('cell')[0]
    except:
        os.chdir("..")
        raise
    init_cell_par = dict([(x, get_elem(cell, x, float)) for x in pname])
    #new_cell_par = dict([(x, get_elem(new_cell, x, float)) for x in pname])
    #new_cell_par = [(get_elem(new_cell, x, float)) for x in pname]
    try:
        zone_list = dom.getElementsByTagName('ZoneScoreList')[0]
    except:
        pass
    # Looking at systematique extinctions
    if zone_list:
        print "\n  Systematique extinctions from pointless:"
        print "  Zone Type            axe len.    #obs     Condition    Prob."
        print "  "+60*"-"
        for node in zone_list.getElementsByTagName('Zone'):
            ztype = get_elem(node, 'ZoneType', str)
            nobs = get_elem(node, 'Nobs', int)
            prob = get_elem(node, 'Prob', float)
            condition = get_elem(node, 'Condition', str)
            axe = init_cell_par[ztype[ztype.index("[")+1:ztype.index("]")]]
            all_dat = (ztype, axe, nobs, condition, prob)
            print "%21s %8.1f Å %6d  %12s  %7.3f" % all_dat
    print "\n  Possible spacegroup from pointless:"
    print "  Symbol      num   TotalProb   SysAbsProb   Reindexing_card"
    print "  "+58*"-"
    # looking for most probable spacegroup
    for node in spg_list.getElementsByTagName('Spacegroup'):
        total_prob = get_elem(node, 'TotalProb', float)
        sys_abs_prob = get_elem(node, 'SysAbsProb', float)
        spg_name = get_elem(node, 'SpacegroupName', str)
        spg_num = get_elem(node, 'SGnumber', int)
        reidx_mat = map(int, get_elem(node, 'ReindexMatrix', str).split())
        all_dat = (spg_name, spg_num, total_prob, sys_abs_prob, reidx_mat)
        prob_max = max(total_prob, prob_max)
        print "%11s   #%d  %9.3f    %9.3f   %s" % all_dat
        if total_prob == prob_max:
            likely_spacegroups.append(all_dat)
    os.chdir("..")
    return likely_spacegroups, new_cell_par

def run_pointless(dir_name, hklinp="XDS_ASCII.HKL"):
    # Run pointless and extract pointgroup and spacegroup determination
    # from its xmlout file.
    cmline = "pointless XDSIN %s XMLOUT XDS_pointless.xml"
    cmline += " HKLOUT XDS_pointless.mtz > XDS_pointless.log"
    os.chdir(dir_name)
    os.system(cmline % hklinp)

def pointless(dir_name, hklinp="XDS_ASCII.HKL"):
    run_pointless(dir_name, hklinp)
    return process_pointless_xml()

if __name__ == '__main__':
    import os
    from xml.dom import minidom
    if not is_pointless_installed():
        print "!!  Error. Pointless program doesn't seems to be properly installed."
    else:
        print process_pointless_xml()
        #pointless(".")
