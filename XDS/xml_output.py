#!/usr/bin/env python2
# -*- coding: utf-8 -*-

""" Module to parse XML file from AIMLESS"""

__version__ = "0.0.4"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "01-10-2018"
__copyright__ = "Copyright (c) 2018 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

import sys
import time
from xml.dom import minidom
PNAME = 'a', 'b', 'c', 'alpha', 'beta', 'gamma'
RNAME = 'Overall', 'Inner', 'Outer'
XML_TEMPL_MAIN = """<?xml version="1.0"?>
<AutoProcContainer>\n%s\n</AutoProcContainer>\n"""

XML_TEMPL_AIMLESS = """  <AutoProc>
    <spaceGroup>%(spgn)s</spaceGroup>
    <wavelength>%(wavelength).5f</wavelength>
    <refinedCell_a>%(a).3f</refinedCell_a>
    <refinedCell_b>%(b).3f</refinedCell_b>
    <refinedCell_c>%(c).3f</refinedCell_c>
    <refinedCell_alpha>%(alpha).3f</refinedCell_alpha>
    <refinedCell_beta>%(beta).3f</refinedCell_beta>
    <refinedCell_gamma>%(gamma).3f</refinedCell_gamma>
  </AutoProc>
  <AutoProcScalingContainer>
    <AutoProcScaling>
      <recordTimeStamp>%(scale_runtime)s</recordTimeStamp>
    </AutoProcScaling>
    <AutoProcScalingStatistics>
      <scalingStatisticsType>overall</scalingStatisticsType>
      <resolutionLimitLow>%(res_low_overall).3f</resolutionLimitLow>
      <resolutionLimitHigh>%(res_high_overall).3f</resolutionLimitHigh>
      <rMerge>%(rmerge_overall_overall).3f</rMerge>
      <rMeasWithinIPlusIMinus>%(rmeas_overall).3f</rMeasWithinIPlusIMinus>
      <rMeasAllIPlusIMinus>%(rmeas_overall_overall).3f</rMeasAllIPlusIMinus>
      <rPimWithinIPlusIMinus>%(rpim_overall).3f</rPimWithinIPlusIMinus>
      <rPimAllIPlusIMinus>%(rpim_overall_overall).3f</rPimAllIPlusIMinus>
      <nTotalObservations>%(numb_obs_overall)d</nTotalObservations>
      <nTotalUniqueObservations>%(numb_refl_overall)d</nTotalUniqueObservations>
      <meanIOverSigI>%(mean_ios_overall).2f</meanIOverSigI>
      <completeness>%(compl_overall).1f</completeness>
      <multiplicity>%(multi_overall).1f</multiplicity>
      <ccHalf>%(cc_half_overall).3f</ccHalf>
      <anomalousCompleteness>%(anom_compl_overall).1f</anomalousCompleteness>
      <anomalousMultiplicity>%(anom_multi_overall).1f</anomalousMultiplicity>
      <ccAnomalous>%(anom_cchalf_overall).3f</ccAnomalous>
      <DanoOverSigDano>.856</DanoOverSigDano>
    </AutoProcScalingStatistics>
    <AutoProcScalingStatistics>
      <scalingStatisticsType>innerShell</scalingStatisticsType>
      <resolutionLimitLow>%(res_low_inner).3f</resolutionLimitLow>
      <resolutionLimitHigh>%(res_high_inner).3f</resolutionLimitHigh>
      <rMerge>%(rmerge_overall_inner).3f</rMerge>
      <rMeasWithinIPlusIMinus>%(rmeas_inner).3f</rMeasWithinIPlusIMinus>
      <rMeasAllIPlusIMinus>%(rmeas_overall_inner).3f</rMeasAllIPlusIMinus>
      <rPimWithinIPlusIMinus>%(rpim_inner).3f</rPimWithinIPlusIMinus>
      <rPimAllIPlusIMinus>%(rpim_overall_inner).3f</rPimAllIPlusIMinus>
      <nTotalObservations>%(numb_obs_inner)d</nTotalObservations>
      <nTotalUniqueObservations>%(numb_refl_inner)d</nTotalUniqueObservations>
      <meanIOverSigI>%(mean_ios_inner).2f</meanIOverSigI>
      <completeness>%(compl_inner).1f</completeness>
      <multiplicity>%(multi_inner).1f</multiplicity>
      <ccHalf>%(cc_half_inner).3f</ccHalf>
      <anomalousCompleteness>%(anom_compl_inner).1f</anomalousCompleteness>
      <anomalousMultiplicity>%(anom_multi_inner).1f</anomalousMultiplicity>
      <ccAnomalous>%(anom_cchalf_inner).3f</ccAnomalous>
      <DanoOverSigDano>1.040</DanoOverSigDano>
    </AutoProcScalingStatistics>
    <AutoProcScalingStatistics>
      <scalingStatisticsType>outerShell</scalingStatisticsType>
      <resolutionLimitLow>%(res_low_outer).3f</resolutionLimitLow>
      <resolutionLimitHigh>%(res_high_outer).3f</resolutionLimitHigh>
      <rMerge>%(rmerge_overall_outer).3f</rMerge>
      <rMeasWithinIPlusIMinus>%(rmeas_outer).3f</rMeasWithinIPlusIMinus>
      <rMeasAllIPlusIMinus>%(rmeas_overall_outer).3f</rMeasAllIPlusIMinus>
      <rPimWithinIPlusIMinus>%(rpim_outer).3f</rPimWithinIPlusIMinus>
      <rPimAllIPlusIMinus>%(rpim_overall_outer).3f</rPimAllIPlusIMinus>
      <nTotalObservations>%(numb_obs_outer)d</nTotalObservations>
      <nTotalUniqueObservations>%(numb_refl_outer)d</nTotalUniqueObservations>
      <meanIOverSigI>%(mean_ios_outer).2f</meanIOverSigI>
      <completeness>%(compl_outer).1f</completeness>
      <multiplicity>%(multi_outer).1f</multiplicity>
      <ccHalf>%(cc_half_outer).3f</ccHalf>
      <anomalousCompleteness>%(anom_compl_outer).1f</anomalousCompleteness>
      <anomalousMultiplicity>%(anom_multi_outer).1f</anomalousMultiplicity>
      <ccAnomalous>%(anom_cchalf_outer).3f</ccAnomalous>
      <DanoOverSigDano>.839</DanoOverSigDano>
   </AutoProcScalingStatistics>"""

XML_INTEGRATION_PART = """  <AutoProcIntegrationContainer>
      <Image>
        <fileName>%(name)s</fileName>
        <fileLocation>%(image_path)s</fileLocation>
      </Image>
      <AutoProcIntegration>
        <cell_a>%(ia).3f</cell_a>
        <cell_b>%(ib).3f</cell_b>
        <cell_c>%(ic).3f</cell_c>
        <cell_alpha>%(ialpha).3f</cell_alpha>
        <cell_beta>%(ibeta).3f</cell_beta>
        <cell_gamma>%(igamma).3f</cell_gamma>
        <startImageNumber>%(image_start)d</startImageNumber>
        <endImageNumber>%(image_last)d</endImageNumber>
        <refinedDetectorDistance>%(det_dist).3f</refinedDetectorDistance>
        <refinedXBeam>%(x_beam).2f</refinedXBeam>
        <refinedYBeam>%(y_beam).2f</refinedYBeam>
        <rotationAxisX>%(rvx).6f</rotationAxisX>
        <rotationAxisY>%(rvy).6f</rotationAxisY>
        <rotationAxisZ>%(rvz).6f</rotationAxisZ>
        <beamVectorX>%(bvx).6f</beamVectorX>
        <beamVectorY>%(bvy).6f</beamVectorY>
        <beamVectorZ>%(bvz).6f</beamVectorZ>
      </AutoProcIntegration>
   </AutoProcIntegrationContainer>
  </AutoProcScalingContainer>"""

XML_PROGRAM_PART = """  <AutoProcProgramContainer>
    <AutoProcProgram>
      <processingCommandLine>%(cmd_line)s</processingCommandLine>
      <processingPrograms>XDSME %(xdsme_version)s with XDS %(xds_version)s</processingPrograms>
      <processingStatus>1</processingStatus>
      <processingMessage></processingMessage>
      <processingStartTime>%(exec_time_start)s</processingStartTime>
      <processingEndTime>%(exec_time_end)s</processingEndTime>
      <processingEnvironment>
  Host      : %(hostname)s
  OS        : %(osname)s
  User      : %(username)s
  Directory : %(run_dir)s
  Date      : %(exec_time_end)s
  xdsme     : /data2/bioxsoft/progs/XDSME/xdsme
      </processingEnvironment>
    </AutoProcProgram>
    <AutoProcProgramAttachment>
      <fileType>Log</fileType>
      <fileName>%(xdsmeout)s</fileName>
      <filePath>%(run_dir_p)s</filePath>
    </AutoProcProgramAttachment>
    <AutoProcProgramAttachment>
      <fileType>Log</fileType>
      <fileName>%(aimlessout)s</fileName>
      <filePath>%(run_dir)s</filePath>
    </AutoProcProgramAttachment>
    <AutoProcProgramAttachment>
      <fileType>Result</fileType>
      <fileName>%(mtz_out)s</fileName>
      <filePath>%(run_dir)s/ccp4if</filePath>
    </AutoProcProgramAttachment>
    <AutoProcProgramAttachment>
      <fileType>Result</fileType>
      <fileName>XSCALE.LP</fileName>
      <filePath>%(run_dir)s</filePath>
    </AutoProcProgramAttachment>
  </AutoProcProgramContainer>"""

TRANSLATE_RESULTS = {
    'ResolutionLow':'res_low', 'ResolutionHigh':'res_high',
    'Rmerge':'rmerge', 'RmergeOverall':'rmerge_overall',
    'Rmeas':'rmeas', 'RmeasOverall':'rmeas_overall',
    'Rpim':'rpim', 'RpimOverall':'rpim_overall',
    'NumberObservations':'numb_obs', 'NumberReflections': 'numb_refl',
    'MeanIoverSD':'mean_ios', 'CChalf':'cc_half', 'Completeness':'compl',
    'Multiplicity':'multi', 'MeanChiSq':'mean_chisq',
    'AnomalousCChalf':'anom_cchalf', 'AnomalousCompleteness':'anom_compl',
    'AnomalousMultiplicity':'anom_multi',
}

get_elem = lambda n, m, f: f(
    n.getElementsByTagName(m)[0].childNodes[0].data.strip())
get_val = lambda m, f: f(m.childNodes[0].data.strip())

def parse_aimless_xml(xml_inp_name):
    "Parse the aimless XML file and return a dictonary with selected values"

    #xml_inp = file_name_id + ""
    extr = {}
    # Raw read of the XML file
    filein = open(xml_inp_name)
    xml_raw = filein.read()
    filein.close()

    # Reading informations from the XML file
    try:
        dom = minidom.parse(xml_inp_name)
        cell = dom.getElementsByTagName('cell')[0]
        #lattice_sym = dom.getElementsByTagName('LatticeSymmetry')[0]
        #new_cell = lattice_sym.getElementsByTagName('cell')[0]
    except:
        raise
    extr.update(dict([(x, get_elem(cell, x, float)) for x in PNAME]))
    extr['wavelength'] = get_val(dom.getElementsByTagName('Wavelength')[0], float)
    extr['spgn'] = get_val(dom.getElementsByTagName('SpacegroupName')[0], str)

    tsi = xml_raw.find('RunTime="')
    #extr['scale_runtime'] = xml_raw[tsi:tsi+40].split('"')[1]
    extr['scale_runtime'] = time.strftime('%F %X')
    results = dom.getElementsByTagName('Result')[0]
    for elem in TRANSLATE_RESULTS:
        xxe = results.getElementsByTagName(elem)[0]
        for shell in RNAME:
            xxk = '%s_%s' % (TRANSLATE_RESULTS[elem].lower(), shell.lower())
            xxv = get_elem(xxe, shell, float)
            extr[xxk] = xxv

    return extr

def xml_aimless_to_autoproc(file_name_id):
    "Convert XML from aimless to autoProc style"
    return XML_TEMPL_AIMLESS % parse_aimless_xml(file_name_id)

def get_process_context():
    return " ".join(sys.argv)

if __name__ == '__main__':
    INP_NAME = "splTsm_2_aimless.xml"
    print get_process_context()
    xmldict = parse_aimless_xml(INP_NAME)
    #print(XML_TEMPL_AIMLESS % xmldict)
    #dim = minidom.parse(inp_name)
    #for e in PNAME:
    #    print e, get_val(dim.getElementsByTagName(e)[0], float)
    #for e in RNAME:
    #    print e, get_val(dim.getElementsByTagName(e)[0], float)
