######################################################################
XDSME: XDS Made Easier
Copyright (C) 2007-2015
          Pierre Legrand.
License:  New BSD License, see LICENSE file or
          http://www.opensource.org/licenses/bsd-license.php
######################################################################

xdsme is a collection of python scripts made to simplify the processing of crystal diffraction images with the XDS Program Package (X-ray Detector Software, http://xds.mpimf-heidelberg.mpg.de/). Provided that the diffraction parameters are well recorded in the diffraction image headers, XDS data processing can be started with a simple command line like:

 $ xdsme pos1_1_???.img

Supported detector image format include: PILATUS, ADSC, MARCCD, MAR345, SATURN. There is an experimental support for RAXIS and MAR555 detectors (limited by lack of test images).

The main scripts are: - xdsme (XDS.py), xscale.py and xdsconv.py for the data processing, scaling andfile conversion. - XOalign.py for the goniometer setting calculation (can be set to work with different type of goniometer including Kappa, mini-Kappa, Euler...). - xds2mos.py or xds2dnz.py ... (for convertion of orientation matrices)

All scripts are pure python code, so the only dependency is Python version >= 2.2. It should work on any linux or mac-osx directly after unpacking by adding the xdsme/bin/noarch dir to your PATH variable.

A typical session will look like that:

 $ xdsme  col1_1_*.img 
 $ xdsme  --O -r 2.1 -s P3121 -c  “59 59 123 90 90 120” col1_1_*.img 
 $ cd  xds_process_col1_1 
 
 $ xscale.py  XDS_ASCII.HKL ../xds_process_lowres/XDS_ASCII.HKL 
 $ xdsconv.py  XSCALE.HKL  8  Se shelx 
 $ xdsconv.py  XSCALE.HKL  8  Se solve 
 $ xdsconv.py  XSCALE.HKL  8  Se ccp4 shelx/XDS_ASCII_F4.hkl 
 $ xdsconv.py  XSCALE.HKL  8  Se phaser ccp4/XDS_ASCII.mtz 
 $ xdsconv.py  XSCALE.HKL  cns 

For most scripts, there is also an in-line help aviable with the argument -h (xdsme -h, xdsconv.py -h). I am happy to read any comments, suggestions or report of bugs or problems for installing or using these scripts.

Included code:
http://cgkit.sourceforge.net
