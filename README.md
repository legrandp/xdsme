# XDSME: XDS Made Easier

The XDSME package contains a collection of python scripts made to simplify the processing of crystal diffraction images with the XDS Program Package (X-ray Detector Software, http://xds.mpimf-heidelberg.mpg.de/). Provided that the diffraction parameters are well recorded in the diffraction image headers, XDS data processing can be started with a simple command line like:

```bash
 $ xdsme pos1_1_???.img
 ```

## Features

Supported detector image format include: PILATUS, EIGER (experimental), ADSC, MARCCD, MAR345, SATURN. There is an experimental support for RAXIS and MAR555 detectors (limited by lack of test images).

The main scripts are:
 - xdsme for data processing,
 - xscale.py for scaling of multiple reflexion files,
 - and xdsconv.py for file format conversion,
 - XOalign.py for goniometer setting calculation (can be set to work with different type of goniometer including Kappa, mini-Kappa, Euler...)
 - xds2mos.py or xds2dnz.py ... (for convertion of orientation matrices)

## Install

All scripts are pure python code, so the only dependency is Python version >= 2.2. It should work on any linux or mac-osx directly after unpacking by adding the xdsme/bin/noarch dir to your PATH variable.

## Examples

A typical session will look like that:
```bash
 $ xdsme  col1_1_*.img
 $ xdsme  --O -r 2.1 -s P3121 -c  “59 59 123 90 90 120” col1_1_*.img
 $ cd  xds_process_col1_1

 $ xscale.py  XDS_ASCII.HKL ../xds_process_lowres/XDS_ASCII.HKL
 $ xdsconv.py  XSCALE.HKL  8  Se shelx
 $ xdsconv.py  XSCALE.HKL  8  Se solve
 $ xdsconv.py  XSCALE.HKL  8  Se ccp4 shelx/XDS_ASCII_F4.hkl
 $ xdsconv.py  XSCALE.HKL  8  Se phaser ccp4/XDS_ASCII.mtz
 $ xdsconv.py  XSCALE.HKL  cns
```
For most scripts, there is also an in-line help aviable with the argument -h (xdsme -h, xdsconv.py -h). I am happy to read any comments, suggestions or report of bugs or problems for installing or using these scripts.

## License

XDSME is available under the New BSD License [Read the copyright statement and license] (/legrandp/xdsme/LICENSE)
          http://www.opensource.org/licenses/bsd-license.php


## Included code

Some code comes from the [CgKit project](http://cgkit.sourceforge.net)
