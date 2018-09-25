# XDSME: XDS Made Easier

The XDSME package contains a collection of python scripts made to simplify the processing of crystal diffraction images with the XDS Program Package (X-ray Detector Software, http://xds.mpimf-heidelberg.mpg.de/). Provided that the diffraction parameters are well recorded in the diffraction image headers, XDS data processing can be started with a simple command line like:

```bash
 $ xdsme pos1_1_???.img
 or
 $ xdsme pos1_1_???.cbf
 or 
 $ xdsme pos1_1_master.h5
 ```

## Features

Supported detector image format include: PILATUS, EIGER, ADSC, MARCCD, MAR345, SATURN. There is an experimental support for RAXIS and MAR555 detectors (limited by lack of test images).

xdsme and xds are able to directly process compressed frames in the .gz or .bz2 formats. 

The main scripts are:
 - xdsme for data processing,
 - xscale.py for scaling of multiple reflexion files,
 - and xdsconv.py for file format conversion,
 - XOalign.py for goniometer setting calculation (can be set to work with different type of goniometer including Kappa, mini-Kappa, Euler...)
 - xds2mos.py or xds2dnz.py ... (for convertion of orientation matrices)

## Install

You can either:
 - *Recommended way:* Download the [lastest released version](https://github.com/legrandp/xdsme/releases/latest) and gunzip it. 
 - or ownload the xdsme-master.zip file and unzip it.
 - or use git to clone the repository. This method makes updating to a new version easier.
 - Then add the xdsme-master/bin/noarch to your PATH variable (export PATH=$PATH:$HOME/progs/xdsme-master/bin/noarch),
```bash
 $ cd $HOME/progs 
 $ git clone https://github.com/legrandp/xdsme.git
 $ export PATH=$PATH:$HOME/progs/xdsme/bin/noarch
```

To update xdsme, you will only need to run the following command:
```bash
 $ git pull origin
```

All scripts are pure python code, so the only dependency is Python version >= 2.5 (or >= 2.7 for processing HDF5 files). It should work on any linux or mac-osx directly after unpacking by adding the xdsme/bin/noarch dir to your PATH variable.

## Examples

A typical session will look like that:
```bash
 $ xdsme  col1_1_*.img
 $ xdsme  --O -r 2.1 -s P3121 -c  “59 59 123 90 90 120” col1_1_*.img
 $ cd  xdsme_col1_1

 $ xscale.py  XDS_ASCII.HKL ../xds_process_lowres/XDS_ASCII.HKL
 $ xdsconv.py  XSCALE.HKL  8  Se shelx -o SePeak
 $ xdsconv.py  XSCALE.HKL  8  Se ccp4 shelx/XDS_ASCII_F4.hkl
 $ xdsconv.py  XSCALE.HKL  8  Se phaser ccp4/XDS_ASCII.mtz
 $ xdsconv.py  XSCALE.HKL  [ccp4if|cns|replace|solve|epmr|amore]
```
For most scripts, there is also an in-line help aviable with the argument -h (xdsme -h, xdsconv.py -h). I am happy to read any comments, suggestions or report of bugs or problems for installing or using these scripts.

## License

XDSME is available under the New BSD License [Read the copyright statement and license] (/legrandp/xdsme/LICENSE)
          http://www.opensource.org/licenses/bsd-license.php

## Citation

If you wish to cite this work, you can use the following reference:

Legrand, P. XDSME: XDS Made Easier (2017) GitHub repository, https://github.com/legrandp/xdsme 
DOI 10.5281/zenodo.837885

## Acknowledgements

Many thanks to the following people how have contributed with pieces of code:
 - Miguel Ortiz Lombardía (CNRS, Marseille, France),
 - Jun Aishima (Australian Synchrotron, Melbourne, Australia), 
 - Clemens Vonrhein (Global Phasing, Cambridge, UK),
 - Ludovic Pecqueur (College de France, Paris, France),
 - Michal Babiak (CEITEC, Brno, Czech Republic)

## Included code

Some code comes from the [CgKit project](http://cgkit.sourceforge.net)
The [pyfive module](https://github.com/jjhelmus/pyfive) is also distributed in the released tarball.
