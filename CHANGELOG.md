Release 0.6.6.0 Changes (from 0.6.2.0):
- Now explicitly invoking python2, instead of python to avoid python version conflicts.
- General script: run_xds2staraniso.sh for runing [staraniso](http://staraniso.globalphasing.org) from an XDS_ASCII.HKL or XSCALE.HKL file.
XDSME
- Fix usage info
- Verifying the mean backgroud level from DEFPIX step to eventually fix scale factors in integrate step if starting after INIT step.
- Removed the limitation on maximum number of proc. Let xds_par decide.
- Set default MAXIMUM_NUMBER_OF_JOBS to 1
- Bug fix for catching critical errors in run_idxref
XIO
- Improved handling of Dectris Eiger HDF5 file format.
- Fixing bugs in the miniCBF beam_XY value extraction with hdf2mini-cbf
- Forcing GAIN= 1. in XIO xds export for all Eiger and Pilatus detector frames.
XUPY
- Completing the SPGlib for some centrosymmetric spacegroups

Release 0.5.9.1 Changes:
XDSME
- new option: --optimized --O[1-3]
- new option: --E --exec
- added XDS_PATH env. variable.
- Fix bug for calling option --ice
- Fix bug for calling option --project
XDSCONV
- force option -a (anomalous difference output) to be the default.
- added -o option to define output file names.
- fixe run_phaser.sh script for default HKLOUT name.
- added XDS_PATH env. variable.
XSCALE2
- added option -z
- cleaning dead code
XIO
- Correct MSC ccd header interpretation for spindle orientation.
XUPY
- Change for XDS_PATH env. variable.
