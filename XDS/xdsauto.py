#!/usr/bin/env python
# TODO.
# Verifications:
#   - verify the input parameter from auto.par (creat a dic of alowed keys)
#       give warning for unalowed keys.
#   - try the xds_par or xds, using which (/...xds or /...xds_par)
#   - Make a test procedure.
# Step III and IV: If more than 1 sector (see lookup_frame())
#               III  -> Profile FIT with each sector with (SPOT RANGE...)
#               III  -> Integrate separetly each sector (INTEGRATE)
#                IV  -> Merge and reindex with XSCALE,,,
# look at statistics of occuring space-group to give hint when systematic (see pdb)
#  abscence are not 

"""
Try to automate XDS data processing
"""
import sys, shutil
import os, re
from xupy import *
import xdsSetupDB

__version__ = "0.6.12"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "21-02-2004"
__copyright__ = "Copyright (c) 2005 Pierre Legrand"
__license__ = "GPL"

_job_angle = 4., 8., 16., 32, 90, 360
_nopt = 2         # Default number of profil optimization steps.
_idxref_fit_limit = 40
_nspot_slice = 8 # number of frames/sequence in COLSPOT and Init. Prof determination.
_wcut_rmeas = 1.8
_nspot_min = 200  # Minimium number of spots for autoindexing.
_nDigits = 3 # Number of digits in frame numbering

# Determine the number of max threads used in the pool
if os.popen("uname -a").read().count("openmosix"):
    _nThreadsMax = len(os.popen("ls /proc/hpc/nodes").readlines()) -1
    _nThreadsMax = 2
else:
    _nThreadsMax = 1


Step_str = ["I","II","III","IV","V","VI"]

def adpParamsParser(inp_str):
    """ return a dict of ADP input parameters. Functional Programming style."""
    lines = inp_str.splitlines()
    nocomm = lambda s: s[0] != '#'
    eq2pair = lambda s: s.split('=')
    strip = lambda s: s.strip()
    def assign(pair):
        try: pair = [strip(pair[0]),int(pair[1])]
        except:
            try: pair = [strip(pair[0]),float(pair[1])]
            except: pair = map(strip,pair)
        return pair
    ass = map(assign,map(eq2pair,filter(nocomm,filter(None,map(strip,lines)))))
    return dict(ass)

def get_adpParams(parf = "auto.par"):                                            
    # Defaults values:
    auto = XParam()
    auto.img_dir = "../img/"
    auto.beamline = "ESRF_BM30A"
    auto.spg = 0
    auto.cell = 0
    auto.anomalous = "no"
    auto.resolution_range = "30 0.0"
    
    try: adpParams_raw = opReadCl(parf)
    except:
        print_o(fmt_err % ("Can't open file %s\n" % parf),f3)
        stop()
    try:
        auto.add(adpParamsParser(adpParams_raw))
    except:
        print_o(fmt_err%("Problems while reading the file %s\n"%parf),f3)
        stop()
    # correct some midification in the automatique auto.par file creation at ESRF by ProDC.
    try:
        auto.phi_init = auto.phi_init_val
    except: pass
    return auto
            
def get_num(seq,num):
    L = []
    for i in num: L.append(seq[i])
    return L
                 
def print_o(s,_files):
    if type(s) == tuple or type(s) == list:
        s = "".join(map(str,s))
    if len(_files) <= 2:
        for _file in _files: _file.write(s)
    elif len(_files) == 3:
        for _file in _files[:2]: _file.write(s)
        _files[2].write(s)
        #files[2].write(s.replace("\n","<br>\n"))

def print_o2(t,v,f):
    print_o("  %-49s%24s\n" % (t,v),f)

def print_o3(t,v1,v2,f):
    print_o("\t   %-22s %14s    %-s\n" % (t,v1,v2),f)
    
def print_stage(s,files):
    print_o(s,files[:2])
    print_o("<h2>"+s+"</h2>\n",(files[-1],))
    
def stop(s=""):
    flog.close()
    fhtml.write("</pre></BODY>")
    fhtml.close()
    os.chdir("../")
    saveLastVersion(["xdsauto.log","xdsauto.html"])
    print "  %s  STOP\n" % s
    sys.exit()
  
if __name__=='__main__':
    step = 0
    if "-test" in sys.argv:
        print "test"
        get_maxResolution()
        sys.exit()
    if len(sys.argv) >= 1:
        flog = open("xdsauto.log","w")
        fhtml = open("xdsauto.html","w")
        fhtml.write("<HTML> <HEAD>\n<TITLE>Automated XDS Data Processing")
        fhtml.write('</TITLE>\n</HEAD>\n<BODY BGCOLOR="#ffffff">\n<pre>')
        
        #dirname = os.path.split(os.getcwd())[1]
        f1, f2, f3 = (sys.stdout), (sys.stdout,flog), (sys.stdout,flog,fhtml)
        print_o("\n",f3)
        
        for arg in sys.argv[1:]:
            if arg[:2] == "-O" and len(arg) == 3 and arg[-1].isdigit():
                opt = int(arg[-1])
                if 6 >= opt >= 1:
                    _nopt = opt
                    print_o2("INFO:  Number of profil optimization steps",_nopt,f3)
                    sys.argv.remove(arg)
                else:
                    print "\n  ERROR: Illegal argument for the -O option\n"
                    sys.exit()
            if arg[:2] == "-T" and len(arg) == 3 and arg[-1].isdigit():
                opt = int(arg[-1])
                if 5 >= opt >= 1:
                    _nThreadsMax = opt
                    print_o2("INFO:  Number of profil optimization steps",_nopt,f3)
                    sys.argv.remove(arg)
                else:
                    print "\n  ERROR: Illegal argument for the -T option\n"
                    sys.exit()
            if arg[:2] == "-S" and len(arg) >= 3 and arg[-1].isdigit():
                opt = int(arg[2:])
                if opt >= 1:
                    _nspot_slice = opt
                    print_o2("INFO:  Number of images read for spot collection",opt*2,f3)
                    sys.argv.remove(arg)
                else:
                    print "\n  ERROR: Illegal argument for the -I option\n"
                    sys.exit()
        for arg in sys.argv[1:]:
            if len(arg) == 2 and arg[-1].isdigit():
                step = int(arg[-1])
                
    if step >= 0:
        fmt_shape = "beam divergence / reflecting range"
        fmt_err = "\n  ERROR!  %s\n\n"
        fmt_stg = "\n\n\n\t Stage %25s\n\n"
        fmt_link_xds = '  <A href="xds/%s">XDS log file</A> for this stage.<br>'
        fmt_cell = "%6.1f%6.1f%6.1f"
        
        BRAVAIS = ""
        
        parf = "auto.par"
        print_o("\n",f3)
        print_o2("Reading Data collection parameters in file", parf,f3)
        print_o("\n\n",f3)
        auto = get_adpParams(parf)
        
        try: shutil.copyfile("XDS.INP","XDS.INP.old")
        except:  pass

        if not "xds" in os.listdir("."): os.mkdir("xds")
        os.chdir("xds")
        
        if auto.img_dir[0] != "/": auto.img_dir = "../" + auto.img_dir
        if auto.img_dir[-1] != "/": auto.img_dir += "/"
	
	dataCol = DataCollectInfo((auto.img_dir,auto.prefix[:-1],
	                                _nDigits, auto.suffix))
        
        keys_show = ["beamline","frame_first",
                     "frame_last","delta_phi","anomalous","wavelength",
                     "distance","beam_x","beam_y","resolution_range",
                     "cell","spg"]
        
        for par in keys_show:
            if auto[par]: print_o("\t%-20s\t%s\n" % (par, auto[par]),f3)
	
        print_o("\n\tLooking for images\t%s\n" % dataCol.getXDSTemplate(),f3)

        xpar = xdsInp2Param(inp_str=xdsinp_base)
        BL = getattr(xdsSetupDB,auto.beamline)
        xpar.mix(BL)
           
        print_o("\n\n",f3)
        if step >= 1: print_o2("INFO:  Jumping directly to step",
	                                          Step_str[step-1],f3)
	
        if not dataCol.lookup_image_ranges():
            print fmt_err%("No frame found with template %s" % 
	                                     dataCol.getXDSTemplate())
            stop()

        print_o2("INFO:  First and last frame number found",
                       "%d - %d" % tuple(dataCol.get_range()),f3)
        if _nThreadsMax > 1:
            print_o2("INFO:  Maximum number of threads used",
                       "%d" % (_nThreadsMax),f3)
        if len(dataCol.image_ranges) >= 2:
            print_o2("WARNING! Discontinuous frame sequence",
	              dataCol.image_ranges,f3)
        minf, maxf = dataCol.get_range(minf=auto.frame_first, maxf=auto.frame_last)
        xpar.DATA_RANGE = [minf, maxf]
        
        nmin_f = 2*_nspot_slice+1
        if maxf - minf + 1 >= nmin_f:
            # use two range ex: i-i+2, f-2,f
            # with f at maximum 90 degre distance of i
            max_frame = min(maxf, minf+int(89./auto.delta_phi+_nspot_slice))
            xpar.SPOT_RANGE=((minf, minf + _nspot_slice - 1),
                             (max_frame - _nspot_slice + 1, max_frame))
        else:
            xpar.SPOT_RANGE = (minf, min(minf + nmin_f - 1,maxf)),
        
        xpar.STARTING_FRAME = dataCol.getClosestImage(auto.frame_first)
        xpar.BACKGROUND_RANGE = [xpar.STARTING_FRAME, 
	                         dataCol.getClosestImage(auto.frame_first + 4)]
        xpar.OSCILLATION_RANGE = auto.delta_phi
        xpar.DELPHI= 4*auto.delta_phi
        xpar.STARTING_ANGLE = auto.phi_init
        xpar.X_RAY_WAVELENGTH = auto.wavelength
        xpar.DETECTOR_DISTANCE = auto.distance
        xpar.NAME_TEMPLATE_OF_DATA_FRAMES = dataCol.getXDSTemplate()
        xpar.ORGX = auto.beam_y/xpar.QX
        xpar.ORGY = auto.beam_x/xpar.QY
        if auto.cell:
            auto.cell = map(float,auto.cell.split())
            if 0. not in auto.cell:
                xpar.UNIT_CELL_CONSTANTS = auto.cell
                xpar.SPACE_GROUP_NUMBER = 1
            else:
                print "WARNING input cell contains zerro parameter!"
                stop()
		
        r = map(float,auto.resolution_range.split())
        xpar.INCLUDE_RESOLUTION_RANGE = max(r), min(r)
        if min(r) != 0.0:
             print_o2("INFO:  Truncating the highest resolution to ",
                                    "%.2f" % min(r),f3)
        xpar.FRIEDELS_LAW = "TRUE" # set to true first for spacegroupe comparaison
        if auto.suffix == "pck":
            xpar.NAME_TEMPLATE_OF_DATA_FRAMES += " CCP4"
        
    if step <= 1:
        Stage = " I -  Peak Search"
        #==================================#
        n = 0
        xdsout = "xds_i_%d.out" % n
        print_stage(fmt_stg % Stage,f3)
        fhtml.write(fmt_link_xds % xdsout)
        
        s = ""
        for r in xpar.SPOT_RANGE: s += "%d to %d and " % tuple(r)
        print_o2("Using images","%s" % s[:-5],f3)
        
        xpar.JOB = "XYCORR", "INIT", "COLSPOT"
        run_xds(xpar, inp_f=None, out_f=xdsout, save=0)
        saveLastVersion(("XDS.INP","XYCORR.LP","INIT.LP","COLSPOT.LP"),"_I")
        
        initlp = opReadCl("INIT.LP")
        colspotlp = opReadCl("COLSPOT.LP")
        ind1 = initlp.find("MEAN GAIN VALUE")
        ind2 = initlp.find("AVERAGE BACKGROUND COUNTS IN A DATA IMAGE PIXEL")
        ind3 = colspotlp.find("NUMBER OF DIFFRACTION SPOTS ACCEPTED")
        nspot = int(colspotlp[ind3:ind3+58].split()[-1])
        m_gain = float(initlp[ind1:ind1+51].split()[-1])
        m_count = float(initlp[ind2:ind2+58].split()[-1])
        
        print_o2("Mean Gain value", "%8.2f" % m_gain ,f3)
        print_o2("Average background counts","%8.1f" % m_count ,f3)
        print_o2("Number of diffraction spots found","%12d" % nspot ,f3)
        
        while nspot <= _nspot_min:
            n += 1
            xpar.STRONG_PIXEL -= 2
            xpar.MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT -= 1
            if xpar.STRONG_PIXEL < 3. or \
               xpar.MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT < 5: pass
            print_o("  WARNING: Number of spots found is too low (%d <= %d)\n"
                                              % (nspot,_nspot_min),f3)
            print_o("  Recollecting spots with minimum I/sigma = %.1f\n\n"
                                               % xpar.STRONG_PIXEL,f3)
            xpar.JOB = "COLSPOT"
            xdsout = "xds_i_%d.out" % n
            fhtml.write(fmt_link_xds % xdsout)
            run_xds(xpar,inp_f=None,out_f=xdsout,save=0)
            saveLastVersion(("XDS.INP","COLSPOT.LP",),"_I%d" % n)
            colspotlp = opReadCl("COLSPOT.LP")
            nspot = int(colspotlp[ind3:ind3+58].split()[-1])
            
        #print_o2("Number of diffraction spots found","%12d" % nspot ,f3)

    else: opWriteCl("XDS.INP","%s" % xpar)
    
    if step <= 2:
        Stage = " II -  Auto-indexing"
        #==================================#
        xdsout = "xds_ii.out"
        print_stage(fmt_stg % Stage,f3)
        fhtml.write(fmt_link_xds % xdsout)
        
        xpar.JOB = "IDXREF"
        run_xds(xpar,out_f=xdsout,save=0)
        saveLastVersion(("XDS.INP","XPARM.XDS","IDXREF.LP",),"_II")
        
        ilp = opReadCl("IDXREF.LP_II")
        ilp_err = ilp.find("!!! ERROR !!!")
        if ilp_err >= 0:
            print_o("\n  WARNING %s\n\n" % ilp[ilp_err:].splitlines()[0],f3)
        ilp_ind = ilp.find("INDEXING OF OBSERVED SPOTS IN SPACE GROUP #   1")+54
        indexed = ilp[ilp_ind:ilp_ind+38].split()
        num_indexed, num_total = int(indexed[0]),  int(indexed[3])
        ilp_xyz = ilp.index("SPOT    POSITION (PIXELS)", ilp_ind)+25
        ilp_xyz = ilp[ilp_xyz:ilp_xyz+67].split()
        stdev_xy, stdev_z = float(ilp_xyz[0]), float(ilp_xyz[-1])
        
        reducedCell = get_xparm_cell("XPARM.XDS_II")
        reducedCell_volum = Lattice(reducedCell,"aP").volum()

        print_o2("Selected reduced cell axes", fmt_cell % reducedCell[:3],f3)
        print_o2("Selected reduced cell angles", fmt_cell % reducedCell[3:],f3)
        print_o2("Fraction of indexed spots","(%d/%d) %7.1f%%" % 
                    (num_indexed, num_total,100.*num_indexed/num_total),f3)
        print_o2("Standard deviation of spot    position (pixels) ",
                                                       "%8.2f" % stdev_xy,f3)
        print_o2("Standard deviation of spindel position (degrees)",
                                                       "%8.2f" % stdev_z,f3)
        print_o("  \n\n  Selected list of possible Bravais lattices\n",f3)
        print_o(("\n  Bravais Lattice Symm   fit"+15*" "+"cell"+17*" "+
                "ORDER\n\t"+26*" "+"a     b     c  alp  bet  gam\n\n"),f3)
                
        for LAT in select_lattices(limit=100, idxref="IDXREF.LP_II"):
             rc_au_vol, au_vol = "", ""
             LAT.idealize()
             brav = BravaisDico[LAT.Bravais_type[0]]
             if LAT.fit <= _idxref_fit_limit:
                  au_vol = LAT.multiplicity/(LAT.volum()/reducedCell_volum)
                  rc_au_vol = int(round(au_vol)) * "*"
                  au_vol = "%.0f" % au_vol
             fmt_l = "  %1s %-13s (%2s) %5.1f "+"%s"+"%3s  %s\n"
             print_o(fmt_l % (LAT.Bravais_type[1],brav,LAT.Bravais_type,LAT.fit,
                              LAT.prt(fmt=3*"%6.0f"+3*"%5.0f"),au_vol,rc_au_vol),f3)
        print_o("\n  ORDER is defined as the ratio: reduce cell volum/asym. unit volum.\n", f3)
        
    if step <= 3:
        Stage = " III - Integrate and Optimize Profiles"
        #==============================================#
        print_stage(fmt_stg % Stage,f3)
        
        xpar.JOB = "DEFPIX", "INTEGRATE"
        saved_file = "XDS.INP","DEFPIX.LP","INTEGRATE.LP"
        xpar.DATA_RANGE = [minf, maxf]
        reducedCell = get_xparm_cell("XPARM.XDS_II")
        print_o2("Selected reduced cell axes",fmt_cell % reducedCell[:3],f3)
        print_o2("Selected reduced cell angles",fmt_cell % reducedCell[3:],f3)
        print_o2("Runing integration with Space Group",SPGlib[1][1],f3)
        print_o2("Anomalous set","OFF",f3)
        print_o2("Number of optimization steps","%s" % _nopt,f3)
        
        for i in range(_nopt):
            dataCol.lookup_image_ranges()
            collect_range = dataCol.get_range(minf=auto.frame_first,
                                              maxf=auto.frame_last)
            print collect_range, auto.frame_last
            nframe = int(_job_angle[i]/auto.delta_phi)-1       
            collect_range[1] = min(collect_range[1],collect_range[0]+nframe)
            print collect_range, auto.frame_last
           
            xpar.DATA_RANGE = collect_range
            print_o("\n",f3)
            print_o2("Step %d/%d  Integrate image range" % (i+1,_nopt),
                                 "%d to %d" % tuple(xpar.DATA_RANGE),f3)
            xdsout = "xds_iii_%d.out"%i
            fhtml.write(fmt_link_xds % xdsout)
            run_xds(xpar, out_f=xdsout, save=0)
            saveLastVersion(saved_file, "_III_%d" % i)
            saveLastVersion(saved_file, "_III_last")
            
            xpar.JOB = "INTEGRATE"
            saved_file = "XDS.INP","INTEGRATE.LP"
            #if i == 0:
            #    ipar = getProfilRefPar(infile="INTEGRATE.LP_III_0",init="yes")
            #    print_o2("Initial fitted "+fmt_shape, "%.3f / %.3f" % 
            #                (ipar.BEAM_DIVERGENCE,ipar.REFLECTING_RANGE),f3)
            xpar.mix(getProfilRefPar(infile="INTEGRATE.LP_III_%d" % i))
            print_o2("New refined    "+fmt_shape, "%.3f / %.3f" % 
                            (xpar.BEAM_DIVERGENCE,xpar.REFLECTING_RANGE),f3)
        xpar.REFINE_INTEGRATE = "ORIENTATION", "BEAM",  "CELL" #"DISTANCE",

    if step <= 4:
        Stage = " IV -  Testing Laue Symmetries"
        #=====================================================#
        print_stage(fmt_stg % Stage,f3)
 
        selected = select_lattices(idxref="IDXREF.LP_II")
        print_o(("\n Lat  fit"+13*" "+"cell"+15*" "+
                "Symm  Compar Compl  Rsym Rmeas I/sigI Misf ORDER\n"+15*" "+
                "a    b    c  al  be  ga"+15*" "+3*"   (%)"+"\n\n"),f3)
        
        stage4= []
        
        reducedCell = get_xparm_cell("XPARM.XDS_II")
        reducedCell_volum = Lattice(reducedCell,"aP").volum()
        
        for newLatt in selected:
             
             newLatt.idealize()
             print_o(("%3s %5.1f" % (newLatt.Bravais_type,newLatt.fit),\
                      " %5.0f%5.0f%5.0f%4.0f%4.0f%4.0f " % newLatt.cell),f3)
             
             n = 0
             for laue in Bravais_to_Laue[newLatt.Bravais_type]:
                 laueLatt = Lattice(newLatt.cell, newLatt.Bravais_type, symmetry=laue[1])
                 laueLatt.multiplicity = SPGlib[laue[1]][2]
                 
                 au_fraction = int(round(laueLatt.multiplicity /
                                    (laueLatt.volum()/reducedCell_volum)))
                 if n != 0: print_o(38*" ",f3)
                 print_o("%6s" % (laue[2]),f3)
                 
                 xpar.JOB = "CORRECT"
                 xpar.SPACE_GROUP_NUMBER = laueLatt.symmetry_num
                 xpar.UNIT_CELL_CONSTANTS = laueLatt.cell
                 xpar.REIDX = newLatt.reindexing
                 run_xds(xpar, inp_f=None, out_f="xds_iv.out", save=0)
                 
                 s = resum_scaling()
                 if not s:
                     print_o(fmt_err % "in CORRECT step.",f3)
                     continue
                     
                 for k in ['compl','compl3','rsym','rmeas','rmeas3','isig','compar']:
                     setattr(laueLatt,k,getattr(s,k))
                 laueLatt.stats = s

                 print_o(" %7d" % (s.compar),f3)
                 for k in ['compl','rsym','rmeas','isig']:
                     print_o( "%6.1f" % (s[k]),f3)

                 print_o("%5d   %d\n" % (s.misfit, au_fraction),f3)
                 saveLastVersion(("XDS.INP","CORRECT.LP"),"_IV_laue_%s"%laue[2])
                 stage4.append(laueLatt)
                 n = 1
        
        import pickle
        pickle.dump(stage4,open(auto.prefix+'LaueSymTest.pickle','w'),bin=1)
        
        listIsig = []
        listRmeas = []
        listCompl = [0.01]
        stage4selected = []
        
        _minNumberScaledRefl = 60
        for L in stage4:
            if L.compar > _minNumberScaledRefl:
                listIsig.append(L.isig)
        
        _isigCutoff = 0.5
        listIsig.sort()
        print listIsig
        for L in stage4:
            if L.isig >= _isigCutoff*listIsig[-1]:
                listRmeas.append([L.rmeas, L.compar])
                listCompl.append(L.compl)
                stage4selected.append(L)
        listCompl.sort()

        print_o("\n  ORDER is defined as the ratio: reduce cell volum/asym. unit volum.\n", f3)
        if len(listRmeas) > 1:
            wMeanRmeas = wMean(listRmeas)
            wRMSRmeas = wStandardDeviation(listRmeas)
            MeanRmeas = mean(listRmeas)
            RMSRmeas = standardDeviation(listRmeas)
        else:
            wMeanRmeas = MeanRmeas = listRmeas[0][0]
            wRMSRmeas = RMSRmeas = 0.
            
        print_o("\n            Mean and RMS for Rmeas =  %7.2f  %7.2f\n" % \
                           (MeanRmeas, RMSRmeas),f3)
        print_o("   Weighted Mean and RMS for Rmeas =  %7.2f  %7.2f\n" % \
                           (wMeanRmeas, wRMSRmeas),f3)
        
        print_o("\n  Chosen Lattice by order of preference:\n\n",f3)
        
        i = 1
        if len(listRmeas) >= 1:
            for a in range(len(stage4selected)):
                if i <= len(stage4selected) and i < 4:
                    for L in stage4selected:
                        #print i, L.symmetry_num, listCompl[-i]
                        if L.compl == listCompl[-i] \
                        and 0 < L.rmeas <= wMeanRmeas+wRMSRmeas*_wcut_rmeas:
                           if not auto.spg:
			       auto.spg = L.symmetry_num
                           if i == 1:
			       print_o("<b>",(fhtml,))
                           print_o("  %d) %s %5s %4d" % \
                                (i, L, L.symmetry_str1, L.symmetry_num),f3)
                           print_o("%5.1f %5.1f\n" % (L.compl,L.rmeas),f3)
                           if i == 1: print_o("</b>",(fhtml,))
                           i += 1
        else:
            print_o("\n  Autoindexing failed.\n\n",f3)
            stop()
               
        
    if step <= 5:
        Stage = " V -  Integration and Scaling"
        #=========================#
        print_stage(fmt_stg % Stage,f3)

        if xpar.has_key("REIDX"): del xpar.REIDX
        if auto.spg:
            BRAVAIS = ""
            spgn = 0
            if type(auto.spg) == str:
                auto.spg = auto.spg.upper()
                for _spg in SPGlib.keys():
                    if auto.spg == SPGlib[_spg][0]:
                        spgn = _spg
            elif type(auto.spg) == int:
                spgn = auto.spg
                if spgn not in SPGlib.keys():
                    print_o(fmt_err % \
                    "Space group not possible for biological macromolecules\n",f3)
                    stop()
            if not spgn:
                print_o(fmt_err % ("Unconsistant spage group: %s\n"%auto.spg),f3)
                stop()
            SPG = SPGlib[spgn]
        
        if not BRAVAIS:
            for bravais in Bravais_to_Laue.keys():
                for laue in Bravais_to_Laue[bravais]:
                    if spgn in laue[3]:
                        BRAVAIS = bravais
        
        #if BRAVAIS == "mI": BRAVAIS = "mC"
        
        selected1 = select_lattices(limit=_idxref_fit_limit, idxref="IDXREF.LP_II")
        selected2 = []
        bestfit = 1000
        for _latt in selected1:
            if _latt.Bravais_type == BRAVAIS:
                if _latt.fit <= bestfit: bestfit = _latt.fit
                selected2.append(_latt)
        nsel = len(selected2)
        if nsel == 0:
            print_o("  No lattice could be selected for Space Group %s\n" % \
                                                SPG[1],f3)
            stop()
        for _latt in selected2:
            if _latt.fit == bestfit: selectedLattice = _latt
        
        selectedLattice.idealize()
        
        xpar.UNIT_CELL_CONSTANTS = selectedLattice.cell
        xpar.SPACE_GROUP_NUMBER = spgn
        if auto.anomalous == "yes" or auto.anomalous == "on":
            xpar.FRIEDELS_LAW = "FALSE"
            ano = "ON"
        else:
            xpar.FRIEDELS_LAW = "TRUE"
            ano = "OFF"
        
	dataCol.lookup_image_ranges()
	xpar.DATA_RANGE = dataCol.get_range(minf=auto.frame_first,\
                                            maxf=auto.frame_last)
        
        print_o2("Selected cell axes",fmt_cell % selectedLattice.cell[:3],f3)
        print_o2("Selected cell angles",fmt_cell % selectedLattice.cell[3:],f3)
        print_o2("Selected Bravais lattice",BRAVAIS,f3)
        print_o2("Runing integration with Space Group", SPG[1],f3)
        print_o2("","(number %d)" % spgn,f3)
        print_o2("Using frames number","%d to %d" % tuple(xpar.DATA_RANGE),f3)
        print_o2("Anomalous set",ano,f3)
        
        xpar.mix(getProfilRefPar("INTEGRATE.LP_III_last"))
        print_o2("Using refined "+fmt_shape, "%.3f / %.3f" %\
                             (xpar.BEAM_DIVERGENCE, xpar.REFLECTING_RANGE),f3)
        
        xpar.JOB = "IDXREF"
        xdsout = "xds_v_1.out"
        fhtml.write(fmt_link_xds % xdsout)
        run_xds(xpar, inp_f=None, out_f=xdsout, save=0)
        saveLastVersion(("XDS.INP","XPARM.XDS","IDXREF.LP",),"_V_1")
        
        xdsout = "xds_v_2.out"
        xpar.JOB = "DEFPIX", "XPLAN", "INTEGRATE", "CORRECT"
        fhtml.write(fmt_link_xds % xdsout)
        run_multi_integrate(xpar,inp_f=None,nThreads=_nThreadsMax,init=1)
        
        #saveLastVersion(("XDS.INP","DEFPIX.LP","INTEGRATE.LP","CORRECT.LP"),"_V_2")
        
        stop()
        
        ## De XSCALE : afficher les correlation 2-2,
        ## De CORRECT: Afficher les Rmeas, Compl,...
        ## De CORRECT: Afficher les extinctions systematiques.
        ## De INTEGRATE: Afficher les profiles individuels et moyens.

        #print_o2("Final refined  "+fmt_shape, "%.3f / %.3f" % \
        #                     (xpar.BEAM_DIVERGENCE, xpar.REFLECTING_RANGE),f3)
