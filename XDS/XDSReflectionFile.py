
import re
import sys
import time

class XDSReflectionFile:
    
    """
    >>> I=XDSReflectionFile("/home/legrand/proj/benchs/xds/BENCH_xds/xds_dir.2007-10-02-10h55/INTEGRATE.HKL")
    >>> I.header["OUTPUT_FILE"]
    'INTEGRATE.HKL'
    >>> I.getData()
    >>> I.calcCorrection()
    
    """
    def __init__(self, fileName):
        
        self.fileName = fileName
        self.raw = open(fileName).read()
        
        headEndMarkup = "!END_OF_HEADER\n"
        dataEndMarkup = "!END_OF_DATA\n"
        indexHeadEnd = self.raw.find(headEndMarkup) + \
                                       len(headEndMarkup)
        self.rawHeader = self.raw[:indexHeadEnd-len(headEndMarkup)]
        self.rawData = self.raw[indexHeadEnd:-len(dataEndMarkup)]
        self.header = self.getHeader()
        self.data = None
        
    def getHeader(self):
        
        _re_xdsPar = r"([^ ]+[=])"
        rec_xdsPar = re.compile(_re_xdsPar)
        
        lpar = []
        for line in self.rawHeader.splitlines():
            l_s = rec_xdsPar.split(line[1:])
            len_s = len(l_s)
            if len_s > 1 and len_s % 2:
                for i in range(1,len_s,2):
                    
                    lpar.append((l_s[i][:-1],l_s[i+1].strip()))
        return dict(lpar)
    
    def getData(self):
        nitem = int(self.header["NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD"])
        t0 = time.time()
        _data = map(float,self.rawData.split())
        if not len(_data) % nitem:
            data = [_data[i:i+nitem] for i in range(0,len(_data),nitem)]
        else:
            print "Error in the interpretation of reflection file"
            sys.exit()
        #print time.time() - t0
        try:
            import Numeric
            t1 = time.time()
            self.data = Numeric.array(data)
        except:
            print "Error the Numeric modul is not installed."
            sys.exit()
        #print time.time() - t1
        #print self.data.shape
    
    def calcCorrection(self):
        """Correction for diamond absorption, to apply on the raw
           integrated intensities.
           Estimated density of diamond: 3.52 g/cm**3
        """
        if self.header["OUTPUT_FILE"] != "INTEGRATE.HKL":
            print "Error. Wrong type of input file."
            print "       It must be an INTEGRATE.HKL type of file!"
        if not self.data: self.getData()
        I = self.data[:,3]
        X = self.data[:,12]*float(self.header["QX"])
        Y = self.data[:,13]*float(self.header["QY"])
        D = float(self.header["DETECTOR_DISTANCE"])
        import pprint
        pprint.pprint(self.header)

        