#!/usr/bin/env python
"""
28/01/02 pierre.legrand@crchul.ulaval.ca
"""
usage   = """
>>>   Usage : 2pck.py FORMAT image1 [image2 ...]\n
        FORMAT = TIFF (marccd)
"""
import  string,sys,os


def find_template(image_name):
    dirname = os.path.split(os.getcwd())[1]
    if image_name[-3:] == ".gz" or image_name[-3:] == ".GZ":
        image_name = image_name[:-3]
    if image_name[-2:] == ".Z" or image_name[-2:] == ".z":
        image_name = image_name[:-2]
    templ = string.join(image_name.split(".")[:-1])
    ext = str(image_name.split(".")[-1])
    n = 0
    while templ[-1*n-1].isdigit():n+=1
    return templ[:-1*n], n, ext
    
#print os.path.split(os.getcwd())
if len(sys.argv) >= 1:
    format = sys.argv[1]
else:
    print usage
    sys.exit()

nx = ny = 2048
templ, n, ext = find_template(sys.argv[2])
templ_in  = templ + n*"?" + "." + ext
templ_out = templ + n*"?" + ".pck"
print ">> input:           %s" % templ_in
print ">> output:          %s" % templ_out

i1, i2 = 1e10, -1
for img in sys.argv[2:]:
    num = img[img.index(templ)+len(templ):img.rindex(ext)-1]
    num = string.atoi(num)
    i1, i2 = min(i1,num), max(i2,num)
print ">> first image: %7d\n>> last  image: %7d" % (i1, i2)

script = """
 NAME_TEMPLATE_OF_DATA_FRAMES= %s DIRECT %s
 DATA_RANGE=   %d %d
 NX= %d NY= %d
 NAME_TEMPLATE_OF_OUTPUT_FRAMES= %s
""" % (templ_in, format,i1,i2,nx,ny,templ_out)

open("2PCK.INP","wb").write(script)
#os.system("2pck")
if os.path.exists("dataframe"):
   os.remove("dataframe")
