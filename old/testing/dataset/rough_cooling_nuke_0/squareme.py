import numpy
import os
import glob
import re

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def driveavg(drivedir):
    dirs=glob.glob(drivedir+"/debugout/conditions.*")
    dirs.sort(key=natural_keys)

    valuelist=[]
    for conddir in dirs:

        formation=numpy.loadtxt(conddir+"/generalized_enthalpy.txt")
        formation2=formation*formation
        formavg=numpy.average(formation)
        form2avg=numpy.average(formation2)

        valuelist.append(form2avg)

    return numpy.array(valuelist)


mudirs=glob.glob("mu-*.*")

for muval in mudirs:
    print muval
    filename=muval+"/debugout/genthalpy2_averages.txt"
    avgvals=driveavg(muval)
    numpy.savetxt(filename, avgvals)

