import glob
import re
import numpy as np

gens=glob.glob("*genthalpy*")
avgs=glob.glob("*averages*")

gens.sort()
avgs.sort()

for gen, avg in zip(gens, avgs):
    muval=re.split("_", gen)[0]
    
    averages=np.loadtxt(avg, comments="#")
    genthalpy2=np.loadtxt(gen, comments="#")
    dummy=np.copy(genthalpy2)
    dummy[:]=np.nan

    columns=[]

    mu1=3+2+3
    mu2=mu1+1
    columns.append(averages[:,mu1])
    columns.append(averages[:,mu2])

    comp1=3
    comp2=comp1+1
    columns.append(averages[:,comp1])
    columns.append(averages[:,comp2])

    temp=3+2+1
    columns.append(averages[:,temp])
    
    beta=temp+1
    columns.append(averages[:,beta])

    gcenergy=3+2
    columns.append(averages[:,gcenergy])

    columns.append(genthalpy2)

    lotexp=gcenergy #lol
    columns.append(averages[:,lotexp])

    #free energy
    columns.append(dummy)

    finaldata=np.vstack(columns)
    finaldata=np.transpose(finaldata)
    np.savetxt("./"+muval+".thin", finaldata)
