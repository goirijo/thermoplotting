import numpy as np
import pandas as pd

headers=["phi", "N0", "N1", "mu0", "mu1", "T"]

numsamples=10

for mu1 in np.arange(0.0,1.0,0.2):
    for mu2 in np.arange(2.0,3.0,0.2):
        T=np.arange(0.0,float(numsamples))

        randdata=np.random.rand(numsamples,len(headers))
        nameddat=pd.DataFrame(randdata, columns=headers)

        nameddat["mu0"]=mu1
        nameddat["mu1"]=mu2
        nameddat["T"]=T
        nameddat["N0"]+=50
        nameddat["N1"]+=70


        filename=format(mu1,".2f")+"-"+format(mu2,".2f")+"-data.txt"
        f=open(filename,'w')
        f.write("#")
        f.write(nameddat.to_string(index=False))
        print nameddat.to_string(index=False)
        f.close()
