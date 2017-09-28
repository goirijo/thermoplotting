import json
import glob
import numpy as np

conditionslist=glob.glob("./conditions.*")

size=8*8*8
numatoms=2

a=np.array([2,0,0])
b=np.array([0,1,1])
o=np.array([1,0,1])

Tmat=np.vstack((a,b,o)).T    #This matrix converts parametric to normal composition (also works for chemical potentials)
d=np.linalg.det(Tmat)


#0->Ni, 1->Va, 2->Al
header="#mu0    mu1    mu2    x0    x1    x2    mua    mub    xa    xb    T    b    u    phi"
print header

for condition in conditionslist:
    datadump=open(condition+"/monte_averages.json").read()
    jsonaverages=json.loads(datadump)


    mu0=jsonaverages["conditions"]["mu"][0][0]
    mu1=jsonaverages["conditions"]["mu"][1][0]
    mu2=jsonaverages["conditions"]["mu"][2][0]

    mu=np.array([mu0,mu1,mu2]).T
    lamb=np.dot(Tmat.T, mu)

    mua=lamb[0]-lamb[-1]
    mub=lamb[1]-lamb[-1]

    x0=jsonaverages["averages"][1][0]/(size*numatoms)
    x1=jsonaverages["averages"][2][0]/(size*numatoms)
    x2=jsonaverages["averages"][3][0]/(size*numatoms)

    n=np.array([x0,x1,x2]).T
    x=np.dot(np.linalg.inv(Tmat),(d*n+o))

    xa=x[0]
    xb=x[1]

    T=jsonaverages["conditions"]["temperature"]
    b=jsonaverages["conditions"]["beta"]

    u=jsonaverages["averages"][0][0]/size
    phi=jsonaverages["averages"][4][0]/size


    row=np.array((mu0,mu1,mu2,x0,x1,x2,mua,mub,xa,xb,T,b,u,phi))
    rowstring="    ".join(map(str,row))
    print rowstring
