import thermoplotting as tp
import glob
import numpy as np
import matplotlib.pyplot as plt

#Specify controlled variables (names of input columns)
controlled=["mua","mub","T"]

#Make list of all averages files
filenames=glob.glob("./datasets/ground0/mu_*.*_*.*/tabulated_averages.txt")

#Create data array (This one has three dimensions: mua, mub and T)
heatdata=tp.ThermoArray(filenames,controlled)
#Identify what each axis corresponds to
muaax=heatdata.axis("mua")
muabx=heatdata.axis("mub")
Tax=heatdata.axis("T")

#To integrate along constant mu we need phi as reference, beta values and mu values
mua=heatdata.data_view("mua")
mub=heatdata.data_view("mub")
b=heatdata.data_view("b")
T=heatdata.data_view("T")

#DO NOT TRUST REPORTED phi VALUES
xa=heatdata.data_view("xa")
xb=heatdata.data_view("xb")
u=heatdata.data_view("u")
phi=u-mua*xa-mub*xb

#The array comes sorted lowest to highest, so take bottom temperature as reference
phiref=tp.bottom_view(phi,Tax)

grandcanon=tp.gc.integrate.beta(b,phi,phiref,Tax)

#Make sure the integration worked
muaind=20
mubind=20
Tind=slice(None)

viewind=[muaind,mubind,Tind]

plt.scatter(T[viewind],grandcanon[viewind])
plt.show()
