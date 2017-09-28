import numpy as np
import xray
import glob
import re
import thermoplotting as tp

import matplotlib.pyplot as plt
from scipy import integrate
from collections import OrderedDict

headerdict=OrderedDict([("formation_energy","U"),
            ("Ni","N0"),
            ("Al","N1"),
            ("generalized_enthalpy","phi"),
            ("temperature","T"),
            ("beta","b"),
            ("mu_Ni","mu0"),
            ("mu_Al","mu1")
            ])

controlledvar=["mu0","mu1","T"]

#Integrate heating run from low T to high T

heatingnames=glob.glob("./dataset/heating_nuke_0/mu-*/tabulated_averages.txt")
heatingdata=tp.ThermoArray(heatingnames, ["mu0","mu1","T"], headerdict)


heatingphidata=heatingdata.data_view("phi")
heatingbetadata=heatingdata.data_view("b")
heatingPHIref=heatingphidata[0,:,:]

heatingPHIdata=tp.grandcanonical.integrate.beta(heatingbetadata, heatingphidata, heatingPHIref, 0)
heatingdata=heatingdata.push_back(heatingPHIdata,"omega")


fig=plt.figure()
ax=fig.add_subplot(311)
y1=heatingPHIdata[:,0,0]
ax.scatter(heatingdata.data_view("T")[:,0,0], y1)
ax=fig.add_subplot(312)
y2=heatingdata.data_view("omega")[:,0,0]
ax.scatter(heatingdata.data_view("T")[:,0,0],y2)
ax=fig.add_subplot(313)
y3=y2-y1
ax.scatter(heatingdata.data_view("T")[:,0,0],y3)
plt.show()
