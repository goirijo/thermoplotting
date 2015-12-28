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

#Generate references for cooling run integration: low mu to higher mu (high T)

lowmunames=glob.glob("./dataset/lowmu_nuke_0/mu-*/tabulated_averages.txt")
lowmudata=tp.ThermoArray(lowmunames, ["mu0","mu1","T"], headerdict)

lowmuphidata=lowmudata.data_view("phi")
lowmumudata=lowmudata.data_view("mu1")


fig=plt.figure()
ax=fig.add_subplot(111)
#ax.scatter(tempdata[:,0,0], phidata[:,0,0])
ax.scatter(heatingdata.data_view("T")[:,0,0], heatingPHIdata[:,0,0])
plt.show()

#coolingnames=glob.glob("./dataset/cooling_nuke_0/mu-*/tabulated_averages.txt")
#coolingdata=tp.ThermoArray(coolingnames, ["mu0","mu1","T"], headerdict)
