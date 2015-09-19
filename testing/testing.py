import numpy as np
import xray
import glob
import re
import thermoplotting as tp

import matplotlib.pyplot as plt
from scipy import integrate


datanames=glob.glob("./dataset/heating/*.txt")

headerdict={"formation_energy":"U",
            "Ni":"N0",
            "Al":"N1",
            "generalized_enthalpy":"phi",
            "temperature":"T",
            "beta":"b",
            "mu_Ni":"mu0",
            "mu_Al":"mu1"
            }

testdata=tp.thermoarray.HyperArray(datanames, ["mu0","mu1","T"], headerdict)
print testdata

testdata.reverse("T")


N1data=testdata.data_view("N1")
phidata=testdata.data_view("phi")
betadata=testdata.data_view("b")
tempdata=testdata.data_view("T")
phiref=phidata[0,:,:]

PHIdata=tp.grandcanonical.integrate.beta(betadata, phidata, phiref, 0)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.scatter(tempdata[:,0,0], phidata[:,0,0])
plt.show()
