import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def doublematch(startstr,endstr,namemap):
    """Return subset of name_map.txt data
    such that only configurations that start
    with startstr (e.g. plain) and end with
    namemap (e.g. _AlAl) are available

    :startstr: str
    :endstr: str
    :namemap: ndarray
    :returns: ndarray

    """
    startre=re.compile(startstr+".*")
    vmatch=np.vectorize(lambda x:bool(startre.match(x)))
    startmatch=vmatch(namemap[:,0])

    endre=re.compile(".*"+endstr+"$")
    vmatch=np.vectorize(lambda x:bool(endre.match(x)))
    endmatch=vmatch(namemap[:,0])

    combo=startmatch*endmatch
    return namemap[combo]

def autoplot(startstr,endstr,namemap,ax):
    """Automatically generate a plot for a particular
    set of defects by extracting data with doublematch.
    Sets energies relative to the largest distance, and
    unnormalizes by multiplying the energy by the
    supercell size (108)

    :startstr: str
    :endstr: str
    :namemap: ndarray
    :ax: matplotlib object
    :returns: ndarray

    """
    #get subset of data
    subset=doublematch(startstr,endstr,datadump)

    #unnormalize and make energies relative to largest distance
    distance=subset[:,6].astype(float)
    maxind=np.argmax(distance)
    formation=subset[:,1].astype(float)
    formation=108*(formation-formation[maxind])
    expansion=subset[:,2].astype(float)
    expansion=108*(expansion-expansion[maxind])


    ax.scatter(distance,formation,s=90,edgecolors='k',color="b", label=r"\textbf{DFT energies}")
    ax.scatter(distance,expansion,s=70,edgecolors='k',color="lightgreen", label=r"\textbf{Expanded energies}")

    ax.legend(loc="lower right",prop={"size":15})
    ax.legend(loc="upper right",prop={"size":15})

    ax.set_xlabel(r"\textbf{Distance [\AA]}")
    ax.set_ylabel(r"\textbf{Relative energy [eV]}")


datadump=np.loadtxt("./name_map.txt",dtype=str)

fig=plt.figure()
ax=fig.add_subplot(111)

titledict={"plain":"Ni","ordered":"L_{12}","delta":"\delta"}

phase="plain"
defect="_AlAlAl"

autoplot(phase,defect,datadump,ax)
ax.set_title(r"\textbf{"+defect[1:3]+"-"+defect[3:5]+" defect energies in $\mathrm{"+titledict[phase]+"}$}")

plt.tight_layout()
fig.savefig("./figures/"+phase+defect+".eps")
#plt.show()
