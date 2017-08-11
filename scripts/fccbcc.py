import matplotlib.pyplot as plt
import numpy as np
import re
import thermoplotting as tp


def E_X(energyfile):
    """Return array with just a column for energy
    and a column for composition based on an energy
    file from casm for a binary system"""

    dump=np.genfromtxt(energyfile,dtype="str")
    names=dump[:,-1]
    E=dump[:,0].astype(float)
    x=dump[:,2].astype(float)

    return names,np.array((E,x)).T

def E_XX(energyfile):
    """Return array with just a column for energy
    and a column for composition based on an energy
    file from casm for a binary system"""

    dump=np.genfromtxt(energyfile,dtype="str")
    names=dump[:,-1]
    E=dump[:,0].astype(float)
    x1=dump[:,2].astype(float)
    x2=dump[:,3].astype(float)

    return names,np.array((E,x1,x2)).T

def flatten_composition(EXX):
    """Convert the ternary composition of B2 into
    a binary by taking out the vacancy composition
    degree of freedom.

    :EXX: ndarray (E, xa, xb)
    :returns: ndarray (E,x)

    """
    E=EXX[:,0]
    a=EXX[:,1]
    b=EXX[:,2]

    xNi=1+a-b
    xAl=1-a
    xVa=b

    x=xNi/(xNi+xAl)
    E=E/(xNi+xAl)

    return np.array((E,x)).T

def config_match(string,names):
    """Find the places where the names match a
    certain pattern

    :string: TODO
    :names: TODO
    :returns: TODO

    """
    r=re.compile(string)
    vmatch=np.vectorize(lambda x:bool(r.match(x)))
    sel=vmatch(names)

    return sel

def relax_bin(scores, xE, goodscore):
    """Split the composition-energy array into two
    bins: stable and unstable configurations.
    Indices of scores and xE must refer to the same
    configuration.

    :scores: tuple of structure scores
    :xE: ndarray with composition and energy (or whatever you want binned)
    :goodscore: int, which of the scores you want to bin by
    :returns: ndarray,ndarray (good,bad)
    
    """
    bestscore=tp.strain.min_scores(scores)
    good=bestscore==goodscore
    bad=bestscore!=goodscore

    return xE[good],xE[bad]

#Read energy files for supercells (deformations.txt does not include _custom)
fccnames,fcc=E_X("./energy-FCC")
bccnames,bcc=E_XX("./energy-BCC")
bcc=flatten_composition(bcc)

_,hull=E_X("./truehull")

sel108=config_match("SCEL108",fccnames)
sel96=config_match("SCEL96",fccnames)

#These are the supercell compositions and energies
fccsuper=np.logical_or(sel96,sel108)
bccsuper=config_match("SCEL64",bccnames)

#For the rest of the configurations load up deformations.txt

############################################################################################
#FCC
#Load everything as a string
stringdump=np.genfromtxt("./deformations-FCC", dtype="S25")
names=stringdump[:,0]

#Recast data as float
datadump=np.array(stringdump[:,2::],dtype=float)

#take views of relevant columns
FCCscore=datadump[:,0]+datadump[:,1]
BCCscore=datadump[:,2]+datadump[:,3]
HCPscore=datadump[:,4]+datadump[:,5]
energy=datadump[:,15]
xNi=datadump[:,12]

fccxE=np.vstack((xNi,energy)).T

fccstable,fccunstable=relax_bin((FCCscore,BCCscore,HCPscore),fccxE,0)

############################################################################################
#BCC
#Load everything as a string
stringdump=np.genfromtxt("./deformations-BCC", dtype="S25")
names=stringdump[:,0]

#Recast data as float
datadump=np.array(stringdump[:,2::],dtype=float)

#take views of relevant columns
BCCscore=datadump[:,0]+datadump[:,1]
FCCscore=datadump[:,2]+datadump[:,3]
HCPscore=datadump[:,4]+datadump[:,5]

#These two add up to 1 (Va not considered)
xAl=datadump[:,16]
xNi=datadump[:,14]

#still need vacancies to know how to normalize energy
xa=datadump[:,12]
xb=datadump[:,13]   #This is the number of vacancies

energy=datadump[:,17]/(2-xb)

bccxE=np.vstack((xNi,energy)).T

bccstable,bccunstable=relax_bin((FCCscore,BCCscore,HCPscore),bccxE,1)


############################################################################################

fig=plt.figure()
ax=fig.add_subplot(111)

ax.set_xlabel(r"$\mathrm{x_{Ni}}$")
ax.set_ylabel(r"\textbf{Formation energy [eV]}")

ax.plot(hull[:,1],hull[:,0],color='k',ls='--',linewidth=4,label=r"\textbf{True Hull}")

#These are FCC after relaxing
ax.scatter(fccstable[:,0],fccstable[:,1],s=70,color='royalblue',edgecolor='k',label=r"\textbf{FCC configurations}")
ax.scatter(bccunstable[:,0],bccunstable[:,1],s=70,color='royalblue',edgecolor='gray')

#These are bCC after relaxing
ax.scatter(bccstable[:,0],bccstable[:,1],s=70,color='red',edgecolor='k',label=r"\textbf{BCC configurations}")
ax.scatter(fccunstable[:,0],fccunstable[:,1],s=70,color='red',edgecolor='gray')

ax.scatter(fcc[:,1][fccsuper],fcc[:,0][fccsuper],s=70,color='yellowgreen',edgecolor='k')
ax.scatter(bcc[:,1][bccsuper],bcc[:,0][bccsuper],s=70,color='yellowgreen',edgecolor='k',label=r"\textbf{Dilute defects}")

ax.scatter(hull[:,1],hull[:,0],facecolor='gray',s=70,edgecolor='k',zorder=5)

ax.set_xlim([0.0,1.0])
ax.set_ylim([-0.75,0.1])

ax.legend(loc="upper center",frameon=True, ncol=2,prop={"size":14},fancybox=True,bbox_to_anchor=(0.5,1.07))

plt.tight_layout()
#plt.savefig("fccbcc.svg",bbox_inches="tight")
plt.show()
