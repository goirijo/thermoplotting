import numpy as np
import thermoplotting as tp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

#Reads deformations.txt, generated via query.sh

#Load everything as a string
stringdump=np.genfromtxt("./deformations.txt", dtype="S25")
names=stringdump[:,0]

#Recast data as float
datadump=np.array(stringdump[:,2::],dtype=float)

#take views of relevant columns
FCCscore=datadump[:,0]+datadump[:,1]
BCCscore=datadump[:,2]+datadump[:,3]
HCPscore=datadump[:,4]+datadump[:,5]
energy=datadump[:,15]
xNi=datadump[:,12]

#Bin structures by type of PRIM
bestscore=tp.strain.min_scores((FCCscore,BCCscore,HCPscore))
truefcc=bestscore==0
notfcc=bestscore!=0

#plot that shit yo
fig=plt.figure()

ax=fig.add_subplot('111')
#plot symmetrized metrics
ax.set_xlim([0.0,1.0])
ax.set_ylim([-0.75,0.1])

ax.scatter(xNi[truefcc],energy[truefcc],s=70,color='blue',edgecolor='k',label=r"\textbf{Enumerated FCC}")
ax.scatter(xNi[notfcc],energy[notfcc],s=70,color='skyblue',edgecolor='k',label=r"\textbf{Unstable}")

ax.set_xlabel(r"$\mathrm{x_{Ni}}$")
ax.set_ylabel(r"\textbf{Formation energy [eV]}")

ax.legend(loc="lower right",prop={"size":10})

plt.tight_layout()
plt.show()
