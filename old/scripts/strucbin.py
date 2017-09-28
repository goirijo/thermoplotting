import numpy as np
import thermoplotting as tp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#Use structure score and lattice score to bin structures as FCC, BCC or HCP.
#Uses struclat_scores.txt, generated through
#casm query -k 'struc_score(PRIM,basis_score)' 'struc_score(PRIM,lattice_score)' 'struc_score(../NiAl-B2/PRIM-true.vasp,basis_score)' 'struc_score(../NiAl-B2/PRIM-true.vasp,lattice_score)' 'struc_score(../NiAl-HCP/PRIM,basis_score)' 'struc_score(../NiAl-HCP/PRIM,lattice_score)' -c selections/enumerated.json -o info/strain/struclat_scores.txt
#Weighs the scores equally, and takes the lowest one to decide the type of structure.

#Load everything as a string
stringdump=np.genfromtxt("./deformations.txt", dtype="S25")
names=stringdump[:,0]

#Recast relevant columns as float
datadump=np.array(stringdump[:,2::],dtype=float)

FCCscore=np.array(datadump[:,0]+datadump[:,1])
BCCscore=np.array(datadump[:,2]+datadump[:,3])
HCPscore=np.array(datadump[:,4]+datadump[:,5])

xAl=datadump[:,13]
xNi=datadump[:,12]

bestscore=tp.strain.min_scores((FCCscore,BCCscore,HCPscore))

sizes=tp.strain.scel_sizes(names)
unique=np.unique(sizes)

fccsizes=sizes[bestscore==0]
fccxNi=xNi[bestscore==0]

fig=plt.figure()
ax=fig.add_subplot(111)

numbins=11

#ax.set_xlim((0.0,1.0))
#ax.set_ylim((0,460))

enumcolor="royalblue"
stablecolor="red"

histkwargs={"log":False,"bins":numbins,"align":"mid","rwidth":0.7}
counts,bins,patches=ax.hist(xNi,color=enumcolor,**histkwargs)
counts1,bins1,patches1=ax.hist(fccxNi,color=stablecolor,**histkwargs)

for p,p1 in zip(patches,patches1):
    p.set_edgecolor('w')
    p1.set_edgecolor('w')
    p.set_linewidth(0)
    p1.set_linewidth(0)

stabilityratio=[]   #Keep track of how many configs relax at each histogram block
histlocations=[]
for count,count1,x in zip(counts,counts1,bins1):
    note=r"$\frac{"+str(int(count1))+"}{"+str(int(count))+"}$"
    ax.annotate(note,xy=(x+0.5/numbins,count+3),xycoords=("data","data"),ha="center",size=22)
    stabilityratio.append(count1/count)
    histlocations.append(x+0.5/numbins)

stabilitydata=np.vstack((histlocations,stabilityratio)).T

enum_patch = mpatches.Patch(color=enumcolor, label=r'\textbf{Enumerated}')
stable_patch = mpatches.Patch(color=stablecolor, label=r'\textbf{Stable}')
ax.legend(handles=[enum_patch,stable_patch],prop={'size':16},loc="upper left")

ax.set_xlabel(tp.plot.texmrm("x_{Ni}"))
ax.set_ylabel(tp.plot.texbf("Number of configs"))

ax2=ax.twinx()
ax2.plot(stabilitydata[:,0],1-stabilitydata[:,1],'--o',c='k',label="Fraction")
ax2.set_ylim(0,1)
ax2.set_ylabel(tp.plot.texbf("Unstable fraction"))

ax.set_xlim(0,1)
ax2.set_xlim(0,1)
fig.tight_layout()
plt.show()
