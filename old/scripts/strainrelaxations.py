import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

#Reads deformations.txt, generated via
#`casm query -k 'struc_score(PRIM,basis_score)' 'struc_score(PRIM,lattice_score)' 'struc_score(../NiAl-B2/PRIM-true.vasp,basis_score)' 'struc_score(../NiAl-B2/PRIM-true.vasp,lattice_score)' 'struc_score(../NiAl-HCP/PRIM,basis_score)' 'struc_score(../NiAl-HCP/PRIM,lattice_score)' relaxation_strain comp atom_frac -c selections/enumerated.json -o deformations.txt`


#Load everything as a string
stringdump=np.genfromtxt("./deformations.txt", dtype="S25")
names=stringdump[:,0]

#Recast data as float
datadump=np.array(stringdump[:,2:-1],dtype=float)

#take views of relevant columns
FCCscore=datadump[:,0]+datadump[:,1]
BCCscore=datadump[:,2]+datadump[:,3]
HCPscore=datadump[:,4]+datadump[:,5]

E11=datadump[:,6]
E22=datadump[:,7]
E33=datadump[:,8]


#Calculate strain order parameters
eta3=-(E11+E22-2*E33)/math.sqrt(6)
eta2=(E11-E22)/math.sqrt(2)

etadump=np.vstack((names,eta2,eta3)).T
np.savetxt("etadump.txt",etadump,fmt="%s")

etamat=np.array([eta2,eta3]).T

#Bin structures by type of PRIM
scoremat=np.vstack((FCCscore,BCCscore,HCPscore)).T
bestscore=np.array([np.argmin(row) for row in scoremat])  #0=FCC. 1=BCC, 2=HCP

#Write the scores out for each structure
strucscores=np.vstack((names,bestscore.astype(str))).T
np.savetxt("./strucscores.txt", strucscores,fmt="%s")

colordict={0:'r',1:'g',2:'b'}
colorscore=[]
for score in bestscore:
    colorscore.append(colordict[score])


#Apply symmetry to strain order parameters
mirrortrans=np.array([[-1,0],[0,1]])
rottrans=np.array([[-0.5, -0.5*math.sqrt(3)],[0.5*math.sqrt(3),-0.5]])

etarotmat1=np.dot(etamat,rottrans)
etarotmat2=np.dot(etarotmat1,rottrans)
etarotmat=np.vstack((etamat,etarotmat1,etarotmat2))
etamirrormat=np.dot(etarotmat,mirrortrans)
redundantcolorscore=np.array(3*colorscore)

#Specify maximum radius in eta space to get configurations from
#maxrad=0.085
#radmat=np.sqrt(eta2*eta2+eta3*eta3)
#centeredidx=[radmat < maxrad]

#plot that shit yo
fig=plt.figure()

ax=fig.add_subplot('111')
#plot symmetrized metrics
ax.scatter(etamirrormat[:,0],etamirrormat[:,1],color=redundantcolorscore,edgecolors='gray',s=50,alpha=0.5)
ax.scatter(etarotmat[:,0],etarotmat[:,1],color=redundantcolorscore,edgecolors='gray',s=50,alpha=0.5)
#plot original data
ax.scatter(eta2,eta3,s=50,edgecolors='k',color=colorscore)

ax.set_aspect('equal')
ax.set_xlim([-0.3,0.3])
ax.set_ylim([-0.3,0.3])

ax.set_xlabel(r"$\mathrm{e_2}$")
ax.set_ylabel(r"$\mathrm{e_3}$")

ax.set_title(r"\textbf{FCC strain relaxations}")
red_patch = mpatches.Patch(color='red', label=r'\textbf{FCC}')
green_patch = mpatches.Patch(color='green', label=r'\textbf{BCC}')
blue_patch = mpatches.Patch(color='blue', label=r'\textbf{HCP}')
ax.legend(handles=[red_patch,green_patch,blue_patch],prop={'size':12},loc="upper left")

plt.tight_layout()
plt.show()


#Save the configurations that are FCC
FCCindx=[(bestscore==0)]
FCCnames=names[FCCindx]
FCCfiltered=np.array(FCCscore[FCCindx],dtype="S25")

FCCdump=np.vstack((FCCnames,FCCfiltered)).T
np.savetxt("FCC_scores.txt",FCCdump,fmt="%s")


#Save the configurations that are BCC
BCCindx=[(bestscore==1)]
BCCnames=names[BCCindx]
BCCfiltered=np.array(BCCscore[BCCindx],dtype="S25")

BCCdump=np.vstack((BCCnames,BCCfiltered)).T
np.savetxt("BCC_scores.txt",BCCdump,fmt="%s")


#Save the configurations that are HCP
HCPindx=[(bestscore==2)]
HCPnames=names[HCPindx]
HCPfiltered=np.array(HCPscore[HCPindx],dtype="S25")

HCPdump=np.vstack((HCPnames,HCPfiltered)).T
np.savetxt("HCP_scores.txt",HCPdump,fmt="%s")


