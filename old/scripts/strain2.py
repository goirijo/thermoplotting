import thermoplotting as tp
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.ticker as ticker
from scipy import ndimage
from scipy import interpolate
from matplotlib import cm
from matplotlib.colors import LogNorm


datadump=np.loadtxt("./datadump.txt",dtype=str)

e1e2E=np.vstack((datadump[:,3].astype(float),datadump[:,4].astype(float),1000*datadump[:,1].astype(float))).T
names=datadump[:,0]

x=e1e2E[:,0]
y=e1e2E[:,1]
z=e1e2E[:,2]

origind=np.where((x==0)*(y==0))
minEind=np.argmin(z)

z-=z[minEind]

e2rad=0.35
e3rad=0.35

xmesh,ymesh=np.meshgrid(np.arange(-e2rad,e2rad,0.01),np.arange(-e3rad,e3rad,0.01))
zmesh=interpolate.griddata((x,y),np.squeeze(z),(xmesh,ymesh))
zmesh=ndimage.gaussian_filter(zmesh,sigma=1.3,order=0)
    
print names[origind]
print e1e2E[origind]
print names[minEind]
print e1e2E[minEind]

fig=plt.figure()

cmap=cm.get_cmap('terrain')

ax=fig.add_subplot(111)

ax.set_xlim([-e2rad,e2rad])
ax.set_ylim([-e3rad,e3rad])

ax.set_xlabel(r"$\mathrm{e_2}$")
ax.set_ylabel(r"$\mathrm{e_3}$")
ax.set_title(r"\textbf{FCC centered Ni}")

#ax.tick_params(axis="both",which="major",labelsize=18)


levels=[0,40,300,700,1500]
label=[10,102,200,500,1000]

ax.set_aspect("equal")

CSnolabel=ax.contourf(xmesh,ymesh,zmesh,levels=sorted((levels+label)),cmap=cmap)
CB=fig.colorbar(CSnolabel,label=r"\textbf{Energy/atom [eV]}")

ax.contour(xmesh,ymesh,zmesh,colors='k',levels=levels,linewidths=3)

fmt=ticker.FormatStrFormatter(r"$\mathbf{%s}$")
CS=ax.contour(xmesh,ymesh,zmesh,colors='k',levels=label,linewidths=3)
ax.clabel(CS,label,inline=1,manual=False,fontsize=13,fmt=fmt)

ax.scatter(x[origind],y[origind],color='b',edgecolor='w',s=50,zorder=2)
ax.scatter(0.0, -0.282976151503,color='r',edgecolor='w',s=50,zorder=2) #BCC variant
ax.scatter(0.24506454,  0.14148808,color='r',edgecolor='w',s=50,zorder=2) #BCC variant
ax.scatter(-0.24506454,  0.14148808,color='r',edgecolor='w',s=50,zorder=2) #BCC variant

#rottrans=np.array([[-0.5, -0.5*math.sqrt(3)],[0.5*math.sqrt(3),-0.5]])


plt.tight_layout()
plt.show()
