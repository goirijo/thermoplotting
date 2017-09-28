import thermoplotting as tp
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from scipy import interpolate
from matplotlib import cm



datadump=np.loadtxt("./datadump.txt",dtype=str)

e1e2E=np.vstack((datadump[:,3].astype(float),datadump[:,4].astype(float),1000*datadump[:,1].astype(float))).T
names=datadump[:,0]

x=e1e2E[:,0]
y=e1e2E[:,1]
z=e1e2E[:,2]

origind=np.where((x==0)*(y==0))
minEind=np.argmin(z)

z-=z[minEind]

xmesh,ymesh=np.meshgrid(np.arange(-0.3,0.3,0.001),np.arange(-0.3,0.3,0.001))
zmesh=interpolate.griddata((x,y),np.squeeze(z),(xmesh,ymesh))

print names[origind]
print e1e2E[origind]
print names[minEind]
print e1e2E[minEind]

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(xmesh,ymesh,zmesh,cmap=cm.nipy_spectral,lw=0, vmin=0, vmax=800)
ax.scatter(x,y,z)
ax.scatter(x[origind],y[origind],z[origind],color='g')
ax.scatter(x[minEind],y[minEind],z[minEind],color='r')

ax.set_xlabel(r"$\mathrm{\eta_2}$")
ax.set_ylabel(r"$\mathrm{\eta_3}$")
ax.set_zlabel(r"$\Delta$\textbf{E [meV]}")
ax.set_title(r"\textbf{FCC centered}")


plt.show()
