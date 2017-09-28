import glob
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#plt.style.use('ggplot')

dirlist=glob.glob("./mu-*.*")

fig=plt.figure()
ax=fig.add_subplot(111)

for run in dirlist:
    rundata=np.loadtxt(run+"/debugout/tabulated_averages.txt")
    x=rundata[:,3]/(rundata[:,3]+rundata[:,4])
    y=rundata[:,6]-273

    if run=="./mu-1.6":
        postx=x
        posty=y

    ax.scatter(x,y,s=70)
        
ax.scatter(x,y,color='red',edgecolor='black', s=70)

ax.set_xlabel(r"\textbf{x$_{\mathrm{Ni}}$}")
ax.set_ylabel(r"\textbf{T [$^{\circ}$C]}")
ax.set_xlim(0.68,1.01)
ax.set_ylim(0.-300,2000)

plt.title(r"\textbf{Cooling runs at constant \mu}")

plt.show()
