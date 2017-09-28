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
    x=rundata[:,4]/(rundata[:,3]+rundata[:,4])
    y=rundata[:,6]

    ax.scatter(x,y)
        
ax.set_xlabel(r"\textbf{x$_{\mathrm{Al}}$}")
ax.set_ylabel(r"\textbf{T [K]}")

plt.show()
