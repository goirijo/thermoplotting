import thermoplotting as tp
import numpy

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

#tp.thermoio.casm_energy3_to_thin("energy")

clobbered=tp.thermoio.clobber()

faces=tp.ternary.pruned_hull_facets(clobbered)

x=tp.ternary.composition(clobbered,0)
y=tp.ternary.composition(clobbered,1)
z=tp.ternary.energy(clobbered)/2

clobbered=tp.ternary.equil_trans(clobbered)
faces=tp.ternary.equil_trans(faces)

x=tp.ternary.composition(clobbered,0)
y=tp.ternary.composition(clobbered,1)
z=tp.ternary.energy(clobbered)

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim(-0.7,0)

#ax.scatter(x,y,z, alpha=0)
ax.add_collection3d(Poly3DCollection(faces, facecolors='w', linewidths=1, alpha=1.0))
ax.scatter(x,y,z)
ax.add_collection3d(Line3DCollection(faces, colors='k', linewidths=0.2, linestyles=':'))

plt.show()
