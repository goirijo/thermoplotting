import thermoplotting as tp
import numpy

import matplotlib.tri as tri
import matplotlib.pyplot as plt

data_list=tp.thermoio.collect()
clobbered=tp.thermoio.clobber()

faces=tp.ternary.pruned_hull_facets(clobbered)

clobbered=tp.ternary.equil_trans(clobbered)
faces=tp.ternary.equil_trans(faces)

x=faces[:,:,0].ravel()
y=faces[:,:,1].ravel()

triangles=tri.Triangulation(x,y)

fig=plt.figure()
plt.triplot(triangles)
ax = fig.add_subplot()

plt.show()
