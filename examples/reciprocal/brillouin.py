import thermoplotting as tp
import numpy as np
import matplotlib.pyplot as plt

def main():
    a,b,c=tp.latmat("./POSCAR").T

    lat=tp.xtals.Lattice(a,b,c)

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')

    lat.draw_brillouin_zone(ax,alpha=0.5)
    lat.draw_real_vectors(ax)
    lat.draw_reciprocal_vectors(ax)
    ax.scatter(0,0,0)

    plt.show()


if __name__ == "__main__":
    main()
