import thermoplotting as tp
import numpy as np

def main():
    a=np.array([5,0,0])
    b=np.array([0,5,0])
    c=np.array([0,0,5])

    a=np.array([-5,5,5])
    b=np.array([5,-5,5])
    c=np.array([5,5,-5])

    a=np.array([0,5,5])
    b=np.array([5,0,5])
    c=np.array([5,5,0])

    a=np.array([0,5,7])
    b=np.array([3,0,7])
    c=np.array([3,5,0])

    lat=tp.xtals.Lattice(a,b,c)
    lat.draw_brillouin_zone()



if __name__ == "__main__":
    main()
