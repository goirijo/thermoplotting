import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import thermoplotting as tp

def main():
    #Make up some data
    size=20
    data=np.random.rand(size,3)
    x,y,z=data.T
    #Rescale data so that compositions can't add to more than 1.0, and make energy negative
    #so that there's actually something worth plotting
    x=0.5*x
    y=0.5*y
    z=z-0.5

    #Create plotting object
    plotter=tp.Energy3()

    #Add the random data
    plotter.add_data(x,y,z)
    #Add edge points as data (Unlikely the random data will have gotten the extrema of the compositions)
    plotter.add_data([1,0,0],[0,1,0],[0,0,0])

    #Plot the convex hull in 3d
    fig1=plt.figure(1)
    ax=fig1.add_subplot(111,projection='3d')
    plotter.draw_convex_hull(ax)
    ax.set_zlim([-0.5,0.1])

    #Plot the convex hull projected onto 2d
    fig2=plt.figure(2)
    ax=fig2.add_subplot(111)
    plotter.draw_projected_convex_hull(ax)

    #Show plots
    plt.show()

if __name__=="__main__":
    main()
