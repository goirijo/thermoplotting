import pandas as pd
import numpy as np
import thermoplotting.ternary
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import plot


class Energy3(object):

    """Handles plotting ternary energy, with
    a convex hull"""

    def __init__(self, xlabel="comp(a)", ylabel="comb(b)", zlabel="Energy [eV]"):
        """Object construction, does pretty much nothing.

        :xlabel: Label for x-axis
        :ylabel: Label for y-axis
        :zlabel: Label for z-axis
        :returns: TODO

        """
        self._xlabel=xlabel
        self._ylabel=ylabel
        self._zlabel=zlabel

        self._running_data=pd.DataFrame(columns=[self._xlabel,self._ylabel,self._zlabel])

    def add_data(self, x, y, z):
        """Apply shear transformation to data to create equilateral triangle and
        save it so that a hull can be constructed later.

        :x: list of composition
        :y: list of composition
        :z: list of energy
        :returns: sheared data (np array)

        """
        #Shear composition
        digestable=np.array((x,y,z)).T
        digested=thermoplotting.ternary.equil_trans(digestable)

        #Save the data
        concatable=pd.DataFrame({self._xlabel:digested[:,0],self._ylabel:digested[:,1],self._zlabel:digested[:,2]})
        self._running_data=pd.concat((self._running_data,concatable))

        return digested

    def scatter(self, ax, x, y, z, *args, **kwargs):
        """Scatter the energy data onto the given matplotlib object,
        applying a shear transformation to it so that it lies on an
        equilateral triangle. Saves the scattered data so that a hull
        can be constructed with it later.

        :ax: matplotlib axis
        :x: list of composition
        :y: list of composition
        :z: list of energy
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        digested=self.add_data(x,y,z)

        #Scatter things
        ax.scatter(digested[:,0],digested[:,1],digested[:,2],*args,**kwargs)

        return ax

    def projected_scatter(self, ax, x, y, z, *args, **kwargs):
        """Scatter the energy data onto the given matplotlib object,
        applying a shear transformation to it so that it lies on an
        equilateral triangle. Saves the scattered data so that a hull
        can be constructed with it later. The z values are not plotted
        here though, only the x and y are scattered onto the
        composition plane

        :ax: matplotlib axis
        :x: list of composition
        :y: list of composition
        :z: list of energy
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        digested=self.add_data(x,y,z)

        #Scatter things
        ax.scatter(digested[:,0],digested[:,1],*args,**kwargs)

        return ax

    def draw_convex_hull(self, ax):
        """Using all the data that has been scattered so far, draw
        the convex hull

        :ax: matplotlib axis
        :returns: matplotlib axis

        """
        digestable=self._running_data.as_matrix(columns=[self._xlabel,self._ylabel,self._zlabel])
        facets=thermoplotting.ternary.pruned_hull_facets(digestable)

        ax.add_collection3d(Poly3DCollection(facets, facecolors='w', linewidths=2))
        ax.add_collection3d(Line3DCollection(facets, colors='k', linewidths=0.2, linestyles=':'))

        return ax

    def _draw_projected_facet(self, ax, facet, kwargs):
        """Draw a single facet onto the composition plane

        :ax: mpl axis
        :facet: from scipy.Hull
        :kwargs: Polygon arguments
        :returns: ax

        """
        coords=facet[:,[0,1]]
        tri=plt.Polygon(coords,**kwargs)
        ax.add_patch(tri)
        return ax

    def draw_projected_convex_hull(self, ax, kwargs={"edgecolor":"black","linewidth":1.5,"fill":False}):
        """Draw the convex hull, but suppress the energy axis, yielding projected facets
        onto the composition plane

        :ax: matplotlib axis
        :kwargs: keyword arguments for the add_patch function
        :returns: matplotlib axis

        """
        digestable=self._running_data.as_matrix(columns=[self._xlabel,self._ylabel,self._zlabel])
        facets=thermoplotting.ternary.pruned_hull_facets(digestable)

        for f in facets:
            ax=self._draw_projected_facet(ax,f,kwargs)

        ax.set_aspect('equal')

        return ax

        
    
class Energy2(object):

    """Handles plotting binary energies, complete with
    a convex hull"""

    def __init__(self, xlabel="x_a", ylabel="Energy [eV]"):
        """Just sets the labels

        :xlabel: str for x-axis
        :ylabel: str for y-axis

        """
        self._xlabel = plot.texmrm(xlabel)
        self._ylabel = plot.texbf(ylabel)
        
        self._running_data={}

    def add_data(self, label, x, y):
        """Add a set of data you want to plot onto a figure

        :label: str
        :x: composition
        :y: energy
        :returns: self

        """
        concatable=pd.DataFrame({self._xlabel:x,self._ylabel:y})
        self._running_data[label]=concatable

        return self
        

    def scatter(self, ax, label, *args, **kwargs):
        """Scatter the energy data onto the given matplotlib object.
        Saves the data so the hull can be constructed later.

        :ax: matplotlib axix
        :label: str
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        x=self._running_data[label][self._xlabel]
        y=self._running_data[label][self._ylabel]
        ax.scatter(x,y,*args,**kwargs)

        return ax

    def draw_convex_hull(self, ax, label, *args, **kwargs):
        """Use the data of all the scattered points to draw the
        convex hull

        :ax: matplotlib axis
        :label: str
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        hullable=self._running_data[label]
        hull=ConvexHull(hullable)

        for simplex in hull.simplices:
            ax.plot(hullable.ix[simplex][self._xlabel],hullable.ix[simplex][self._ylabel],*args,**kwargs)
            # ax.scatter(hullable.ix[simplex][self._xlabel],hullable.ix[simplex][self._ylabel],s=20,c='k')

        return ax
