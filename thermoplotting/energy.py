from __future__ import absolute_import
from __future__ import division
from builtins import object

from . import ternary
from . import plot
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


def _digest_data(x, y, z):
    """Given 3D data, shear the x and y values so that the
    space is an equilateral triangle

    :x: list of composition
    :y: list of composition
    :z: list of energy
    :returns: np array

    """
    digestable = np.array((np.ravel(x), np.ravel(y), np.ravel(z))).T
    digested = ternary.equil_trans(digestable)

    return digested


def scatter3(ax, x, y, z, *args, **kwargs):
    """Scatter the given data onto the given matplotlib object,
    applying a shear transformation to it so that it lies on an
    equilateral triangle.

    :ax: matplotlib axis
    :x: list of composition
    :y: list of composition
    :z: list of energy
    :*args: stuff to pass to matplotlib
    :**kwargs: more stuff to pass to matplotlib
    :returns: matplotlib axis

    """
    digested = _digest_data(x, y, z)

    #Scatter things
    return ax.scatter(digested[:, 0], digested[:, 1], digested[:, 2], *args, **kwargs)

    return ax


def projected_scatter(ax, x, y, *args, **kwargs):
    """Scatter the energy data onto the given matplotlib object,
    applying a shear transformation to it so that it lies on an
    equilateral triangle.

    :ax: matplotlib axis
    :x: list of composition
    :y: list of composition
    :*args: stuff to pass to matplotlib
    :**kwargs: more stuff to pass to matplotlib
    :returns: matplotlib axis

    """
    #z values don't matter, since we're projecting
    digested = _digest_data(x, y, y)

    ax.set_xlim([-0.1, 1.1])
    ax.set_ylim([-0.1, 1.1])
    ax.set_aspect('equal')

    #Scatter things
    return ax.scatter(digested[:, 0], digested[:, 1], *args, **kwargs)

    return ax


def draw_convex_hull3(ax, x, y, z):
    """Using all the given data, plot a 3D convex hull
    the convex hull

    :ax: matplotlib axis
    :x: list of composition
    :y: list of composition
    :z: list of energy
    :returns: matplotlib axis

    """
    digested = _digest_data(x, y, z)
    facets = ternary.pruned_hull_facets(digested)

    ax.add_collection3d(Poly3DCollection(facets, facecolors='w', linewidths=2))
    ax.add_collection3d(
        Line3DCollection(facets, colors='k', linewidths=0.2, linestyles=':'))

    return ax


def _draw_projected_facet(ax, facet, kwargs):
    """Draw a single facet onto the composition plane

    :ax: mpl axis
    :facet: from scipy.Hull
    :kwargs: Polygon arguments
    :returns: ax

    """
    coords = facet[:, [0, 1]]
    tri = plt.Polygon(coords, **kwargs)
    ax.add_patch(tri)
    return ax

def _scatter_projected_facet(ax, facet, kwargs):
    """Scatter the corners of a facet onto the composition plane

    Parameters
    ----------
    :ax: mpl axis
    :facet: from scipy.Hull
    :kwargs: Polygon arguments
    :returns: ax


    Returns
    -------
    mpl axis

    """
    coords = facet[:, [0, 1]]
    ax.scatter(facet[:,0],facet[:,1],**kwargs)
    return ax

def draw_projected_convex_hull(ax,
                               x,
                               y,
                               z,
                               kwargs={
                                   "edgecolor": "black",
                                   "linewidth": 1.5,
                                   "fill": False
                               }):
    """Draw the convex hull, but suppress the energy axis, yielding projected facets
    onto the composition plane

    :ax: matplotlib axis
    :kwargs: keyword arguments for the add_patch function
    :returns: matplotlib axis

    """
    digested = _digest_data(x, y, z)
    facets = ternary.pruned_hull_facets(digested)

    for f in facets:
        ax = _draw_projected_facet(ax, f, kwargs)

    ax.set_xlim([-0.1, 1.1])
    ax.set_ylim([-0.1, 1.1])
    ax.set_aspect('equal')

    return ax

def scatter_projected_convex_hull(ax, x, y, z, kwargs):
    """Draw the points that make up the convex hull, but plot
    them projected onto the composition plane   

    Parameters
    ----------
    ax : mpl axis
    x : composition data
    y : composition data
    z : composition data
    :kwargs: keyword arguments for the add_patch function

    Returns
    -------
    mpl axis

    """
    digested = _digest_data(x, y, z)
    facets = ternary.pruned_hull_facets(digested)

    for f in facets:
        ax = _scatter_projected_facet(ax, f, kwargs)

    #Set projected view? axis, aspect, etc?
    return ax


class Energy3(object):
    """Handles plotting ternary energy, with
    a convex hull"""

    def __init__(self, xlabel="comp(a)", ylabel="comb(b)",
                 zlabel="Energy [eV]"):
        """Object construction, does pretty much nothing.

        :xlabel: Label for x-axis
        :ylabel: Label for y-axis
        :zlabel: Label for z-axis
        :returns: TODO

        """
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._zlabel = zlabel

        self._running_data = pd.DataFrame(
            columns=[self._xlabel, self._ylabel, self._zlabel])

    def add_data(self, x, y, z):
        """Apply shear transformation to data to create equilateral triangle and
        save it so that a hull can be constructed later.

        :x: list of composition
        :y: list of composition
        :z: list of energy
        :returns: sheared data (np array)

        """
        #Save the data
        concatable = pd.DataFrame({
            self._xlabel: np.ravel(x),
            self._ylabel: np.ravel(y),
            self._zlabel: np.ravel(z)
        })
        self._running_data = pd.concat((self._running_data, concatable))

        #shearing happens at plot time

        return self._running_data

    def scatter(self, ax, *args, **kwargs):
        """Scatter the saved data onto the given matplotlib object,
        applying a shear transformation to it so that it lies on an
        equilateral triangle.

        :ax: matplotlib axis
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        ax = scatter3(ax, self._running_data[self._xlabel],
                      self._running_data[self._ylabel],
                      self._running_data[self._zlabel])
        return ax

    def projected_scatter(self, ax, *args, **kwargs):
        """Scatter the stored composition data onto a plane,
        ignoring the energy data.

        :ax: matplotlib axis
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        ax = projected_scatter(ax, self._running_data[self._xlabel],
                               self._running_data[self._ylabel],
                               self._running_data[self._zlabel])
        return ax

    def draw_convex_hull(self, ax):
        """Using all the data that has been scattered so far, draw
        the convex hull

        :ax: matplotlib axis
        :returns: matplotlib axis

        """
        ax = draw_convex_hull3(ax, self._running_data[self._xlabel],
                               self._running_data[self._ylabel],
                               self._running_data[self._zlabel])
        return ax

    def draw_projected_convex_hull(
            self,
            ax,
            kwargs={"edgecolor": "black",
                    "linewidth": 1.5,
                    "fill": False}):
        """Draw the convex hull, but suppress the energy axis, yielding projected facets
        onto the composition plane

        :ax: matplotlib axis
        :kwargs: keyword arguments for the add_patch function
        :returns: matplotlib axis

        """
        ax = draw_projected_convex_hull(ax, self._running_data[self._xlabel],
                                        self._running_data[self._ylabel],
                                        self._running_data[self._zlabel],
                                        kwargs)

        ax.set_aspect('equal')

        return ax

    def scatter_projected_convex_hull(self, ax, kwargs={}):
        """Scatter the points of the convex hull, but project all points onto the composition
        space.

        Parameters
        ----------
        :ax: matplotlib axis
        :kwargs: keyword arguments for the add_patch function
        :returns: matplotlib axis

        Returns
        -------
        mpl axis

        """
        ax = scatter_projected_convex_hull(ax, self._running_data[self._xlabel],
                                        self._running_data[self._ylabel],
                                        self._running_data[self._zlabel],
                                        kwargs)

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

        self._running_data = pd.DataFrame(columns=[self._xlabel, self._ylabel])

    def add_data(self, x, y):
        """Add a set of data you want to plot onto a figure

        :x: composition
        :y: energy
        :returns: np array

        """

        concatable = pd.DataFrame({
            self._xlabel: np.ravel(x),
            self._ylabel: np.ravel(y)
        })

        self._running_data = pd.concat((self._running_data, concatable))

        return self._running_data

    def scatter(self, ax, *args, **kwargs):
        """Scatter the energy data onto the given matplotlib object.
        Saves the data so the hull can be constructed later.

        :ax: matplotlib axix
        :label: str
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        x = self._running_data[self._xlabel]
        y = self._running_data[self._ylabel]
        ax.scatter(x, y, *args, **kwargs)

        return ax

    def draw_convex_hull(self, ax, *args, **kwargs):
        """Use the data of all the scattered points to draw the
        convex hull

        :ax: matplotlib axis
        :label: str
        :*args: stuff to pass to matplotlib
        :**kwargs: more stuff to pass to matplotlib
        :returns: matplotlib axis

        """
        hullable = self._running_data.as_matrix()
        hull = ConvexHull(hullable)

        for simplex in hull.simplices:
            ax.plot(hullable[simplex, 0], hullable[simplex, 1], *args, **kwargs)
            # ax.scatter(hullable.ix[simplex][self._xlabel],hullable.ix[simplex][self._ylabel],s=20,c='k')

        return ax
