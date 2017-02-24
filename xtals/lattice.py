from thermoplotting.ternary import normal
from thermoplotting.misc import *
import itertools
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
from mpl_toolkits.mplot3d.art3d import Line3DCollection 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from scipy.spatial import Voronoi, ConvexHull

def simplex_bin(hull):
    """Given a convex hull, check the equations of the
    hyperplanes and bin them, returning a set of simplex
    groups with a common equation (i.e. coplanar simplices)

    :hull: convex hull
    :returns: array of array of simplex

    """
    equations=np.vstack({tuple(q) for q in hull.equations})

    binned=[[] for q in equations]

    for q,s in zip(hull.equations,hull.simplices):
        #whichever row is zero has the same equation as the current simplex
        single_zero_row=equations-q
        index=np.where((single_zero_row==0).all(axis=1))[0]

        assert(index.shape==(1,))
        index=index[0]

        binned[index].append(s)

    return [np.unique(a) for a in binned]

def signed_angle_3d(v0,v1,vn):
    """Get the signed angle for two vectors in 3d space.

    :v0: np vector
    :v1: np vector
    :vn: np vector (normal vector)
    :returns: float (rad)

    """
    v0n=v0/np.linalg.norm(v0)
    v1n=v1/np.linalg.norm(v1)

    #Avoid float point pain with 8 decimal places. Close enough
    angle=np.arccos(round(np.dot(v0n,v1n),8))

    cross=np.cross(v0,v1)
    
    if np.dot(vn, cross) < 0:
        angle=-angle

    return angle


def polygonal_sort(points):
    """Given a set of points that define a polygon,
    sort them so
    that they all go around the center in order.

    :points: np array
    :returns: np array

    """
    n=normal(points[0:3])
    c=np.sum(points,axis=0)/len(points)

    ref=points[0]-c

    angles=np.array([signed_angle_3d(c-p,ref,n) for p in points])
    sortinds=np.argsort(angles)

    return points[sortinds]


def reciprocal_lattice(latmat):
    """Cross vectors and multiply by 2 pi to get the
    reciprocal of the given lattice

    :latmat: np 3x3 (vectors in columns)
    :returns: np 3x3

    """
    a,b,c=latmat.T

    vol=np.dot(a,np.cross(b,c))

    astar=2*np.pi*np.cross(b,c)/vol
    bstar=2*np.pi*np.cross(c,a)/vol
    cstar=2*np.pi*np.cross(a,b)/vol

    return np.array([astar,bstar,cstar]).T

def wigner_seitz_points(latmat):
    """Determine the edges of the Wigner Seitz cell, given the lattice.
    Generates just enough lattice points to generate a single WS cell,
    then selects points from the only full region.
    If the reciprocal lattice is given, then the points define the first Brillouin
    zone.

    :latmat: 3x3 vectors as columns
    :returns: np list of points (as rows)

    """
    a,b,c=latmat.T
    #Range of lattice points that will be enough to enclose the Weigner Seitz cell
    radpoints=range(-1,2)
    counterpoints=[(x,y,z) for x in radpoints for y in radpoints for z in radpoints]
    gridpoints=np.array([x*a+y*b+z*c for x,y,z in counterpoints])

    #Construct Voronoi cell
    vor=Voronoi(gridpoints,furthest_site=False)
    vorpoints=vor.vertices
    vorregions=vor.regions

    #Only one full Voronoi cell should have been constructed
    goodregions=[x for x in vorregions if len(x) > 0 and x[0] is not -1]
    assert(len(goodregions)==1)

    return vorpoints[goodregions[0]]

def wigner_seitz_facets(latmat):
    """Returns a list of polygons corresponding to the Weigner Seitz cell
    :returns: Poly3DCollection

    """
    vorpoints=wigner_seitz_points(latmat)

    ch=ConvexHull(vorpoints)
    binned=simplex_bin(ch)

    polygons=[polygonal_sort(ch.points[b]) for b in binned]

    return polygons


class Lattice(object):

    """Simple class to hold the lattice vectors of a lattice, with
    a few routines to do things in reciprocal space"""

    def __init__(self, a, b, c):
        """Define the lattice with three lattice vectors, stored
        vertically in a matrix

        :a: 3x1
        :b: 3x1
        :c: 3x1

        """
        self._latmat=np.array([a,b,c]).T
        self._recipmat=reciprocal_lattice(self._latmat)

    def _draw_voronoi_cell(self,vectormat,ax):
        """Plot the Voronoi cell using the given lattice

        :vectormat: Either the real or reciprocal lattice
        :ax: matplotlib subplot
        :returns: ax

        """
        norms=np.linalg.norm(vectormat,axis=0)
        maxrange=np.amax(norms)

        polygons=wigner_seitz_facets(vectormat)
        ax.add_collection(Poly3DCollection(polygons,facecolors='w',linewidth=2,alpha=1,zorder=0))
        ax.add_collection(Line3DCollection(polygons,colors='k',linewidth=0.8, linestyles=':'))

        ax.set_xlim([-maxrange,maxrange])
        ax.set_ylim([-maxrange,maxrange])
        ax.set_zlim([-maxrange,maxrange])

        return ax

    def draw_wigner_seitz_cell(self,ax):
        """Plot the Wigner Seitz cell of the lattice
        (Voronoi of real lattice)

        :ax: matplotlib subplot
        :returns: ax

        """
        return self._draw_voronoi_cell(self._latmat,ax)

    def draw_brillouin_zone(self,ax):
        """Plot the first Brillouin zone in reciprocal space
        (Voronoi of reciprocal lattice)

        :ax: matplotlib subplot
        :returns: ax

        """
        return self._draw_voronoi_cell(self._recipmat,ax)


if __name__ == "__main__":
    main()
