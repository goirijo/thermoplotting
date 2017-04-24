from thermoplotting.ternary import normal
from thermoplotting.misc import *
import itertools
from mpl_toolkits.mplot3d.art3d import Poly3DCollection 
from mpl_toolkits.mplot3d.art3d import Line3DCollection 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from scipy.spatial import Voronoi, ConvexHull
import warnings

def to_cartesian(vecmat,points):
    """Given lattice vectors and points in fractional coordinates,
    convert the points to Cartesian

    :vecmat: np array (lattice vectors as columns)
    :points: np array (vertically stacked coordinates)
    :returns: np array

    """
    return np.dot(vecmat,points.T).T

def to_fractional(vecmat,points):
    """Given lattice vectors and points in Cartesian coordinates,
    convert the points to fractional

    :vecmat: np array (lattice vectors as columns)
    :points: np array (vertically stacked coordinates)
    :returns: np array

    """
    return np.dot(np.linalg.inv(vecmat),points.T).T

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

def polygon_facet_center(points):
    """Given a set of points that define a polygon,
    find the center of the polygon

    :points: np array
    :returns: np array

    """
    center=np.average(points,axis=0)
    return center

def polygon_edge_centers(points):
    """Given a set of points that define a polygon,
    find the centers of the edges.

    :points: np array
    :returns: np array

    """
    rolled=np.roll(points,1,axis=0)
    centers=(rolled+points)/2
    return centers


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

    if len(goodregions)!=1:
        warnings.warn("Could not isolate a single Voronoi cell! Results may be wonky.")


    return vorpoints[goodregions[-1]]

def wigner_seitz_facets(latmat):
    """Returns a list of polygons corresponding to the Weigner Seitz cell
    :returns: Poly3DCollection

    """
    vorpoints=wigner_seitz_points(latmat)

    ch=ConvexHull(vorpoints)
    binned=simplex_bin(ch)

    polygons=[polygonal_sort(ch.points[b]) for b in binned]

    return polygons


def draw_voronoi_cell(vectormat,ax,alpha):
    """Plot the Voronoi cell using the given lattice

    :vectormat: Either the real or reciprocal lattice
    :ax: matplotlib subplot
    :returns: ax

    """
    norms=np.linalg.norm(vectormat,axis=0)
    maxrange=np.amax(norms)

    polygons=wigner_seitz_facets(vectormat)
    ax.add_collection(Poly3DCollection(polygons,facecolors='w',linewidth=2,alpha=alpha,zorder=0))
    ax.add_collection(Line3DCollection(polygons,colors='k',linewidth=0.8, linestyles=':'))

    ax.set_xlim([-maxrange,maxrange])
    ax.set_ylim([-maxrange,maxrange])
    ax.set_zlim([-maxrange,maxrange])

    return ax

def voronoi_facet_centers(vectormat, fractional=True):
    """Calculate the centers of facets of either the brillouin zone,
    or Wigner Seitz cell, depending on the given vectormat

    :vectormat: Either the real or reciprocal lattice
    :fractional: bool
    :returns: np array

    """
    polygons=wigner_seitz_facets(vectormat)
    centers=np.stack([polygon_facet_center(p) for p in polygons])

    if fractional:
        centers=to_fractional(vectormat,centers)

    return centers

def voronoi_edge_centers(vectormat, fractional=True):
    """Calculate the centers of the edges of either the brillouin zone,
    or Wigner Seitz cell, depending on the given vectormat

    :vectormat: Either the real or reciprocal lattice
    :fractional: bool
    :returns: np array

    """
    polygons=wigner_seitz_facets(vectormat)
    for p in polygons:
        print polygon_edge_centers(p)
    centers=np.concatenate([polygon_edge_centers(p) for p in polygons],axis=0)


    if fractional:
        centers=to_fractional(vectormat,centers)

    return np.vstack({tuple(row) for row in centers})

def voronoi_vertexes(vectormat, fractional=True):
    """Get the coordinates of the corners/vertexes of the brillouin zone

    :vectormat: Either the real or reciprocal lattice
    :fractional: bool
    :returns: np array

    """
    polygons=wigner_seitz_facets(vectormat)
    points=np.concatenate(polygons,axis=0)

    return np.vstack({tuple(row) for row in points})


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

    def real_to_cartesian(self, points):
        """Convert a list of fractional coordinates into Cartesian
        for the real lattice

        :points: np array (vertically stacked coordinates)
        :returns: np array

        """
        return to_cartesian(self._latmat,points)

    def real_to_fractional(self, points):
        """Convert a list of Cartesian coordinates into fractional
        for the real lattice

        :points: np array (vertically stacked coordinates)
        :returns: np array

        """
        return to_fractional(self._latmat,points)

    def reciprocal_to_cartesian(self, points):
        """Convert a list of fractional coordinates into Cartesian
        for the reciprocal lattice

        :points: np array (vertically stacked coordinates)
        :returns: np array

        """
        return to_cartesian(self._recipmat,points)

    def reciprocal_to_fractional(self, points):
        """Convert a list of Cartesian coordinates into fractional
        for the reciprocal lattice

        :points: np array (vertically stacked coordinates)
        :returns: np array

        """
        return to_fractional(self._recipmat,points)

    def draw_wigner_seitz_cell(self,ax,alpha=1):
        """Plot the Wigner Seitz cell of the lattice
        (Voronoi of real lattice)

        :ax: matplotlib subplot
        :returns: ax

        """
        return self._draw_voronoi_cell(self._latmat,ax,alpha)

    def draw_brillouin_zone(self,ax,alpha=1):
        """Plot the first Brillouin zone in reciprocal space
        (Voronoi of reciprocal lattice)

        :ax: matplotlib subplot
        :returns: ax

        """
        return draw_voronoi_cell(self._recipmat,ax,alpha)

    def brillouin_facet_centers(self,fractional=True):
        """Calculate the center of all facets of the brillouin zone

        :returns: np array
        """
        return voronoi_facet_centers(self._recipmat,fractional)

    def brillouin_edge_centers(self,fractional=True):
        """Calculate the center of all facets of the brillouin zone

        :returns: np array
        """
        return voronoi_edge_centers(self._recipmat,fractional)

    def brillouin_vertexes(self,fractional=True):
        """Get the coordinates of the vertexes of the brillouin zone

        :returns: np array
        """
        return voronoi_vertexes(self._recipmat,fractional)

    def draw_real_vectors(self, ax):
        """Draw the real lattice vectors

        :ax: matplotlib subplot
        :returns: ax

        """
        for v,color in zip(self._latmat.T,['r','g','b']):
            arr=Arrow3D([0,v[0]],[0,v[1]],[0,v[2]],lw=3,arrowstyle="-|>",mutation_scale=20,color=color,linestyle="-")
            ax.add_artist(arr)

        return ax

    def draw_reciprocal_vectors(self, ax):
        """Draw the reciprocal lattice vectors

        :ax: matplotlib subplot
        :returns: ax

        """
        for v,color in zip(self._recipmat.T,['r','g','b']):
            arr=Arrow3D([0,v[0]],[0,v[1]],[0,v[2]],lw=3,arrowstyle="-|>",mutation_scale=20,color=color,linestyle="--")
            ax.add_artist(arr)

        return ax

    def angles(self, rad=True, reciprocal=False):
        """Return the value of alpha, beta and gamma, i.e. the angles
        between the lattice vectors.
        :returns: (float,float,float)

        """
        if not reciprocal:
            a,b,c=self._latmat.T
        else:
            a,b,c=self._recipmat.T

        alpha=angle_between(b,c)
        beta=angle_between(c,a)
        gamma=angle_between(a,b)

        if not rad:
            alpha=alpha*180/np.pi
            beta=beta*180/np.pi
            gamma=gamma*180/np.pi

        return alpha,beta,gamma

    def lengths(self,reciprocal=False):
        """Return the length of each lattice vector
        :returns: TODO

        """
        if not reciprocal:
            a,b,c=self._latmat.T
        else:
            a,b,c=self._recipmat.T

        al=np.linalg.norm(a)
        bl=np.linalg.norm(b)
        cl=np.linalg.norm(c)

        return al,bl,cl


if __name__ == "__main__":
    main()
