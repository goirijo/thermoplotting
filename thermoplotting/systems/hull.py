from __future__ import division
# from past.utils import old_div
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

def pruned_hull_facet_inds(hull,tol=0.00001):
    """Given a ConvexHull, return indexes to the facets
    that make up the bottom of the hull

    Parameters
    ----------
    hull : scipy.spatial ConvexHull

    Returns
    -------
    list of int

    """
    inds=[ix for ix, eq in enumerate(hull.equations) if eq[-2]<0-tol]
    return inds

def project_point_on_facet(point,equation):
    """Given a point and an equation for a facet, project
    the point onto the facet and return the coordinates
    of the point projected onto the facet

    Parameters
    ----------
    point : np array
    equation : np array (has one more element than point)

    Returns
    -------
    np array

    """
    if len(point) + 1 != len(equation):
        raise ValueError(
            "Dimensions of point ({}) not consistent with dimensions of equation ({})".
            format(len(point), len(equation)))

    norm_energy_component = equation[-2]
    offset = equation[-1]

    #The equation comes in form a*x+b*y+c*z+d=0
    #so dot a,b (from equation) with x,y (from point)
    non_energy_dot=np.dot(point[0:-1],equation[0:-2])
    projected_facet_value=(-non_energy_dot-offset)/norm_energy_component

    projected_point=np.copy(point)
    projected_point[-1]=projected_facet_value

    return projected_point

def distance_from_facet(point, equation):
    """Using the given facet equation, determine
    what the projection of the point onto the
    facet is, and return the difference between
    the given point and it's projection, yielding
    the distance to the facet

    Parameters
    ----------
    point : np array
    equation : np array (has one more element than point)

    Returns
    -------
    float

    """
    proj_point=project_point_on_facet(point,equation)
    return point[-1]-proj_point[-1]


def distance_from_hull(point, hull):
    """Calculate the projected distance from the
    pruned convex hull by checking the equations
    of all the facets

    Parameters
    ----------
    point : np array
    hull : scipy.spatial ConvexHull

    Returns
    -------
    float

    """
    bottom_inds=pruned_hull_facet_inds(hull)
    bottom_eqs=[hull.equations[ix] for ix in bottom_inds]
    distances = [distance_from_facet(point, eq) for eq in bottom_eqs]
    return min(d for d in distances)

def distances_from_hull(points,hull):
    """Call distance_from_hull on every point given,
    returning an array of all the distances to the hull

    Parameters
    ----------
    points : np array
    hull : scipy.spatial ConvexHull

    Returns
    -------
    np array

    """
    distances=[distance_from_hull(p,hull) for p in points]
    return np.array(distances)

