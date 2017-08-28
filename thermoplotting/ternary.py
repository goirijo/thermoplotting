import scipy.spatial as spa
import scipy.linalg as lin
import numpy as np

def prepare_input(file_name):       #pylint: disable=unused-argument
    """Read in from list of files and output them as new .thin
    files that are compatible with the expected format:
    composition0, composition1, energy

    :file_name: name of file, presumably casm output
    :returns: void, but you'll end up with new a file "oldfile.thin"

    """
    pass

def hull_facets(data_list):
    """Creates convex hull from your data set and returns
    a list of all the facets that make up the hull. Meant to be used with a
    data list that has two compositions and energy.

    :data_list: numpy array of numpy array
    :returns: array of facets (collection of 3 points)

    """
    # scipy.spatial.Delauny is numerically unstable, we use scipy.sptial.ConvexHull instead.
    new_tri = spa.ConvexHull(data_list) #pylint: disable=no-member
    return new_tri.points[new_tri.simplices]

def normal(three_points):
    """Makes two vectors out of three points,
    crosses them and returns its normalized value.

    :three_points: numpy array, presumably out of a facet
    :returns: numpy array of size 3

    """

    # get our vector differences w.r.t. the first end-state and roll
    vecs = (three_points[1:] - three_points[0])
    base_list = np.arange(three_points.shape[1])

    crossvec = np.array([(-1)**i * np.linalg.det(vecs[:, base_list[base_list != i]])
                            for i in range(three_points.shape[1])])

    return crossvec/lin.norm(crossvec)

def endstate(data_list, component, tolerance=0.0001):
    """Returns point with lowest energy for the edge of
    the phase diagram

    :data_list: full list of data, presumably from clobber()
    :component: int specifying component (0, 1, 2=origin)
    :returns: numpy array of size 3

    """
    # Grab the max energy for reference
    returnpoint = data_list[data_list[:, -1].argmax()]
    # Mask to entries along our selected axis, with lower energy
    mask = np.logical_and(data_list[:, -1] <= data_list[:, -1].max(),
                             abs(composition(data_list, component) - 1.0) < tolerance)

    # Only move forward if mask has at least one 'true' value
    if np.any(mask):
        returnpoint = data_list[mask][data_list[mask][:, -1].argmin()]
    # Otherwise use our default value
    else:
        pass

    return returnpoint

def endstates_normal(data_list):
    """Return normal vector to the plane spanned by
    the three endstates, pointing in a positive z
    direction

    :data_list: array of numpy array. Presumably from clobber()
    :returns: numpy array

    """
    # Use a list comprehension to grab the end states for every composition axes
    endstates = np.array([endstate(data_list, i) for i in range(data_list.shape[1])])
    statenormal = normal(endstates)

    if statenormal[-1] < 0.:
        statenormal = -statenormal

    return statenormal

def truncated_data(data_list):
    """Eliminates all data points with energy above the plane
    that the end states make up. Useful for creating an initial
    convex hull that doesn't include the top.

    :data_list: full list of data, presumably from clobber()
    :returns: truncated version of input

    """
    # # Grab the normal vector to the endstates
    statenormal = endstates_normal(data_list)
    # # Re-orient by the first end state
    refendstate = endstate(data_list, 0)

    # Return only entries that are below the plane defined by the 3 end-states
    return data_list[np.dot(data_list - refendstate, statenormal) <= 0]

def sliced_facets(facet_list, mynorm, refstate, tolerace=0.0001):
    """Eliminates any facets that contain points above
    above the specified plane

    :facet_list: full convex hull facets
    :normal: normal vector to plane to slice through
    :refstate: reference point for normal
    :returns: subset of the given facets

    """
    # Find only the faces where all 3 points are below the reference normal
    mask = np.all((np.dot((facet_list - refstate), mynorm) <= tolerace), axis=1)
    # Return the masked array
    return facet_list[mask]

def pruned_facets(facet_list, normalvec, tolerace=0.0001):
    """Goes through a list of facets and removes any facet
    whose three points are coplanar with a normal vector
    parallel to the given normalvec

    :facet_list: list of facets in hull
    :normalvec: reference normal vector to know which facets to remove
    :returns: truncated version of input

    """

    normalvec = normalvec/lin.norm(normalvec)

    returnlist = []

    for facet in facet_list:
        facetnorm = normal(facet)
        if abs(np.dot(normalvec, facetnorm)) < (1.0-tolerace):
            returnlist.append(facet)

    return np.array(returnlist)

def composition(data_list, component):
    """Returns one composition column from the given data, which
    was presumably made using the standard clobber() function.

    :data_list: double numpy array
    :component: int specifying which composition to return (2=origin)
    :returns: list of first composition values

    """
    if component == data_list.shape[1]-1:
        origincomp = 1.0 - data_list[:, :-1].sum(axis=1)
        return origincomp

    else:
        return data_list[:, component]

def energy(data_list):
    """Returns energy column from the given data, which
    was presumably made using the standard clobber() function

    :data_list: double numpy array
    :returns: list of second composition values

    """
    return data_list[:, -1]

def pruned_hull_facets(data_list):
    """Calls facets() to get all facets of the convex hull, but then
    removes all facets on binary subspace, as well as any resulting
    facets above the reference endpoints.

    :data_list: numpy array of numpy array
    :returns: array of facets (collection of 3 points)

    """
    new_tri = spa.ConvexHull(data_list) #pylint: disable=no-member
    good_simplex=[]

    for eq,sx in zip(new_tri.equations,new_tri.simplices):
        if eq[2]<0:
            good_simplex.append(sx)
    good_simplex=np.vstack(good_simplex)
    return new_tri.points[good_simplex]

    # facet_list = hull_facets(data_list)

    # statenormal = endstates_normal(data_list)
    # refendstate = endstate(data_list, 0)
    # facet_list = sliced_facets(facet_list, statenormal, refendstate)


    # norm = np.array([0, 1, 0])
    # facet_list = pruned_facets(facet_list, norm)

    # norm = np.array([1, 1, 0])
    # facet_list = pruned_facets(facet_list, norm)

    # norm = np.array([1, 0, 0])
    # facet_list = pruned_facets(facet_list, norm)

    # #norm=endstates_normal(data_list)
    # #facet_list=pruned_facets(facet_list, norm);

    # return facet_list

def equil_trans(data_list):
    """Recursively goes through all the data points and applies shear matrix
    so that the composition space looks like an equilateral triangle

    :data_list: arbitrarily dimensioned array, can be clobber() or faces()
    :returns: double numpy array with transformed composition values

    """
    transmat = np.array([[1, 0.5, 0], [0, 3**(0.5)/2, 0], [0, 0, 1]])

    # If we have 2 or 1 dimnesions, we can directly dot
    if data_list.ndim < 3:
        return np.dot(transmat, data_list.T).T
    # Otherwise use recursion to dot each entry in the data_list
    else:
        return np.array([equil_trans(dlist) for dlist in data_list])


def point_is_over_facet(point,facet):
    """Given a point, and a facet of a convex hull, determine
    if the point lies above the given facet along the z (energy)
    axis. This is meant to identify which facet of the convex hull
    should be used to calculate the hull distance.

    Use this to understand algorithm:
    http://blackpawn.com/texts/pointinpoly/

    Basically in the 2D space, for the point p, and three vertices
    of the triangle r1, r2 and r3, redefine p with new basis vectors
    in terms of r1, r2 and r3.

    p-r0=u*(r1-r0)+v*(r2-r0)
    v2=u*v0+v*v1

    v2v0=u*v0v0+v*v1v0
    v2v1=u*v0v1+v*v1v1

    Solve for u and v, which tells you if p is within the triangle.
    It's outside if u or v is negative.
    It's outside if u or v is greater than 1.
    It's outside if u+v is greater than 1.

    :point: single coordinate as np
    :facet: three coordinates as np
    :returns: bool

    """
    p=point[0:2]
    r0=facet[0,0:2]
    r1=facet[1,0:2]
    r2=facet[2,0:2]


    v0=r1-r0
    v1=r2-r0
    v2=p-r0

    u=(np.dot(v2,v1)*np.dot(v1,v0)-np.dot(v2,v0)*np.dot(v1,v1))/(np.dot(v0,v1)*np.dot(v1,v0)-np.dot(v0,v0)*np.dot(v1,v1))
    v=(np.dot(v2,v0)*np.dot(v0,v1)-np.dot(v2,v1)*np.dot(v0,v0))/(np.dot(v1,v0)*np.dot(v0,v1)-np.dot(v1,v1)*np.dot(v0,v0))

    #Round to some decimal points to avoid floating noise
    u=round(u,10)
    v=round(v,10)


    if v < 0 or u < 0 or v > 1 or u > 1 or v+u > 1:
        return False
    else:
        return True

def find_facet_under_point(point, pruned_facets):
    """Given a list of facets that make up the convex hull (assumes that
    the top of the hull and the sides were pruned off, not sure if that makes
    a difference), find which facet lies directly below the given point.

    :point: single coordinate as np
    :pruned_facets: list of three coordinates as np
    :returns: three coordinates as np

    """
    for f in pruned_facets:
        if point_is_over_facet(point,f):
            return f

    return None


def distance_from_hull(point, pruned_facets):
    """Determine the distance of a point from the convex hull. The given
    convex hull must have had the top and sides pruned off!

    :point: single coordinate as np
    :pruned_facets: list of three coordinates as np
    :returns: float

    """
    facet=find_facet_under_point(point,pruned_facets)
    assert(point_is_over_facet(point,facet))

    #translate everything so the point is at the origin
    transfacet=facet-point

    #Three point define a plane now
    p1,p2,p3=transfacet

    #Two vectors spawn the plane
    v1=p3-p1
    v2=p2-p1

    #And the normal vector is used to define ax+by+cz=d
    n=np.cross(v1,v2)
    a,b,c=n
    d=np.dot(n,p3)

    #The intercept is defined by z when x and y are 0
    z=d/c

    #The distance from the hull is measured in the other direction
    return -z

    
