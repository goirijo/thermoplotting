"""A Docstring"""

import scipy.spatial as spa
import scipy.linalg as lin
import numpy

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
    # tri = spa.Delaunay(data_list)     #pylint: disable=no-member
    # simplex = tri.find_simplex
    # facets = []
    # for i, j, k in new_tri.simplices:
    #     facets.append(data_list[[i, j, k]])

    # return facets

    # scipy.spatial.Delauny is numerically unstable, we use scipy.sptial.ConvexHull instead.
    new_tri = spa.ConvexHull(data_list) #pylint: disable=no-member
    # the ConvexHull object contains members indexing into the simplices, which can be used
    #   on the (referenced) set of points from data_list
    return new_tri.points[new_tri.simplices]

def normal(three_points):
    """Makes two vectors out of three points,
    crosses them and returns its normalized value.

    :three_points: numpy array, presumably out of a facet
    :returns: numpy array of size 3

    """

    # get our vector differences w.r.t. the first end-state and roll
    # vecs = numpy.roll(three_points[1:] - three_points[0], -1, 0)
    vecs = (three_points[1:] - three_points[0])
    base_list = numpy.arange(three_points.shape[1])

    crossvec = numpy.array([(-1)**i * numpy.linalg.det(vecs[:, base_list[base_list != i]])
                            for i in range(three_points.shape[1])])


    # vec1 = three_points[0]-three_points[1]
    # vec2 = three_points[0]-three_points[2]

    # crossvec = numpy.cross(vec1, vec2)
    return crossvec/lin.norm(crossvec)

def endstate(data_list, component, tolerance=0.0001):
    """Returns point with lowest energy for the edge of
    the phase diagram

    :data_list: full list of data, presumably from clobber()
    :component: int specifying component (0, 1, 2=origin)
    :returns: numpy array of size 3

    """
    # maxindex=numpy.argmax(energy(data_list))
    # returnpoint=numpy.array([data_list[maxindex]])
    # for point in data_list:
    #     point=numpy.array([point])
    #     comp=composition(point,component)

    #     if(abs(1.-comp)<tolerance and energy(point)<=energy(returnpoint)):
    #         returnpoint=point
    # return returnpoint[0]
    # Grab the max energy for reference
    returnpoint = data_list[data_list[:, -1].argmax()]
    # Mask to entries along our selected axis, with lower energy
    mask = numpy.logical_and(data_list[:, -1] <= data_list[:, -1].max(),
                             abs(composition(data_list, component) - 1.0) < tolerance)

    # Only move forward if mask has at least one 'true' value
    if numpy.any(mask):
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
    # endstates=[]
    # endstates.append(endstate(data_list,0))
    # endstates.append(endstate(data_list,1))
    # endstates.append(endstate(data_list,2))

    # Use a list comprehension to grab the end states for every composition axes
    endstates = numpy.array([endstate(data_list, i) for i in range(data_list.shape[1])])
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
    return data_list[numpy.dot(data_list - refendstate, statenormal) <= 0]

    # belowceiling=[]
    # for point in data_list:
    #     tvec=point-refendstate
    #     if(numpy.dot(tvec,statenormal)<=0):
    #         belowceiling.append(point)

    # return numpy.array(belowceiling)

def sliced_facets(facet_list, mynorm, refstate, tolerace=0.0001):
    """Eliminates any facets that contain points above
    above the specified plane

    :facet_list: full convex hull facets
    :normal: normal vector to plane to slice through
    :refstate: reference point for normal
    :returns: subset of the given facets

    """
    # Find only the faces where all 3 points are below the reference normal
    mask = numpy.all((numpy.dot((facet_list - refstate), mynorm) <= tolerace), axis=1)
    # Return the masked array
    return facet_list[mask]

    # belowplane=[]
    # for facet in facet_list:
    #     keep=True
    #     for point in facet:
    #         tvec=point-refstate
    #         if(numpy.dot(tvec,normal)>0):
    #             keep=False
    #             break
    #     if(keep==True):
    #         belowplane.append(facet)

    # return numpy.array(belowplane)

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
        if abs(numpy.dot(normalvec, facetnorm)) < (1.0-tolerace):
            returnlist.append(facet)

    return numpy.array(returnlist)

def composition(data_list, component):
    """Returns one composition column from the given data, which
    was presumably made using the standard clobber() function.

    :data_list: double numpy array
    :component: int specifying which composition to return (2=origin)
    :returns: list of first composition values

    """
    if component == data_list.shape[1]-1:
    # if component == 2:
        # oneslist = numpy.ones(data_list.shape[0])
        # origincomp = oneslist-data_list[:, 0]-data_list[:, 1]
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
    facet_list = hull_facets(data_list)

    statenormal = endstates_normal(data_list)
    refendstate = endstate(data_list, 0)
    facet_list = sliced_facets(facet_list, statenormal, refendstate)


    norm = numpy.array([0, 1, 0])
    facet_list = pruned_facets(facet_list, norm)

    norm = numpy.array([1, 1, 0])
    facet_list = pruned_facets(facet_list, norm)

    norm = numpy.array([1, 0, 0])
    facet_list = pruned_facets(facet_list, norm)

    #norm=endstates_normal(data_list)
    #facet_list=pruned_facets(facet_list, norm);

    return facet_list

def equil_trans(data_list):
    """Recursively goes through all the data points and applies shear matrix
    so that the composition space looks like an equilateral triangle

    :data_list: arbitrarily dimensioned array, can be clobber() or faces()
    :returns: double numpy array with transformed composition values

    """
    transmat = numpy.array([[1, 0.5, 0], [0, 3**(0.5)/2, 0], [0, 0, 1]])

    # If we have 2 or 1 dimnesions, we can directly dot
    if data_list.ndim < 3:
        return numpy.dot(transmat, data_list.T).T
    # Otherwise use recursion to dot each entry in the data_list
    else:
        return numpy.array([equil_trans(dlist) for dlist in data_list])

    # returndata=[]
    # for point in data_list:
    #     if point.ndim==1:
    #         point=numpy.dot(transmat,point)
    #     else:
    #         point=equil_trans(point)

    #     returndata.append(point)

    # return numpy.array(returndata)
