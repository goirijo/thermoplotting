import scipy.spatial as spa
import scipy.linalg as lin
import numpy

def faces(data_list):
    """Creates convex hull from your data set and returns
    a list of all the faces that make up the hull. Meant to be used with a
    data list that has two compositions and energy.

    :data_list: numpy array of numpy array
    :returns: array of faces (collection of 3 points)

    """
    tri=spa.Delaunay(data_list)
    simplex=tri.find_simplex
    faces=[]
    for ia, ib, ic in tri.convex_hull:
        faces.append(data_list[[ia,ib,ic]])

    return faces


def normal(three_points):
    """Makes two vectors out of three points,
    crosses them and returns its normalized value.

    :three_points: numpy array, presumably out of a facet
    :returns: numpy array of size 3

    """
    vec1=three_points[0]-three_points[1]
    vec2=three_points[0]-three_points[2]

    crossvec=numpy.cross(vec1,vec2)
    return crossvec/lin.norm(crossvec)


def endstate(data_list, component, tolerace=0.0001):
    """Returns point with lowest energy for the edge of
    the phase diagram

    :data_list: full list of data, presumably from clobber()
    :component: int specifying component (0, 1, 2=origin)
    :returns: numpy array of size 3

    """
    maxindex=numpy.argmax(energy(data_list))
    returnpoint=numpy.array([data_list[maxindex]])
    for point in data_list:
        point=numpy.array([point])
        comp=composition(point,component)

        if(abs(1-comp)<tolerace and energy(point)<energy(returnpoint)):
            returnpoint=point

    return returnpoint[0]


def endstates_normal(data_list):
    """Return normal vector to the plane spanned by
    the three endstates, pointing in a positive z
    direction

    :data_list: array of numpy array. Presumably from clobber()
    :returns: numpy array

    """
    endstates=[]
    endstates.append(endstate(data_list,0))
    endstates.append(endstate(data_list,1))
    endstates.append(endstate(data_list,2))
    statenormal=normal(endstates)

    if statenormal[2]<0:
        statenormal=-statenormal

    return statenormal


def truncated_data(data_list):
    """Eliminates all data points with energy above the plane
    that the end states make up. Useful for creating an initial
    convex hull that doesn't include the top.

    :data_list: full list of data, presumably from clobber()
    :returns: truncated version of input

    """
    statenormal=endstates_normal(data_list)
    refendstate=endstate(data_list,0)

    belowceiling=[]
    for point in data_list:
        tvec=point-refendstate
        if(numpy.dot(tvec,statenormal)<=0):
            belowceiling.append(point)

    return numpy.array(belowceiling)


def pruned_facets(facet_list, normalvec, tolerace=0.0001):
    """Goes through a list of facets and removes any facet
    whose three points are coplanar with a normal vector
    parallel to the given normalvec

    :face_list: list of facets in hull
    :normalvec: reference normal vector to know which facets to remove
    :returns: truncated version of input

    """

    normalvec=normalvec/lin.norm(normalvec)

    returnlist=[]

    for facet in facet_list:
        crossvec=numpy.cross(normal(facet),normalvec)
        if abs(abs(lin.norm(crossvec))-1) > tolerace:
            returnlist.append(facet)

    return numpy.array(returnlist)


def composition(data_list, component):
    """Returns one composition column from the given data, which
    was presumably made using the standard clobber() function.

    :data_list: double numpy array
    :component: int specifying which composition to return (2=origin)
    :returns: list of first composition values

    """
    if component==2:
        oneslist=numpy.ones(data_list.shape[0])
        origincomp=oneslist-data_list[:,0]-data_list[:,1]
        return origincomp

    else:
        return data_list[:,component]


def energy(data_list):
    """Returns energy column from the given data, which
    was presumably made using the standard clobber() function

    :data_list: double numpy array
    :returns: list of second composition values

    """
    return data_list[:,2]
