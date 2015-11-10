import numpy as np

def facet_intercept(facet):
    """Determine to intercept of the plane defined
    by a facet with the axis that define the composition
    space.

    :facet: ndarray
    :returns: array

    """
    pass

def indexed_view(indx, ndarray, axis):
    """Return last view along specified axis
    E.g. ndarray[:,:,indx,:] if axis=2

    :indx:
    :ndarray: 
    :axis:
    :returns: ndarray

    """
    n=ndarray.ndim
    pretuple=[slice(None)]*ndarray.ndim
    pretuple[axis]=indx
    return ndarray[tuple(pretuple)]

def top_view(ndarray, axis):
    """Return last view along specified axis
    E.g. ndarray[:,:,-1,:] if axis=2

    :ndarray: 
    :axis:
    :returns: ndarray

    """
    return indexed_view(-1,ndarray,axis)

def bottom_view(ndarray, axis):
    """Return 0th view along specified axis
    E.g. ndarray[:,:,0,:] if axis=2

    :ndarray: 
    :axis:
    :returns: ndarray

    """
    return indexed_view(0,ndarray,axis)

def reverse(ndarray, axis, inplace=False):
    """Invert order of N-dimensional array
    along the given axis
    
    :ndarray:
    :axis:
    :inplace: bool
    :returns: ndarray
    """
    resorted=np.swapaxes(ndarray,axis,0)
    resorted=resorted[::-1,...]
    resorted=np.swapaxes(resorted,0,axis)

    #views into ndarray will be updated
    if inplace:
        ndarray[...]=resorted[...]
        return ndarray

    #views into ndarray are not views into resorted
    else:
        return resorted

def minimize(arraylist):
    """Given a list of arrays of matching shapes, return
    a new array of the same shape that has all the minimum
    values from the list.

    :arraylist: array of ndarray
    :returns: ndarray

    """
    mininds=np.argmin(arraylist,axis=0)
    return np.array([arraylist[midx][idx] for idx, midx in np.ndenumerate(mininds)]).reshape(mininds.shape)
