import numpy as np

def bottom_view(ndarray, axis):
    """Return 0th view along specified axis
    E.g. ndarray[:,:,0,:] if axis=2

    :ndarray: 
    :axis:
    :returns: ndarray

    """
    n=ndarray.ndim
    pretuple=[slice(None)]*ndarray.ndim
    pretuple[axis]=0
    return ndarray[tuple(pretuple)]

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
