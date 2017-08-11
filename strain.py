import numpy as np
import scipy.linalg

#https://en.wikipedia.org/wiki/Finite_strain_theory#Seth-Hill_family_of_generalized_strain_tensors

def _F(startmat,endmat):
    """Calculate the deformation tensor
    to go from start to end

    :startmat: ndarray
    :endmat: ndarray
    :returns: ndarray

    """
    F=np.dot(endmat,np.linalg.inv(startmat))
    return F

def _C(startmat,endmat):
    """Calculate right Cauchy-Green deformation
    tensor to go from start to end

    :startmat: ndarray
    :endmat: ndarray
    :returns: ndarray

    """
    F=_F(startmat,endmat)
    C=np.dot(F.T,F)
    return C

def hencky(startmat,endmat):
    """Calculate the Hencky strain to go from
    start to end. Expects column vectors (use misc.latmat)

    :startmat: ndarray
    :endmat: ndarray
    :returns: ndarray

    """
    C=_C(startmat,endmat)
    return 0.5*scipy.linalg.logm(C)

def _matsplit(mat):
    """Grab the 6 metrics out of a matrix
    """

    xx=mat[0,0]
    yy=mat[1,1]
    zz=mat[2,2]
    yz=mat[1,2]
    xz=mat[0,2]
    xy=mat[0,1]

    return xx,yy,zz,yz,xz,xy

def _e1(xx,yy,zz):
    return (xx+yy+zz)/np.sqrt(3)

def _e2(xx,yy):
    return (xx-yy)/np.sqrt(2)

def _e3(xx,yy,zz):
    return (2*zz-xx-yy)/np.sqrt(6)

def _e4(yz):
    return np.sqrt(2)*yz

def _e5(xz):
    return np.sqrt(2)*xz

def _e6(xy):
    return np.sqrt(2)*xy

def parameters(mat):
    """Calculate the six strain order parameters
    for the given strain metric

    :mat: ndarray
    :returns: float,float,float,float,float,float

    """

    xx,yy,zz,yz,xz,xy=_matsplit(mat)

    e1=_e1(xx,yy,zz)
    e2=_e2(xx,yy)
    e3=_e3(xx,yy,zz)
    e4=_e4(yz)
    e5=_e5(xz)
    e6=_e6(xy)

    return e1,e2,e3,e4,e5,e6

def min_scores(scores):
    """Given a tuple of scores for different types
    of structures, return an array of indices to the
    best (lowest) scores.

    :scores: tuple of array of float. A list of scores for
    each structure you try to map to
    :returns: list of int

    """
    scoremat=np.vstack(scores).T
    bestscore=np.array([np.argmin(row) for row in scoremat])
    return bestscore

def scel_sizes(namelist):
    """Run through the list of SCEL names and return an
    array of the same size that notes the size of each
    supercell with an integer

    :namelist: ndarray(dtype=str)
    :returns: ndarray(dtype=int)

    """
    getsize = lambda s : s[4:s.index("_")]

    sizes=np.array([getsize(name) for name in namelist])
    return sizes.astype(int)
    
