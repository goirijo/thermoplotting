from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import zip
from builtins import range
from builtins import object

import numpy as np

from . import crystal

def bijk_to_coord(prim,bijk,sd=(True,True,True)):
    """Convert the bijk index into an atomic coordinate,
    given the primitive structure

    Parameters
    ----------
    prim : Crystal
    bijk : (int,int,int,int)

    Returns
    -------
    SelectiveAtomCoord

    """
    cart=np.dot(prim.lattice().column_lattice(),bijk[1::])+prim.cart()[bijk[0]]
    name=prim.basis()[0].name()
    return crystal.SelectiveAtomCoord(cart[0],cart[1],cart[2],name,sd)

def stamp_site(jumbo,site,tol):
    """
    Find the closest matching site in the superstructure, if the dot
    product between the distance is within the tolerance, replace
    the site with the new one.

    Parameters
    ----------
    jumbo : Crystal
    site : SelectiveAtomCoord
    tol : float

    Returns
    -------
    Crystal

    """
    ix=crystal.argsort_periodic_coord_match(jumbo.basis(),[site],jumbo.lattice())[0]
    dot=crystal.shortest_periodic_coord_distance(site,jumbo.basis()[ix],jumbo.lattice())

    if dot>tol:
        raise ValueError("You're coordinates don't match. What do.")

    jumbo._basis[ix]=site
    return jumbo

def stamp_bijk(jumbo,prim,bijk,name,trans=(0,0,0),tol=0.00001):
    """
    Find the closest matching site in the superstructure, if the dot
    product between the distance is within the tolerance, replace
    the site with the new one. The site is specified in terms of
    primitive vectors.

    Parameters
    ----------
    jumbo : Crystal
    prim : Crystal
    bijk : [int,int,int,int]
    names : str
    tol : float, optional

    Returns
    -------
    Crystal

    """
    bijk=[i for i in bijk]
    bijk[1]+=trans[0]
    bijk[2]+=trans[1]
    bijk[3]+=trans[2]
    stamp=bijk_to_coord(prim,bijk)
    stamp._name=name
    stamp_site(jumbo,stamp,0.00001)

    return jumbo
