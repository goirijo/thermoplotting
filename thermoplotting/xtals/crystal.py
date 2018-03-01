from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import Lattice
import numpy as np
import pandas as pd
import itertools
import copy

def bring_fractional_within(frac):
    """Drops the leading integer from all the 

    Parameters
    ----------
    frac : np.array 1x3

    Returns
    -------
    np.array 1x3

    """
    within_frac=frac.copy()
    for ix,elem in enumerate(within_frac):
        within_frac[ix]-=int(elem)
    return within_frac

class AtomCoord(object):

    """Contains the Cartesian coordinates of a crystal site and a label
    associated with the species residing there"""

    def __init__(self,x,y,z,name):
        """Initialize with the x, y and z components of the coordinates as
        well as a string to identify the species

        Parameters
        ----------
        x : float
        y : float
        z : float
        name : string

        """
        self._cartesian=np.array((x,y,z))
        self._name=name

    def cart(self):
        """Get the Cartesian coordinates
        Returns
        -------
        np.array

        """
        return self._cartesian

    def frac(self, lat):
        """Get the fractional coordinates
        for the given lattice

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        np.array 

        """
        return lat.real_to_fractional(self.cart())

    def name(self):
        """Returns the name of the species
        Returns
        -------
        str

        """
        return self._name

    def moved_to(self, x, y, z):
        """Return a copy of self that has been moved to a new location

        Parameters
        ----------
        x : float
        y : float
        z : float

        Returns
        -------
        AtomCoord

        """
        moved=copy.copy(self)
        moved._cartesian=np.array([x,y,z])
        return moved

    def moved_within(self, lat):
        """Return a copy of self that has been moved within the given lattice

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        AtomCoord

        """
        raw_frac=self.frac(lat)
        within_frac=bring_fractional_within(raw_frac)
        within_cart=lat.real_to_cartesian(within_frac)

        return self.moved_to(*within_cart)
        

class SelectiveAtomCoord(AtomCoord):

    """Contains Cartesian coordinates of a crystal site along with 
    the label defining what species it is, in addition to also having
    three bool values defining the selective dynamics relaxation"""

    def __init__(self,x,y,z,name,sd=(True,True,True)):
        """Initialize with the Cartesian coordinates and the species name. Default
        for the selective dynamics is to relax in all directions

        Parameters
        ----------
        x : float
        y : float
        z : float
        name : string
        sd : (bool,bool,bool), optional


        """
        AtomCoord.__init__(self,x,y,z,name)
        self._sd = sd

    def is_selective_dynamics(self):
        """Returns true if any of the selective dynamics flags are true
        Return
        -------
        bool

        """
        return (False in self._sd)
        

class Crystal(object):

    """Holds a Lattice and a basis in the form of a list of AtomCoords or
    SelectiveAtomCoords"""

    def _initialize(self):
        """Sets up the Cartesian and Fractional coordinates
        for the current basis
        Returns
        -------
        void

        """
        self._cart_coords=np.array([c.cart() for c in self._basis])
        self._frac_coords=self._lattice.real_to_fractional(self.cart())
        return

    def _is_selective_dynamics(self):
        """Returns true if any of the selective dynamics flags are true
        Returns
        -------
        bool

        """
        for b in self._basis:
            if b.is_selective_dynamics():
                return True
        return False

    def __init__(self,lat,basis=[],title="TP",scaling=1.0):
        """Initialize with Lattice and list of AtomCoords"""
        if scaling!=1.0:
            raise NotImplemented("Cannot deal with non unity scaling for crystal structures")

        self._lattice=lat
        self._basis=basis

        self._title=title
        self._scaling=scaling

        self._initialize()

    def lattice(self):
        """Returns the lattice of the structures
        Returns
        -------
        Lattice

        """
        return self._lattice

    def cart(self):
        """Return all the Cartesian coordinates of the basis sites
        Returns
        -------
        np.array nx3

        """
        return self._cart_coords

    def frac(self):
        """Return all the fractional coordinates of the basis sites
        Returns
        -------
        np.array nx3

        """
        return self._frac_coords

    def append_to_basis(self, new_site):
        """Add a new site to the current list of basis sites

        Parameters
        ----------
        new_site : SelectiveAtomCoord

        Returns
        -------
        void

        """
        self._basis.append(new_site)
        self._initialize()

    def _species_counts(self):
        """Returns a dict of species names with how many times that
        type appears in the basis
        Returns
        -------
        dict

        """
        sdict={}

        for s in self._basis:
            name=s.name()
            if name in sdict:
                sdict[name]+=1
            else:
                sdict[name]=1

        return sdict

    def to_vasp5(self, filename):
        """Dump the contents in vasp5 format. No fancy options
        here. Will print selective dynamics if any of the flags
        are turned on

        Parameters
        ----------
        filename : string

        Returns
        -------
        void

        """
        sdict=self._species_counts()
        with open(filename, 'w') as posdump:
            posdump.write(self._title+'\n')
            posdump.write(str(self._scaling)+'\n')

            np.savetxt(posdump,self._lattice.row_lattice())

            for s in sdict:
                posdump.write(s+"    ")
            posdump.write('\n')
            for s in sdict:
                posdump.write(str(sdict[s])+"    ")
            posdump.write('\n')

            if self._is_selective_dynamics():
                posdump.write("Selective dynamics\n")

            #Only fractional allowed for the moment
            posdump.write("Direct\n")

            #Print in groups of specie type
            for s in sdict:
                for b,f in zip(self._basis, self.frac()):
                    if b.name()==s:
                        np.savetxt(posdump,np.expand_dims(f,axis=0))
            return

def _vasp5_lines_to_lattice(poscar_lines):
    """Extract the lattice from the line dump of a poscar file

    Parameters
    ----------
    poscar : list of str

    Returns
    -------
    Lattice

    """
    abc=np.array([l.split() for l in poscar_lines[2:5]])
    abc=abc.astype(float)

    return Lattice(*abc)

def _vasp5_lines_species(poscar_lines):
    """Extract the number of each species of the structre
    and return the values as tuples

    Parameters
    ----------
    poscar : list of str

    Returns
    -------
    list of (str, int)

    """
    species=poscar_lines[5].split()
    num_species=poscar_lines[6].split()
    return zip(species,[int(n) for n in num_species])

def _vasp5_lines_is_cartesian(poscar_lines):
    """Returns true if the coordinates have been set for Cartesian mode in
    the poscar files

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    bool

    """
    slot=7
    if _vasp5_lines_is_selective_dynamics(poscar_lines):
        slot+=1
    tag=poscar_lines[slot].split()[0][0]
    if tag in "cCkK":
        return True
    else:
        return False
        
def _vasp5_lines_is_selective_dynamics(poscar_lines):
    """Returns true if selective dynamics has been turned on for
    the poscar files

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    bool

    """
    tag=poscar_lines[7].split()[0][0]
    if tag in "sS":
        return True
    else:
        return False

def _vasp5_lines_raw_coordinates(poscar_lines):
    """Return the coordinates of the sites and the corresponding
    selective dynamics. Agnostic to Direct/Cartesian, this just converts
    the lines into the appropriate type.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    (np.array nx3, list of (bool,bool,bool))

    """
    slot=8
    is_selective_dynamics=_vasp5_lines_is_selective_dynamics(poscar_lines)
    if is_selective_dynamics:
        slot+=1

    numsites=0
    species=_vasp5_lines_species(poscar_lines)
    for s,n in species:
        numsites+=n

    basis_lines=poscar_lines[slot:slot+numsites]

    points=np.array([l.split()[0:3] for l in basis_lines],dtype=float)

    sd=[(True,True,True) for p in points]
    if is_selective_dynamics:
        print("_vasp5_lines_raw_coordinates is untested for selective dynamics")
        sd=[(x=='T',y=='T',z=='T') for x,y,z in l.split()[3:6] for l in basis_lines]

    return points,sd
        
def vasp5_to_lattice(poscar):
    """Extract the lattice from a poscar file

    Parameters
    ----------
    poscar : filename

    Returns
    -------
    Lattice

    """
    with open(poscar,'r') as posf:
        lines=posf.readlines()

    return _vasp5_lines_to_lattice(lines)

def vasp5_to_crystal(poscar):
    """Read a poscar file and create a crystal from it

    Parameters
    ----------
    poscar : filename

    Returns
    -------
    Crystal

    """
    with open(poscar,'r') as posf:
        lines=posf.readlines()

    name=lines[0]
    scaling=lines[1]
    lat=_vasp5_lines_to_lattice(lines)
    species=_vasp5_lines_species(lines)
    is_selective_dynamics=_vasp5_lines_is_selective_dynamics(lines)
    is_cartesian=_vasp5_lines_is_cartesian(lines)

    returner=lambda points : points
    coord_generator=returner
    if(not is_cartesian):
        coord_generator=lat.real_to_cartesian

    raw_coords,sds=_vasp5_lines_raw_coordinates(lines)
    cart_coords=coord_generator(raw_coords)
    coord_species=[name  for name,num in species for i in xrange(num)]

    selective_basis=[SelectiveAtomCoord(*cart,name=name,sd=sd) for cart,name,sd in zip(cart_coords,coord_species,sds)]

    return Crystal(lat,selective_basis)

def blind_cart_match(coord, to_basis):
    """Returns the smallest dot product between the coordinate and all
    of the coordinates in the to basis. Only compares the Cartesian
    coordinates and has no concept of periodic images.

    Parameters
    ----------
    coord : np.array 1x3
    to_basis : np.array mx3

    Returns
    -------
    float

    """
    distances=to_basis-coord
    norms=np.linalg.norm(distances,axis=1)
    return np.amin(norms)

def neighboring_basis_images(accord, lat):
    """Given a coordinate and the corresponding lattice, return
    a list of the same coordinate within the neighboring lattice
    cells (27 total, including the original coordinate)

    Parameters
    ----------
    accord : AtomCoord
    lat : Lattice

    Returns
    -------
    list[AtomCoord] (27 entries)

    """
    #This creates fractional coordinates of the form [-1,-1,-1],[-1,-1,0]...[0,0,0],[0,0,1]...[1,1,0],[1,1,1]
    frac_image_origins=np.array(list(itertools.product([-1,0,1],repeat=3)))
    frac_images=frac_image_origins+accord.frac(lat)

    image_neighbors=[accord.moved_to(*lat.real_to_cartesian(im)) for im in frac_images]
    return image_neighbors


def shortest_periodic_coord_distance(acoord0, acoord1, lat):
    """Given two AtomCoords, get the shortest distance between the two
    basis sites after taking periodic images of each other into
    account

    Parameters
    ----------
    acoord0 : AtomCoord
    acoord1 : AtomCoord
    lat : Lattice

    Returns
    -------
    float

    """
    within0=acoord0.moved_within(lat)
    within1=acoord1.moved_within(lat)

    images0=neighboring_basis_images(within0,lat)
    image0_carts=np.array([im.cart() for im in images0])
    return blind_cart_match(within1.cart(),image0_carts)

def argsort_periodic_coord_match(from_basis, to_basis, lat):
    """Rank the from basis by how well it can map onto the
    to basis taking periodic images into account. Best matches
    are at the top. Ignores species.

    Parameters
    ----------
    from_basis : list of AtomCoord
    to_basis : list of AtomCoord

    Returns
    -------
    list of int

    """
    shortest=[np.amin([shortest_periodic_coord_distance(f,t,lat) for t in to_basis]) for f in from_basis]
    return np.argsort(shortest)
    

