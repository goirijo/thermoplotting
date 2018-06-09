from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import object

from .. import xtals
from ..misc import *

def fake_properties_calc_json(poscar_file):
    """Create a fake properties.calc.json file from a poscar.
    Fills all the fields with ideal coordinates, no displacements,
    and no energy.

    Parameters
    ----------
    poscar_file : path to file

    Returns
    -------
    dict

    """
    with open(poscar_file,'r') as posf:
        pos_lines=posf.readlines()

    faker={}
    species_tups=xtals.crystal._vasp5_lines_species(pos_lines)
    types,numbers=zip(*species_tups)
    faker["atom_type"]=types
    faker["atoms_per_type"]=numbers

    if xtals.crystal._vasp5_lines_is_cartesian(pos_lines):
        faker["coord_mode"]="cartesian"
    else:
        faker["coord_mode"]="direct"

    coords,selectives=xtals.crystal._vasp5_lines_raw_coordinates(pos_lines)

    faker["relaxed_basis"]=coords.tolist()
    faker["relaxed_energy"]=0.0

    blanks=(coords*0).tolist()

    faker["relaxed_forces"]=blanks
    faker["relaxed_forces"]=blanks

    abc=xtals.crystal._vasp5_lines_to_abc(pos_lines)

    faker["relaxed_lattice"]=abc.tolist()

    return faker
    

class Properties(object):
    """A simple wrapper around a properties.calc.json dictionary
    with routines to access the members. This is for a *single*
    calculation, not something like a barrier, where there's a
    collection of images."""

    def __init__(self, properties,ignore=[]):
        """Initialize with the json dictionary

        Parameters
        ----------
        properties : dict
        ignore : list of keys that you don't care about


        """
        for p in self._required_keys():
            if p not in properties.keys() and p not in ignore:
                raise ValueError("Missing key '{}' in dictionary!".format(p))
        self._properties = properties
        return

    @staticmethod
    def _required_keys():
        """Returns all the keys that should exist in the
        properties.calc.json files

        Returns
        -------
        list of str

        """
        return [
            "atom_type", "atoms_per_type", "coord_mode", "relaxed_basis",
            "relaxed_forces", "relaxed_lattice"
        ]

    @classmethod
    def from_json(cls, filename, ignore=[]):
        """Initialize instance with a filename

        Parameters
        ----------
        filename : string or path
        ignore : list of keys you don't care about

        Returns
        -------
        Properties

        """
        filedict = json_from_file(filename, ignore)
        return cls(filedict)

    @classmethod
    def fake(cls, posfile):
        """Construct a fake set of properties by assigning
        a value of zero everywhere possible, only paying
        attention to the atoms of the specified poscar
        file

        Parameters
        ----------
        posfile : str or path

        Returns
        -------
        Properties

        """
        faked=fake_properties_calc_json(posfile)
        return cls(faked)

    def to_json(self, filename):
        """Save the properties to a file

        Parameters
        ----------
        filename : string or path

        Returns
        -------
        void

        """
        json_to_file(self._properties, filename)

    def species(self):
        """Returns the names of the atom types as
        a list

        Returns
        -------
        list of str

        """
        return self._properties["atom_type"]

    def num_species(self):
        """Returns list of how many of each atom
        there are

        Returns
        -------
        np.array int

        """
        return np.array(self._properties["atoms_per_type"])

    def total_atoms(self):
        """Returns the total number of atoms in the simulation

        Returns
        -------
        int

        """
        return np.sum(self.num_species())

    def compositions(self):
        """Returns list of atomic compositions

        Returns
        -------
        list of float

        """
        return self.num_species() / self.total_atoms()

    def composition(self, specie):
        """Return the composition of the specified
        specie

        Parameters
        ----------
        specie : str

        Returns
        -------
        float

        """
        return self.compositions()[self.species().index(specie)]
        pass

    def energy(self):
        """Return the calculated vasp energy

        Returns
        -------
        float

        """
        return self._properties["relaxed_energy"]

    def as_dict(self):
        """Return all the data as a dictionary

        Returns
        -------
        dict

        """
        return self._properties


class BarrierProperties(object):
    """Combines a list of Properties into a single object
    so that it can calculate kra values."""

    def _atomic_sanity_throw(self):
        """Ensure that all images have the same number of
        each type of atom
        Returns
        -------
        void

        """
        for p in ["atom_type", "atoms_per_type"]:
            value = self._image_properties[0].as_dict()[p]
            for props in self._image_properties[1::]:
                compare = props.as_dict()[p]
                if value != compare:
                    raise ValueError(
                        "The property '{}' is inconsistent throughout the images!"
                    )
        return

    def __init__(self, image_properties):
        """Initialize with a list of Properties, in the correct
        interpolation order

        Parameters
        ----------
        image_properties : TODO


        """
        self._image_properties = image_properties
        self._atomic_sanity_throw()

    @classmethod
    def from_json(cls, filename, ignore=[]):
        """Initialize from a json file

        Parameters
        ----------
        filename : str or path
        ignore : a list of keys you don't care about

        Returns
        -------
        BarrierProperties

        """
        filedict = json_from_file(filename)
        keys = filedict.keys()
        keys.sort()
        props = [Properties(filedict[k],ignore) for k in keys if k.isdigit()]
        return cls(props)

    def species(self):
        return self._image_properties[0].species()

    def num_species(self):
        return self._image_properties[0].num_species()

    def total_atoms(self):
        return self._image_properties[0].total_atoms()

    def compositions(self):
        return self._image_properties[0].compositions()

    def composition(self,specie):
        return self._image_properties[0].composition(specie)

    def kra(self):
        """Calculate the KRA value. At the moment the only way
        this is done is by assuming that the middle image has
        the highest barrier.

        Returns
        -------
        float

        """
        # if len(self._image_properties) % 2 == 0:
        #     raise ValueError(
        #         "Under the current implementation, an odd number of images is required"
        #         "to calculate the KRA value.")

        midpoint=(self._image_properties[0].energy()+self._image_properties[-1].energy())/2
        energies=np.array([im.energy() for im in self._image_properties])
        maxen=np.max(energies)
        return maxen-midpoint


