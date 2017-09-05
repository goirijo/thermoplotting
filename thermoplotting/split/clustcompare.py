import json
import pandas as pd
import numpy as np

def species_from_basis(basis):
    """Returns a list of the species that make up the basis of the
    cluster expansion
    :returns: list of str

    """
    species_dump=[]
    for sf in basis["site_functions"]:
        for phi in sf["basis"]:
            for elem in sf["basis"][phi]:
                species_dump.append(elem)

    return set(species_dump)

def compare_cluster(clust0, clust1):
    """Given two clusters, compare the multiplicity
    and prototype to determine whether they describe
    the same set of sites. This should work fine, but
    there's no guarantee from casm that the prototype
    will match across different basis sets.

    :clust0: json
    :clust1: json
    :returns: bool

    """
    is_equal=True

    for prop in ("mult","prototype"):
        if clust0[prop]!=clust1[prop]:
            is_equal=False

    return is_equal

def species_to_basis_dict(basis):
    """Run through the basis functions and return a map that tells
    you which species are associated with each basis function. This
    is meant to be used for occupation basis only.
    E.g.: Ni->phi_0_0

    :basis: json
    :returns: dict

    """
    species_to_basis={}
    for unit in basis["site_functions"]:
        for basis in unit["basis"]:
            for b in unit["basis"][basis]:
                if unit["basis"][basis][b]==1:
                    if b in species_to_basis:
                        species_to_basis[b].append(basis)
                    else:
                        #Strip the \\ from the string so you're just left with phi_x_y
                        species_to_basis[b]=[basis[1::]]

    return species_to_basis


def formula_has_bfunc(formula, bfunc):
    """Checks if bfunc is contained in formula string

    :formula: str
    :bfunc: str
    :returns: bool
    """
    return formula.count(bfunc)>0


def formula_has_multi_bfunc(formula, bfunc):
    """Checks for repetitions of a particular basis function
    (e.g. phi_0_0 representing Cr) in a given cluster function
    formula

    :formula: str
    :bfunc: str
    :returns: bool
    """
    return formula.count(bfunc)>1

def formula_has_any_bfunc(formula,bfunc_set):
    """Checks if any of the provided basis functions appear in
    the formula

    :formula: str
    :bfunc_set: list of str
    :returns: bool

    """
    for b in bfunc_set:
        if formula_has_bfunc(formula,b):
            return True
    return False

def formula_has_multi_any_bfunc(formula, bfunc_set):
    """Returns true if the total times any of the provided
    basis functions appear in the formula is more than unity.

    :formula: str
    :bfunc_set: list of str
    :returns: bool

    """
    instances=0
    for b in bfunc_set:
        instances+=formula.count(b)

    return instances>1
    


class FitBlob(object):

    """Combines a set of correlations, energies and ECI values to give you
    the cluster expanded energies
    
    TODO: Still unimplemented. Might be dead.
    """

    def _partial_check(self):
        """Checks to make sure the eci values aren't greater than the correlations.
        This is not enough to ensure a proper combination of eci and correlations
        was given, but it's something

        :returns: bool

        """
        pass

    def _vectorized_eci(self):
        """Extract the eci values from the json format into a vector of mostly zeros with
        the appropriate values in the right positions

        :returns: pd DataFrame

        """
        zeros=np.zeros(len(self._corr.columns.values))

        for ix,eci in self._ecidump["fit"]["eci"]:
            zeros[ix]=eci

        #This should properly take care of the index-eci relationship
        return pd.DataFrame({"eci":zeros})

    def _pure_corr_data(self):
        """From the DataFrame containing the correlations, extract the subset of columns
        that have the form "corr(x)"

        :returns: pd DataFrame

        """
        cols=self._corrdump.columns.values
        corrcols=[c for c in cols if "corr(" in c]
        corrdata=self._corrdump[corrcols]
        return corrdata

    def __init__(self, corr, eci):
        """Initialize with correlations in the form corr(0), corr(1), etc,
        and a casm formatted eci.json object.
        The correlation table might also need to include the configname.

        :corr: pd DataFrame
        :eci: json (eci.json)

        """
        self._corrdump = corr
        self._ecidump = eci

        self._corr=self._pure_corr_data()
        self._eci=self._vectorized_eci()

