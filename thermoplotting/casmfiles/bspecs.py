from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import zip
from builtins import object
from builtins import range

import json

#TODO: Consider a module for this type of stuff, like raise.py or something
def _raise_invalid_choice(choice, choices):
    """Test to see if the given choice is part of the allowed list. If not,
    throw an error, otherwise return.

    :choice: obj
    :choices: list of obj
    :returns: void

    """
    if choice not in choices:
        raise ValueError("The specified basis function type '{}' is not in the allowed list.".format(choice))
    else:
        return

def custom_cluster(sites,mode="Integral",subclust=True):
    _raise_invalid_choice(mode,["Integral","Direct"])

    cluster={"coordinate_mode":mode,"include_subclusters":subclust,"prototype":sites}
    return cluster

def _insert_basis_functions(bspecs, basis_functions):
    """Insert field regarding which type of basis function is desired
    (occupation or chebyshev)

    :bspecs: dict
    :basis_functions: str
    :returns: dict

    """
    allowed_types=["occupation", "chebyshev"]
    _raise_invalid_choice(basis_functions, allowed_types)

    bspecs["basis_functions"]={"site_basis_functions":basis_functions}
    return bspecs

def _insert_branch_specs(bspecs, lengths):
    """Given a list of lengths, sequentially assign values to the max length of
    each cluster size, beginning with pairs.

    :bspecs: dict
    :lengths: list of float
    :returns: dict

    """
    clust_sizes=[i+2 for i in range(len(lengths))]
    bspecs["orbit_branch_specs"]={str(s):{"max_length":l} for s,l in zip(clust_sizes,lengths)}
    return bspecs

def _insert_cluster(bspecs, cluster):
    """Add a single custom cluster in direct or integral mode to the bspecs by specifying the prototype.

    :bspecs: dict
    :mode: str
    :prototype: list of list of int or float (depending on mode)
    :include_subclusters: bool
    :returns: dict

    """
    if "orbit_specs" not in bspecs:
        bspecs["orbit_specs"]=[]

    bspecs["orbit_specs"].append(cluster)

    return bspecs

def _set_custom_clusters(bspecs, clusters):
    """Given a list of custom clusters, set them to the custom
    clusters to add to bspecs 

    Parameters
    ----------
    bspecs : dict
    clusters : list of properly formatted clusters

    Returns
    -------
    dict

    """
    bspecs["orbit_specs"]=clusters
    return bspecs

class Bspecs(object):

    """A class for generating a bspecs file, including adding custom clusters
    to the enumeration"""

    def __init__(self, basis_functions, lengths):
        """Initialize with values for type of basis functions and the max length
        of the branches

        :basis_functions: str ("occupation" or "chebyshev")
        :lengths: list of float

        """
        self._bspecsdict=_insert_basis_functions({},basis_functions)
        self._bspecsdict=_insert_branch_specs(self._bspecsdict,lengths)

    def insert_custom_integral_cluster(self, prototype, subclusters=True):
        """Insert a single custom cluster in integral mode

        :prototype: list of list of int
        :returns: dict

        """
        intcluster=custom_cluster(prototype,"Integral",subclusters)
        self._bspecsdict=_insert_cluster(self._bspecsdict, intcluster)
        return self._bspecsdict

    def miniaturize_basis_set(self, eci_source):
        """Given a eci.json fit, create a set of custom clusters from
        the basis functions with active eci, and set that list to be
        the custom clusters

        Parameters
        ----------
        eci_source : json

        Returns
        -------
        dict

        """
        custom_clusters=[custom_cluster(cf["prototype"]["sites"],"Integral",False) for cf in eci_source["cluster_functions"] if "eci" in cf]
        _set_custom_clusters(self._bspecsdict, custom_clusters)
        return self._bspecsdict

    def dump(self):
        """Mostly for debugging. Print the bspecs dictionary to the screen.
        :returns: void

        """
        print(json.dumps(self._bspecsdict, sort_keys=True, indent=4))
        return


