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
    clust_sizes=[i+2 for i in xrange(len(lengths))]
    bspecs["orbit_branch_specs"]={str(s):{"max_length":l} for s,l in zip(clust_sizes,lengths)}
    return bspecs

def _insert_cluster(bspecs, mode, prototype, include_subclusters):
    """Add a single custom cluster in direct or integral mode to the bspecs by specifying the prototype.

    :bspecs: dict
    :mode: str
    :prototype: list of list of int or float (depending on mode)
    :include_subclusters: bool
    :returns: dict

    """
    _raise_invalid_choice(mode,["Integral","Direct"])

    if "orbit_specs" not in bspecs:
        bspecs["orbit_specs"]=[]

    intclust={"coordinate_mode":mode, "prototype":prototype,"include_subclusters":include_subclusters}
    bspecs["orbit_specs"].append(intclust)

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
        :returns: TODO

        """
        self._bspecsdict=_insert_cluster(self._bspecsdict, "Integral", prototype, subclusters)
        return

    def dump(self):
        """Mostly for debugging. Print the bspecs dictionary to the screen.
        :returns: void

        """
        print json.dumps(self._bspecsdict, sort_keys=True, indent=4)
        return


