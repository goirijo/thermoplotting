from clustcompare import *

def _non_sequential_raise(basis):
    """Ensure that the linear function index of the cluster functions
    is sane. That is, all indexes are simply 0, 1, 2... n

    :basis: json
    :returns: void or raises

    """
    expected_ix=range(len(basis["cluster_functions"]))
    given_ix=[cf["linear_function_index"] for cf in basis["cluster_functions"]]

    try:
        assert(not np.any(np.array(expected_ix)-np.array(given_ix)))
    except:
        raise ValueError("The provided indexes of the detector are not sequential! If you are sure you're not missing"
                "basis functions, reindex the basis functions and try again.")
    return


class Detector(object):

    """Holds the basis functions for a particular system
    and can determine which basis functions have which species.
    Your basis must be occupation for this to make any sense."""

    def __init__(self, basis):
        """Initialize with basis.json or similar

        :basis: json

        """
        _non_sequential_raise(basis)

        self._basis = basis
        self._species_to_basis=species_to_basis_dict(self._basis)

    def _raise_if_invalid_species(self, specie):
        """raise error informing that the requested species is not
        part of the allowed set.

        :specie: str
        :returns: void

        """
        if specie not in self._species_to_basis:
            raise ValueError("The required specie "+str(specie)+" doesn't participate in the basis set")
        else:
            return

    def _formula_has_specie(self, formula, specie):
        """Checks if the formula's cluster has at least one instance
        of the given site

        :specie: str
        :formula: str
        :returns: bool

        """
        self._raise_if_invalid_species(specie)
        return formula_has_any_bfunc(formula,self._species_to_basis[specie])
        

    def _formula_has_multi_specie(self, formula, specie):
        """Checks if the formula's cluster has multiple instances
        of the given site

        :specie: str
        :formula: str
        :returns: bool

        """
        self._raise_if_invalid_species(specie)
        return formula_has_multi_any_bfunc(formula,self._species_to_basis[specie])

    def _formula_has_exclusively(self, formula, species):
        """Returns true if the basis functions in the formula all correspond to one
        of the specified species

        :species: list of str
        :returns: bool

        """
        for k in self._species_to_basis:
            if k not in species:
                if self._formula_has_specie(formula,k):
                    return False
        return True

    def exclusive_indexes(self, species):
        """Return the indexes for all the basis functions that are made up exclusively
        from basis functions corresponding to the specified species. Use this function
        to extract a subset of basis functions, such as a binary subspace from a ternary.

        :species: list of str
        :returns: list of int

        """
        clust_funcs=self._basis["cluster_functions"]
        indexes=[cf["linear_function_index"] for cf in clust_funcs 
                if self._formula_has_exclusively(cf["prototype_function"],species)]
        return indexes

    def detect_clusters(self, detector,expect_all=True):
        """Map the indexes of shared clusters from the given detector to 
        self. As a default, the expectation is that every
        single cluster of the given detector should map
        somewhere on self.

        :detector: Detector
        :expect_all: bool
        :returns: list of (int,int)

        """
        ix_map=[]
        for subfunc in detector._basis["cluster_functions"]:
            found=False
            for func in self._basis["cluster_functions"]:
                if compare_cluster(subfunc,func):
                    found=True
                    ix_map.append((subfunc["linear_function_index"],func["linear_function_index"]))
            if found==False and expect_all:
                raise ValueError("Could not map cluster "+str(subfunc["linear_function_index"])+" onto this basis!")
        return ix_map

    def basis_species(self):
        """Return list of species that make up the basis

        :returns: list of str

        """
        return [k for k in self._species_to_basis]


    def species(self):
        """Return list of all species, including background

        :returns: list of str

        """
        return list(species_from_basis(self._basis))

    def index_map(self,detector,expect_all=True):
        return self.detect_clusters(detector,expect_all)

    def active_indexes(self,eps=0.0005):
        """Return list of indexes corresponding to non-zero ECI values
        :returns: list of int

        """
        active_ix=[]
        for cf in self._basis["cluster_functions"]:
            try:
                eci=cf["eci"]
            except:
                continue

            if abs(eci)>eps:
                active_ix.append(cf["linear_function_index"])

        return active_ix

    def detect_indexes(self, detector, expect_all=True):
        """First detect where the clusters of the given detector lie in self,
        then from that set of basis functions, return only those that contain
        the same species as the basis set of the given detector.

        :detector: Detector
        :expect_all: bool
        :returns: list of (int,int)

        """
        #Find which clusters we're dealing with
        clust_map=self.detect_clusters(detector)

        #Find which species we're dealing with
        subspecies=detector.species()

        #Find basis functions that deal only with species of the given detector
        subspecies_indexes=self.exclusive_indexes(subspecies)

        #Get intersection between the clusters and the basis functions with the appropriate species
        indexes=[p for p in clust_map if p[1] in subspecies_indexes]
        return indexes

    def detect_active_indexes(self, detector, expect_all=True, eps=0.0005):
        """Map the indexes from the given detector onto self, but return
        only the ones that have an active ECI

        :detector: Detector
        :expect_all: bool
        :returns: list of (int,int)
        :eps: float tolerance

        """
        all_indexes=self.detect_indexes(detector,expect_all)
        subactive=detector.active_indexes(eps)
        active_indexes=[p for p in all_indexes if p[0] in subactive]

        return active_indexes

    def detect_inactive_indexes(self, detector, expect_all=True, eps=0.0005):
        """Map the indexes from the given detector onto self, but return
        only the ones that DO NOT have an active ECI

        :detector: Detector
        :expect_all: bool
        :returns: list of (int,int)
        :eps: float tolerance

        """
        all_indexes=self.detect_indexes(detector,expect_all)
        subactive=detector.active_indexes(eps)
        inactive_indexes=[p for p in all_indexes if p[0] not in subactive]

        return inactive_indexes

    def vectorized_eci(self):
        """Extract the eci values from the json format into a vector of mostly zeros with
        the appropriate values in the right positions

        :returns: pd DataFrame

        """
        zeros=np.zeros(len(self._basis["cluster_functions"]))

        for cf in self._basis["cluster_functions"]:
            if "eci" in cf:
                ix=cf["linear_function_index"]
                val=cf["eci"]
                zeros[ix]=val


        # if "fit" in self._basis:
        #     print "NAILED IT"
        #     for ix,eci in self._basis["fit"]["eci"]:
        #         zeros[ix]=eci

        #This should properly take care of the index-eci relationship
        return pd.DataFrame({"eci":zeros})

    def basis(self):
        """Return copy of the basis json used at construction
        :returns: json

        """
        return self._basis.copy()

    def traced_eci_basis(self, ecivals):
        """Go through the given eci values and add entries into the basis json object.
        The eci values must have an index associated with each value, which corresponds
        to the linear function index. Any preexisting eci values will be removed.

        :ecivals: pd Series
        :returns: json

        """
        traced_basis=self.basis()

        for cf in traced_basis["cluster_functions"]:
            if "eci" in cf:
                cf.pop("eci")
            if cf["linear_function_index"] in ecivals.index:
                cf["eci"]=ecivals.loc[cf["linear_function_index"]]

        return traced_basis

        

