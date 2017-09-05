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

        


class FitBlob(object):

    """Combines a set of correlations, energies and ECI values to give you
    the cluster expanded energies"""

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


class Subtracter(object):

    """Holds correlations, energies, and eci values for subspace fits you already have.
    As subspace fits are provided, correlations are dropped and energies are subtracted,
    so that the only basis functions left are those you haven't fit to already.
    
    For example, a Ni-Al-Cr system can be reduced with Ni-Al and Ni-Cr fits, so that only
    difference in energies are used to fit the ternary basis functions."""

    def _pure_corr_data(self):
        """From the DataFrame containing the correlations, extract the subset of columns
        that have the form "corr(x)"

        :returns: pd DataFrame

        """
        cols=self._datadump.columns.values
        corrcols=[c for c in cols if "corr(" in c]
        corrdata=self._datadump[corrcols]
        return corrdata

    def __init__(self, fitdata, basis):
        """Initialize with a pandas DataFrame that contains correlation values, confignames
        and formation_energy values

        :fitdata: pd DataFrame
        :basis: json

        """
        self._detector=Detector(basis)
        self._datadump=fitdata
    
        self._datadump["subtracted_energy"]=0

        self._corr=self._pure_corr_data()
        self._eci=self._detector.vectorized_eci()

        self._dropped_corr=[]

    def _indexed_eci(self, detector):
        """Extract the eci values from the given detector and reassign index values so that the
        slots match appropriately with self

        :detector: Detector
        :returns: pd DataFrame

        """
        indexes=self._detector.detect_indexes(detector)
        sub_eci=detector.vectorized_eci()

        from_ix=np.array(indexes)[:,0]
        to_ix=np.array(indexes)[:,1]

        #Ensure that the from indexes are simply values of 0,1,2,3,4...n
        #This part is now redundant due to check at detector construction,
        #but I'll leave it just in case
        try:
            assert(not np.any(from_ix-range(len(from_ix))))
        except:
            raise ValueError("The provided indexes of the detector are not sequential! If you are sure you're not missing\
                    basis functions, reindex the basis functions and try again.")

        sub_eci.index=to_ix
        return sub_eci

    def _trace_eci(self, detector):
        """Take the eci values from the given detector and place them into the appropriate
        slots of self

        :detector: Detector
        :returns: void

        """
        sub_eci=self._indexed_eci(detector)
        self._eci.loc[sub_eci.index]+=sub_eci

        return

    def _expanded_energy(self):
        """Dot the full correlations with all the active ECI to determine the energy contribution
        of the active basis functions
        :returns: pd Series

        """
        corrmat=self._corr.as_matrix()
        ecivec=self._eci.as_matrix()

        exp_energy=np.dot(corrmat,ecivec)
        return pd.Series(exp_energy[:,0], index=self._corr.index)


    def _sub_correlation_matrix(self):
        """Run through the correlations and return a view where all columns corresponding
        to ECI values that are non-zero have been dropped
        :returns: pd DataFrame

        """
        inactive_corr_cols=[col for col in self._corr if col not in self._dropped_corr]

        return self._corr[inactive_corr_cols]

    def _drop_correlations(self, indexes):
        """Add the provided indexes to the list of correlations that should be removed when
        looking at the sub-correlation matrix
        :indexes: list of int
        :returns: void

        """
        new_dropped_columns=["corr("+str(i)+")" for i in indexes]
        self._dropped_corr=self._dropped_corr+new_dropped_columns
        self._subcorr=self._sub_correlation_matrix()

        return


    def subtract(self, detector, drop_correlations=True):
        """Map basis functions of the given detector onto self and subtract the active
        eci values out, effectively eliminating a subspace from the fit you need to make.

        :detector: Detector
        :returns: void

        """
        #insert ECI values
        self._trace_eci(detector)

        #Subtract energy contribution of active ECI
        self._datadump["subtracted_energy"]=self._datadump["formation_energy"]-self._expanded_energy()

        if drop_correlations==True:
            #Eliminate basis functions from detector that have active ECI
            indexes=self._detector.detect_indexes(detector)
            to_ix=np.array(indexes)[:,1]
            self._drop_correlations(to_ix)

        return

    def merged_eci(self, hallcandidate):
        """Given a candidate from the hall of fame, add the eci to the current values.
        Return a Detector object with the updated eci values, which can be used as a final
        eci.json file.

        The idea behind this routine is that an instance of a Subtracter has been created,
        had values subtracted, then the difference has been fit. The eci of the fit need to be
        combined back with the initial eci that were subtracted.

        :hallcandidate: Deap object, one of the fits from your hall of fame
        :returns: json

        """
        eci_dict={ix:e for ix,e in hallcandidate.eci}

        merged_eci=self._detector.basis()
        for cf in merged_eci["cluster_functions"]:
            lfix=cf["linear_function_index"]
            eci_val=self._eci.loc[lfix]["eci"]

            if lfix in eci_dict:
                eci_val+=eci_dict[lfix]

            if eci_val!=0:
                cf["eci"]=eci_val

        return merged_eci

    def _non_dropped_ix(self):
        """Return the indexes corresponding to correlations that have not been dropped during the
        subtraction

        :returns: np array
        """
        full_cols=self._corr.columns.values
        non_dropped=[c for c in full_cols if c not in self._dropped_corr]
        #extract the integer n out of corr(n)
        non_dropped=np.array([c[5:-1] for c in non_dropped],dtype=int)

        return non_dropped

    def sub_halloffame_trace(self, candidate):
        """Given a halloffame candidate from the subtracted fit, trace the eci values back onto
        self, adding the eci back onto the existing ones.

        :candidate: Deap object, one of the fits from your hall of fame
        :returns: pd DataFrame

        """
        eci_dict=[(ix,e) for ix,e in candidate.eci]
        ix,ecivals=zip(*eci_dict)

        subeci=pd.Series(ecivals,ix)
        return self.sub_eci_trace(subeci)

    def sub_eci_trace(self, subeci):
        """Given a vector of eci corresponding to the current subtracted view
        of self, insert the eci values into the global eci vector.

        For example, after subtracting out binary basis functions, the resulting
        correlations and subtracted energies are used to determine a some eci. This
        eci vector will be shorter than the global one, but can be traced back in
        with this routine.

        eci values are *added* back in, not replaced.

        :subeci: pd Series
        :returns: pandas DataFrame

        """
        non_dropped=self._non_dropped_ix()

        non_dropped_eci=self._eci.loc[non_dropped]
        non_dropped_eci.reset_index(inplace=True)
        non_dropped_eci.loc[subeci.index,"eci"]+=subeci
        non_dropped_eci.set_index("index",inplace=True)

        self._eci.loc[non_dropped_eci.index,"eci"]=non_dropped_eci["eci"]

        return self.eci()

    def eci(self):
        """Return copy of all ECI values that have been traced onto self
        :returns: pd DataFrame

        """
        return self._eci.copy()

    def sub_correlation_matrix(self):
        """Return a copy of the correlation matrix after having dropped whatever columns.
        :returns: pd DataFrame

        """
        return self._sub_correlation_matrix().copy()

    def sub_energy(self):
        """Return copy of the energies after you've subtracted out certain ECI contributions.
        :returns: pd.DataFrame

        """
        return self._datadump["subtracted_energy"].copy()

    def sub_data(self):
        """Return a copy the DataFrame used to initialize self, but add a column for subtracted
        energy, and drop the appropriate correlation matrices.
        :returns: pd DataFrame

        """
        extra_data_cols=[c for c in self._datadump.columns if "corr" not in c]
        extra_data=self._datadump[extra_data_cols]

        corrdata=self._sub_correlation_matrix()

        datacombo=pd.concat((extra_data,corrdata),axis=1)
        return datacombo.copy()

    def as_json_fit(self):
        """Trace all the stored eci values onto the basis json object and return it so that
        it can be used as a casm fit
        :returns: json

        """
        non_zero_eci=self._eci.loc[self._eci["eci"]!=0]
        return self._detector.traced_eci_basis(non_zero_eci["eci"])
