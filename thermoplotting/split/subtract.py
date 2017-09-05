from clustcompare import *
from detect import Detector

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
        # eci_dict=[(ix,e) for ix,e in candidate.eci]
        # ix,ecivals=zip(*eci_dict)
        ix,ecivals=zip(*candidate.eci)

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
