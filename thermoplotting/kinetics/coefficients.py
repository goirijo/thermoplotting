from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import pandas as pd
import numpy as np

class OnsagerCalculator(object):

    """This class will calculate the Onsager coefficients for the different
    species using the information about temperature and unit cell size
    given at construction. The number of species is determined from the
    trajectory data provided.

    The input for the trajectories must be in the form of KineticData.
    Units are in K, eV, A3

    """

    def __init__(self,T,Omega,d=3):
        """Initialize with the kinetic data of a particular KMC run.

        Parameters
        ----------
        T : float (temperature)
        Omega : float (the volume of your unit cell
        d : int (dimensionality of your simulation)

        """
        self._T=T
        self._Omega=Omega
        self._d=d
        self._kb=8.6173303e-5 #eV/K

    def __stack_time_values_for_specie(self, kmc_datas, time, specie):
        """Given a list of KineticData, extract the trajectory vectors
        at the particular time for the specified specie

        Parameters
        ----------
        kmc_datas : list of KineticData
        time : float
        specie : str

        Returns
        -------
        np.array (num_data, num_atoms, 3)

        """
        return np.array([kd.specie_values_at_time(time,specie).stack("atomic")[["x","y","z"]].values for kd in kmc_datas])

    def _stacked_values_for_specie(self, kmc_data, specie):
        """Given a list of KineticData, extract the trajectory vectors
        the specified specie and return the values.

        Parameters
        ----------
        kmc_datas : list of KineticData
        specie : str

        Returns
        -------
        np.array (num_data, num_atoms, 3)

        """
        data_view=kmc_data.specie_data(specie).stack("atomic")[["x","y","z"]].values
        num_atoms=len(kmc_data.specie_cols(specie))
        num_data=kmc_data.values().shape[0]
        return data_view.reshape(num_data,num_atoms,3)

    def _calculate_single(self, kmc_data, specie0, specie1):
        """Use all the data points in the data to determine the average
        Onsager coefficient.

        Parameters
        ----------
        kmc_data : KineticData
        specie0 : str
        specie1 : str

        Returns
        -------
        float or nan

        """
        if specie0 not in kmc_data.specie_names() or specie1 not in kmc_data.specie_names():
            return np.nan

        vals0=self._stacked_values_for_specie(kmc_data, specie0)
        vals1=self._stacked_values_for_specie(kmc_data, specie1)

        sum0=np.sum(vals0,axis=1)
        sum1=np.sum(vals1,axis=1)

        dots=np.sum(sum0*sum1,axis=1)
        # ensemble=np.mean(dots)

        Lsquiggle=dots/(2*self._d*kmc_data.t()*kmc_data.num_atoms())
        return np.mean(Lsquiggle)/(self._Omega*self._T*self._kb)

    def calculate(self, kmc_datas, specie0, specie1):
        """Given a list of KineticData, calculate the Onsager coefficient
        for the specified species (i.e. Lij). For each KineticData provided,
        take the averaging over all the points in that run.

        Parameters
        ----------
        kmc_datas : list of KineticData
        specie0 : str
        specie1 : str

        Returns
        -------
        float

        """
        return np.array([self._calculate_single(kd,specie0,specie1) for kd in kmc_datas])

    # def calculate(self, kmc_datas, time, specie0, specie1):
    #     """Given a list of KineticData, calculate the Onsager coefficient
    #     for the specified species (i.e. Lij). The time used
    #     for the calculation is passed as a value, so all the segments of
    #     data must have at least that amount of time elapsed

    #     Parameters
    #     ----------
    #     kmc_datas : list of KineticData
    #     time : float
    #     specie0 : str
    #     specie1 : str

    #     Returns
    #     -------
    #     float

    #     """
    #     vals0=self._stack_time_values_for_specie(kmc_datas, time, specie0)
    #     vals1=self._stack_time_values_for_specie(kmc_datas, time, specie1)

    #     sum0=np.sum(vals0,axis=1)
    #     sum1=np.sum(vals1,axis=1)

    #     dots=np.sum(sum0*sum1,axis=1)
    #     ensemble=np.mean(dots)

    #     Lsquiggle=ensemble/(2*self._d*time*kmc_datas[0].num_atoms())
    #     return Lsquiggle/(self._Omega*self._T*self._kb)
