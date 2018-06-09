from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import pandas as pd
import numpy as np


class KineticTrajectory(object):
    """A trajectory is a list of x,y,z and time coordinates for a single
    atom in a kinetic Monte Carlo simulation, which has the values of that
    atom after every hop that happens in the simulation. When dealing with data
    for several atoms, do not use this class. Instead use KineticData."""

    def __init__(self, x, y, z, t, copy=False):
        """Initialize with a list of coordinates

        Parameters
        ----------
        x : x component of coordinate
        y : y component of coordinate
        z : z component of coordinate
        t : time elapsed for the current coordinate
        copy : if True, creates copy of the data passed


        """
        #Not convinced this is managing the memory the way you think, but "not copying" appears to be faster
        if copy:
            self._data = pd.DataFrame(data={"x": x.copy(), "y": y.copy(), "z": z.copy(), "t": t.copy()})
        else:
            self._data = pd.DataFrame(data={"x": x, "y": y, "z": z, "t": t})

        #Add the norm of the distances
        self._data["r"]=np.sqrt(np.square(self._data[["x","y","z"]]).sum(axis=1))

    def x(self):
        return self._data["x"]

    def y(self):
        return self._data["y"]

    def z(self):
        return self._data["z"]

    def t(self):
        return self._data["t"]

    def r(self):
        return self._data["r"]

    def data(self):
        return self._data

    def size(self):
        return len(self.t())

    def as_matrix(self):
        return self._data[["x","y","z","t"]].as_matrix()

    def segment(self, n):
        """Split the trajectory into n independent looking
        trajectories. If the number of samples is not divisible
        by n, the remainder will be discarded.

        Parameters
        ----------
        n : int

        Returns
        -------
        list[KineticTrajectory]

        """
        block_size=self.size()//n
        data_blocks=[self._data.loc[i*block_size:(i+1)*block_size,["x","y","z","t"]] for i in xrange(n)]

        for ix,d in enumerate(data_blocks[1::]):
            d-=self._data.loc[block_size*(ix+1)-1]

        return [KineticTrajectory(**d) for d in data_blocks]




class KineticData(object):
    """Store and retrieve kinetic Monte Carlo data by type of
    species, and other conveniences. This is meant to store a single
    KMC simulation from start to finish"""

    def _input_sanity_raise(self,trajectories, time, occ_species):
        if(trajectories.shape[0]!=len(occ_species)):
            raise ValueError("There must be an xyz trajectory for each species to name")

        if(trajectories.shape[1]!=len(time)):
            raise ValueError("There must be as many time data points as there are coordinates for each atom")

        if(trajectories.shape[2]!=3):
            raise ValueError("The trajectories arrays must hold only values for the x, y, and z coordinates")

        return

    def _master_dataframe(self, trajectories, time, occ_species):
        """Given the constructor data, create the master DataFrame that holds
        all the information about the trajectories of each atom, including what
        species each one is and where it was sitting at the beginning of the
        KMC simulation cell.

        Parameters
        ----------
        trajectories : list of tx3 arrays of length s as np.array
        time : array of float of length t
        occ_species list of str of length s

        Returns
        -------
        pd.DataFrame

        """
        #Create the labels for each atom, with the species name and the index into the starting configdof
        occ_labels=[o+"({})".format(ix) for o,ix in zip(occ_species,xrange(len(occ_species)))]

        #Calculate the norm of the displacements for every atom at every time step
        norms=np.linalg.norm(trajectories,axis=2)

        assert(len(occ_labels)==len(trajectories))
        assert(len(norms)==len(trajectories))
        #The concatenated numpy array now has shape[2]==4 with the norm travelled as a new value
        full_trajectory_data=np.concatenate((trajectories,np.expand_dims(norms,2)),axis=2)
        assert(full_trajectory_data.shape[2]==4)

        #Create MultiIndex for columns, which will group x,y,z,r by atom doing the trajectory
        labels0=[ix for ix,_ in enumerate(occ_labels) for i in xrange(4)]
        assert(labels0[0]==labels0[3] and labels0[-1]==labels0[-4])
        labels1=[i for ix,_ in enumerate(occ_labels) for i in xrange(4)]
        assert(labels0[1]==labels1[-4])
        col_mix=pd.MultiIndex(levels=[occ_labels,["x","y","z","r"]],labels=[labels0,labels1],names=["atomic","cart"])

        #Reshape the trajectory data so that it's 2 dimensional, with the xyzr columns side by side
        nats,ntime,ndim=full_trajectory_data.shape
        data_digest=full_trajectory_data.transpose(0,2,1).reshape(nats*ndim,ntime).T

        #Include the time into the set of data as an additional Index
        time_ix=np.arange(ntime)
        timed_mix=pd.MultiIndex(levels=[time_ix,time],labels=[time_ix,time_ix],names=["index","time"])

        #Create the master DataFrame, this has all the things and has columns at two levels:
        #by species and by trajectory. There are two index levels, sample index and time
        master_frame=pd.DataFrame(data_digest,index=timed_mix,columns=col_mix)
        return master_frame
        

    def __init__(self, trajectories, time, occ_species,direct=None):
        """Initialize with a list of trajectories, the elapsed time per step,
        and a list of the occupation name for each atom. Assumes all data
        comes in incremental time (will not sort anything).

        Internally this is a multi-indexed Pandas array, where one level
        deals with the atoms, naming each "column" things like "Ni(0)", "Al(1)",
        etc, to indicate the species and the index into the unrolled configuration
        of the starting config, as well as the elapsed time, which is common across
        every atom. The other level deals with columns of type "x", "y",
        or "z" to keep track of the trajectory of each atom. The master data
        should always remain in a state where Level 0 refers to the atom labels
        and Level 1 refers to the trajectories

        Parameters
        ----------
        trajectories : list of 3xt arrays of length s as np.array
        time : array of float of length t
        occ_species : list of str of length s
        direct : pd.DataFrame, bypasses the normal construction


        """
        if(direct is None):
            self._input_sanity_raise(trajectories, time, occ_species)
            self._master=self._master_dataframe(trajectories,time,occ_species)
        
        else:
            self._master=direct

        return

    def atom_cols(self, va_as_specie=False):
        """Return array of the column names for every atom.
        If specified, include the vacancies as a specie.

        Parameters
        ----------
        va_as_specie : bool

        Returns
        -------
        list

        """
        everything=self._master.columns.get_level_values("atomic").unique()
        if va_as_specie:
            return everything
        else:
            return [x for x in everything if "Va" not in x]
        
    def specie_cols(self, specie):
        """Return an array of column names that can be used to index into
        every trajectory of a particular specie

        Parameters
        ----------
        specie : str

        Returns
        -------
        list of str

        """
        return [s for s in self.atom_cols() if specie in s]

    def num_atoms(self,va_as_specie=False):
        """Returns total number of sites that there is data for
        If specified, include the vacancies as a specie.

        Parameters
        ----------
        va_as_specie : bool

        Returns
        -------
        int

        """
        return len(self.atom_cols(va_as_specie))

    def composition(self, specie, va_as_specie=False):
        """Returns the ratio of number of specie to total number of atoms
        (not including vacancies unless specified)

        Parameters
        ----------
        specie : str

        Returns
        -------
        float

        """
        return len(self.specie_cols(specie))/self.num_atoms(va_as_specie)

    def index_trajectory(self, index):
        """Return the x, y, z, and t values of a particular atom throughout
        the simulation, specifying only the index and not the specie

        Parameters
        ----------
        atom : int

        Returns
        -------
        pd.DataFrame with x,y,z columns and t as secondary index

        """
        for a in self.atom_cols():
            if "({})".format(index) in a:
                return self.atomic_trajectory(a)

    def atomic_trajectory(self, atom):
        """Return the x, y, z, and t values of a particular atom throughout
        the simulation

        Parameters
        ----------
        atom : str (e.g. Ni(9))

        Returns
        -------
        pd.DataFrame with x,y,z columns and t as secondary index

        """
        return self._master[atom]

    def specie_data(self, specie):
        """Return only the data for a particular species

        Parameters
        ----------
        specie : str

        Returns
        -------
        pd.DataFrame

        """
        return self._master[self.specie_cols(specie)]

    def specie_names(self):
        """Returns the names of all species present
        Returns
        -------
        set of str

        """
        all_cols=self.atom_cols(va_as_specie=True)
        return set([col.split("(")[0] for col in all_cols])

    def _column_swap(self):
        """return the master data with cart over atomic
        Returns
        -------
        DataFrame

        """
        return self._master.swaplevel("atomic","cart",axis=1)

    def x(self):
        return self._column_swap["x"]

    def y(self):
        return self._column_swap["y"]

    def z(self):
        return self._column_swap["z"]

    def r(self):
        return self._column_swap["r"]

    def t(self):
        return self._master.index.get_level_values("time").values

    def _index_at_time(self, time):
        """Return the index (row) corresponding to the data
        for the instant just after (or equal to) the specified time

        Parameters
        ----------
        time : float

        Returns
        -------
        int

        """
        return self._master[self.t()>=time].index.get_level_values("index")[0]

    def values_at_time(self, time):
        """Return the values of everything just below the value of
        the time specified.

        Parameters
        ----------
        time : float

        Returns
        -------
        pd.DataFrame 

        """
        return self._master.loc[self._index_at_time(time)]

    def specie_values_at_time(self, time, specie):
        """Return the values of everything just below the value of
        the time specified, but only for the desired specie

        Parameters
        ----------
        time : float
        specie : str

        Returns
        -------
        pd.DataFrame

        """
        specie_dump=self.specie_data(specie)
        return specie_dump.loc[self._index_at_time(time)]

    def independized_measurements(self):
        """Similar to segmenting the data into multiple apparently independent run,
        this routine will make every point appear to have started at t=0 and r=0.
        This can be useful for data you collect where you don't sample every step,
        and you'd like to keep all the "final" data points in the same array.

        Returns
        -------
        KineticData
        

        """
        #create copy of data and subtract out values
        indep=self._master.copy()
        indep.iloc[1::]=indep.iloc[1::].values-indep.iloc[0:-1]

        #fix the distance
        stacked=indep.stack("atomic")
        stacked["r"]=np.linalg.norm(stacked[["x","y","z"]],axis=1)
        indep=stacked.unstack("atomic").stack("cart").unstack("cart")

        #set the time
        reset_time=self._master.index.get_level_values("time").values
        reset_time[1::]=reset_time[1::]-reset_time[0:-1]
        indep.index.set_levels(reset_time,"time",inplace=True)

        return KineticData(None,None,None,direct=indep)


    def _indexed_segmentation(self, end_inds):
        """Given indexes into the sampled data, split
        the master DataFrame into the specified chunks,
        and reset the elapsed time and coordinates
        such that each segment appears to be an independent
        run


        Parameters
        ----------
        end_inds : list of int, each int is the "up to" index of each segment

        Returns
        -------
        list of KineticData

        """
        start_inds=[0]+end_inds[0:-1]
        raw_segments=[self._master.iloc[ix:nx] for ix,nx in zip(start_inds,end_inds)]
        # raw_segments=[self._master.iloc[seg_length*s:seg_length*(s+1)] for s in xrange(n)]
        n=len(raw_segments)

        #We will subtract the values of the "previous simulation", starting with
        #the final segment
        #These are indexes in reverse that exclude zero
        rev_seg_ix=np.arange(n-1)[::-1]+1
        for rix in rev_seg_ix:
            raw_segments[rix]=raw_segments[rix]-raw_segments[rix-1].iloc[-1]
            #The norm (r) needs to be recalculated
            raw_segments[rix]=raw_segments[rix].stack("atomic")
            raw_segments[rix]["r"]=np.linalg.norm(raw_segments[rix][["x","y","z"]],axis=1)
            raw_segments[rix]=raw_segments[rix].unstack("atomic").stack("cart").unstack("cart")
            #The time also needs to be reset
            reset_time=self._master.index.get_level_values("time")-raw_segments[rix-1].index.get_level_values("time")[-1]
            raw_segments[rix].index.set_levels(reset_time,"time",inplace=True)

        return [KineticData(None,None,None,direct=raw) for raw in raw_segments]

    def sampled_segmentation(self, n):
        """Split the data into n KineticData as if the data
        had been run independently, subtracting out time and
        coordinates so that they start at zero. Remainder data
        is discarded.

        Parameters
        ----------
        n : int

        Returns
        -------
        list of KineticData

        """
        seg_length=len(self._master)//n
        seg_inds=[seg_length*(i+1) for i in xrange(n)]
        return self._indexed_segmentation(seg_inds)

    def timed_segmentation(self, n):
        """Return segments of data in which equal sets
        of time have elapsed

        Parameters
        ----------
        time : int

        Returns
        -------
        list of KineticData

        """
        time_length=self.total_time()/n
        time_inds=[self._index_at_time(time_length*(i+1)) for i in xrange(n)]
        return self._indexed_segmentation(time_inds)

    def values(self):
        """Return all the data ever
        Returns
        -------
        pd.DataFrame

        """
        return self._master

    def total_time(self):
        """Returns the most amount of time elapsed
        Returns
        -------
        float

        """
        return self._master.index.get_level_values("time")[-1]




