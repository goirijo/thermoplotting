from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object

import numpy as np
import pandas as pd
import copy
from .misc import *

class ThermoArray(object):

    """Organize all your Monte Carlo data into a big N-dimensional
    array for easy sorting, viewing and integration. Expects headers
    in the input files to know what each column represents.
    You must provide a list of strings of the header fields that correspond to
    the controlled variables of your run (e.g. mu, T).

    In addition, you can provide a dictionary to set new strings for the
    field names in case you don't like the ones that your input files have or
    want to shorten them to something else when referencing them.

    """
    def _raise_bad_input_data(self, datalist):
        """Ensure that the given data makes sense. Controlled variables must be
        regularly gridded

        :datalist: input set of data at construction
        :returns: void or raises

        """
        length=len(datalist[0])
        for d in datalist:
            if len(d)!=length:
                raise ValueError("The list of data objects must all have the same lengh!")
                
        return


    def _raise_not_controlled(self, parameter):
        """If the given parameter is not a controlled variable,
        raise an error

        Parameters
        ----------
        parameter : str

        Returns
        -------
        void or raises

        """
        if parameter not in self._controlled_var:
            raise KeyError("The field {} was is not controlled".format(field))
        else:
            return

    def _prepare_data(self, datalist, headerdict, decimals):
        """Put all the data sets together and rename the columns
        to match the headerdict

        :datalist: list of pd DataFrame
        :headerdict: dict
        :decimals: int
        :returns: pandas DataFrame

        """
        #At this point the data is assumed to be good, all columns match, etc
        combined_data=pd.concat(datalist,ignore_index=True)
        new_col_names=[headerdict[c] if c in headerdict else c for c in combined_data.columns]
        combined_data.columns=new_col_names

        #rounding does not remove -0.0, thought it shouldn't matter
        return combined_data.round(decimals)

    def __init__(self, datalist, controlled_var, headerdict={}, decimals=8):
        """Stack DataFrames together that all have the same controlled variables
        it all into an N-dimensional array with axis T, mu0, mu1...

        :datalist: List of data to read (presumably from all the results.json)
        :controlled_var: List of strings, specify controlled parameters of data
        :headerdict: Translate the file headers to the internal standard, drop all other fields
        :decimals: int Specify the floating point tolerance. Values will be rounded. Default is 1E-8 (value 8)

        """
        self._raise_bad_input_data(datalist)

        # self._readfilelist = readfilelist
        self._headerdict=headerdict

        unrolled_data=self._prepare_data(datalist,headerdict,decimals)

        #store dependent and independent variables in separate lists
        colnames=unrolled_data.columns
        self._controlled_var=[c for c in colnames if c in controlled_var]
        self._dependent_var=[c for c in colnames if c not in controlled_var]

        controlled_data=unrolled_data[self._controlled_var].as_matrix()
        dependent_data=unrolled_data[self._dependent_var].as_matrix()

        #At this point there is:
        #_datalist, which is a list of all the data we started with
        #_controlled_var, which is a list such as ["mu0","mu1","T"]
        #controlled_data, which is the data corresponding to the list above
        #_dependent_var, which is a list such as ["phi","N0","N1"]
        #dependent_data, which is the data corresponding to the list above

        #reshape the array with one dimension per controlled variable

        #Find the number of unique values for each controlled variable
        self._params_shape=tuple(np.unique(col).shape[0] for col in controlled_data.T)

        #get the set of row indexes that sort the controlled variables in ascending order, starting with last column
        idx=np.lexsort(controlled_data[:,::-1].T)

        #sort and reshape the controlled variables
        self._controlled_params=controlled_data[idx].T.reshape((len(self._controlled_var),)+self._params_shape)

        #sort and reshape the dependent variables
        self._dependent_params=dependent_data[idx].T.reshape((len(self._dependent_var),)+self._params_shape)

    def _tuple_index(self, field):
        """Returns a tuple to access all values along an axis, such as
        [1,:,:,:]
        Used indirectly for getting indexes of different fields
        :returns: tuple

        """
        #if field in self._headerdict:
        #    field=self._headerdict[field]

        pretuple=[slice(None)]*(len(self._params_shape)+1)

        if field in self._controlled_var:
            pretuple[0]=self._controlled_var.index(field)
        elif field in self._dependent_var:
            pretuple[0]=self._dependent_var.index(field)
        else:
            print(self._controlled_var)
            print(self._dependent_var)
            # raise KeyError("The field "+str(field)+" was neither in the controlled or dependent variables lists")
            raise KeyError("The field {} was neither in the controlled or dependent variables lists".format(field))

        return tuple(pretuple)

    def data_view(self, field):
        """Returns view (!!) of the data requested by field. Dimensions
        will match self (implementation of _tuple_index)
        """
        tupindx=self._tuple_index(field)
        
        if field in self._controlled_var:
            return self._controlled_params[tupindx]
        else:
            return self._dependent_params[tupindx]

    def __getitem__(self, field):
        """Calls data_view

        :field: str
        :returns: np view

        """
        return self.data_view(field)

    def top_data_view(self, datafield, axisfield):
        """Return the first value along the axisfield dimension
        of the datafield data

        :datafield: string
        :axisfield: string
        :returns: ndarray

        """
        data=self.data_view(datafield)
        axis=self.axis(axisfield)
        return top_view(data,axis)
    
    def bottom_data_view(self, datafield, axisfield):
        """Return the last value along the axisfield dimension
        of the datafield data

        :datafield: string
        :axisfield: string
        :returns: ndarray

        """
        data=self.data_view(datafield)
        axis=self.axis(axisfield)
        return bottom_view(data,axis)
    
    def nan_fields(self):
        """Make array of NaN with appropriate shape for dealing with self data
        :returns: ndarray

        """
        emptynp=np.empty(self._params_shape)
        emptynp.fill(np.nan)
        return emptynp

    def axis(self, parameter):
        """Returns which axis the given parameter is aligned.
        The given parameter must be a controlled variable.
        :returns: int

        """
        self._raise_not_controlled(parameter)
        axis=self._controlled_var.index(parameter)
        return axis

    def reverse(self, stringaxis, inplace=False):
        """Invert the order along a given axis
        :stringaxis: string label for axis
        :inplace: bool, make true if you want views into data updated
        :returns: ThermoArray

        """
        axis=self.axis(stringaxis)+1
        self._controlled_params=reverse(self._controlled_params, axis,inplace)
        self._dependent_params=reverse(self._dependent_params, axis,inplace)
        return self

    def push_back(self, new_param, param_name):
        """Append a new dependent parameter to the end
        of the list.

        :new_param: ndarray (should have "data_view" dimensions)
        :param_name: string
        :returns: ThermoArray

        """
        if param_name in self._dependent_var:
            raise ValueError(param_name+" is already a dependent parameter!")

        self._dependent_var.append(param_name)
        self._dependent_params=np.append(self._dependent_params,[new_param],axis=0)

        return self

    def _template(self):
        """Return copy of self with all dependent parameter values as nan
        :returns: TODO

        """
        duplicate=copy.deepcopy(self)
        duplicate._dependent_params[...]=np.NaN
        return duplicate

    def controlled_parameter_values(self, parameter):
        """Get a list of the values of one of the controlled parameters

        :parameter: str
        :returns: 1D np

        """
        self._raise_not_controlled(parameter)

        pretuple=[0]*(len(self._params_shape)+1)
        pretuple[0]=self._controlled_var.index(parameter)
        axis=self.axis(parameter)
        pretuple[axis+1]=slice(None)
        
        return self._controlled_params[pretuple]
        

    def mask_list_along(self, parameter):
        """Create a list on indexes that correspond to slices along a particular
        parameter

        :parameter: str
        :returns: list of np indexes

        """
        return [np.where(self.data_view(parameter)==v) for v in self.controlled_parameter_values(parameter) ]

    def as_dataframe(self):
        """Unroll the entire array and return the values as a pandas
        DataFrame

        Returns
        -------
        pd DataFrame

        """
        datadict={}
        for d in self._dependent_var:
            datadict[d]=self[d].ravel()

        for c in self._controlled_var:
            datadict[c]=self[c].ravel()

        return pd.DataFrame(datadict)

    def axis_range(self, parameter):
        """Get list of sampled values along a controlled variable
        axis

        Parameters
        ----------
        parameter : str

        Returns
        -------
        np array

        """
        self._raise_not_controlled(parameter)

        paramix=self._controlled_var.index(parameter)
        sliceix=[0]*(len(self._params_shape)+1)
        sliceix[0]=paramix
        sliceix[paramix+1]=slice(None)

        return self._controlled_params[tuple(sliceix)]
        

    @staticmethod
    def minimize(arraylist,axistag):
        """Create a new ThermoArray by taking the values corresponding
        to the minimum values of a certain parameter (e.g. free energy)

        :arraylist: [ThermoArray], must match dimensions
        :axistag: string
        :returns: ThermoArray

        """
        mininds=argmin_stack([data.data_view(axistag) for data in arraylist])
        minimizedarray=arraylist[0]._template()

        paramstack=stack([data._dependent_params for data in arraylist])
        minimizedarray._dependent_params=np.choose(mininds,paramstack)

        return minimizedarray
