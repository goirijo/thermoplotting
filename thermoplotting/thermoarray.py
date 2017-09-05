import numpy as np
import thermoio
import copy
from misc import *

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

    def __init__(self, readfilelist, controlled_var, headerdict={}, decimals=0.000000001):
        """Read a bunch of files, stack them together, and reshape
        it all into an N-dimensional array with axis T, mu0, mu1...

        :readfilelist: List of files to read
        :controlled_var: List of strings, specify controlled parameters of data
        :headerdict: Translate the file headers to the internal standard, drop all other fields
        :decimals: Specify the floating point tolerance. Values will be rounded. Default is 1E-8

        """

        self._readfilelist = readfilelist
        #self._headerdict=headerdict

        self._headerwords,dataclob=thermoio.safe_clobber(readfilelist,headerdict)

        try:
            self._headerwords,dataclob=thermoio.safe_clobber(readfilelist,headerdict)
        except:
            print "Bad input data! All input must be in the same order."
            print self._headerwords
            print headerdict
            raise

        #self._headerwords=thermoio.header_split(readfilelist[0])
        #self._standardheaderwords=[headerdict[field] if field in headerdict else field for field in self._headerwords]


        #store dependent and independent variables in separate lists
        self._controlled_var=[]
        controlled_data_list=[]
        self._dependent_var=[]
        dependent_data_list=[]
        for index, var in enumerate(self._headerwords):
            if(var in controlled_var):
                self._controlled_var.append(var)
                controlled_data_list.append(dataclob[:,index])
            else:
                self._dependent_var.append(var)
                dependent_data_list.append(dataclob[:,index])


        controlled_data=np.vstack(controlled_data_list).T
        dependent_data=np.vstack(dependent_data_list).T

        #squash floating point errors and +0/-0 issues
        controlled_data[controlled_data==0.]=0.
        dependent_data[dependent_data==0.]=0.


        #At this point there is:
        #_readfilelist, which is a list of all the files we read from
        #_controlled_var, which is a list such as ["mu0","mu1","T"]
        #controlled_data, which is the data corresponding to the list above
        #_dependent_var, which is a list such as ["phi","N0","N1"]
        #dependent_data, which is the data corresponding to the list above

        #reshape the array with one dimension per controlled variable

        #Find the number of unique values for each controlled variable
        self._params_shape=tuple(np.unique(col).shape[0] for col in controlled_data.T)   #consider using col.round(decimals=5))

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
            print self._controlled_var
            print self._dependent_var
            raise KeyError("The field "+str(field)+" was neither in the controlled or dependent variables lists")

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
            raise RuntimeError(param_name+" is already a dependent parameter!")

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
        if parameter not in self._controlled_var:
            raise KeyError("The field "+str(field)+" was is not controlled")

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
