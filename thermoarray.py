import numpy as np
import thermoio
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

        try:
            self._headerwords,dataclob=thermoio.safe_clobber(readfilelist,headerdict)
        except:
            print "Bad input data! All input must be in the same order."
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
        :returns: void

        """
        axis=self.axis(stringaxis)+1
        self._controlled_params=reverse(self._controlled_params, axis)
        self._dependent_params=reverse(self._dependent_params, axis)
        return

