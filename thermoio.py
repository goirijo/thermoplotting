import glob
import numpy

def collect(directory='./'):
    """Collect all the data in *.thdat files and return values as numpy array
    :returns: array of double numpy array

    """
    globable=directory+'*.thin'
    input_data_files=glob.glob(globable)

    input_data_list=[]
    for filename in input_data_files:
        datadump=numpy.loadtxt(filename, comments='#')
        input_data_list.append(datadump)

    return input_data_list

def clobber():
    """Returns a double numpy array with all input data given by collect()

    :returns: double numpy array

    """

    clobbered_array=numpy.vstack(collect())
    return clobbered_array
