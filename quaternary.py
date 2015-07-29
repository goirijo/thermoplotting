import scipy.spatial as spa
import scipy.linalg as lin
import numpy

def prepare_input(file_name):
    """Read in from list of files and output them as new .thin
    files that are compatible with the expected format:
    composition0, composition1, energy

    :file_name: name of file, presumably casm output
    :returns: void, but you'll end up with new a file "oldfile.thin"

    """
    pass


def composition(data_list, component):
    """Returns one composition column from the given data, which
    was presumably made using the standard clobber() function.

    :data_list: double numpy array
    :component: int specifying which composition to return (2=origin)
    :returns: list of first composition values

    """
    if component==3:
        oneslist=numpy.ones(data_list.shape[0])
        origincomp=oneslist-data_list[:,0]-data_list[:,1]-data_list[:,2]
        return origincomp

    else:
        return data_list[:,component]


def energy(data_list):
    """Returns energy column from the given data, which
    was presumably made using the standard clobber() function

    :data_list: double numpy array
    :returns: list of second composition values

    """
    return data_list[:,3]
