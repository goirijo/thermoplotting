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

def casm_energy3_to_np(energyfile):
    """Reads in a ternary casm energy file and strips columns down
    to return a numpy array as if you'd just read a .thin file

    :energyfile: output from `casm energy` command
    :returns: numpy matrix

    """
    energy_data=numpy.genfromtxt(energyfile,usecols=(2,3,0))
    return energy_data

def casm_energy4_to_np(energyfile):
    """Reads in a quaternary casm energy file and strips columns down
    to return a numpy array as if you'd just read a .thin file

    :energyfile: output from `casm energy` command
    :returns: numpy matrix

    """
    energy_data=numpy.genfromtxt(energyfile,usecols=(2,3,4,0))
    return energy_data

def casm_energy3_to_thin(energyfile):
    """Reads in a ternary casm energy file and strips columns down
    to create a new .thin file that with only composition
    and energy columns

    :energyfile: output from `casm energy` command
    :returns: void

    """
    energy_data=casm_energy3_to_np(energyfile)
    numpy.savetxt(energyfile+".thin", energy_data)

def casm_energy4_to_thin(energyfile):
    """Reads in a quaternary casm energy file and strips columns down
    to create a new .thin file that with only composition
    and energy columns

    :energyfile: output from `casm energy` command
    :returns: void

    """
    energy_data=casm_energy4_to_np(energyfile)
    numpy.savetxt(energyfile+".thin", energy_data)

