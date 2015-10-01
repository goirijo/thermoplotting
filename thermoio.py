import glob
import numpy as np

def collect(directory='./'):
    """Collect all the data in *.thdat files and return values as np array
    :returns: array of double np array

    """
    globable=directory+'*.thin'
    input_data_files=glob.glob(globable)

    input_data_list=[]
    for filename in input_data_files:
        datadump=np.loadtxt(filename, comments='#')
        input_data_list.append(datadump)

    return input_data_list

def clobber():
    """Returns a double np array with all input data given by collect()

    :returns: double np array

    """

    clobbered_array=np.vstack(collect())
    return clobbered_array

def header_split(filename):
    """Strip '#' from first line and get string list of column

    :filename: file to get header from
    :returns: list of strings

    """
    with open(filename,'r') as f:
        header=f.readline()
        #strip header and split into words
        header=header.translate(None,'#')
        headernames=header.split()

    return headernames

def safe_clobber(readfilelist):
    """Strip header from list of files and stack data onto
    a np array. Checks to make sure all headers match.

    :readfilelist: List of files to read
    :returns: ndarray

    """
    finalheader=header_split(readfilelist[-1])

    datalist=[]

    for filename in readfilelist:
        currentheader=header_split(filename)
        
        if finalheader!=currentheader:
            print "Final header:"
            print finalheader
            print "Current header:"
            print currentheader
            print "Working on file "+str(filename)
            raise AssertionError("Header mismatch while loading files!")

        npdata=np.loadtxt(filename)
        datalist.append(npdata)
    
    dataclob=np.vstack(datalist)
    return dataclob

def casm_energy3_to_np(energyfile):
    """Reads in a ternary casm energy file and strips columns down
    to return a np array as if you'd just read a .thin file

    :energyfile: output from `casm energy` command
    :returns: np matrix

    """
    energy_data=np.genfromtxt(energyfile,usecols=(2,3,0))
    return energy_data

def casm_energy4_to_np(energyfile):
    """Reads in a quaternary casm energy file and strips columns down
    to return a np array as if you'd just read a .thin file

    :energyfile: output from `casm energy` command
    :returns: np matrix

    """
    energy_data=np.genfromtxt(energyfile,usecols=(2,3,4,0))
    return energy_data

def casm_energy3_to_thin(energyfile):
    """Reads in a ternary casm energy file and strips columns down
    to create a new .thin file that with only composition
    and energy columns

    :energyfile: output from `casm energy` command
    :returns: void

    """
    energy_data=casm_energy3_to_np(energyfile)
    np.savetxt(energyfile+".thin", energy_data)

def casm_energy4_to_thin(energyfile):
    """Reads in a quaternary casm energy file and strips columns down
    to create a new .thin file that with only composition
    and energy columns

    :energyfile: output from `casm energy` command
    :returns: void

    """
    energy_data=casm_energy4_to_np(energyfile)
    np.savetxt(energyfile+".thin", energy_data)

