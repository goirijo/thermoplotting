import glob
import numpy as np
import json

#This file may be completely obsolete at this point

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

def txt_split(filename):
    """Strip '#' from first line and get string list of column

    :filename: file to get header from
    :returns: [string],ndarray

    """
    with open(filename,'r') as f:
        header=f.readline()
        #strip header and split into words
        header=header.translate(None,'#')
        headernames=header.split()

    npdata=np.genfromtxt(filename)
    return headernames,npdata

def json_split(filename):
    """Return the keys from a Monte Carlo results.json as
    a list of strings. Functionally the same as header_split
    but for json files.

    :filename: string
    :returns: [string],ndarray

    """
    datadump=open(filename).read()
    jsonresults=json.loads(datadump)

    headernames=jsonresults.keys()
    datacolumns=[np.array(jsonresults[key],dtype=float) for key in headernames]

    return headernames,np.vstack(datacolumns).T

def truncate((headerstrings,datacolumns),headerdict):
    """After reading in a file, drop all the columns
    not in headerdict, and rename the ones that are there
    appropriately.

    :headerstrings: [string]
    :datacolumns: ndarray
    :headerdict: dictionary
    :returns: [string],ndarray

    """
    if len(headerdict)==0:
        return headerstrings,datacolumns

    truncheaders=[]
    trunccols=[]
    for key in headerdict:
        idx=headerstrings.index(key)

        truncheaders.append(headerdict[key])
        trunccols.append(datacolumns[:,idx])

    return truncheaders,np.vstack(trunccols).T

def header_split(filename,headerdict):
    """Get data in columns with a corresponding set of strings
    for each one. Works on txt and json.

    :filename: string
    :returns: [string],ndarray

    """
    try:
        if filename[-4::]==".txt" or filename[-4::]==".csv":
            return truncate(txt_split(filename),headerdict)
        elif filename[-5::]==".json":
                return truncate(json_split(filename),headerdict)
        else:
            raise RuntimeError("File name "+filename+" was neither json or txt")

    except:
        raise ValueError("Something went horribly wrong with your input file "+filename+". Check for non equilibrated conditions.")



def safe_clobber(readfilelist,headerdict={}):
    """Strip header from list of files and stack data onto
    a np array. Checks to make sure all headers match.
    Expects txt files (header with columns).

    :readfilelist: List of files to read
    :returns: [string],ndarray

    """
    finalheader,_=header_split(readfilelist[-1],headerdict)

    datalist=[]

    for filename in readfilelist:
        currentheader,npdata=header_split(filename,headerdict)
        
        if finalheader!=currentheader:
            print "Final header:"
            print finalheader
            print "Current header:"
            print currentheader
            print "Working on file "+str(filename)
            raise AssertionError("Header mismatch while loading files!")

        datalist.append(npdata)
    
    dataclob=np.vstack(datalist)

    return currentheader,dataclob

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

