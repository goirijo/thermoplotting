import numpy

def prepare_input(filename, num_components):
    """Read in from file of averages and output them as new .thin
    files that are compatible with the expected format:
    mu_i, ... , mu_N, T, beta, gc energy, low temp expansion


    :filename: file name, presumably casm output
    :num_components: integer of how many components are in your system (how many chemical potentials)
    :returns: void, but you'll end up with a new file

    """ 
    datadump=numpy.loadtxt(filename, comments='#')

    desiredcolumns=[]
    for compind in range(0,num_components):
        mucolumnind=3+num_components+3+compind
        desiredcolumns.append(datadump[:, mucolumnind])

    tempcolumnind=3+num_components+1
    desiredcolumns.append(datadump[:, tempcolumnind])

    betacolumnind=tempcolumnind+1
    desiredcolumns.append(datadump[:, betacolumnind])

    gcecolumnind=3+num_components
    desiredcolumns.append(datadump[:, gcecolumnind])

    lowtexpcolumnind=0  #not really
    desiredcolumns.append(datadump[:, lowtexpcolumnind])
    desiredcolumns[-1]=numpy.zeros(numpy.shape(desiredcolumns[-1]))

    desireddata=numpy.array(desiredcolumns)
    desireddata=numpy.transpose(desireddata)

    outputname=filename+".thin"
    numpy.savetxt(outputname, desireddata)

    return

def mu_ind(component, num_components):
    """Return index into chemical potential column of specified component

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :num_components: integer of total number of components
    :returns: integer

    """
    return component


def temperature_ind(num_components):
    """Return index into temperature column

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: integer

    """
    return num_components


def beta_ind(num_components):
    """Return index into beta column

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: integer

    """
    return num_components+1

def energy_ind(num_components):
    """Return index into energy column

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: integer

    """
    return num_components+2

def low_T_expansion_ind(num_components):
    """Return index into low temperature expansion column

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: integer

    """
    return num_components+3

def mu(data, component, num_components):
    """Fetch values for the specified component

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :num_components: integer of total number of components
    :returns: numpy array of chemical potentials

    """
    return data[:,mu_ind(component, num_components)]


def temperature(data, num_components):
    """Fetch values for temperature

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of temperature

    """
    return data[:, temperature_ind(num_components)]


def beta(data, num_components):
    """Fetch values for beta

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of beta

    """
    return data[:, beta_ind(num_components)]


def energy(data, num_components):
    """Fetch values for grand canonical energy

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of grand canonical energy

    """
    return data[:, energy_ind(num_components)]


def low_T_expansion(data, num_components):
    """Fetch values for low temperature expansion

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of low temperature expansion

    """
    return data[:, low_T_expansion_ind(num_components)]
