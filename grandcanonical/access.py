import numpy

#def prepare_input(filename, num_components):
#    """Read in from file of averages and output them as new .thin
#    files that are compatible with the expected format:
#    mu_i, ... , mu_N, T, beta, gc energy, low temp expansion
#
#
#    :filename: file name, presumably casm output
#    :num_components: integer of how many components are in your system (how many chemical potentials)
#    :returns: void, but you'll end up with a new file
#
#    """ 
#    datadump=numpy.loadtxt(filename, comments='#')
#
#    desiredcolumns=[]
#    for compind in range(0,num_components):
#        mucolumnind=3+num_components+3+compind
#        desiredcolumns.append(datadump[:, mucolumnind])
#
#    for compind in range(0,num_components):
#        speciescolumnind=3+compind
#        desiredcolumns.append(datadump[:, speciescolumnind])
#
#    tempcolumnind=3+num_components+1
#    desiredcolumns.append(datadump[:, tempcolumnind])
#
#    betacolumnind=tempcolumnind+1
#    desiredcolumns.append(datadump[:, betacolumnind])
#
#    gcecolumnind=3+num_components
#    desiredcolumns.append(datadump[:, gcecolumnind])
#
#    lowtexpcolumnind=0  #not really
#    desiredcolumns.append(datadump[:, gcecolumnind])    #make this a real thing to actually account for low T expansion
#
#
#
#    #column for free energy
#    emptyfreeenergy=numpy.empty(numpy.shape(desiredcolumns[-1]))
#    emptyfreeenergy.fill(numpy.nan)
#    desiredcolumns.append(emptyfreeenergy)
#
#    desireddata=numpy.array(desiredcolumns)
#    desireddata=numpy.transpose(desireddata)
#
#    outputname=filename+".thin"
#    numpy.savetxt(outputname, desireddata)
#
#    return

def mu_ind(component, num_components):
    """Return index into chemical potential column of specified component

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :num_components: integer of total number of components
    :returns: integer

    """
    return component


def species_ind(component, num_components):
    """Return index into column of number of species
    for the specified component

    :component: ind
    :num_components: ind
    """
    return num_components+component


def temperature_ind(num_components):
    """Return index into temperature column

    :num_components: integer of total number of components
    :returns: integer

    """
    return 2*num_components


def beta_ind(num_components):
    """Return index into beta column

    :num_components: integer of total number of components
    :returns: integer

    """
    return temperature_ind(num_components)+1

def energy_ind(num_components):
    """Return index into energy column

    :num_components: integer of total number of components
    :returns: integer

    """
    return beta_ind(num_components)+1

def energy2_ind(num_components):
    """Return index into average of square energy column

    :num_components: integer of total number of components
    :returns: integer

    """
    return energy_ind(num_components)+1

def low_T_expansion_ind(num_components):
    """Return index into low temperature expansion column

    :num_components: integer of total number of components
    :returns: integer

    """
    return energy2_ind(num_components)+1


def free_energy_ind(num_components):
    """Return index into free energy column

    :num_components: integer of total number of components
    :returns: integer

    """
    return low_T_expansion_ind(num_components)+1


def mu(data, component, num_components):
    """Fetch values for the specified component chemical potential

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :num_components: integer of total number of components
    :returns: numpy array of chemical potentials

    """
    return data[:,mu_ind(component, num_components)]


def species(data, component, num_components):
    """Fetch values for the specified component number of species

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :num_components: integer of total number of components
    :returns: numpy array of chemical potentials

    """
    return data[:,species_ind(component, num_components)]


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

def energy2(data, num_components):
    """Fetch average values for the square of the grand
    canonical energy

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of grand canonical energy

    """
    return data[:, energy2_ind(num_components)]


def low_T_expansion(data, num_components):
    """Fetch values for low temperature expansion

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of low temperature expansion

    """
    return data[:, low_T_expansion_ind(num_components)]


def free_energy(data, num_components):
    """Fetch values for grand canonical free energy

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array of grand canonical energy

    """
    return data[:, free_energy_ind(num_components)]

def backwards(data):
    """Invert the order of your data and return a view.


    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: numpy array in new order

    """
    data=data[::-1]
    return data
