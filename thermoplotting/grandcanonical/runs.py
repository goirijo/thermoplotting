import numpy
import access

def prepare_input(filenamelist, num_components):
    """Read in from file list of averages and outputs them as new .thin
    files that are compatible with the expected format:
    mu_i, ... , mu_N, T, beta, gc energy, low temp expansion


    :filenamelist: list of file names, presumably casm output
    :num_components: integer of how many components are in your system (how many chemical potentials)
    :returns: void, but you'll end up with new files

    """ 
    for filename in filenamelist:
        access.prepare_input(filename, num_components)

    return
    
def ceil_copy_mu(data, component, maxmu, num_components):
    """Return slice without any entries that have a chemical potential
    larger than the specified one.

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :maxmu: maximum allowed value for chemical potential
    :num_components: integer of total number of components
    :returns: COPY of double numpy array without certain entries

    """
    return data[data[:,access.mu_ind(component,num_components)]<=maxmu]



def floor_copy_mu(data, component, minmu, num_components):
    """Return slice without any entries that have a chemical potential
    lower than the specified one.

    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :minmu: minimum allowed value for chemical potential
    :num_components: integer of total number of components
    :returns: COPY of double numpy array without certain entries

    """
    return data[data[:,access.mu_ind(component,num_components)]>=minmu]

def slice_copy_mu(data, component, muval, num_components, tolerance=0.000001):
    """Return a slice that only has chemical potential values
    equal to the specified one (within a tolerance. Stupid floats...)
    
    :data: double numpy array of all monte values
    :component: index into component of interest. Start counting at 0.
    :muval: value of chemical potential you want a slice of
    :num_components: integer of total number of components
    :tolerance: float comparison tolerance, defaults to 0.000001
    :returns: COPY of double numpy array without certain entries
    """

    return data[abs(data[:,access.mu_ind(component,num_components)]-muval)<tolerance]


    
def ceil_copy_T(data, maxtemp, num_components):
    """Return slice without any entries that have a temperature
    larger than the specified one.

    :data: triple numpy array of all monte values, for all runs.
    :maxtemp: maximum allowed value for temperature
    :num_components: integer of total number of components
    :returns: COPY of double numpy array without certain entries

    """
    return data[data[:,access.temperature_ind(num_components)]<=maxtemp]


def floor_copy_T(data, mintemp, num_components):
    """Return slice without any entries that have a temperature
    larger than the specified one.

    :data: double numpy array of all monte values
    :mintemp: minimum allowed value for temperature
    :num_components: integer of total number of components
    :returns: COPY of double numpy array without certain entries

    """
    return data[data[:,access.temperature_ind(num_components)]>=mintemp]

def slice_copy_T(data, tempval, num_components, tolerance=0.000001):
    """Return a slice that only has temperature values
    equal to the specified one (within a tolerance. Stupid floats...)
    
    :data: double numpy array of all monte values
    :tempval: value of temperature you want a slice of
    :num_components: integer of total number of components
    :tolerance: float comparison tolerance, defaults to 0.000001
    :returns: COPY of double numpy array without certain entries
    """

    return data[abs(data[:,access.temperature_ind(num_components)]-tempval)<tolerance]

def sort_mu_copy(data, component, num_components):
    """Return copy of your data set, sorted by the chemical potential
    of component

    :data: double numpy array of all monte values
    :component: which component to sort by
    :num_components: integer of total number of components
    :returns: COPY of double numpy array in new order
    """

    reorder=numpy.lexsort((access.mu(data,component,num_components),))
    return data[reorder]

def sort_T_copy(data, num_components):
    """Return copy of your data set, sorted by temperature

    :data: double numpy array of all monte values
    :num_components: integer of total number of components
    :returns: COPY of double numpy array in new order
    """

    reorder=numpy.lexsort((access.temperature(data, num_components),))
    return data[reorder]

def heat_capacity(data, num_components):
    """Use fluctuations in grand canonical energy to
    calculate the heat capacity. Returns new array.

    :type data: ndarray
    :type num_components: int
    :returs: ndarray
    """

    genavg2=access.energy(data, num_components)*access.energy(data,num_components)
    gen2avg=access.energy2(data, num_components)
    oneoverT=(access.temperature(data,num_components))**(-1)
    betavals=access.beta(data,num_components)

    heatcap=oneoverT*((gen2avg-genavg2)*betavals*oneoverT)
    return heatcap

def gibbs(data, num_components):
    """Take Legendre transform of grand canonical
    free energy and return new array with Gibbs
    free energy

    :type data: ndarray
    :type num_components: int
    :returs: ndarray
    """

    gcfree=access.free_energy(data,num_components)
    gibbsfree=numpy.copy(gcfree)
    
    for component in xrange(0,num_components):
        currmu=access.mu(data,component,num_components)
        curspecies=access.species(data,component,num_components)
        gibbsfree+=currmu*curspecies

    return gibbsfree
    

