import access
import numpy
from scipy import integrate

import matplotlib.pyplot as plt

def smartrapz(y,x,reference):
    """Use the usual cumtrapz, but add the reference
    to every evaluation and also insert the reference
    as the first value

    :x: ndarray 1D integration path
    :y: ndarray 1D values for each x
    :reference: starting value of your integration
    :returns: ndarray 1D

    """
    integrated=integrate.cumtrapz(y,x,initial=0)
    integrated[:]+=reference

    return integrated

def mu(mupath, component, freeenergyref, num_components):
    """Given a 1D path through varying chemical potential,
    integrate the Grand Canonical free energy, filling
    the Grand Canonical energy column as you go. Expects
    data to be in an order that you can integrate in!

    Derivation:
    Phi=U-mu*N
    PHI=U-T*S-mu*N
    b=1/k*T

    b*PHI=b*U-S/k-b*mu*N
    d(PHI*b)=U*db-b*N*dmu-mu*Ndb

    For constant b and all other mu:
    dPHI=-b*N*dmu
    PHI(mu_f)-PHI(mu_i)=-b*INT(N*dmu) from mu_i to mu_f
    we want PHI(mu_f), and pass PHI(mu_i) as freeenergyref

    :mupath: ndarray
    :component: int, which component you're integrating over
    :freeenergyref: starting value for free energy
    :num_components: int, total number of components
    :returns: void

    """
    muvals=access.mu(mupath, component, num_components)
    speciesvals=access.species(mupath, component, num_components)
    betaval=access.beta(mupath,num_components)

    intchunk=integrate.cumtrapz(speciesvals,muvals,initial=0)
    freeenergy=-betaval*intchunk+freeenergyref

    mupath[:,access.free_energy_ind(num_components)]=freeenergy

    return
    

def beta(betapath, freeenergyref, num_components):
    """Given a 1D path through varying temperature, integrate
    the Grand Canonical free energy, filling the Grand Canonical
    free energy column as you go. EXPECTS THE DATA TO BE IN
    AN ORDER THAT YOU CAN INTEGRATE IN!!

    Derivation:
    Phi=U-mu*N
    PHI=U-T*S-mu*N
    b=1/k*T

    b*PHI=b*U-S/k-b*mu*N
    d(PHI*b)=U*db-b*N*dmu-mu*Ndb

    For constant mu:

    d(PHI*b)=Phi*db

    PHI(b_f)*b_f-PHI(b_i)*b_i=INT(Phi*db) from b_i to b_f
    we want PHI(b_f), and pass PHI(b_i) as freeenergyref

    :betapath: numpy.ndarray
    :freeenergyref: initial value for free energy needed as reference
    :returns: the given betapath with new values for the free energy

    """
    
    betavals=access.beta(betapath,num_components)
    tempvals=access.temperature(betapath,num_components)
    grandcanon=access.energy(betapath, num_components)
    
    beta_i=betavals[0]

    intchunk=integrate.cumtrapz(grandcanon,betavals,initial=0)
    freeenergy=1/betavals*(intchunk+freeenergyref*beta_i)

    #Ti=tempvals[0]
    #intchunk=integrate.cumtrapz(grandcanon/(tempvals*tempvals), tempvals, initial=0)
    #freeenergy=(-intchunk+freeenergyref/Ti)*tempvals

    betapath[:,access.free_energy_ind(num_components)]=freeenergy

    return
    
