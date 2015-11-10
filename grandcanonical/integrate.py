import numpy as np
from scipy import integrate
from thermoplotting.misc import *

def beta(betavals, phivals, Omegaref, axis):
    """Integrate grand canonical energy through beta axis
    in N-dimensional array and return N-dimensional array
    of free energy values. Integrates in given order!

    Derivation:
    Phi=U-mu*N
    Omega=U-T*S-mu*N
    b=1/k*T

    b*Omega=b*U-S/k-b*mu*N
    d(Omega*b)=U*db-b*N*dmu-mu*Ndb

    For constant mu:

    d(Omega*b)=Phi*db

    Omega(b_f)*b_f-Omega(b_i)*b_i=INT(Phi*db) from b_i to b_f
    we want Omega(b_f), and pass Omega(b_i) as Omegaref

    :betavals: ndarray, (rank=N) 1/kT
    :phivals: ndarray, (rank=N) U-mu*N
    :Omegaref: ndarray, (rank=N-1) U-mu*N-T*S reference values
    :axis: int, axis corresponding to beta in ndarrays
    :returns: ndarray, (rank=N) U-mu*N-T*S

    """
    intchunk=integrate.cumtrapz(phivals, betavals, initial=0, axis=axis)
    beta_i=bottom_view(betavals,axis)
    Omegadata=1/betavals*(intchunk+np.expand_dims(Omegaref*beta_i,axis))

    return Omegadata

def mu(muvals, xvals, Omegaref, axis):
    """Integrate grand canonical energy through a mu axis
    in N-dimensional array and return N-dimensional array
    of free energy values. Integrates in given order!

    Derivation:
    Phi=U-mu*N
    Omega=U-T*S-mu*N
    b=1/k*T

    b*Omega=b*U-S/k-b*mu*N
    d(Omega*b)=U*db-b*N*dmu-mu*Ndb

    For constant b and all other mu:
    dOmega=-b*N*dmu
    Omega(mu_f)-Omega(mu_i)=-b*INT(N*dmu) from mu_i to mu_f
    we want Omega(mu_f), and pass Omega(mu_i) as Omegaref

    :muvals: ndarray, (rank=N) mu
    :phivals: ndarray, (rank=N) U-mu*N
    :Omegaref: ndarray, (rank=N-1) U-mu*N-T*S reference values
    :axis: int, axis corresponding to beta in ndarrays
    :returns: ndarray, (rank=N) U-mu*N-T*S

    """
    intchunk=integrate.cumtrapz(xvals,muvals,initial=0,axis=axis)
    Omegadata=np.expand_dims(Omegaref,axis)-intchunk

    return Omegadata
    
