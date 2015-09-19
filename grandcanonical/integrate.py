import numpy as np
from scipy import integrate
from thermoplotting.misc import *

def beta(betavals, phivals, PHIref, axis):
    """Integrate grand canonical energy through beta axis
    in N-dimensional array and return N-dimensional array
    of free energy values. Integrates in given order!

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

    :betavals: ndarray, (rank=N) 1/kT
    :phivals: ndarray, (rank=N) U-mu*N
    :PHIref: ndarray, (rank=N-1) U-mu*N-T*S reference values
    :axis: int, axis corresponding to beta in ndarrays
    :returns: ndarray, (rank=N) U-mu*N-T*S

    """
    intchunk=integrate.cumtrapz(phivals, betavals, initial=0, axis=axis)
    beta_i=bottom_view(betavals,axis)
    PHIdata=1/betavals*(intchunk+PHIref*beta_i)

    return PHIdata
