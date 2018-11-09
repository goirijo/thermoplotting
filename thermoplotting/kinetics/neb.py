from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def simple_plot(ax, data):
    """Given a few neb image energies, interpolate
    the data and plot it.
    Uses cubic interpolation, and adds extra flat datapoints
    at the edges to get a minima at the ends.

    Parameters
    ----------
    ax : plt plottable canvas
    data : list of NEB image energies, in order

    Returns
    -------
    ax

    """
    data=np.array(data)-(data[0]+data[-1])/2
    x=xrange(len(data))
    mid_ix=len(data)//2
    print mid_ix
    paddata=np.array([data[2],data[1]]+list(data)+[data[-2],data[-3]])
    padx=np.arange(len(paddata))-2
    f=interp1d(padx,paddata,kind='cubic')
    ipodense=20
    xipol=np.linspace(padx[0], padx[-1], num=len(padx)*ipodense, endpoint=True)

    # print len(f(xipol))
    ax.plot(xipol[2*ipodense:-2*ipodense],f(xipol)[2*ipodense:-2*ipodense],c='r',ls='--')
    ax.scatter(x,data,s=70,c='royalblue')
    ax.set_ylabel(r"$\mathrm{Energy[eV]}$")
    ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off

    return ax
