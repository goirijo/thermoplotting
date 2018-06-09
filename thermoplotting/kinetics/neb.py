import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def simple_plot(ax, data):
    """Given a few neb image energies, interpolate
    the data and plot it

    Parameters
    ----------
    ax : TODO
    data : TODO

    Returns
    -------
    TODO

    """
    x=xrange(len(data))
    f=interp1d(x,data,kind='cubic')
    xipol=np.linspace(x[0], x[-1], num=91, endpoint=True)

    ax.plot(xipol,f(xipol),c='r',ls='--')
    ax.scatter(x,data,s=70,c='royalblue')
    ax.set_ylabel(r"$\mathrm{Energy[eV]}$")
    ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off

    return ax
