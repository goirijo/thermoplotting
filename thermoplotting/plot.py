def texbf(string):
    """Format the given string so that
    it becomes bold when plotting:

    r"\\textbf{string}"

    :string: str
    :returns: str

    """
    return r"\textbf{"+string+"}"

def texmrm(string):
    """Format given string into mathmode
    with mathrm

    r"\\$\mathrm{string}$"

    :string: str
    :returns: str

    """
    return r"$\mathrm{"+string+"}$"

def plot_label_for_cluster_size(clustersize):
    """Return string for the legend of clusters of the
    given size

    :clustersize: int
    :returns: str

    """
    if clustersize==0:
        return texbf("Empty")

    elif clustersize==1:
        return texbf("Point")

    elif clustersize==2:
        return texbf("Pairs")

    else:
        return texbf(str(clustersize)+" body")
    
