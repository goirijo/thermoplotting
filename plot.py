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

