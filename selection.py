import pandas as pd
import numpy as np

def generate_table(confignames):
    """Create a DataFrame with a "selected"
    column

    :confignames: list of str
    :returns: pd DataFrame

    """
    selectrow=np.ones(len(confignames))

    selectpd=pd.DataFrame({"configname":confignames,"selected":selectrow})
    return selectpd

class Selection(object):

    """Holds a list of confignames and knows
    how to write them to a selection file"""

    def __init__(self, confignames):
        """Saves an initial list of configuration
        names."""
        
        self._confignames=confignames
        self._pdrep=generate_table(self._confignames)

    def _set_hashmode(self, hash_on):
        """Name the "configname" column with a leading
        '#' or not

        :hash_on: bool
        :returns: pd DataFrame

        """
        cols=self._pdrep.columns

        if hash_on and "configname" in cols:
            self._pdrep=self._pdrep.rename(columns={"configname":"#configname"})

        elif not hash_on and "#configname" in cols:
            self._pdrep=self._pdrep.rename(columns={"#configname":"configname"})

        return self._pdrep

    def to_csv(self, filename, dumbhash=True):
        """Prints the selection in csv format

        :dumbhash: bool (should header include '#'?)
        :filename: str (path to output file)
        :returns: void

        """
        self._pdrep=self._set_hashmode(dumbhash)
        self._pdrep.to_csv(filename,sep=' ',index=False)

        return
