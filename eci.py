import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import plot

def json_eci_to_pandas(jsoneci):
    """Go through an eci.json file and extract all
    the relevant information into a pandas data frame
    with

    eci
    index
    size
    max_length
    min_length
    prototype_function
    multiplicity

    :jsoneci: json of eci.json
    :returns: Pandas DataFrame

    """
    data=jsoneci

    eci_df=pd.DataFrame()
    for func in data["cluster_functions"]:
        eci=0.0
        if "eci" in func:
            eci=func["eci"]

        index=func["linear_function_index"]
        size=len(func["prototype"]["sites"])
        min_length=func["prototype"]["min_length"]
        max_length=func["prototype"]["max_length"]
        prototype_function=func["prototype_function"]
        multiplicity=func["mult"]

        data_dict={"eci":eci,
            "size":size,
            "min_length":min_length,
            "max_length":max_length,
            "prototype_function":prototype_function,
            "index":index,
            "multiplicity":multiplicity}

        eci_df=eci_df.append(data_dict,ignore_index=True)

    return eci_df

def squash_eci(pdeci, max_value=0.0):
    """Run through the values of the eci and return
    a view with eci values that are of greater magnitude
    that the provided maximum value

    :pdeci: pandas DataFrame of eci values
    :max_value: float
    :returns: pandas DataFrame

    """
    return pdeci.loc[pdeci["eci"].abs()>max_value]


class ECI(object):

    """Container for the eci.json file created by CASM.
    Stores values as pandas DataFrame in addition to the
    json format"""

    def __init__(self, filename):
        """Constructs the object, saving all the data from the json
        file, making a pandas representation, and setting different
        parameters for plotting, etc.

        :filename: eci.json path

        """
        self._filename = filename

        datadump=open(filename).read()
        data=json.loads(datadump)

        self._jsoneci=data

        #Store as pandas, but discard eci with values of 0.0
        self._pdeci=squash_eci(json_eci_to_pandas(self._jsoneci))
        self._pdeci=self._pdeci.reset_index()

    def cluster_size(self, numsites):
        """Return view into data for clusters with the
        given number of sites

        :numsites: int
        :returns: pandas DataFrame

        """
        return self._pdeci.loc[self._pdeci["size"]==numsites]

    def plot_values(self, ax, clustersizes, normalized=False):
        """Make a bar graph for the clusters specified in
        the list of sizes. Uses colors from the current color
        cycler of your mpl axes.

        :clustersizes: int
        :ax: matplotlib axes
        :normalized: bool, divides by multiplicity if true
        :returns: matplotlib axes

        """
        prop=ax._get_lines.prop_cycler

        for param,s in zip(prop,clustersizes):
            subdata=self.cluster_size(s)
            if not normalized:
                plottable=subdata["eci"]
            else:
                plottable=subdata["eci"]/subdata["multiplicity"]

            ax.bar(subdata.index,plottable,color=param['color'],label=plot.plot_label_for_cluster_size(s))

        ax.legend(fontsize=14)

        if not normalized:
            ax.set_ylabel(plot.texbf("ECI"))
        else:
            ax.set_ylabel(plot.texbf("ECI/multiplicity"))

        return ax


def main():
    eci_data=query_eci("./eci.json")
    eci_data=eci_data.loc[eci_data["eci"]!=0.0]
    eci_data=eci_data.reset_index()


    point=eci_data.loc[eci_data["size"]==1]
    pairs=eci_data.loc[eci_data["size"]==2]
    triple=eci_data.loc[eci_data["size"]==3]
    quadruple=eci_data.loc[eci_data["size"]==4]

    eci=eci_data["eci"]
    index=eci_data["index"]
    ind=np.arange(len(point)+len(pairs)+len(triple)+len(quadruple))

    fig,ax=plt.subplots()

    ax.bar(point.index,point["eci"],color='#e66101',label=r"\textbf{Point}")
    ax.bar(pairs.index,pairs["eci"],color='#5e3c99',label=r"\textbf{Pairs}")
    ax.bar(triple.index,triple["eci"],color='#b2abd2',label=r"\textbf{3 body}")
    ax.bar(quadruple.index,quadruple["eci"],color='#fdb863',label=r"\textbf{4 body}")

    ax.legend(fontsize=14)
    ax.set_ylabel(r"\textbf{ECI [eV/atom]}")


    plt.tight_layout()
    plt.show()

