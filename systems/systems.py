import pandas as pd
import pickle
import os
from thermoplotting.misc import casm_query
import thermoplotting.misc
import casm

class GenericSystem(object):

    """Class for a container that holds information about
    all the groundstates of a particular system so that you can
    get the full convex hull incorporated into your casm project."""

    def __init__(self, components):
        """Specify the compositions that span your system,
        e.g. Ni, Al, Cr

        :components: list of str
        """
        self._components=list(set(components))
        self._numcomp=len(self._components)

        return


class QueryData(object):

    """Holds numpy DataFrame of different casm
    queries for a single casm project"""

    def current_selection(self):
        """Returns the list of configuration names that are currently
        selected
        :returns: pd DataFrame

        """
        return self._curr_selection

    def current_project(self):
        """Returns the list of configuration names that are currently
        selected
        :returns: casm Project

        """
        return self._casmproj

    def set_selection(self, key):
        """Select one of the available selections for the
        current view.

        :key: str
        :returns: pd DataFrame

        """
        self._curr_selection=self._selections[key]
        return self.current_view()

    def add_selection(self, new_config_list, selection_name, set_selection=True):
        """Add a new set of configurations to choose from

        :new_config_list: list of str
        :selection_name: str
        :set_selection: bool, sets the new selection if True
        :returns: list of str

        """
        if selection_name in self._selections:
            raise ValueError("Name already in list of available selections")

        self._selections[selection_name]=new_config_list

        if set_selection:
            self.set_selection(selection_name)

        return self.current_selection()

    def select_by_sizes(self, sizes):
        """Set the current selection to focus only on supercell sizes
        of that match the given list

        :sizes: list of int
        :returns: list of str

        """
        sized=thermoplotting.misc.confignames_of_size(sizes, proj=self.current_project())
        self._curr_selection=sized

        return self.current_selection()

    def current_view(self, cols=None):
        """Return a view into the current data, but only for
        the currently selected configurations

        :cols: list of str
        :returns: pd DataFrame

        """
        subconfigs=self._running_data[self._running_data["configname"].isin(self.current_selection())]

        if cols is None:
            cols=subconfigs.columns


        return subconfigs[cols]
        

    def __init__(self, project_dir):
        """Saves root directory of casm project and creates
        an initial list of selections

        :project_dir: str

        """
        self._project_dir = os.path.abspath(project_dir)
        self._casmproj=casm.project.Project(path=self._project_dir)

        all_configs=thermoplotting.misc.all_confignames(proj=self._casmproj)
        calc_configs=thermoplotting.misc.calculated_confignames(proj=self._casmproj)

        self._selections={"all":all_configs,"calculated":calc_configs}
        self._curr_selection=self._selections["all"]

        self._running_data=casm_query(self._selections["all"],["configname"],proj=self._casmproj)


    def query(self, query_args):
        """Given a list of query-able arguments, query for them
        and add the entries to the stored data.

        Currently undefined behavior if you query for a value that you already
        queried for, since I'm not totally sure how the merge routine is going
        to work with columns of float.

        :query_args: list of str
        :returns: pd DataFrame

        """
        if "configname" not in query_args:
            query_args.append("configname")


        new_data=casm_query(self.current_selection(), query_args, proj=self.current_project())

        non_conf_cols=[c for c in new_data.columns if c != "configname"]
        new_data[non_conf_cols]=new_data[non_conf_cols].apply(pd.to_numeric, errors='coerce')
        self._running_data=pd.merge(self._running_data,new_data,how='outer')

        ret_cols=non_conf_cols

        if "configname" in query_args:
            ret_cols.append("configname")

        return self.current_view(ret_cols)

    def determine_structure(self, struc_paths_dict, struc_col_name="struc",weight=0.5):
        """Given a dict of paths to some prototype structures, query for the total
        structure score, and then add a new column to the queried data that identifies
        each configuration by name. The name of each prototype structure are the keys
        for each prototype structure path.

        :struc_paths_dict: str to path dictionary
        :struc_col_name: name to give the column identifying structures by type
        :weight: value to weigh basis vs lattice scores
        :returns: pd DataFrame

        """
        def q_arg(path):
            return "struc_score("+os.path.abspath(path)+",total_score,"+str(weight)+")"

        query_args=[q_arg(struc_paths_dict[k]) for k in struc_paths_dict]

        self.query(query_args)

        labelled=pd.DataFrame({k:self._running_data[q_arg(struc_paths_dict[k])] for k in struc_paths_dict})
        self._running_data[struc_col_name]=labelled.idxmin(axis=1)

        return self.current_view(["configname",struc_col_name])

    def to_pickle(self, filename):
        """Dump data to pickle file

        :filename: str
        :returns: void

        """
        pickle.dump(self, open(filename,'wb'))
        return

    @staticmethod
    def from_pickle(filename):
        """Load pickle file to create instance

        :filename: str
        :returns: QueryData

        """
        return pickle.load(open(filename,'rb'))

    def __getstate__(self):
        """When saving the pickle, you gotta destroy
        the casm project, because casm is stupid.
        :returns: dict

        """
        self._casmproj=None
        return self.__dict__

    def __setstate__(self, d):
        """When reading the pickle, reconstruct the casm
        project, because casm is dumb.

        :d: dict
        :returns:

        """
        self.__dict__=d
        self._casmproj=casm.project.Project(path=self._project_dir)
        return


