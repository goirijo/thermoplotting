from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import object

import pandas as pd
import numpy as np
import json
from ..misc import all_confignames, json_to_file
import string
import copy
import os
import glob


def generate_condition_from_dict(cond_dict):
    """Given a dictionary of settings values, generate a condition
    for a Monte Carlo run. The keys associated with the parametric
    chemical potential are single characters, while temperature is
    "T" and tolerance is "tol"

    :cond_dict: dict
    :returns: dict

    """
    pass


def generate_condition(pots, temperature, tolerance=0.001):
    """Create a json object that contains the information for one
    of the custom conditions, so that you can stick it into a list

    :pots: list of float (length should match number of compositions)
    :temperature: float
    :tolerance: float
    :returns: dict
    """
    settings = dict()
    settings["param_chem_pot"] = dict()
    for ix, mu in enumerate(pots):
        label = string.ascii_lowercase[ix]
        settings["param_chem_pot"][label] = mu
    settings["temperature"] = temperature
    settings["tolerance"] = tolerance

    return settings


def generate_incremental_condition_set(inits, finals, incs):
    """Create a dict of three conditions: initial, incremental, and
    final, by passing three sets of arguments you would pass to
    generate_condition

    :inits: pots,temperature,tolerance
    :incs: pots,temperature,tolerance
    :finals: pots,temperature,tolerance
    :returns: dict

    """
    settings_combo = {}
    slots = ("initial_conditions", "final_conditions", "incremental_conditions")

    for slot, sett in zip(slots, (inits, finals, incs)):
        settings[slot] = generate_condition(*sett)

    return settings_combo


def dirname_from_condition(condition,
                           ignore=[],
                           mu_format="{0:.2f}",
                           T_format="{0:.1f}"):
    """Given a particular Monte Carlo settings condition, generate
    a name for a directory with format mua_$mua__mub_$mub__T_$T

    :condition: dict
    :ignore: list of str (suppress these conditions from the directory name)
    :returns: str

    """
    dirname = ""
    pots = condition["param_chem_pot"]
    for c in pots:
        if c in ignore:
            continue
        else:
            mustr = "mu{}".format(c)
            dirname += mustr + "_" + mu_format.format(pots[c]) + "__"

    if "T" not in ignore:
        dirname += "T_" + T_format.format(condition["temperature"])

    return dirname

def streak_dirs(streakname, groundname, jsonfile):
    pathname=os.path.join(".",streakname,groundname,"mu*_*__T_*",jsonfile)
    return glob.glob(pathname)


# def direname_from_incremental_driver(driver_settings):
#     """Create a directory name as specified by dirname_from_condition, but suppress
#     the settings that have increment

#     :driver_settings: dict
#     :returns: str

#     """
#     ignore=[]

#     initials=driver_settings["initial_conditions"]
#     incs=driver_settings["incremental_conditions"]

#     inctol=incs["tolerance"]

#     pots=incs["param_chem_pot"]
#     for c in pots:
#         if abs(pots[c])<inctol:
#             ignore.append(c)

#     if abs(incs["temperature"])<inctol:
#         ignore.append("T")

#     return dirname_from_condition(initials,ignore)


def atom_frac_cols(fracs):
    """Return columns corresponding to the average atom fraction

    :fracs: list str
    :returns: list str

    """
    return ["<atom_frac({})>".format(f) for f in fracs]


def comp_cols(comps):
    """Return columns corresponding to the average composition

    :comps: str
    :returns: list str

    """
    return ["<comp({})>".format(c) for c in comps]


def chem_pot_cols(chems):
    """Return columns corresponding to the chemical potential

    :chems: str
    :returns: list str

    """
    return ["param_chem_pot({})".format(c) for c in chems]


def singular_streak_along_axis(settings, fixed_cond, axis, axis_range, tols):
    """Edit the conditions for initial, final and incremental conditions
    at a single set of fixed values, except for the axis

    :settings: thermoplotting.casmfiles.monte.Monte
    :fixed_cond: dict of the fixed conditions
    :axis: str, either a letter for parametric chemical potential or "T"
    :axis_range: initial,final,increment values for the axis
    :tols: tolerances for initial,final,increment
    :returns: dict

    """
    fixed_cond[axis] = axis_range[0]
    settings.set_initial_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[0])

    fixed_cond[axis] = axis_range[1]
    settings.set_final_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[1])

    for k in fixed_cond:
        if k in string.ascii_lowercase:
            fixed_cond[k] = 0.0
    fixed_cond[axis] = axis_range[2]
    settings.set_incremental_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[2])

    return settings.settings()


def streak_along_axis(settings, domain, axis, tols=(0.001, 0.001, 0.001)):
    """Choose a direction in thermodynamic space that the Monte Carlo simulation
    is going to run in, and construct a grid of values for all the other variables.
    Return settings that keep every condition on the grid constant, but vary
    along the specified axis.

    For example, for a ternary system, a settings object could be created for every
    possible combination of T and mu_a, keeping those values fixed, but having the
    conditions of mu_b vary from the initial to the final values

    :settings: thermoplotting.casmfiles.monte.Monte
    :domain: thermoplotting.casmfiles.monte.SettingsDomain
    :axis: str ('a', 'b', ...., 'T')
    :returns: list of dict

    """
    if settings.driver_mode() != "incremental":
        raise UserWarning(
            "The settings passed for creating a streak are NOT incremental!")
    if axis not in domain._axes:
        raise ValueError(
            "The specified axis '{}' is not in bounded by the condition domain.".
            format(axis))

    non_streak_axes = [x for x in domain._axes if x != axis]
    fixed_cond_vals = domain.grid(non_streak_axes)

    dict_fix = fixed_cond_vals.to_dict(orient='records')
    streak_range = domain._ranges[axis]

    settings_streak = [
        singular_streak_along_axis(settings, fixed, axis, streak_range, tols)
        for fixed in dict_fix
    ]
    return settings_streak


class Monte(object):
    """Construct a json file for the `casm monte` command. Currently requires
    a preexisting json object to start with that has all the fieds, letting you
    edit them."""

    def _default_base(self):
        """Create a default json file with all the required fields filled up
        :returns: json

        """
        raise NotImplementedError(
            "Default json values for monte are not currently implemented")
        return {}

    def __init__(self, base):
        """Initialize with a starting point

        :base: dict

        """
        self._json = base

    def set_motif(self, configname, check_value=True):
        """Set the initial configuration of the simulation

        :configname: str
        :returns: json

        """
        if check_value and configname!="restricted_auto":
            all_configs = all_confignames()
            for c in all_configs:
                if not any(all_configs == configname):
                    raise ValueError(
                        "The requested motif is not an existing configuration!")

        self._json["driver"]["motif"]["configname"] = configname
        return self._json

    def set_driver_mode(self, mode):
        """Set the driver mode ("custom" or "incremental")

        :mode: str
        :returns: json

        """
        if mode not in ("custom", "incremental"):
            raise ValueError(
                "The requested mode is neither 'custom' or 'incremental'")

        self._json["driver"]["mode"] = mode
        return self._json

    def driver_mode(self):
        """Set the driver mode ("custom" or "incremental")

        :mode: str
        :returns: json

        """
        return self._json["driver"]["mode"]

    def _set_condition_type(self, typeid, pots, temperature, tol):
        """Use to set either the initial, incremental or final conditions

        :typeid: str
        :pots: list of float (length should match number of compositions)
        :temperature: float
        :tolerance: float
        :returns: json

        """
        assert (typeid in ["initial", "incremental", "final"])
        cond = generate_condition(pots, temperature, tol)
        self._json["driver"][typeid + "_conditions"] = cond

        return self._json

    def set_initial_conditions(self, pots, temperature, tol=0.001):
        """Given a list of chemical potentials and a temperature, set
        the initial conditions for the simulation

        :pots: list of float (length should match number of compositions)
        :temperature: float
        :tolerance: float
        :returns: json

        """
        self._json = self._set_condition_type("initial", pots, temperature, tol)
        return self._json

    def set_final_conditions(self, pots, temperature, tol=0.001):
        """Given a list of chemical potentials and a temperature, set
        the final conditions for the simulation

        :pots: list of float (length should match number of compositions)
        :temperature: float
        :tolerance: float
        :returns: json

        """
        self._json = self._set_condition_type("final", pots, temperature, tol)
        return self._json

    def set_incremental_conditions(self, pots, temperature, tol=0.001):
        """Given a list of chemical potentials and a temperature, set
        the incremental conditions for the simulation

        :pots: list of float (length should match number of compositions)
        :temperature: float
        :tolerance: float
        :returns: json

        """
        self._json = self._set_condition_type("incremental", pots, temperature,
                                              tol)
        return self._json

    def set_simulation_cell(self, transf_mat):
        """Set the size and shape of the simulation cell via a transformation matrix
        relative to the primitive cell

        :transf_mat: np array
        :returns: json

        """
        serializable = transf_mat.tolist()
        # serializable=transf_mat.T.tolist()        #WHICH IS IT????
        self._json["supercell"] = serializable
        return self._json

    def to_json(self, filename):
        """Write the settings out to a file

        :filename: str
        :returns: void

        """
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4, sort_keys=True)
        return

    def settings(self):
        """Return copy of the current state of the settings
        :returns: json

        """
        return copy.deepcopy(self._json)


class ConditionsDomain(object):
    """Given initial, final and incremental values, store
    the domain of the different chemical potentials"""

    def _arange(self, initial, final, increment):
        """Exacly like np.arange, but if the increment
        is zero, return a single value, and always include
        the final value in the range of increments.

        :initial: float
        :final: float
        :increment: float
        :returns: list of float

        """
        #Because the MC runs will go from initial *into* the final conditions
        #(final conditions are calculated), you must include the end of the
        #specified ranges in the grid
        if increment == 0.0:
            return [initial]
        else:
            return np.arange(initial, final + increment, increment)

    def __init__(self, mu_domains, T_domain):
        """Initialize by specifying the edges of each chemical
        potential domain, and the desired grid resolution

        :mu_domains: list of tuple (initial,final,increment), in same order as composition axis (a,b,c...)
        :T_domain: initial,final,increment of temperature

        """
        labels = string.ascii_lowercase[0:len(mu_domains)] + "T"
        all_domains = mu_domains + [T_domain]

        self._ranges = {labels[ix]: r for ix, r in enumerate(all_domains)}

        initials, finals, increments = list(zip(*all_domains))
        self._initials = list(initials)
        self._finals = list(finals)
        self._increments = list(increments)

        self._domain = {}
        for ix, (ini, fin, inc) in enumerate(
                zip(self._initials, self._finals, self._increments)):
            l = labels[ix]
            self._domain[l] = {}
            self._domain[l] = self._arange(ini, fin, inc)

        self._axes = [k for k in self._domain]

    def grid(self, axes):
        """Grid up the conditions for all the specified axes and return
        all combined values

        :axes: list of str
        :returns: pandas DataFrame

        """
        mesh = np.meshgrid(*(self._domain[k] for k in axes))
        cond_combos = np.reshape(mesh, (len(axes), -1))
        return pd.DataFrame({k: cond_combos[ix] for ix, k in enumerate(axes)})

    def axes(self):
        return self._axes[:]


def singular_streak_along_axis(settings, fixed_cond, axis, axis_range, tols):
    """Edit the conditions for initial, final and incremental conditions
    at a single set of fixed values, except for the axis

    :settings: thermoplotting.casmfiles.monte.Monte
    :fixed_cond: dict of the fixed conditions
    :axis: str, either a letter for parametric chemical potential or "T"
    :axis_range: initial,final,increment values for the axis
    :tols: tolerances for initial,final,increment
    :returns: dict

    """
    fixed_cond[axis] = axis_range[0]
    settings.set_initial_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[0])

    fixed_cond[axis] = axis_range[1]
    settings.set_final_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[1])

    #All increments except the axis should be zero
    for k in fixed_cond:
        fixed_cond[k] = 0.0
    fixed_cond[axis] = axis_range[2]
    settings.set_incremental_conditions(
        [fixed_cond[k] for k in string.ascii_lowercase
         if k in fixed_cond], fixed_cond["T"], tols[2])

    return settings.settings()


def streak_along_axis(settings, domain, axis, tols=(0.001, 0.001, 0.001)):
    """Choose a direction in thermodynamic space that the Monte Carlo simulation
    is going to run in, and construct a grid of values for all the other variables.
    Return settings that keep every condition on the grid constant, but vary
    along the specified axis.

    For example, for a ternary system, a settings object could be created for every
    possible combination of T and mu_a, keeping those values fixed, but having the
    conditions of mu_b vary from the initial to the final values

    :settings: thermoplotting.casmfiles.monte.Monte
    :domain: thermoplotting.casmfiles.monte.SettingsDomain
    :axis: str ('a', 'b', ...., 'T')
    :returns: list of dict

    """
    if settings.driver_mode() != "incremental":
        raise UserWarning(
            "The settings passed for creating a streak are NOT incremental!")
    if axis not in domain._axes:
        raise ValueError(
            "The specified axis '{}' is not in bounded by the condition domain.".
            format(axis))

    non_streak_axes = [x for x in domain._axes if x != axis]
    fixed_cond_vals = domain.grid(non_streak_axes)

    dict_fix = fixed_cond_vals.to_dict(orient='records')
    streak_range = domain._ranges[axis]

    settings_streak = [
        singular_streak_along_axis(settings, fixed, axis, streak_range, tols)
        for fixed in dict_fix
    ]
    return settings_streak


def prepare_streak_directories(settings,
                               domain,
                               streak_axis,
                               motif,
                               settings_filename="monte_settings.json",
                               streak_dirname=None,
                               motif_dirname=None,
                               target=".",
                               **kwargs):
    """
    Given starting settings and a domain for the conditions, select a thermodynamic variable as an axis
    through which to create the incremental conditions. A settings file is generated for every gridpoint for a
    run that holds those grid points constant, but increments along the specified axis

    :settings: Monte
    :domain: ConditionsDomain
    :streak_axis: str ('a', 'b', ... or 'T')
    :motif: str (configname)
    :settings_filename: str (name for the final json settings file)
    :streak_dirname: str (name for the top level directory, default streak_{}.format(streak_axis))
    :motif_dirname: str (name for directory inside streak_direname, default is the configname)
    :target: str (name of directory to work from)
    :returns: void

    """
    if streak_dirname is None:
        streak_dirname = "streak_{}".format(streak_axis)

    if motif_dirname is None:
        motif_dirname = motif.replace('/', '.')

    basedir = os.path.join(target,streak_dirname, motif_dirname)

    settings.set_motif(motif)
    streaks = streak_along_axis(settings, domain, streak_axis)

    for s in streaks:
        streakdir = dirname_from_condition(s["driver"]["initial_conditions"],
                                           [streak_axis],**kwargs)
        targetdir = os.path.join(basedir, streakdir)
        os.makedirs(targetdir)
        json_to_file(s, os.path.join(targetdir, settings_filename))

    return
