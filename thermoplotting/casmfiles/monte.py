import pandas as pd
import json

def generate_condition(a,b,temperature,tolerance=0.001):
    """Create a json object that contains the information for one
    of the custom conditions, so that you can stick it into a list

    :a: float
    :b: float
    :temperature: float
    :tolerance: float
    :returns: json
    """
    settings=dict()
    settings["param_chem_pot"]=dict()
    settings["param_chem_pot"]["a"]=a
    settings["param_chem_pot"]["b"]=b
    settings["temperature"]=temperature
    settings["tolerance"]=tolerance

    return settings

