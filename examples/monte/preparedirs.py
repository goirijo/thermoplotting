from __future__ import print_function
import pandas as pd
import json
import os
import os.path
import numpy as np
import thermoplotting as tp
import string
from thermoplotting.casmfiles.monte import ConditionsDomain, prepare_streak_directories


def main():
    #Round my shit yo
    from json import encoder
    encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    mua_range = (-2.5, -0.25, 0.05)
    mub_range = (-2.5, 0.5, 0.05)
    T_range = (1173, 1173, 0)

    domain = ConditionsDomain([mua_range, mub_range], T_range)

    #print domain._domain

    motifdict = {
        "a": ['SCEL1_1_1_1_0_0_0/2', 'SCEL4_2_2_1_1_1_0/6'],
        "b": ['SCEL1_1_1_1_0_0_0/2']
    }

    settingsdump = tp.json_from_file("./monte_settings_template.json")
    settings = tp.casmfiles.Monte(settingsdump)

    for streak in motifdict:
        for motif in motifdict[streak]:
            print("Working on streak {} with motif {}".format(streak, motif))
            prepare_streak_directories(settings, domain, streak, motif)


if __name__ == "__main__":
    main()
