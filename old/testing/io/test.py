import thermoplotting as tp
import numpy as np
import json


names,arr=tp.thermoio.safe_clobber(["./results.json"])

coldict={"Beta":"b","mu_a":"mu_a","mu_b":"mu_b"}

print names
print arr

names,arr=tp.thermoio.truncate(names,arr,coldict)

print names
print arr
