import pandas as pd
import numpy as np
import thermoplotting as tp
import matplotlib.pyplot as plt
import casm
from scipy.spatial import ConvexHull

def data_from_project(project_path, queryargs):
    # Create the casm project
    system=casm.project.Project(project_path)

    #Get a list of all the calculated configurations for that project
    confignames=tp.calculated_confignames(proj=system)

    #Query the project for data
    data=tp.casm_query(confignames, queryargs, proj=system)

    #insert path into queried data
    data["project_path"]=project_path

    return data

def main():
    #Paths to different casm projects
    paths=["../../../NiAlCr","../../../NiAlCr-B2-all","../../../CrAlNi-BCC"]
    #Values to query for
    queryargs=["configname", "atom_frac", "formation_energy_per_atom","hull_dist(CALCULATED,atom_frac)"]
    #store queries for each casm project
    datalist=[data_from_project(p,queryargs) for p in paths]

    #Put all the data together into a single pandas DataFrame
    alldata=pd.concat(datalist)
    #Store data that we want to create convex hull out of in numpy array
    hulldata=alldata[["atom_frac(Ni)","atom_frac(Cr)","formation_energy_per_atom"]].as_matrix()
    #Create convex hull out of all the projects
    global_hull=ConvexHull(hulldata)

    #Get distances from hull for all the points we made the convex hull with
    distances=tp.systems.hull.distances_from_hull(hulldata,global_hull)
    #Put the calculated distances back with the rest of the data
    alldata["global_hull_dist"]=distances

    #The casm hull distance will match only some of the global hull distance, due to the
    #introductions of additional stable structures that a single casm project is not aware of

    return
    

if __name__=="__main__":
    main()
