import numpy as np
import pandas as pd
import casm.project
import os
import hashlib
import json
import re
from scipy.spatial.distance import squareform, pdist
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def latmat(poscar):
    """Extract the lattice from a POSCAR file
    and return as column vector matrix

    :poscar: string
    :returns: ndarray

    """
    with open(poscar) as f:
        name=f.readline()
        scale=f.readline()

        vector=f.readline()
        a=vector.split()
        vector=f.readline()
        b=vector.split()
        vector=f.readline()
        c=vector.split()
        
        lat=np.array([a,b,c]).T

    return float(scale)*lat.astype(float)

def latvol(lattice):
    """Compute volume of lattice

    :lattice: ndarray
    :returns: float

    """
    a=lattice[:,0]
    b=lattice[:,1]
    c=lattice[:,2]

    return np.dot(a,np.cross(b,c))


def facet_intercept(facet):
    """Determine to intercept of the plane defined
    by a facet with the axis that define the composition
    space.

    :facet: ndarray
    :returns: array

    """
    pass

def indexed_view(indx, ndarray, axis):
    """Return last view along specified axis
    E.g. ndarray[:,:,indx,:] if axis=2

    :indx:
    :ndarray: 
    :axis:
    :returns: ndarray

    """
    n=ndarray.ndim
    pretuple=[slice(None)]*ndarray.ndim
    pretuple[axis]=indx
    return ndarray[tuple(pretuple)]

def top_view(ndarray, axis):
    """Return last view along specified axis
    E.g. ndarray[:,:,-1,:] if axis=2

    :ndarray: 
    :axis:
    :returns: ndarray

    """
    return indexed_view(-1,ndarray,axis)

def bottom_view(ndarray, axis):
    """Return 0th view along specified axis
    E.g. ndarray[:,:,0,:] if axis=2

    :ndarray: 
    :axis:
    :returns: ndarray

    """
    return indexed_view(0,ndarray,axis)

def reverse(ndarray, axis, inplace=False):
    """Invert order of N-dimensional array
    along the given axis
    
    :ndarray:
    :axis: int
    :inplace: bool
    :returns: ndarray
    """
    resorted=np.swapaxes(ndarray,axis,0)
    resorted=resorted[::-1,...]
    resorted=np.swapaxes(resorted,0,axis)

    #views into ndarray will be updated
    if inplace:
        ndarray[...]=resorted[...]
        return ndarray

    #views into ndarray are not views into resorted
    else:
        return resorted

def stack(arraylist,axis=0):
    """Take a bunch of arrays with matching dimensions
    and stack them along a new axis.

    :arraylist: [ndarray]
    :axis: int
    :returns: ndarray

    """
    candidates=[np.expand_dims(data,axis=0) for data in arraylist]
    stack=np.concatenate(candidates,axis=0)
    return stack
    
def amin_stack(arraylist):
    """Given a list of arrays of matching shapes, return
    a new array of the same shape that has all the minimum
    values from the list.

    :arraylist: array of ndarray
    :returns: ndarray

    """
    candidates=[np.expand_dims(data,axis=0) for data in arraylist]
    stack=np.concatenate(candidates,axis=0)
    return np.amin(stack,axis=0)

def argmin_stack(arraylist):
    """Given a list of arrays of matching shapes, return
    a new array of the same shape that has all the minimum
    indexes from the list.

    :arraylist: array of ndarray
    :returns: ndarray

    """
    candidates=[np.expand_dims(data,axis=0) for data in arraylist]
    stack=np.concatenate(candidates,axis=0)
    #return np.argmin(stack,axis=0)
    return np.nanargmin(stack,axis=0)

def nearest_entry(data,values,columns=None):
    """From an array, find the entry that most closely
    matches the given values. E.g. try to find the entry
    for z with the closest given x,y values.

    Do not attempt on anything more complicated than a 2D table.


    :data: ndarray
    :values: tuple(float), entry values
    :columns: tuple(int), entry columns
    :axis: int
    :returns: int

    """
    if columns is None:
        columns=np.arange(len(values))

    if len(columns)!=len(values):
        raise ValueError("There must be one column index per column value!")
    
    zerosum=np.zeros(len(data))

    for val,ind in zip(values,columns):
        zerosum+=abs(data[:,ind]-val)

    return np.argmin(zerosum)

def configlist_path(proj=None):
    if proj==None:
        proj=casm.project.Project()

    casmroot=proj.path
    return os.path.join(casmroot,".casm","config_list.json")

def shatter_path(pathname):
    """Run os.path.split until there's nothing left to split
    and return all the components as a list

    :pathname: str
    :returns: list of str

    """
    pathname=os.path.normpath(pathname)
    return pathname.split(os.sep)


def casm_query(configlist,queryargs,proj=None):
    """Inefficiently query casm by specifying a list of confignames
    rather than a selection path. Assumes you want the project you're
    sitting in.
    Does not guarantee keeping the order of the configlist!

    :configlist: list of str
    :queryargs: list of str (casm query --column arguments)
    :returns: Pandas DataFrame

    """
    if proj==None:
        proj=casm.project.Project()

    requestednames=("configname" in queryargs)

    if not requestednames:
        queryargs.append("configname")
    
    testselection=simulated_selection(proj,configlist)
    subquery=casm.project.query(proj,queryargs,selection=simulated_selection(proj,configlist))

    if not requestednames:
        subquery=subquery.drop("configname",1)
    
    #return subquery.convert_objects(convert_numeric=True)   #Things get stupidly stored as strings without this
    return subquery.apply(pd.to_numeric, errors="ignore") #Things get stupidly stored as strings without this

def simulated_selection(proj,configlist):
    """Simulate a selection file by specifying a list of configurations
    that you want selected

    :configlist: list of str
    :returns: StringIO

    """
    selectionstr="#configname selected\n"
    for c in configlist:
        selectionstr+=c
        selectionstr+=" 1\n"

    filename=hashlib.md5(selectionstr).hexdigest()+"_selection.txt"
    filepath=os.path.join(os.sep, "tmp",filename)

    selectionfile=open(filepath, 'w')
    selectionfile.write(selectionstr)

    return casm.project.Selection(proj, filepath)

def scrub_query(queried_data):
    """Drop all rows that contain an "unknown" value

    :queried_data: Pandas DataFrame
    :returns: Pandas DataFra

    """
    subset=queried_data
    for col in queried_data:
        if subset[col].dtype==object:
            subset=subset[queried_data[col]!="unknown"]                              
    return subset.apply(pd.to_numeric, errors="ignore")

def all_confignames(proj=None):
    """Get list of all the confignames
    
    :returns" Pandas list of str
    """
    if proj==None:
        proj=casm.project.Project()

    configdump=casm.project.query(proj, ["configname"], all=True)

    return configdump["configname"]

def selected_confignames(proj=None):
    """Ger list of the confignames that are currently selected
    (i.e. "MASTER")

    :proj: casm.project.Project
    :returns: pd list of str

    """
    if proj==None:
        proj=casm.project.Project()

    master=casm.project.Selection(proj,"MASTER",False)
    return casm.project.query(proj,["configname"],master)["configname"]
    

def calculated_confignames(invert=False,proj=None):
    """Get list of all the confignames that have is_calculated==1
    
    :invert: bool
    :returns: list of str
    """
    if proj==None:
        proj=casm.project.Project()

    configdump=casm.project.query(proj, ["configname", "is_calculated"], all=True)

    if not invert:
        subset=configdump.loc[configdump["is_calculated"]==1]
    else:
        subset=configdump.loc[configdump["is_calculated"]!=1]

    return subset["configname"]

def confignames_of_size(sizes,proj=None):
    """Get list of the confignames that have a supercell
    of the specified sizes

    :sizes: list of int
    :returns: list of str

    """
    if proj==None:
        proj=casm.project.Project()

    configdump=casm.project.query(proj,["configname"],all=True)

    sizes=["SCEL"+str(s)+"_" for s in sizes]
    subset=[config for config in configdump["configname"] if any(size in config for size in sizes)]
    return subset


def unique_rows(arr, thresh=0.0, metric='euclidean'):
    """Returns subset of rows that are unique, in terms of Euclidean distance
    http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    """
    distances = squareform(pdist(arr, metric=metric))
    idxset = {tuple(np.nonzero(v)[0]) for v in distances <= thresh}

    return arr[[x[0] for x in idxset]]

def unit_vector(v):
    """Normalize a vector so that it has unit length.

    :v: np array
    :returns: np array

    """
    return v/np.linalg.norm(v)

def angle_between(v1, v2):
    """Compute the angle between two vectors. Not smart for signed angles.

    :returns: float
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def json_from_file(filename):
    """Load the given path to file and return json object

    :filename: path
    :returns: json

    """
    with open(filename) as data_file:
        data=json.load(data_file)
    return data

def json_to_file(data,filename):
    """Save the json data to a file

    :data: json
    :filename: path
    :returns: void

    """
    with open(filename, 'w') as f:
        json.dump(data,f,indent=4,sort_keys=True)

    return


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

def set_np_print_format():
    """Make numpy print 8 siginificant figures

    :returns: np module

    """
    np.set_printoptions(precision=8,suppress=True,formatter={'float': '{: 0.8f}'.format})
    return np

def natural_sort_key(l):
    """Sort list numerically by the integers in the string, so that
    values like corr(0), corr(1), corr(2), etc come out in the
    right order.

    Use this function as they key when sorting, e.g.:
    mylist.sort(key=natural_sort_key)

    :l: list of str
    :returns: list of str

    """
    _nsre = re.compile('([0-9]+)')
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, l)]   

def sequential_columns(df, colname):
    """Extract and sort the columns from a pandas dataframe
    containing a sequence of columns, such as
    corr(0), corr(1), corr(2)... etc

    :df: pd DataFrame
    :colname: str
    :returns: list of str

    """
    allcols=df.columns.values
    cols=[c for c in allcols if re.compile(colname+"\([0-9]+\)").search(c)]

    return sorted(cols,key=natural_sort_key)
    
def sequential_dataframe(array, colname):
    """Construct a pandas DataFrame from a numpy array
    with a sequential naming of columns, such as
    corr(0), corr(1), corr(2)... etc

    :array: np
    :colname: str
    :returns: pd DataFrame

    """
    df=pd.DataFrame({"{}({})".format(colname,ix):col for ix,col in enumerate(array.T)})
    return df[sequential_columns(df, colname)]

def dataframe_from_numpy(npdata, colnames):
    """Create a pandas DataFrame from a set of numpy data that
    has the specified column names. Only for 2D data at the
    moment.

    :npdata: np array
    :colnames: list of str
    :returns: pd DataFrame

    """
    if len(colnames)!=len(npdata.T):
        raise ValueError("Number of columns does not match the number of column names!")

    return pd.DataFrame({name:col for name,col in zip(colnames,npdata.T)})

def configmatch(match_to_data, match_from_row, match_values, eps=0.0000001):
    """Find a particular configuration in a list of configurations by checking particular
    values, such as composition, formation energy, etc. Returns the index in the configuration
    list or None

    :match_to_data: pd DataFrame, must have match_values as columns
    :match_from_row: row from pd DataFrame, must have match_values as columns
    :match_values: columns that should be compared to determine if the configuration matches
    :eps: tolerance value
    :returns: list of int

    """
    tmp_col="_____________________________tmp"

    #If your data frame has a columns with such a ridiculous name, I can't deal with it right now...
    if tmp_col in match_to_data.columns.values:
        raise ValueError("Why in the world do you have a column named {}???".format(tmp_col))

    match_to_data[tmp_col]=0.0
    diffsum=match_to_data[tmp_col]

    for m in match_values:
        diffsum+=abs(match_to_data[m]-match_from_row[m])

    ixmatch=diffsum.loc[diffsum<eps].index

    # if len(ixmatch)!=1:
    #     print ixmatch
    #     print '----------------------'
    #     print match_from_row
    #     print '----------------------'
    #     print match_to_data.sort_values(tmp_col)
    #     exit()

    del match_to_data[tmp_col]
    return ixmatch

