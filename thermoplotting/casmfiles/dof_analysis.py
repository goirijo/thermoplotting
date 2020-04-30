import numpy as np
import json
import re
import glob
import os

def angle(v1, v2):
    """Return unsigned angle between two 3d vectors

    Parameters
    ----------
    v1 : 1x3 array
    v1 : 1x3 array

    Returns
    -------
    float

    """
    angle = np.arccos(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle * 180 / np.pi


#TODO: Make sure the vectors are formatted the way you think
def make_column_vector_matrix(sym_json):
    """Extract the lattice from the symmetry analysis
    json dictionary, and retun it as a 3x3 matrix,
    with each lattice vector stored in a column

    Parameters
    ----------
    sym_json : json

    Returns
    -------
    3x3 array

    """
    vectors = sym_json["initial_configuration"]["lattice_vectors"]
    return np.array(vectors).T

def make_lengths(column_vector_matrix):
    """Returns the length of each column vector

    Parameters
    ----------
    column_vector_matrix : 3x3 array

    Returns
    -------
    1x3 array

    """
    return [np.linalg.norm(v) for v in column_vector_matrix.T]


def make_angles(column_vector_matrix):
    """Retun the alpha, beta, and gamma angles for the
    specified lattice.

    Parameters
    ----------
    column_vector_matrix : 3x3 array

    Returns
    -------
    1x3 array

    """
    alpha = angle(column_vector_matrix[1], column_vector_matrix[2])
    beta = angle(column_vector_matrix[0], column_vector_matrix[2])
    gamma = angle(column_vector_matrix[0], column_vector_matrix[1])

    return [alpha, beta, gamma]

def make_all_sites_dict(sym_json):
    """Return a dictionary of all sites of the structure, where each
    each key is the site label (for glossary), and the value is a tuple of
    string(species) and 1x3 array (Cartesian position)

    Parameters
    ----------
    sym_json : json

    Returns
    -------
    dict(str:(str,1x3 array))

    """
    sites_dict = sym_json["initial_configuration"]["sites"]
    sites = {entry: (s, np.array(sites_dict[entry][s]))
             for entry in sites_dict
             for s in sites_dict[entry]}
    return sites

def make_symmetry_adapted_axis_names(sym_json):
    """Run through all the symmetry axes and return a list of
    the name of each one

    Parameters
    ----------
    sym_json : json

    Returns
    -------
    list(str)

    """
    return [axis for axis in sym_json["irreducible_representations"]["adapted_axes"]]

def make_selected_sites_dict(sym_json):
    """Return a dictionary of the selected sites of the structure, where each
    each key is the site label (for glossary), and the value is a tuple of
    string(species) and 1x3 array (Cartesian position)

    Parameters
    ----------
    sym_json : json

    Returns
    -------
    dict(str:(str,1x3 array))

    """
    all_sites=make_all_sites_dict(sym_json)
    selected_sites = sym_json["initial_configuration"]["selected_sites"]
    return {s:all_sites[s] for s in all_sites if int(s[0]) in selected_sites}

#TODO: I'm just assuming magnetic spin, so the glossary is like "sx[12]"
#I don't know what other degrees of freedom might look like
def make_glossary_indexer(sym_json):
    """Reads the glossary from the symmetry report, and returns a
    list of pairs that describe the label of the site, and the index
    into the coordinate of the site (x,y, or z).

    The indexer is used to determine the meaning of each component
    of the symmetry adapted axis. For example, if the nth entry of
    the returned indexer is ("5",1), this means that the nth entry
    of the symmetry adapted axes are describing the y component of
    the site labelled "5".

    Parameters
    ----------
    sym_json : json

    Returns
    -------
    list((str,int))

    """
    glossary=sym_json["glossary"]
    xyz_map={"x":0,"y":1,"z":2}
    coord_components=[xyz_map[g[1]] for g in glossary]
    site_labels=[re.search('\[(.*)\]',g).group(1) for g in glossary]
    return [z for z in zip(site_labels,coord_components)]

def make_site_index_map(sites):
    """Returns a map from site label to index.

    Parameters
    ----------
    sites : dict

    Returns
    -------
    dict{str : int}

    """
    labels=[l for l in sites]
    labels.sort()
    return {l:ix for ix,l in enumerate(labels)}
    
def make_adapted_axis_matrix(sym_json,axis,sites=None):
    """Return a reshaped spin symmetry adapted axis, where each column
    has been ordered to represent the x,y,z values for each of the
    given sites.

    Parameters
    ----------
    sym_json : json
    axis : str
    sites : list((str,1x3 array))

    Returns
    -------
    3xlen(sites) array

    """
    if sites is None:
        sites=make_selected_sites_dict(sym_json)

    glossary_indexer=make_glossary_indexer(sym_json)
    axis=sym_json["irreducible_representations"]["adapted_axes"][axis]

    site_index_map=make_site_index_map(sites)
    axis_matrix=np.zeros((3,len(sites)))
    for ix,(label,xyz) in enumerate(glossary_indexer):
        axis_matrix[xyz,site_index_map[label]]=axis[ix]

    return axis_matrix

################################################# Move to class or separate module:

def make_unit_cell_cif_block(lengths,angles):
    """Create the block of cif file that deals with
    the unit cell of the structure. Symmetry is
    hard coded to P1
    
    Example:

    _symmetry_space_group_name_H-M   'P 1'
    _cell_length_a   5.37820988
    _cell_length_b   5.37820988
    _cell_length_c   14.23794427
    _cell_angle_alpha   90.00000000
    _cell_angle_beta   90.00000000
    _cell_angle_gamma   90.00000000

    Parameters
    ----------
    legths : [float,float,float]
    angles : [float,float,float]

    Returns
    -------
    str

    """
    return f"""_symmetry_space_group_name_H-M   'P 1'
    _cell_length_a    {lengths[0]}
    _cell_length_b    {lengths[1]}
    _cell_length_c    {lengths[2]}
    _cell_angle_alpha {angles[0]}
    _cell_angle_beta  {angles[1]}
    _cell_angle_gamma {angles[2]}"""

def make_occupancy_cif_block(column_vector_lattice,sites):
    """Create a block of cif file that deals with the existing
    sites in the structures. Sites should be provided in
    Cartesian coordinates, which will then be converted to
    fractional using the column vector matrix of the lattice.

    Example:

    loop_
     _atom_site_type_symbol
     _atom_site_label
     _atom_site_symmetry_multiplicity
     _atom_site_fract_x
     _atom_site_fract_y
     _atom_site_fract_z
     _atom_site_occupancy
      Ni0+  Ni1  1  0.000000  0.000000  0.000000  1
      Al0+  Al2  1  0.500000  0.500000  0.500000  1

    Parameters
    ----------
    column_vector_lattice : 3x3 array
    sites : dict(str:(str,1x3 array))

    Returns
    -------
    str

    """
    labels=[l for l in sites]
    labels.sort()
    species=[sites[l][0] for l in labels]
    cartesian=np.array([sites[l][1] for l in labels]).T
    fractional=np.dot(np.linalg.inv(column_vector_lattice),cartesian)

    header_block="""loop_
        _atom_site_type_symbol
        _atom_site_label
        _atom_site_symmetry_multiplicity
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy\n"""

    site_lines=[f"{s}0+ {s}{l} 1 {f[0]} {f[1]} {f[2]} 1" for s,l,f in zip(species,labels,fractional.T)]
    sites_block='\n'.join(site_lines)

    return header_block+sites_block

def make_magnetic_spin_cif_block(column_vector_lattice,axis_matrix,sites,cutoff):
    """Create a block of cif file that assigns spin to each site.
    Converts each column of the axis matrix into a magnetic spin,
    and assign it to the appropriate site. Directions with a magnintude smaller
    than the cutoff will be ignored.

    Example:

    loop_
     _atom_site_moment_label
     _atom_site_moment_crystalaxis_x
     _atom_site_moment_crystalaxis_y
     _atom_site_moment_crystalaxis_z
      Mn9  1.03185  -2.49179  -0.00458
      Mn10  -2.49280  -1.03754  -0.00494

    Parameters
    ----------
    column_vector_lattice : 3x3 array
    axis_matrix : 3xlen(sites) array
    sites : dict(str:(str,1x3 array))
    normalization : cutoff

    Returns
    -------
    str

    """
    header_block="""loop_
        _atom_site_moment_label
        _atom_site_moment_crystalaxis_x
        _atom_site_moment_crystalaxis_y
        _atom_site_moment_crystalaxis_z\n"""

    labels=[l for l in sites]
    labels.sort()
    species=[sites[l][0] for l in labels]

    index_map=make_site_index_map(sites)
    fractional_axis_matrix=np.dot(np.linalg.linalg.inv(column_vector_lattice),axis_matrix)
    norms=np.linalg.norm(axis_matrix,axis=0)
    spins=[fractional_axis_matrix.T[index_map[l]] for l in labels]

    spin_lines=[f"{s}{l} {x[0]} {x[1]} {x[2]}" for s,l,x in zip(species,labels,spins)
            if norms[index_map[l]]>cutoff]
    spins_block='\n'.join(spin_lines)

    return header_block+spins_block

def make_magnetic_spin_symmetry_axes_as_cif(sym_json,axis,cutoff):
    """Read the symmetry analysis for magnetic spin dof,
    and return a cif string where the symmetry axis directions
    have been applied as spin values on each atom. The applied spins
    are normalized by the specified normalization (default 1), and
    any directions smaller than the cutoff are ignored (default 1e-5).

    Parameters
    ----------
    symfile : json
    symfile : str
    symfile : float
    symfile : float

    Returns
    -------
    str

    """
    lat = make_column_vector_matrix(sym_json)
    lengths = make_lengths(lat)
    angles = make_angles(lat)
    cif_lattice_block=make_unit_cell_cif_block(lengths,angles)

    all_sites = make_all_sites_dict(sym_json)
    cif_occupancy_block=make_occupancy_cif_block(lat,all_sites)

    selected_sites = make_selected_sites_dict(sym_json)
    axis_matrix=make_adapted_axis_matrix(sym_json,axis)
    cif_spin_block=make_magnetic_spin_cif_block(lat,axis_matrix,selected_sites,cutoff)

    return '\n'.join((cif_lattice_block,cif_occupancy_block,cif_spin_block))

#######################################################################

def read_json(json_file):
    with open(json_file) as json_data:
        d = json.load(json_data)
    return d

def main():
    scel_dirs=glob.glob("./test_input/SCEL*")

    for scel_dir in scel_dirs:
        dof_spin_sym_file=os.path.join(scel_dir,"dof_analysis_magspin.json")
        sym_dump=read_json(dof_spin_sym_file)

        axis_names=make_symmetry_adapted_axis_names(sym_dump)
        magnetic_cifs=[make_magnetic_spin_symmetry_axes_as_cif(sym_dump,x,1e-5) for x in axis_names]

        for x,cif in zip(axis_names,magnetic_cifs):
            output_path=os.path.join(scel_dir,f"spin_{x}.cif")
            print(f"Writing to {output_path}...")

            output_file=open(output_path,"w")
            output_file.write(cif)
            output_file.close()


if __name__ == "__main__":
    main()
