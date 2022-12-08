"""Generated file from the analysis.org file. Do not edit directly."""

# standard library
import os
import os.path as osp
import shutil as sh
import pickle
from pathlib import Path
from copy import copy, deepcopy
import itertools as it
import gc
import time

# de facto standard library
import numpy as np
import pandas as pd
import sqlalchemy as sqla
import matplotlib.pyplot as plt

# extra non-domain specific
import joblib

# 3rd party domain specific
import mdtraj as mdj
import simtk.unit as tkunit

from tabulate import tabulate

# auxiliary supporting repos
import geomm
import wepy


# some utilities
from wepy.util.mdtraj import traj_fields_to_mdtraj


# project specific auxiliary
import seh_prep

# truly ad hoc

TREE_PATH = Path(osp.expandvars("$TREE"))

PROJECT_PATH = TREE_PATH / "lab/projects/seh.pathway_hopping"

## Paths

# for localizing paths to very commonly used resources and resrouces
# which may change schema. The directory structure for the rest is the
# schema, so just use osp.join(project_path(), 'subpath/to/resource')
# for the rest so a lot of work is reduced in specifying all of them

def data_path():
    return Path(osp.join(PROJECT_PATH, 'data'))

def hpcc_path():
    return Path(osp.join(PROJECT_PATH, 'hpcc'))

def db_path():
    return Path(osp.join(PROJECT_PATH, 'db'))

def media_path():
    return Path(osp.join(PROJECT_PATH, 'media'))

def scratch_path():
    return Path(osp.join(PROJECT_PATH, 'scratch'))

def scripts_path():
    return Path(osp.join(PROJECT_PATH, 'scripts'))

def src_path():
    return Path(osp.join(PROJECT_PATH, 'src'))

def tmp_path():
    return Path(osp.join(PROJECT_PATH, 'tmp'))

def troubleshoot_path():
    return Path(osp.join(PROJECT_PATH, 'troubleshoot'))


# specific things
def sqlite_path():
    return Path(osp.join(PROJECT_PATH, 'db/db.sqlite'))

def joblib_cache_path():
    return Path(osp.join(PROJECT_PATH, 'cache/joblib'))


# other projects paths

def projects_path():
    return TREE_PATH / "lab/projects"


# inhibitors
def inhibitors_proj_path():

    return projects_path() / "seh.inhibitors"

## Setup

# create the sqlite database

# set up the joblib cache
jlmem = joblib.Memory(joblib_cache_path())



# set this when you want to do some recursion stuff with contigtrees
def set_recursion_limit():
    recursion_limit = 5000
    import sys; sys.setrecursionlimit(recursion_limit)
    print("Setting recursion limit to {}".format(recursion_limit))

# set the recursion depth since it is always needing to be increased
set_recursion_limit()

def load_obj(filepath):

    import os.path as osp
    import pickle

    import joblib

    fname = osp.basename(filepath)

    # use the file extension for how to load it
    if fname.endswith('jl.pkl'):
        # it is a joblib object so use joblib to load it
        with open(filepath, 'rb') as rf:
            obj = joblib.load(rf)

    elif fname.endswith('pkl'):
        # it is a pickle object so use joblib to load it
        with open(filepath, 'rb') as rf:
            obj = pickle.load(rf)


    return obj


def save_obj(obj_path, obj, overwrite=False, ext='jl.pkl'):

    import os
    import os.path as osp
    import pickle
    import joblib

    if ext == 'jl.pkl':
        pickler_dump = joblib.dump
    elif ext == 'pkl':
        pickler_dump = pickle.dump
    else:
        raise ValueError("Must choose an extension for format selection")

    # if we are not overwriting check if it exists
    if not overwrite:
        if osp.exists(obj_path):
            raise OSError("File exists ({}), not overwriting".format(obj_path))

    # otherwise make sure the directory exists
    os.makedirs(osp.dirname(obj_path), exist_ok=True)

    # it is a joblib object so use joblib to load it
    with open(obj_path, 'wb') as wf:
        pickler_dump(obj, wf)


def load_table(filepath):

    import os.path as osp

    import pandas as pd

    fname = osp.basename(filepath)

    # use the file extension for how to load it
    if fname.endswith('csv'):

        df = pd.read_csv(filepath, index_col=0)

    elif fname.endswith('pkl'):

        df = pd.read_pickle(filepath)

    else:
        raise ValueError("extension not supported")



    return df

def save_table(table_path, df, overwrite=False, ext='csv'):

    import os
    import os.path as osp
    import pickle

    import pandas as pd

    # if we are not overwriting check if it exists
    if not overwrite:
        if osp.exists(table_path):
            raise OSError("File exists ({}), not overwriting".format(table_path))

    # otherwise make sure the directory exists for this observable
    os.makedirs(osp.dirname(table_path), exist_ok=True)

    if ext == 'csv':

        df.to_csv(table_path)

    elif ext == 'pkl':

        df.to_pickle(table_path)

    else:
        raise ValueError("extension not supported")

## Parameters

# most constants we should take from seh_prep so we don't have
# duplication and potential error
import seh_prep.parameters as seh_params

# TODO: separate the ligand ids from the gexps, since now we have two
# gexps for the same ligand (TPPU/17)
GEXPS = ('3', '10', '17', '18', '20', 'TPPU-legacy',)

LEGACY_GEXPS = ('TPPU-legacy',)

# TODO: currently the LIG_IDS is used in a lot of places that the
# GEXPS should be used. These need to be refactored
LIG_IDS = ('3', '10', '17', '18', '20')

# maps the gexps to the ligand they use, since this is reused
GEXP_LIG_IDS = (
    ('3', '3'),
    ('10', '10'),
    ('17', '17'),
    ('18', '18'),
    ('20', '20'),
    ('TPPU-legacy', '17'),
)


# number of individual steps per cycle in simulations
NUM_CYCLE_STEPS = 20000

CYCLE_TIME = NUM_CYCLE_STEPS * seh_params.STEP_TIME

GEXP_NUM_WALKERS = {
    '3' : 48,
    '10' : 48,
    '17' : 48,
    '18' : 48,
    '20' : 48,
    'TPPU-legacy' : 48,

}

# load the inhibitor data from the csv table

def inhibitors_df():
    import os.path as osp
    import pandas as pd

    return pd.read_csv(inhibitors_proj_path() / 'inhibitors.csv')

# the key for the indexes we use in this analysis
INDEX_KEY = 'inhibitor_index'

# TODO: refactor the MFPT to HALFTIME

# its not the raw rate since it is ln(2) / koff

# TODO: refactor the RATE to KOFF

# these will help with the ambiguity and confusion between the units

# keys for different values. We have a mapping since in the table
# there is a bunch of complexity in the experimental conditions etc.
HALFTIME_KEY = 'half-time_human'
KOFF_KEY = 'koff'

HALFTIME_ERROR_KEY = 'half-time_error_human'
KOFF_ERROR_KEY = 'koff_error'


EXPERIMENTAL_HALFTIME_UNIT = tkunit.minute
EXPERIMENTAL_HALFTIME_FACTOR = 1

EXPERIMENTAL_KOFF_UNIT = (1/tkunit.second).unit
EXPERIMENTAL_KOFF_FACTOR = 1e-4

# e.g.:
# koff_q = koff_value * EXPERIMENTAL_KOFF_FACTOR * EXPERIMENTAL_KOFF_UNIT


def convert_rate_to_halftime(rate_q, time_base=tkunit.minute):

    return (np.log(2) / rate_q).value_in_unit(time_base)

def convert_halftime_to_rate(halftime_q):

    # lol its the same

    return np.log(2) / halftime_q

def convert_halftime_to_rt(halftime_q):

    return 1 / (np.log(2) / halftime_q)

def convert_rt_to_halftime(rt_q):

    # lol its the same

    return np.log(2) / (1 / halftime_q)


def experimental_halftimes():

    inh_df = inhibitors_df()
    values = [(
        lig_id,
        inh_df[inh_df[INDEX_KEY] == int(lig_id)][HALFTIME_KEY].values[0],
    )
              for lig_id in LIG_IDS]

    # convert to the proper units/factors
    quantities = [(
        lig_id,
        val * EXPERIMENTAL_HALFTIME_FACTOR * EXPERIMENTAL_HALFTIME_UNIT,
    )
                  for lig_id, val in values]

    return tuple(quantities)

def experimental_koffs():

    inh_df = inhibitors_df()
    values = [(
        lig_id,
        inh_df[inh_df[INDEX_KEY] == int(lig_id)][KOFF_KEY].values[0],
    )
              for lig_id in LIG_IDS]

    # convert to the proper units/factors
    quantities = [(
        lig_id,
        val * EXPERIMENTAL_KOFF_FACTOR * EXPERIMENTAL_KOFF_UNIT,
    )
                  for lig_id, val in values]


    return tuple(quantities)

def experimental_halftime_errors():


    inh_df = inhibitors_df()
    values = [(
        lig_id,
        inh_df[inh_df[INDEX_KEY] == int(lig_id)][HALFTIME_ERROR_KEY].values[0],
    )
              for lig_id in LIG_IDS]

    # convert to the proper units/factors
    quantities = [(
        lig_id,
        val *  EXPERIMENTAL_HALFTIME_FACTOR * EXPERIMENTAL_HALFTIME_UNIT,
    )
                  for lig_id, val in values]


    return tuple(quantities)



def experimental_koff_errors():



    inh_df = inhibitors_df()
    values = [(
        lig_id,
        inh_df[inh_df[INDEX_KEY] == int(lig_id)][KOFF_ERROR_KEY].values[0],
    )
              for lig_id in LIG_IDS]

    # convert to the proper units/factors
    quantities = [(
        lig_id,
        val * EXPERIMENTAL_KOFF_FACTOR * EXPERIMENTAL_KOFF_UNIT,
    )
                  for lig_id, val in values]


    return tuple(quantities)

## Details on how to save figures and media

# all the extensions of figures to save from matplotlib
FIG_EXTENSIONS = (
    'pdf',
)

## Categorical color scheme

CATEGORY_COLORS = (

    ('blue', (
        ('base', "#C6DBEF"),
        ('pair-light', "#1F77B4"),
    )),

    ('orange', (
        ('base', "#FF7F0E"),
        ('pair-light', "#FFBB78"),
    )),

    ('green', (
        ('base', "#2CA02C"),
        ('pair-light', "#98DF8A"),
    )),

    ('red', (
        ('base', "#D62728"),
        ('pair-light', "#FF9896"),
    )),

    ('purple', (
        ('base', "#9467BD"),
        ('pair-light', "#C5B0D5"),
    )),

    ('brown', (
        ('base', "#8C564B"),
        ('pair-light', "#C49C94"),
    )),

    ('pink', (
        ('base', "#E377C2"),
        ('pair-light', "#F7B6D2"),
    )),

    ('grey', (
        ('base', "#7F7F7F"),
        ('pair-light', "#C7C7C7"),
    )),

    ('yellow', (
        ('base', "#BCBD22"),
        ('pair-light', "#DBDB8D"),
    )),

    ('cyan', (
        ('base', "#17BECF"),
        ('pair-light', "#9EDAE5"),
    )),

)





# Each GEXP gets a color assigned to it to be used for all things if
# possible

GEXP_COLORS = (
      ('3', 'blue'),
      ('10', 'orange'),
      ('17', 'green'),
      ('18', 'red'),
      ('20', 'purple'),
      ('TPPU-legacy', 'pink'),
)

REPLICATE_COLORS = (

    (0, (
        ('base', "#1B9E77"),
    )),

    (1, (
        ('base', "#D95F02"),
    )),

    (2, (
        ('base', "#7570B3"),
    )),

    (3, (
        ('base', "#E7298A"),
    )),

    (4, (
        ('base', "#66A61E"),
    )),

    (5, (
        ('base', "#E6AB02"),
    )),

    (6, (
        ('base', "#A6761D"),
    )),

    # for aggregates of replicates
    ('agg', (
        ('base', "#666666"),
        ('light', "#CCCCCC")
    )),
)

# plot colors

# plots with different lines for different runs with averages with
# errors need to have consistent colors.

def legacy_project_path():

    return projects_path() / 'seh.tppu_unbinding'

### PDB topology stuff

def get_legacy_real_rep_src_top():
    """Load the two source files necessary to get a complete system. Note
    that the box will not be perfect but it will do for our purposes
    since the ligand at this point is not through the boundaries

    """

    src_pdb_path = legacy_project_path() / 'nowater.pdb'
    src_box_vectors_path = legacy_project_path() / 'image_frame.pdb'

    src_traj = mdj.load_pdb(str(src_pdb_path))
    src_bv_traj = mdj.load_pdb(str(src_box_vectors_path))

    # combine the box vectors with the src traj
    src_traj.unitcell_vectors = src_bv_traj.unitcell_vectors

    return src_traj

def make_real_rep_top():
    """Using the reference file make a processed version compatible with
    our modern naming schemes.

    """
    import json

    from wepy.util.mdtraj import (
        mdtraj_to_json_topology,
        json_to_mdtraj_topology,
    )

    src_traj = get_legacy_real_rep_src_top()

    # convert to json_top
    src_json_top_str = mdtraj_to_json_topology(src_traj.top)

    src_json_top = json.loads(src_json_top_str)


    # make the new json top starting with some top level stuff
    real_rep_json_top = {
        'bonds' : src_json_top['bonds'],
        # there is only one chain so hardcode
        'chains' : [
            {
                'index' : src_json_top['chains'][0]['index'],
                'residues' : [],
            },
        ],
    }

    # process this to rename stuff
    for residue in src_json_top['chains'][0]['residues']:

        # initialize with the same
        new_residue = residue

        # match it against some cases and alter if necessary

        if residue['segmentID'] == 'SML':

            new_residue['segmentID'] = 'HETA'

        if residue['name'] == '2RV':

            new_residue['name'] = 'UNL'

        if (
                residue['segmentID'] == 'HETA' and
                residue['name'] == 'SOD'
        ):

            new_residue['segmentID'] = 'SOD'

        # add the new residue to the chain
        real_rep_json_top['chains'][0]['residues'].append(new_residue)

    # make the real rep traj with this topology

    real_rep_mdj_top = json_to_mdtraj_topology(json.dumps(real_rep_json_top))

    real_rep_traj = mdj.Trajectory(
        src_traj.xyz,
        topology=real_rep_mdj_top,
        unitcell_lengths=src_traj.unitcell_lengths,
        unitcell_angles=src_traj.unitcell_angles,
    )

    return real_rep_traj, real_rep_json_top

def save_real_rep_top():

    import json
    import os
    os.makedirs(
        data_path() / 'legacy/top',
        exist_ok=True
    )

    # make the top from the original source
    real_rep_traj, real_rep_json_top = make_real_rep_top()

    # where to write the mdtraj normalized top
    target_pdb_path = data_path() / 'legacy/top/real_rep.pdb'
    target_json_path = data_path() / 'legacy/top/real_rep.top.json'

    real_rep_traj.save_pdb(str(target_pdb_path))

    with open(target_json_path, 'w') as wf:

        wf.write(json.dumps(real_rep_json_top))

    return real_rep_traj, real_rep_json_top

def get_real_rep_tops():
    """Get the 'real_rep' for the legacy TPPU system as both mdtraj Traj
    from the PDB and the JSON top.

    This doesn't remake it, just loads from a generated file.
    """

    target_pdb_path = data_path() / 'legacy/top/real_rep.pdb'
    target_json_path = data_path() / 'legacy/top/real_rep.top.json'

    real_rep_traj = mdj.load_pdb(str(target_pdb_path))

    with open(target_json_path, 'r') as rf:

        real_rep_json_top = rf.read()


    return real_rep_traj, real_rep_json_top


### Writing out centered and transformed stuff for the reference topology

def real_rep_selection_idxs_tops(main_sel_idxs, main_sel_tops):
    """Get the special 'real_rep' selection idxs for the legacy top

    Warning
    -------

    Does not give the binding site stuff selection see: real_rep_bs_selection_idxs_tops

"""

    from seh_prep.modules import (
        ligand_idxs,
        protein_idxs,
    )

    from seh_prep.modules import (
        binding_site_idxs,
    )
    from seh_prep.parameters import BINDING_SITE_CUTOFF

    sel_idxs = {}
    sel_tops = {}


    # the initial state
    init_state = get_ref_state('TPPU-legacy')

    # get some things from this
    init_box_vectors = init_state['box_vectors'] # * init_state.box_vectors_unit

    # get the init state for another known gexp
    tmp_state = get_ref_state('17')

    # use this for the units
    positions_unit = tmp_state.positions_unit
    box_vectors_unit = tmp_state.box_vectors_unit

    # load the uncentered reference for bootstrapping purposes
    real_rep_centered_path = data_path() / 'legacy/top/real_rep.pdb'

    real_ref_traj, real_rep_json_top = get_real_rep_tops()

    sel_tops['real_rep'] = real_rep_json_top

    real_rep_positions = real_ref_traj.xyz[0]

    # get the ligand, protein, and binding site indices for this
    real_ligand_idxs = ligand_idxs(real_rep_json_top)
    real_protein_idxs = protein_idxs(real_rep_json_top)


    sel_idxs['real_rep/ligand'] = real_ligand_idxs
    sel_idxs['real_rep/protein'] = real_protein_idxs

    sel_idxs['real_rep/ligand/homology'] = \
        sel_idxs['real_rep/ligand'][
            main_sel_idxs['ligand/homology']
        ]


    # binding site stuff
    real_bs_idxs = binding_site_idxs(
        real_rep_json_top,
        real_rep_positions * positions_unit,
        init_box_vectors * box_vectors_unit,
        BINDING_SITE_CUTOFF,
    )

    sel_idxs['real_rep/binding_site'] = real_bs_idxs


    # get them relative to the all_atoms selections, which is all of
    # them, so we just set them the same
    sel_idxs['all_atoms/real_rep'] = range(real_ref_traj.n_atoms)
    sel_idxs['all_atoms/real_rep/ligand'] = sel_idxs['real_rep/ligand']
    sel_idxs['all_atoms/real_rep/ligand'] = sel_idxs['real_rep/ligand']


    return sel_idxs, sel_tops

# NOTE: don't need it

# def make_centered_real_rep_top():

#     from geomm.grouping import group_pair
#     from geomm.centering import center_around

#     traj, _ = get_real_rep_tops()

#     sel_idxs, _ = real_rep_selection_idxs_tops()

#     protein_idxs = sel_idxs['real_rep/protein']
#     lig_idxs = sel_idxs['real_rep/ligand']

#     # group
#     grouped_positions = group_pair(
#         traj.xyz[0],
#         traj.unitcell_lengths,
#         protein_idxs,
#         lig_idxs,
#     )

#     # center
#     centered_positions = center_around(
#         grouped_positions,
#         protein_idxs,
#     )


#     centered_traj = traj_fields_to_mdtraj(centered_traj_fields)

#     return centered_traj

# def save_centered_real_rep_top():
#     """Save the centered version of the topology for the legaxy
#     real_rep."""

#     target_path = data_path() / 'legacy/top/real_rep_center_ref.pdb'

#     centered_traj = make_centered_real_rep_top()

#     centered_traj.save_pdb(str(target_path))

#     return centered_traj





### legacy HDF5 stuff
def legacy_new_template_path():

    return data_path() / 'legacy/linker_template.wepy.h5'

def legacy_results_file_paths():

    f = osp.join(legacy_project_path(), 'all_trajs.wepy.h5')
    f_ext = osp.join(legacy_project_path(), 'all_trajs_ext.wepy.h5')

    return [f, f_ext]


def legacy_linker_file_path():

    return osp.join(data_path(), 'legacy/17.wepy.h5')

def make_legacy_linker_header_template():

    import h5py

    # use ligand 17 as the template
    wepy_h5 = get_gexp_wepy_h5('17')

    # get the reference topologies
    real_rep_traj, real_rep_json = get_real_rep_tops()

    n_atoms = real_rep_traj.n_atoms

    with wepy_h5:
        # we will make the base one the templat
        legacy_template_path = legacy_new_template_path()

        hdf = h5py.File(legacy_template_path, mode='w')

        with hdf:

            # topology
            hdf.create_dataset(
                'topology',
                data=real_rep_json,
            )

            # units, the old simulations units are the same as the new ones
            _ = hdf.create_group('units/units')
            for unit_key in list(wepy_h5.h5['units/units']):

                hdf['units/units'].create_dataset(
                    unit_key,
                    data=wepy_h5.h5[f'units/units/{unit_key}'][()]
                )

            # settings
            hdf.create_group('_settings')

            hdf.create_group('_settings/alt_reps_idxs')

            hdf.create_dataset('_settings/continuations',
                               shape=(0,2),
                               dtype=wepy_h5.h5['_settings/continuations'].dtype,
                               maxshape=wepy_h5.h5['_settings/continuations'].maxshape)

            hdf.create_group('_settings/field_feature_dtypes')
            hdf.create_group('_settings/field_feature_shapes')

            # TODO: field_feature shapes and types actually set, don't think these
            # are really used outside of data generation

            hdf.create_dataset("_settings/n_atoms", data=n_atoms)

            n_dims = 3
            hdf.create_dataset("_settings/n_dims", data=n_dims)

            # record fields
            hdf.create_group('_settings/record_fields')
            # for rec_field_key in list(wepy_h5.h5['_settings/record_fields']):

            #     p = f'_settings/record_fields/{rec_field_key}'
            #     hdf.create_dataset(
            #         p,
            #         data=wepy_h5.h5[p],
            #     )

            # sparse fields
            hdf.create_dataset(
                '_settings/sparse_fields',
                data=[],
            )

            main_rep_idxs = range(n_atoms)
            hdf.create_dataset('_settings/main_rep_idxs', data=main_rep_idxs)

    return hdf


def make_legacy_results_linker_h5():

    from wepy.hdf5 import WepyHDF5

    # first make sure the linker header template exists
    _ = make_legacy_linker_header_template()

    template_path = legacy_new_template_path()
    component_paths = legacy_results_file_paths()
    linker_path = legacy_linker_file_path()

    template_wepy = WepyHDF5(template_path, mode='r')

    with template_wepy:
        linker_wepy_h5 = template_wepy.clone(linker_path, mode='w')

    with linker_wepy_h5:

        for component_path in component_paths:
            with WepyHDF5(component_path, mode='r') as other_wepy:
                run_idxs = other_wepy.run_idxs

            for run_idx in run_idxs:
                new_run_idx = linker_wepy_h5.link_run(component_path, run_idx)

def get_legacy_h5(
        mode='r',
):

    from wepy.hdf5 import WepyHDF5

    path = legacy_linker_file_path()

    # don't open in swmr mode if we have write intent
    if mode in ('r+', ):
        wepy_h5 = WepyHDF5(
            path,
            mode=mode,
        )

    # if its just read use swmr mode
    elif mode in ('r',):
        wepy_h5 = WepyHDF5(
            path,
            mode=mode,
            swmr_mode=True,
        )
    else:
        raise ValueError(f"Invalid mode {mode}")


    return wepy_h5

### Data specific to old WExplore
def get_hist_df():

    import pandas as pd

    path = osp.join(legacy_project_path(), 'hist.csv')

    return pd.read_csv(path, index_col=0)

@jlmem.cache
def legacy_walker_traj_map():
    """For a multi-run hist file returns (run_idx, walker_idx)->[(traj_idx, frame_idx),... ]

    Maps the (run_idx, traj_idx) for the HDF5 to the (legacy_traj_idx,
    frame_idx) of the raw set of trajectories.

    """

    from collections import defaultdict

    hist_df = get_hist_df()

    first_frame_map = defaultdict(list)
    second_frame_map = defaultdict(list)

    for run_idx, run_df in hist_df.groupby('run_number'):
        for walker_idx, walker_df in run_df.groupby('walker_idx'):
            for row_idx, row in walker_df.iterrows():

                first_frame_map[(run_idx, walker_idx)].append(
                    ( int(row['traj_idx']), int(row['first_frame']) )
                )

                second_frame_map[(run_idx, walker_idx)].append(
                    ( int(row['traj_idx']), int(row['second_frame']) )
                )

    return first_frame_map, second_frame_map

@jlmem.cache
def legacy_second_frame_map():
    """(leg_traj_idx, leg_frame_idx) -> (run_idx, traj_idx, frame_idx) """

    _, wt_map = legacy_walker_traj_map()

    frame_map = {}
    for rec, leg_trace in wt_map.items():

        run_idx, traj_idx = rec

        for frame_idx, leg_rec in enumerate(leg_trace):

            leg_traj_idx, leg_frame_idx = leg_rec

            frame_rec = (run_idx, traj_idx, frame_idx)

            frame_map[(leg_traj_idx, leg_frame_idx)] = frame_rec

    return frame_map

@jlmem.cache
def legacy_first_frame_map():
    """(leg_traj_idx, leg_frame_idx) -> (run_idx, traj_idx, frame_idx) """

    wt_map, _ = legacy_walker_traj_map()

    frame_map = {}
    for rec, leg_trace in wt_map.items():

        run_idx, traj_idx = rec

        for frame_idx, leg_rec in enumerate(leg_trace):

            leg_traj_idx, leg_frame_idx = leg_rec

            frame_rec = (run_idx, traj_idx, frame_idx)

            frame_map[(leg_traj_idx, leg_frame_idx)] = frame_rec

    return frame_map


def legacy_frame_map():

    return legacy_second_frame_map()

def get_legacy_cluster_legacy_traj_assignments():

    import pickle

    assg_path = osp.join(legacy_project_path(), 'cluster_frames_map.pkl')
    with open(assg_path, 'rb') as rf:
        traj_assgs = pickle.load(rf)

    return traj_assgs

def legacy_trace_update(legacy_trace):
    """convert trace of legacy idxs (legacy_traj_idx, legacy_frame_idx) to
    a trace of (run_idx, traj_idx, frame_idx).

    This will not ignore any data at all. If a frame cannot be found a
    None will be used instead.

    """

    ff_map = legacy_first_frame_map()
    sf_map = legacy_second_frame_map()

    new_trace = []
    for leg_rec in legacy_trace:
        if leg_rec in ff_map:
            new_rec = ff_map[leg_rec]
            new_trace.append(new_rec)

        elif leg_rec in sf_map:
            new_rec = sf_map[leg_rec]
            new_trace.append(new_rec)

        else:
            new_trace.append(None)

    return new_trace


def legacy_trace_canonical_update(legacy_trace):
    """convert trace of legacy idxs (legacy_traj_idx, legacy_frame_idx) to
    a trace of (run_idx, traj_idx, frame_idx).

    This ignores any data that is not a 'second_frame' since the
    firs_frame is redundant. ALso ignores frames that can't be found.

    Puts Nones in their place so you can do indexing if you need.

    Use 'legacy_trace_canonical_filtered_update' to just get the valid
    ones.

    """

    sf_map = legacy_second_frame_map()

    new_trace = []
    for leg_rec in legacy_trace:
        if leg_rec in sf_map:
            new_rec = sf_map[leg_rec]
            new_trace.append(new_rec)

        else:
            new_trace.append(None)

    return new_trace


def legacy_trace_canonical_filtered_update(legacy_trace):

    sf_map = legacy_second_frame_map()

    new_trace = []
    for leg_rec in legacy_trace:
        if leg_rec in sf_map:
            new_rec = sf_map[leg_rec]
            new_trace.append(new_rec)

    return new_trace

### Things related to the old clustering

@jlmem.cache
def get_legacy_cluster_canonical_traj_assignments():

    leg_assgs = get_legacy_cluster_legacy_traj_assignments()

    # convert these to the canonical trace filtereing out the junk
    canon_assgs = {}
    for clust_id, leg_trace in leg_assgs.items():

        print("Working on cluster {}".format(clust_id))

        canon_assgs[clust_id] = legacy_trace_canonical_filtered_update(leg_trace)

    return canon_assgs


def get_legacy_ts_cluster_idxs():

    path = osp.join(legacy_project_path(), 'TSE_idxs.dat')

    with open(path, 'r') as rf:

        cluster_idxs = []
        for line in rf.readlines():
            cluster_idxs.append(int(line.strip()))

    return cluster_idxs

def get_gexp_sim_dir(gexp):

    return hpcc_path() / f'simulations/{gexp}_simulations'

def get_gexp_orch_path(gexp):

    return get_gexp_sim_dir(gexp) / f'orchs/master_sEH_lig-{gexp}.orch.sqlite'

def get_gexp_master_orch(gexp):

    from wepy.orchestration.orchestrator import Orchestrator

    orch_path = get_gexp_orch_path(gexp)

    orch = Orchestrator(str(orch_path), mode='r')

    return orch

def all_gexp_wepy_h5s():
    d = {}
    for gexp in GEXPS:
        d[gexp] = get_gexp_wepy_h5(gexp)
    return d

def gexp_wepy_h5_path(gexp):

    return data_path() / f"results/{gexp}/all_results.wepy.h5"

def get_gexp_wepy_h5(
        gexp,
        mode='r',
):

    # short circuit on special case for legacy
    if gexp in LEGACY_GEXPS:
        return get_legacy_h5(mode=mode)

    from wepy.hdf5 import WepyHDF5

    path = gexp_wepy_h5_path(gexp)

    # don't open in swmr mode if we have write intent
    if mode in ('r+', ):
        wepy_h5 = WepyHDF5(
            path,
            mode=mode,
        )

    # if its just read use swmr mode
    elif mode in ('r',):
        wepy_h5 = WepyHDF5(
            path,
            mode=mode,
            swmr_mode=True,
        )
    else:
        raise ValueError(f"Invalid mode {mode}")


    return wepy_h5


def delete_wepy_h5_observable(
        gexp,
        obs_name,
):

    wepy_h5 = get_gexp_wepy_h5(
        gexp,
        mode='r+',
    )

    with wepy_h5:

        for traj in wepy_h5.iter_trajs():

            try:
                del traj[f'observables/{obs_name}']
            except KeyError:
                print("already deleted")

# the parameters used for each test group
from wepy.resampling.resamplers.resampler import NoResampler
from wepy.resampling.resamplers.revo import REVOResampler
from wepy.resampling.resamplers.wexplore import WExploreResampler

from wepy.resampling.decisions.decision import NoDecision
from wepy.resampling.decisions.clone_merge import MultiCloneMergeDecision


GEXP_METADATA = {
    '3' : {
        'resampler' : WExploreResampler,
        'decision' : MultiCloneMergeDecision,
        'num_walkers' : 48,
    },
    '10' : {
        'resampler' : WExploreResampler,
        'decision' : MultiCloneMergeDecision,
        'num_walkers' : 48,
    },
    '17' : {
        'resampler' : WExploreResampler,
        'decision' : MultiCloneMergeDecision,
        'num_walkers' : 48,
    },
    '18' : {
        'resampler' : WExploreResampler,
        'decision' : MultiCloneMergeDecision,
        'num_walkers' : 48,
    },
    '20' : {
        'resampler' : WExploreResampler,
        'decision' : MultiCloneMergeDecision,
        'num_walkers' : 48,
    },
    'TPPU-legacy' : {
        'num_walkers' : 48,
    },

}

# we use the job management tables for getting the span IDs and
# segments so we can reference them in a sane way

# Datatypes for the ones we care about parsing
JOBS_TABLE_COLUMNS = (
    ('job_ID', 'string'),
    ('contig_ID', 'string'),
    ('start_hash', 'string'),
    ('patched_start_hash', 'string'),
    ('end_hash', 'string'),
)

def get_gexp_jobs_df(gexp):


    jobs_path = data_path() / f'sim_management/lig_{gexp}.csv'

    print("WARNING: Make sure this is a recent table export")
    print(f"reading: {jobs_path}")

    jobs_df = pd.read_csv(jobs_path,
                          dtype=dict(JOBS_TABLE_COLUMNS)
    )

    # add in the gexp as a column
    jobs_df['gexp'] = gexp

    # parse the contig_ID column and get the span_id and segment_idx
    span_id_col = []
    segment_idx_col = []
    for contig_id in jobs_df['contig_ID']:

        span_id, segment_idx = [int(e) for e in contig_id.split('-')]

        span_id_col.append(span_id)
        segment_idx_col.append(segment_idx)

    jobs_df['span_id'] = span_id_col
    jobs_df['segment_idx'] = segment_idx_col

    return jobs_df

def get_master_jobs_df():

    dfs = []

    for gexp in GEXPS:
        dfs.append(get_gexp_jobs_df(gexp))

    return pd.concat(dfs)

def get_gexp_jobs():

    gexp_jobs = {}
    for gexp in GEXPS:

        jobs_df = get_gexp_jobs_df(gexp)

        gexp_jobs[gexp] = list(grp_df['job_ID'].values)

    return gexp_jobs


def get_gexp_paths():
    """Get a dictionary of all of the paths for the jobs dirs for each
    experimental group."""

    raise NotImplementedError

    # gexp_jobs = get_gexp_jobs()

    # # TODO:
    # sims_dir = data_path() / 'sims/wepy/jobs'

    # gexp_paths = defaultdict(list)
    # for grp_name, jobids in gexp_jobs.items():

    #     for job_id in jobids:
    #         path = sims_dir / str(job_id)
    #         gexp_paths[grp_name].append(path)

    # return gexp_paths


def get_gexp_jobs_to_runs_map(gexp):
    """Get the runs in the HDF5 for the gexp that are mapped to by the
    jobids."""

    wepy_h5 = get_gexp_wepy_h5(gexp)

    jobs_df = get_gexp_jobs_df(gexp)

    jobs_runs_d = {}
    with wepy_h5:

        for run_idx in wepy_h5.run_idxs:
            run_grp = wepy_h5.run(run_idx)

            # this is an attribute as well
            run_idx = run_grp.attrs['run_idx']

            # get the start and end hashes from the HDF5 runs and use
            # this to get the job id for this run idx. Just use the
            # end hash for this

            end_hash = run_grp.attrs['end_snapshot_hash']

            # get the row from the dataframe for this
            job_row = jobs_df[jobs_df['end_hash'] == end_hash]

            jobid = job_row['job_ID'].values[0]

            jobs_runs_d[jobid] = run_idx

    return jobs_runs_d


def get_gexp_span_ids(gexp):

    gexp_jobs_df = get_gexp_jobs_df(gexp)

    span_ids = sorted(list(set(gexp_jobs_df['span_id'].values)))

    return span_ids

def get_gexp_span_ids_run_idxs(gexp):
    """Get a mapping of the span_ids (names given to them) to the run idxs
    (in order of the segments in the span via 'segment_idxs').

    The span_ids and segment_idxs come from the jobs_df table.

    """

    gexp_jobs_df = get_gexp_jobs_df(gexp)

    # get the mapping of jobids -> h5 run idxs
    jobs_runs_d = get_gexp_jobs_to_runs_map(gexp)

    # then go through each span and get the run idxs for it
    spans_runs = {}
    for span_id, span_df in gexp_jobs_df.groupby('span_id'):

        # make records for each segment so that we can sort them
        span_segs = []
        for idx, row in span_df.iterrows():

            seg_rec = (row['segment_idx'], row['job_ID'])
            span_segs.append(seg_rec)

        # sort them, dereference the run_idx
        span_run_idxs = [jobs_runs_d[jobid]
                         for seg_idx, jobid in sorted(span_segs)]

        spans_runs[span_id] = span_run_idxs

    return spans_runs

def gexps_contigtree():
    d = {}
    for gexp in GEXPS:
        d[gexp] = get_contigtree(gexp)

    return d


# the base contigtree is the one that can be cached not the one with
# the HDF5 file in it, so we cache this, you can use the other
# function to get the contigtree with the HDF5, but the openening of
# file handles is annoying and this might be the best way to do this
# anyhow
@jlmem.cache
def get_base_contigtree(gexp,
                        runs=Ellipsis,
):

    # uncache

    from wepy.analysis.contig_tree import BaseContigTree

    from wepy.resampling.decisions.clone_merge import MultiCloneMergeDecision
    from wepy.boundary_conditions.receptor import UnbindingBC

    wepy_h5 = get_gexp_wepy_h5(gexp)

    set_recursion_limit()
    base_contigtree = BaseContigTree(
        wepy_h5,
        runs=runs,
        decision_class=MultiCloneMergeDecision,
        boundary_condition_class=UnbindingBC)

    return base_contigtree

def get_contigtree(gexp,
                   runs=Ellipsis):

    from wepy.analysis.contig_tree import ContigTree

    wepy_h5 = get_gexp_wepy_h5(gexp)
    base_contigtree = get_base_contigtree(gexp,
                                          runs=runs
    )

    return ContigTree(wepy_h5,
                      base_contigtree=base_contigtree,
    )

@jlmem.cache
def get_gexp_span_ids_span_idxs_map(gexp):
    """Get a mapping of the span_ids to the span_idxs in the
    contigtree."""

    contigtree = get_base_contigtree(gexp)

    # get the first run idx -> span_id map
    run_idxs_to_span_ids = {run_idxs[0] : span_id
                            for span_id, run_idxs
                            in get_gexp_span_ids_run_idxs(gexp).items()}

    span_ids_to_span_idxs_d = {}
    for span_idx, trace in contigtree.span_traces.items():
        run_idx = trace[0][0]

        span_id = run_idxs_to_span_ids[run_idx]

        span_ids_to_span_idxs_d[span_id] = span_idx

    return span_ids_to_span_idxs_d


def get_gexp_span_contig(gexp, span_id):

    # get the span idx in the contigtree from the span_id in the jobs
    # table
    span_ids_to_span_idxs_d = get_gexp_span_ids_span_idxs_map(gexp)

    span_idx = span_ids_to_span_idxs_d[span_id]

    # get the contig for this span
    contigtree = get_contigtree(gexp)
    contig = contigtree.span_contig(span_idx)

    return contig

def ref_states():

    ref_states = {}

    for gexp in GEXPS:
        ref_states[gexp] = get_ref_state(gexp)

    return ref_states


def get_ref_state(gexp):

    from wepy.walker import WalkerState

    if gexp in LEGACY_GEXPS:

        # get the PDB mdtraj traj
        ref_traj, _ = get_real_rep_tops()

        # convert to a state object
        kwargs = {
            'positions' : ref_traj.xyz[0],
            'box_vectors' : ref_traj.unitcell_vectors[0],
        }

        init_state = WalkerState(**kwargs)

    else:

        init_state_path = data_path() / \
              f"md_systems/{gexp}/sEH_lig-{gexp}_equilibrated.state.pkl"

        with open(init_state_path, 'rb') as rf:
            init_state = pickle.load(rf)

    return init_state


def get_ref_state_traj_fields(gexp, rep_key='main_rep'):

    import numpy as np

    selection_idxs = lig_selection_idxs(gexp)
    tops = lig_selection_tops(gexp)
    ref_state = get_ref_state(gexp)

    ref_state_traj_fields = {}

    if rep_key is None or rep_key == 'all_atoms':
        positions = ref_state['positions']
    else:
        positions = ref_state['positions'][
                         selection_idxs['all_atoms/{}'.format(rep_key)]]

    ref_state_traj_fields['positions'] = np.array([positions])
    ref_state_traj_fields['box_vectors'] = np.array([ref_state['box_vectors']])

    return ref_state_traj_fields

def get_centered_ref_state_traj_fields(gexp, rep_key='main_rep'):
    """Get the reference coordinates for this gexp and rep.

    This is typically used for grouping, aligning, and recentering
    other structures to a common reference point.

    """

    import numpy as np
    from wepy.util.util import box_vectors_to_lengths_angles

    from geomm.grouping import group_pair
    from geomm.centering import center_around

    # get the selection idxs
    sel_idxs = lig_selection_idxs(gexp)

    # get the required selections for the grouping
    lig_idxs = sel_idxs['{}/ligand'.format(rep_key)]
    bs_idxs = sel_idxs['{}/binding_site'.format(rep_key)]

    # get the native state as a traj field from the all atoms
    # reference state
    ref_fields = get_ref_state_traj_fields(gexp, rep_key=rep_key)

    # then group and center the reference state
    ref_box_lengths, _ = box_vectors_to_lengths_angles(ref_fields['box_vectors'][0])

    grouped_ref_positions = group_pair(ref_fields['positions'][0], ref_box_lengths,
                                       bs_idxs, lig_idxs)

    centered_ref_positions = center_around(grouped_ref_positions, bs_idxs)

    ref_fields['positions'] = centered_ref_positions

    return ref_fields



def ligs_ref_selection_trajs():

    d = {}
    for gexp in GEXPS:
        d[gexp] = lig_ref_selection_trajs(gexp)

    return d

@jlmem.cache
def lig_ref_selection_trajs(gexp):

    import numpy as np

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    selection_idxs = lig_selection_idxs(gexp)
    tops = lig_selection_tops(gexp)
    ref_state = get_ref_state(gexp)

    ref_state_fields = {key : np.array([value]) for key, value in ref_state.dict().items()}

    # make fields for the different representations

    # DEBUG: we have issues with the main rep here so we just don't do this one
    if gexp != 17:
        ref_state_fields['main_rep'] = ref_state['positions'][selection_idxs['all_atoms/main_rep']]
    ref_state_fields['image'] = ref_state['positions'][selection_idxs['all_atoms/image']]
    ref_state_fields['correct_rep'] = ref_state['positions'][selection_idxs['all_atoms/correct_rep']]
    ref_state_fields['ligand'] = ref_state['positions'][selection_idxs['all_atoms/ligand']]

    # for all atoms
    all_atoms_traj = traj_fields_to_mdtraj(ref_state_fields, tops['all_atoms'],
                                           rep_key='positions')

    # main rep
    # DEBUG
    if gexp != 17:
        main_rep_traj = traj_fields_to_mdtraj(ref_state_fields, tops['main_rep'],
                                               rep_key='main_rep')

    # image
    image_traj = traj_fields_to_mdtraj(ref_state_fields, tops['image'],
                                           rep_key='image')

    # correct rep
    correct_rep_traj = traj_fields_to_mdtraj(ref_state_fields, tops['correct_rep'],
                                           rep_key='correct_rep')

    # ligand
    ligand_traj = traj_fields_to_mdtraj(ref_state_fields, tops['ligand'],
                                           rep_key='ligand')

    ref_selections = {}
    ref_selections['all_atoms'] = all_atoms_traj
    if gexp != 17:
        ref_selections['main_rep'] = main_rep_traj
    ref_selections['image'] = image_traj
    ref_selections['correct_rep'] = correct_rep_traj
    ref_selections['ligand'] = ligand_traj

    return ref_selections



def lig_centered_ref_selection_trajs(gexp):

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    # DEBUG
    if gexp == 17:
        rep_keys = ('all_atoms', 'correct_rep',)
    else:
        rep_keys = ('all_atoms', 'main_rep', 'correct_rep',)

    tops = lig_selection_tops(gexp)

    ref_selections = {}
    for rep_key in rep_keys:
        top = tops[rep_key]
        ref_state_fields = get_centered_ref_state_traj_fields(gexp, rep_key=rep_key)
        ref_selections[rep_key] = traj_fields_to_mdtraj(ref_state_fields, top,
                                                        rep_key='positions')

    return ref_selections

def get_topologies():
    tops = {}
    for gexp in GEXPS:

        tops[gexp] = get_topology(gexp)

    return tops


def get_topology(gexp):


    if gexp in LEGACY_GEXPS:

        _, json_top = get_real_rep_tops()

    else:
        top_path = data_path() / f"md_systems/{gexp}/sEH_lig-{gexp}_system.top.json"

        with open(top_path, 'rb') as rf:
            json_top = rf.read()

    return json_top

# was previously named get_lig_num_walkers
def get_gexp_span_num_walkers(gexp, span_id):

    # all are the same
    return GEXP_NUM_WALKERS[gexp]

def write_ligs_ref_selection_pdbs():

    for lig_id in LIG_IDS:
        write_lig_ref_selection_pdbs(lig_id)

def write_lig_ref_selection_pdbs(lig_id):

    import os
    import os.path as osp

    ref_selections = lig_ref_selection_trajs(lig_id)

    lig_top_dir = osp.join(data_path(), 'top/{}'.format(lig_id))

    # make sure the directory exists
    os.makedirs(lig_top_dir,
                exist_ok=True)

    for sel_name, sel_ref_traj in ref_selections.items():

        pdb_path = osp.join(lig_top_dir, "{}_ref.pdb".format(sel_name))
        sel_ref_traj.save_pdb(pdb_path)

def write_lig_centered_ref_selection_pdbs(lig_id):

    import os
    import os.path as osp

    ref_selections = lig_centered_ref_selection_trajs(lig_id)

    lig_top_dir = osp.join(data_path(), 'top/{}'.format(lig_id))

    # make sure the directory exists
    os.makedirs(lig_top_dir,
                exist_ok=True)

    for sel_name, sel_ref_traj in ref_selections.items():

        pdb_path = osp.join(lig_top_dir, "{}_center_ref.pdb".format(sel_name))
        sel_ref_traj.save_pdb(pdb_path)

def save_span_warp_table(gexp, span_idx, overwrite=True):

    table_path = data_path() / f'warp-table/spans/{gexp}_{span_idx}__warps.csv'

    warps_df = get_span_warps_df(gexp, span_idx)

    # only write if there is data
    if len(warps_df) > 0:
        save_table(table_path, warps_df, overwrite=overwrite)

def save_gexp_warp_table(gexp, overwrite=True):

    table_path = data_path() / f'warp-table/gexps/{gexp}_warps.csv'

    warps_df = get_gexp_warps_df(gexp)

    # only write if there is data
    if len(warps_df) > 0:
        save_table(table_path, warps_df, overwrite=overwrite)

def save_span_stats_table(gexp, overwrite=True):

    table_path = data_path() / f'span-stats/{gexp}_span-stats.csv'

    spans_df = get_span_stats_df(gexp)

    save_table(table_path, spans_df, overwrite=overwrite)

def traj_fields_to_correct_rep(
        traj_fields,
        gexp,
):

    # lookup the ligand
    lig_id = dict(GEXP_LIG_IDS)[gexp]

    return traj_fields_to_correct_rep_lig(
        traj_fields,
        lig_id
    )

# based on ligand
def traj_fields_to_correct_rep_lig(
        traj_fields,
        lig_id,
):
    """Given a traj_fields struct with both the 'positions'
    (i.e. 'main_rep') and the 'alt_reps/missing' fields returns a
    traj_fields struct where 'positions' is the 'correct_rep'.

    Will just copy all the other fields over.

    """

    # there is no correct rep data for ligand 3 so skip it
    if lig_id == '3':
        raise ValueError("Cannot make correct rep for ligand 3.")

    assert 'alt_reps/missing' in traj_fields, \
        "The field 'alt_reps/missing' must be in the traj_fields"

    assert 'positions' in traj_fields, \
        "The field 'positions' must be in the traj_fields"

    # get the number of frames, we know positions will be here
    n_frames = traj_fields['positions'].shape[0]

    # get the selection idxs for this ligand
    sel_idxs = lig_selection_idxs(lig_id)

    correct_idxs = sel_idxs['all_atoms/correct_rep']

    # the indices of the missing atoms in the correct rep
    correct_missing_idxs = sel_idxs['correct_rep/missing']
    correct_main_rep_idxs = sel_idxs['correct_rep/main_rep']

    correct_rep_positions = np.zeros((n_frames, len(correct_idxs), 3))

    correct_rep_positions[:, correct_main_rep_idxs, :] = traj_fields['positions'][...]
    correct_rep_positions[:, correct_missing_idxs, :] = traj_fields['alt_reps/missing'][...]

    # make the new struct, but leave out the old stuff
    correct_traj_fields = {
        key : value for key, value in traj_fields.items()
        if key not in ('positions', 'alt_reps/missing')
    }

    # now set the correct one in
    correct_traj_fields['positions'] = correct_rep_positions

    return correct_traj_fields

def get_model(
        model_name,
        gexp,
        ext="jl.pkl",
):

    model_path = data_path() / f'models/{model_name}/{gexp}.{ext}'

    return load_obj(model_path)

def ls_h5_observables(gexp):
    """List the observables that are in the HDF5 file."""

    wepy_h5 = get_gexp_wepy_h5(gexp)
    with wepy_h5:
        obs_names = wepy_h5.observable_field_names

    print('\n'.join(obs_names))

    return obs_names


def ls_dir_observables(gexp):
    """List the observables that are in the file system as a pickle
    file."""

    import os
    import os.path as osp

    observables_dir = osp.join(data_path(), 'observables')

    obs_names = []
    for f in os.listdir(observables_dir):
        obs_names.append(f.split('.')[0])

    print('\n'.join(obs_names))

    return obs_names

def get_observable_h5(
        observable_name,
        gexp,
):

    wepy_h5 = get_gexp_wepy_h5(gexp)

    observable_key = f'observables/{observable_name}'

    with wepy_h5:

        # the reshaped run data initialized for the runs and trajs
        data = [[[] for traj_idx in range(wepy_h5.num_run_trajs(run_idx))]
                for run_idx in range(wepy_h5.num_runs)]

        for traj_id, traj_fields in wepy_h5.iter_trajs_fields(
                [observable_key],
                idxs=True,
        ):

            run_idx = traj_id[0]
            traj_idx = traj_id[1]

            # add the observable to the data, we destructure since
            # traj_fields can have more than one field
            data[run_idx][traj_idx] = traj_fields[observable_key]

    return data

def get_observable_fs(
        observable_name,
        gexp,
):

    obs_dir = data_path() / f'observables/{observable_name}'

    # check the folder for the file you want to load depending on the
    # ligand ID, we will use the extension to figure out how to load
    # it
    matches = []
    for f in os.listdir(obs_dir):
        if f.startswith(f'gexp-{gexp}'):
            matches.append(f)

    if len(matches) > 1:
        raise ValueError("Multiple matches")
    elif len(matches) < 1:
        raise ValueError("No match")
    else:
        fname = matches[0]

        obs_path = obs_dir / fname

        return load_obj(obs_path)

def get_observable(
        observable_name,
        gexp,
        source="h5",
):
    """Get an observable in a standard format. Specify 'h5' or 'fs' to get
    it from the HDF5 file or the filesystem."""

    if source == 'h5':
        return get_observable_h5(
            observable_name,
            gexp,
        )

    elif source == 'fs':
        return get_observable_fs(
            observable_name,
            gexp,
        )

def ls_clustering_models():
    """List the observables that are in the file system as a pickle
    file."""

    observables_dir = data_path() / 'clustering/models'

    obs_names = []
    for f in os.listdir(observables_dir):
        obs_names.append(f)

    print('\n'.join(obs_names))

    return obs_names



def get_clustering_model(
        clf_id,
        gexp,
        ext='jl.pkl',
):

    model_path = data_path() / f'clustering/models/clfid-{clf_id}/gexp-{gexp}.{ext}'

    return load_obj(model_path)


def get_clustering_traintest(
        clf_id,
        gexp,
        ext='jl.pkl',
):

    traintest_dir = data_path() / f'clustering/train_test/clfid-{clf_id}'
    traintest_path = traintest_dir / f'gexp-{gexp}.{ext}'

    # WARN: for some reason its always saving this as a tuple, so we
    # just strip it off
    return load_obj(traintest_path)[0]

def ls_networks():
    """List the observables that are in the file system as a pickle
    file."""

    import os
    import os.path as osp

    observables_dir = osp.join(data_path(), 'csn/network')

    obs_names = []
    for model in os.listdir(observables_dir):
        model_names = []
        for lag_time in os.listdir(osp.join(observables_dir, model)):
            model_spec = "{}/{}".format(model, lag_time)
            obs_names.append(model_spec)

    print('\n'.join(obs_names))

    return obs_names



def get_msn(
        csn_id,
        gexp,
        tag=None,
        ext='jl.pkl'
):

    directory = data_path() / f'csn/network/csnid-{csn_id}'
    path = directory / f'gexp-{gexp}_tag-{tag}.{ext}'

    # SNIPPET: I don't think I need this
    # # check the folder for the file you want to load depending on the
    # # ligand ID, we will use the extension to figure out how to load
    # # it
    # matches = []
    # for f in os.listdir(model_dir):
    #     if f.startswith(str(lig_id)):
    #         matches.append(f)

    # if len(matches) > 1:
    #     raise ValueError("Multiple matches")
    # elif len(matches) < 1:
    #     raise ValueError("No match")
    # else:
    #     fname = matches[0]

    # model_path = osp.join(model_dir, fname)

    return load_obj(path)

def get_h5_msn(
        csn_id,
        gexp,
        tag=None,
        ext='jl.pkl'
):
    """Get a full MSN backed by the HDF5 so you have access to microstate data."""

    from wepy.analysis.network import MacroStateNetwork

    # get the base MSN
    base_msn = get_msn(
        csn_id,
        gexp,
        tag=tag,
        ext=ext,
    )

    # get the HDF5
    contigtree = get_contigtree(gexp)

    msn = MacroStateNetwork(
        contigtree,
        base_network=base_msn,
    )

    return msn

def load_gephi_gexf(path):
    """Load and clean a gexf file generated by gephi."""

    import io

    from networkx.readwrite.gexf import GEXFReader


    # load the gexf and then clean it, and use a temporary buffer for
    # it
    with open(path) as rf:
        gexf_buf = io.StringIO(clean_gephi_gexf(rf.read()))

    return GEXFReader(node_type=int)(gexf_buf)


def load_gephi_graph(
        csn_id,
        gexp,
        tag=None,
        layout_id='main',
):
    """Load the networkx graph from the standard gephi generated gexf
    file."""

    directory = data_path() / f'csn/gexf/csnid-{csn_id}'
    model_path = directory / f'gexp-{gexp}_tag-{tag}_layout-{layout_id}.gexf'

    # NOTE: this was how I did it before, thats alright
    # fname = '{}_{}.gephi.gexf'.format(lig_id, name_pattern)

    if not osp.exists(model_path):
        raise OSError("Gexf does not exist")

    return load_gephi_gexf(model_path)

# SNIPPET: this was doing too much so removing it
# def get_gephi_network(
#         msm_id,
#         gexp,
#         csn_id='main',
#         layout_id='main',
# ):
#     """Load the gephi (CSN) network and combine with the network we
#     produce, so that the data is combined (incorporating changes from
#     gephi editing)

#     """

#     graph = load_gephi_graph(
#         msm_id,
#         gexp,
#         csn_id=csn_id,
#         layout_id=layout_id,
#     )

#     net = get_msn(
#         msm_id,
#         gexp,
#         csn_id=csn_id,
#     )

#     # merge them
#     return update_network_from_gephi(
#         net,
#         graph,
#         layout_name='main',
#     )

def update_network_layout_from_gephi(
        net,
        graph,
        layout_name=None,
):

    node_values = {node : graph.nodes[node]['viz']
                   for node in graph.nodes}

    net.set_nodes_layout(layout_name, node_values)

    return net


def update_network_from_gephi(
        net,
        graph,
        fields=None,
        layout_name=None,
):
    """This takes a gexf networkx object and updates the net/msn with
    the special data from this file.

    This is for if you want to add data from the gexf to the network.

    If fields is specified import only those fields.

    WARNING
    -------

    This will overwrite all data. If you just want the layout use the
    other function.

    """

    ignore_fields = ('label',)

    if fields is None:
        node = list(graph.nodes.keys())[0]
        fields = graph.nodes[node].keys()
        del node


    # assume all nodes have the same fields
    for field in fields:

        # don't update columns which are the same name, we assume they are
        # immutable and if there is new data it has a new name
        if (field in net.graph.nodes[node] or
            field in ignore_fields or
            (field == 'viz' and layout_name is None)
        ):
            continue

        # get all values dict for this attribute
        node_values = {node : graph.nodes[node][field]
                       for node in graph.nodes}


        # we treat the layout data differently, and save them in
        # layouts structure
        if field == 'viz':

            net.set_nodes_layout(layout_name, node_values)

        elif field.startswith('_groups'):

            # if it is a node grouping we also need to get the indices for
            # which the group was true
            group_nodes = [node_id for node_id, node_value in node_values.items()
                           if node_value is True]

            group_name = '/'.join(field.split('/')[1:])

            net.set_node_group(group_name, group_nodes)

        elif field.startswith('_observables'):

            observable_name = '/'.join(field.split('/')[1:])

            net.set_nodes_observable(observable_name, node_values)

        else:
            net.set_nodes_attribute(field, node_values)

    return net

def get_network_nodes(model_name, lig_id, lag_time):

    import os
    import os.path as osp

    model_dir = osp.join(data_path(),
                         'csn/nodes/{}/{}'.format(model_name, lag_time))

    # check the folder for the file you want to load depending on the
    # ligand ID, we will use the extension to figure out how to load
    # it
    matches = []
    for f in os.listdir(model_dir):
        if f.startswith(str(lig_id)):
            matches.append(f)

    if len(matches) > 1:
        raise ValueError("Multiple matches")
    elif len(matches) < 1:
        raise ValueError("No match")
    else:
        fname = matches[0]

    model_path = osp.join(model_dir, fname)

    return load_table(model_path)



def get_network_edges(model_name, lig_id, lag_time):

    import os
    import os.path as osp

    model_dir = osp.join(data_path(),
                         'csn/edges/{}/{}'.format(model_name, lag_time))

    # check the folder for the file you want to load depending on the
    # ligand ID, we will use the extension to figure out how to load
    # it
    matches = []
    for f in os.listdir(model_dir):
        if f.startswith(str(lig_id)):
            matches.append(f)

    if len(matches) > 1:
        raise ValueError("Multiple matches")
    elif len(matches) < 1:
        raise ValueError("No match")
    else:
        fname = matches[0]

    model_path = osp.join(model_dir, fname)

    return load_table(model_path)

def load_msm(
        msm_id,
        gexp,
        ext="jl.pkl",
):

    model_dir = data_path() / f'msm/models/msmid-{msm_id}'
    path = model_dir / f'gexp-{gexp}.{ext}'

    pyemma_msm, trimming_mapping = load_obj(path)

    return pyemma_msm, trimming_mapping

def load_tspca_model(
        tspca_id,
        ext="jl.pkl",
):

    model_dir = data_path() / f'ts_pca/models'
    path = model_dir / f'tspca-id-{tspca_id}.{ext}'

    pca_model, model_score, mode_scores, test_size, gexps = load_obj(path)

    return pca_model, model_score, mode_scores, test_size, gexps

def save_model(
        model_name,
        gexp,
        model,
        overwrite=True,
        ext="jl.pkl",
):

    model_dir = data_path() / f'models/{model_name}'
    model_path = model_dir / '{gexp}.{ext}'

    os.makedirs(model_dir, exist_ok=True)

    save_obj(
        model_path,
        model,
        overwrite=overwrite,
        ext=ext,
    )

    return model_path

def get_obs_path(
        observable_name,
        gexp,
        ext="jl.pkl",
):

    all_obs_dir = data_path() / 'observables'

    obs_dir = all_obs_dir / observable_name

    obs_path = obs_dir / f'gexp-{gexp}.{ext}'

    return obs_path

def save_observable(
        observable_name,
        gexp,
        obs_data,
        overwrite=True,
        ext="jl.pkl",
):

    obs_path = get_obs_path(observable_name,
                            gexp,
                            ext=ext)

    save_obj(obs_path,
             obs_data,
             overwrite=overwrite,
             ext=ext)

    return obs_path


def attach_observable(
        observable_name,
        gexp,
):
    """Save one of the saved observables from the filesystem to the
    corresponding WepyHDF5.
    """

    with get_gexp_wepy_h5(gexp, mode='r+') as wepy_h5:

        # get the observable (only from the filesystem), and add it to
        # the wepyHDF5

        wepy_h5.add_observable(
            observable_name,
            get_observable(
                observable_name,
                gexp,
                source='fs',
            )
        )

def save_external_observable(observable_name,
                             gexp,
                             overwrite=True,
                             ext="jl.pkl",
):

    """Save one of the observables from the HDF5 to the filesystem."""

    with get_gexp_wepy_h5(gexp) as wepy_h5:

        save_observable(
            observable_name,
            gexp,
            get_observable(
                observable_name,
                gexp,
                source='h5',
            ),
            overwrite=overwrite,
            ext=ext,
        )

def save_clustering_model(
        clf_id,
        gexp,
        model,
        overwrite=True,
        ext="jl.pkl",
):

    import os
    import os.path as osp

    model_dir = data_path() / f'clustering/models/clfid-{clf_id}'
    model_path = model_dir / f'gexp-{gexp}.{ext}'

    os.makedirs(model_dir, exist_ok=True)

    save_obj(
        model_path,
        model,
        overwrite=overwrite,
        ext=ext,
    )

    return model_path

def save_clustering_traintest(
        clf_id,
        gexp,
        train_idxs,
        test_idxs,
        overwrite=True,
        ext='jl.pkl',
):

    traintest_dir = data_path() / f'clustering/train_test/clfid-{clf_id}'
    traintest_path = traintest_dir / f'gexp-{gexp}.{ext}'

    os.makedirs(traintest_dir, exist_ok=True)

    train_test_d = {
        'training_idxs' : train_idxs,
        'testing_idxs' : test_idxs,
    },


    save_obj(
        traintest_path,
        train_test_d,
        overwrite=overwrite,
        ext=ext,
    )

    return traintest_path

def save_cluster_centers(
        msm_id,
        gexp,
        traj,
        overwrite=True,
        ext='dcd',
):

    traj_dir = data_path() / f'clustering/cluster_center_trajs/msm-{msm_id}'
    traj_path = traj_dir / f'gexp-{gexp}.{ext}'

    os.makedirs(traj_dir, exist_ok=True)

    traj.save(str(traj_path), force_overwrite=overwrite)

    return traj_path

def save_span_warp_lineages_dcds(gexp, span_id):
    """Save the warps as individual trajectory files named by the index of
    the warp event.

    WARNING: probably you will want to use `save_gexp_warp_lineages`
    so you can match that to the the table of warping.

    """

    lineages_dir_path = data_path() / f'warp-lineages_by-span/gexp-{gexp}/span-{span_id}_warp-lineages'

    lineages_dir_path.mkdir(exist_ok=True, parents=True)

    contig = get_gexp_span_contig(gexp, span_id)

    lineage_trajs_gen = get_warp_lineage_trajs(gexp, span_id)

    total_frames = 0
    for warp_idx, lineage_traj in enumerate(lineage_trajs_gen):

        n_frames = len(lineage_traj)
        total_frames += n_frames

        print(f"warp: {warp_idx} has {n_frames} frames")

        path = lineages_dir_path / f'warp-{warp_idx}.dcd'

        print("saving")
        lineage_traj.save_dcd(str(path))

    print(f"total frames: {total_frames}")


def save_gexp_warp_lineages_dcds(gexp):
    """Save the warps as individual trajectory files named by the index of
    the warp event.

    WARNING: probably you will want to use `save_gexp_warp_lineages`
    so you can match that to the the table of warping.

    """

    lineages_dir_path = data_path() / f'warp-lineages_by-gexp/gexp-{gexp}'

    # /span-{span_id}_warp-lineages
    filename_tmpl = "{warp_idx}_span-{span_id}_{span_warp_idx}"

    lineages_dir_path.mkdir(exist_ok=True, parents=True)

    # collect all of the generators for each span
    lineages_gens = []
    for span_id in get_gexp_span_ids(gexp):

        contig = get_gexp_span_contig(gexp, span_id)

        lineage_trajs_gen = get_warp_lineage_trajs(gexp, span_id)

        lineages_gens.append(lineage_trajs_gen)

    # as we go along increment the warp index to keep track
    warp_idx = 0
    for span_id, lineage_gen in enumerate(lineages_gens):
        for span_warp_idx, lineage_traj in enumerate(lineage_gen):

            n_frames = len(lineage_traj)

            print(f"warp: {warp_idx} span: {span_id}:{span_warp_idx} has frames: {n_frames}")

            filename_root = filename_tmpl.format(
                warp_idx=warp_idx,
                span_id=span_id,
                span_warp_idx=span_warp_idx
            )

            path = lineages_dir_path / f'{filename_root}.dcd'

            print("saving")
            lineage_traj.save_dcd(str(path))

            # increment the warp index
            warp_idx += 1



# UGLY,SNIPPET: this did them all in one file which just doesn't work
# unfortunately since its just too big

# def save_warp_lineages_dcd(gexp, span_id):
#     """Save all the warps to the same DCD concatenated together"""

#     lineages_path = data_path() / f'warp-lineages/gexp-{gexp}/span-{span_id}_warp-lineages.dcd'

#     lineages_path.parent.mkdir(exist_ok=True, parents=True)

#     lineages_trace = it.chain(*get_warp_lineage_traces(gexp, span_id))

#     print(f"total frames: {len(lineage_trace)}")

#     contig = get_gexp_span_contig(gexp, span_id)

#     print("Loading")
#     with contig:
#         traj = contig.wepy_h5.traj_fields_to_mdtraj(
#             contig.wepy_h5.get_trace_fields(lineages_trace,
#                                             ['positions', 'box_vectors']
#             ))

#     print("saving DCD")
#     traj.save_dcd(str(lineages_path))

def save_gexp_final_lineages_dcds(gexp):
    """Save the final walkers as individual trajectory files.

    """

    lineages_dir_path = data_path() / f'final-walker-lineages_by-gexp/gexp-{gexp}'

    # /span-{span_id}_warp-lineages
    filename_tmpl = "span-{span_id}_{span_final_idx}"

    lineages_dir_path.mkdir(exist_ok=True, parents=True)

    # collect all of the generators for each span
    lineages_gens = []
    for span_id in get_gexp_span_ids(gexp):

        contig = get_gexp_span_contig(gexp, span_id)

        lineage_trajs_gen = get_final_lineage_trajs(gexp, span_id)

        lineages_gens.append(lineage_trajs_gen)

    for span_id, lineage_gen in enumerate(lineages_gens):
        for span_final_idx, lineage_traj in enumerate(lineage_gen):

            n_frames = len(lineage_traj)

            print(f"span: {span_id}:{span_final_idx} has frames: {n_frames}")

            filename_root = filename_tmpl.format(
                span_id=span_id,
                span_final_idx=span_final_idx
            )

            path = lineages_dir_path / f'{filename_root}.dcd'

            print("saving")
            lineage_traj.save_dcd(str(path))

def save_gexp_high_progress_walkers(
        gexp,
        top_n,
        progress_key='min_distances',
  ):
    """Save the highest progress walkers to trajectory files for:

    - each span
    - across each gexp

    This will get the highest 'top_n' from each span. The full gexp
    one will be all of those.

    """

    from wepy.util.util import concat_traj_fields
    from wepy.util.mdtraj import traj_fields_to_mdtraj

    dir_path = data_path() / f'high-progress_walkers/gexp-{gexp}'

    dir_path.mkdir(exist_ok=True, parents=True)

    suffix = f"top-{top_n}_prog-{progress_key}"

    span_filename_tmpl = f"span-{{span_id}}_{suffix}"

    gexp_filename = f"all_{suffix}"

    span_traj_fields = {}
    for span_id in get_gexp_span_ids(gexp):

        traj_fields = get_high_progress_traj_fields(
                    gexp,
                    span_id,
                    top_n,
                    progress_key,
                )

        ## then save it to the mapping for writing them all out later
        span_traj_fields[span_id] = traj_fields

        # make a traj for this one to save individually

        with contig:
            span_traj = contig.wepy_h5.traj_fields_to_mdtraj(traj_fields)

        del traj_fields
        gc.collect()

        ## write the span traj
        span_filename = span_filename_tmpl.format(
            span_id=span_id,
        )
        span_path = dir_path / f"{span_filename}.dcd"

        print(f"saving for span {span_id}")

        span_traj.save_dcd(str(span_path))


        del span_traj
        gc.collect()

    contigtree = get_contigtree(gexp)

    with contigtree:
        # write out all of them
        all_traj = contigtree.wepy_h5.traj_fields_to_mdtraj(
            concat_traj_fields(
                [span_traj_fields[span_id]
                 for span_id in  get_gexp_span_ids(gexp)]
            ))

    del span_traj_fields
    gc.collect()

    print(f"saving all traj")

    all_path = dir_path / f"{gexp_filename}.dcd"

    all_traj.save_dcd(str(all_path))

def save_gexp_high_progress_walkers_lineages_dcds(
        gexp,
        top_n,
        progress_key='min_distances',
):
    """Save lineages of the high progress walkers to dcds

    """

    lineages_dir_path = data_path() / f'high-progress_lineages_by-gexp/gexp-{gexp}'

    filename_tmpl = "span-{span_id}_{span_lineage_idx}"

    lineages_dir_path.mkdir(exist_ok=True, parents=True)

    # collect all of the generators for each span
    lineages_gens = []
    for span_id in get_gexp_span_ids(gexp):


        contig = get_gexp_span_contig(gexp, span_id)

        lineage_trajs_gen = get_high_progress_lineages_trajs(
            gexp,
            span_id,
            top_n,
            progress_key)

        lineages_gens.append(lineage_trajs_gen)

    for span_id, lineage_gen in enumerate(lineages_gens):
        for span_lineage_idx, lineage_traj in enumerate(lineage_gen):

            n_frames = len(lineage_traj)

            print(f"span: {span_id}:{span_lineage_idx} has frames: {n_frames}")

            filename_root = filename_tmpl.format(
                span_id=span_id,
                span_lineage_idx=span_lineage_idx
            )

            path = lineages_dir_path / f'{filename_root}.dcd'

            print("saving")
            lineage_traj.save_dcd(str(path))

def save_fig(
        group,
        name_root,
        fig,
        tags=None,
):
    """Save a figure of any kind.

    Useful for when all gexps are present or some other snowflake plot.

    See Also
    --------

    save_gexp_fig : helper to do automatic and consistent name
                    formatting for gexp specific plots

    """

    # if there are tags render to string, and add to the root name

    if tags is None:
        name_tagged_root = name_root
    else:
        tag_str = "_".join(
            [f'{key}-{val}' for key, val in tags.items()])

        name_tagged_root = f"{name_root}__{tag_str}"


    # make the path
    fig_stem = media_path() / group / name_tagged_root

    # ensure the directory exists for the group
    fig_stem.parent.mkdir(exist_ok=True, parents=True)

    # save a copy for each extension
    for ext in FIG_EXTENSIONS:

        fig_path = f"{fig_stem}.{ext}"

        print(f"PLOTTING: {fig_path}")

        fig.savefig(fig_path)

    return fig_stem


def save_gexp_fig(gexp,
                  group,
                  fig,
                  tags=None):
    """Save a figure to media group

    group is a grouping of figures that are the same kind. For
    instance if you have plots of rates, you might call this 'rates'
    or 'fe/spans'.

    Each gexp will get it's own figure in this group in the different
    formats etc.

    If there are any tags, these will be added to the end of the
    filename stem as '{key}-{value}' separated by '_' characters.

    """

    name_root = f"gexp-{gexp}"

    return save_fig(
        group,
        name_root,
        fig,
        tags=tags,
    )

def save_csn_gexf(
            csn_id,
            gexp,
            msn,
            tag=None,
            layout_id='main',
            overwrite=True,
):

    directory = data_path() / f'csn/gexf/csnid-{csn_id}'
    path = directory / f'gexp-{gexp}_tag-{tag}_layout-{layout_id}.gexf'

    os.makedirs(directory, exist_ok=True)

    # don't short circuit to the None layout, this is confusing
    if layout_id is None or layout_id == "None":
        msn.write_gexf(path, layout=None)

    else:
        msn.write_gexf(path, layout=layout_id)

    return path

def save_csn(
            csn_id,
            gexp,
            msn,
            tag=None,
            overwrite=True,
            ext='jl.pkl'):

      directory = data_path() / f'csn/network/csnid-{csn_id}'
      path = directory / f'gexp-{gexp}_tag-{tag}.{ext}'

      os.makedirs(directory, exist_ok=True)

      save_obj(
            path,
            msn,
            overwrite=overwrite,
            ext=ext
      )

      return path

def save_csn_nodes(
            csn_id,
            gexp,
            df,
            tag=None,
            overwrite=True,
            ext='csv',
):

      directory = data_path() / f'csn/nodes/csnid-{csn_id}'
      path = directory / f'gexp-{gexp}_tag-{tag}.{ext}'

      os.makedirs(directory, exist_ok=True)

      save_table(
            path,
            df,
            overwrite=overwrite,
            ext=ext,
      )

      return path

def save_csn_edges(
            csn_id,
            gexp,
            df,
            tag=None,
            overwrite=True,
            ext='csv',
):

      directory = data_path() / f'csn/edges/csnid-{csn_id}'
      path = directory / f'gexp-{gexp}_tag-{tag}.{ext}'

      os.makedirs(directory, exist_ok=True)

      save_table(
            path,
            df,
            overwrite=overwrite,
            ext=ext,
      )

      return path

def save_csn_tables(
            csn_id,
            gexp,
            msn,
            tag=None,
            overwrite=True,
            ext='csv',
):

      nodes_df = msn.nodes_to_dataframe()
      edges_df = msn.edges_to_dataframe()

      save_csn_nodes(
            csn_id,
            gexp,
            nodes_df,
            tag=tag,
            overwrite=overwrite,
            ext=ext,
      )

      save_csn_edges(
            csn_id,
            gexp,
            edges_df,
            tag=tag,
            overwrite=overwrite,
            ext=ext,
      )

def save_all_csn_stuff(
            csn_id,
            gexp,
            msn,
            tag=None,
            layout_id=Ellipsis,
            overwrite=True,
):
    """Save all output files for the relevant MSM.

    There is an additional 'csn_id' which can be used to save
    different CSNs with different attached data, so you don't have
    to overwrite other ones.

    For the 'tag' the layout_id if the value is 'None' then it will
    be saved to a default 'nil' file.

    If the values are Ellipsis then all of the available IDs will be
    saved.

    """


    # save the MSN as a pickle object, this will have all the data
    save_csn(
        csn_id,
        gexp,
        msn,
        tag=tag,
    )

    # save the data from the network as both an edge and a node
    # table CSVs
    save_csn_tables(
        csn_id,
        gexp,
        msn,
        tag=tag,
    )

    # save the gexf XML file with a layout for visualization
    save_csn_gexf(
        csn_id,
        gexp,
        msn,
        tag=tag,
        layout_id=layout_id,
    )

def save_csn_group_traj(model_name, lig_id,
                        net, group_name,
                        downsample=None,
                        overwrite=True, ext='dcd'):

    import os
    import os.path as osp

    # transform the name to a file name
    grp_filename = group_name.replace('/', '-')

    traj_dir = osp.join(data_path(), 'csn/group_trajs/{}'.format(model_name))
    traj_path = osp.join(traj_dir, '{}_{}.{}'.format(lig_id, grp_filename, ext))

    os.makedirs(traj_dir, exist_ok=True)

    # make the traj, downsampling if requested
    traj = None
    for node_id in net.node_groups[group_name]:
        print("working on node {}".format(node_id))
        if traj is None:
            traj = get_cluster_micro_traj(model_name, lig_id, node_id,
                                          downsample=downsample)
        else:
            traj += get_cluster_micro_traj(model_name, lig_id, node_id,
                                           downsample=downsample)

        print("{} frames".format(traj.n_frames))

    print("saving traj")
    traj.save(traj_path, force_overwrite=overwrite)

    return traj_path

def save_msm(
        msm_id,
        gexp,
        msm,
        trimming_mapping,
        overwrite=True,
        ext="jl.pkl",
):

    import os

    model_dir = data_path() / f'msm/models/msmid-{msm_id}'
    model_path = model_dir / f'gexp-{gexp}.{ext}'

    os.makedirs(model_dir, exist_ok=True)


    # save the MSM and the mapping together
    obj = (msm, trimming_mapping)

    save_obj(
        model_path,
        obj,
        overwrite=overwrite,
        ext=ext,
    )

    return model_path

def save_tspca_model(
        tspca_id,
        pca_model,
        model_score,
        mode_scores,
        test_size,
        gexps,
        overwrite=True,
        ext="jl.pkl",
):

    model_dir = data_path() / f'ts_pca/models'
    model_path = model_dir / f'tspca-id-{tspca_id}.{ext}'

    os.makedirs(model_dir, exist_ok=True)

    # save the MSM and the mapping together
    obj = (pca_model, model_score, mode_scores, test_size, gexps,)

    save_obj(
        model_path,
        obj,
        overwrite=overwrite,
        ext=ext,
    )

    return model_path


def save_tspca_score_table(
        score_df,
):

    table_path = data_path() / f'ts_pca/model_scores/scores.csv'

    save_table(
        table_path,
        score_df,
        overwrite=True,
    )

def save_com_trajs(
        ts_id,
        gexp,
        com_traj,
        pc_projections,
        overwrite=True,
):

    dir_path = data_path() / f'ts_pca/com_trajs/tsid-{ts_id}'

    os.makedirs(dir_path, exist_ok=True)

    # write these trajectories out. Make one for each PC and put these
    # as the bfactors so we can color them in a visualizer

    # ensure the directory
    os.makedirs(
        f'data/ts_pca/gexp-{gexp}/ts-{ts_id}',
        exist_ok=True,
    )

    for pc_idx, pc_proj in enumerate(pc_projections):

        fname = f'coms_pc-{pc_idx}.pdb'

        fpath = dir_path / fname

        com_traj.save_pdb(
            str(fpath),
            bfactors=pc_proj
        )


    return dir_path

def single_atom_json_rec(atom_idx):
    atom_rec = {'index' : atom_idx,
                'name' : 'H',
                'element' : 'He'}

    return atom_rec

def n_atom_top(n_atoms):

    import json

    n_atom_top_d = {'bonds' : [],
                         'chains' : [
                             {'index' : 0,
                              'residues' : [
                                  {'index' : 0,
                                   'name' : 'DUM',
                                   'resSeq' : 0,
                                   'segmentID' : 'DUMA',
                                   'atoms' :
                                   [single_atom_json_rec(i) for i in range(n_atoms)]}
                              ]}
                         ]}
    n_atom_top = json.dumps(n_atom_top_d)

    return n_atom_top

def n_atom_mdj_top(n_atoms):

    from wepy.util.mdtraj import json_to_mdtraj_topology

    return json_to_mdtraj_topology(n_atom_top(n_atoms))


def save_n_atom_top_pdb(n_atoms):

      import os
      import os.path as osp

      import numpy as np

      import mdtraj as mdj

      top_dir = osp.join(data_path(), 'top/util')

      # make sure the directory exists
      os.makedirs(top_dir,
                  exist_ok=True)

      pdb_path = osp.join(top_dir, "{}_atom_ref.pdb".format(n_atoms))

      ref_traj = mdj.Trajectory(np.zeros((1,1,3)),
                                n_atom_mdj_top(n_atoms))

      ref_traj.save_pdb(pdb_path)

def convert_rate_to_rt(rate_q, time_base=tkunit.minute):

    return (np.float64(1.0) / rate_q.value_in_unit(rate_q.unit)) * (1/rate_q.unit).in_units_of(time_base)


def convert_rt_to_rate(rt_q, time_base=tkunit.second):

    time_unit = (1/time_base).unit

    rate = (1 / rt_q).in_units_of(time_unit)

    return rate

### Plotting Functions

def bin_centers_from_edges(bin_edges):
    return np.array([(bin_edges[i] + (bin_edges[i + 1] - bin_edges[i]))
                            for i in range(bin_edges.shape[0] - 1)])


def move_axes(ax, fig, subplot_spec=111):
    """Move an Axes object from a figure to a new pyplot managed Figure in
    the specified subplot."""

    # get a reference to the old figure context so we can release it
    old_fig = ax.figure

    # remove the Axes from it's original Figure context
    ax.remove()

    # set the pointer from the Axes to the new figure
    ax.figure = fig

    # add the Axes to the registry of axes for the figure
    fig.axes.append(ax)
    # twice, I don't know why...
    fig.add_axes(ax)

    # then to actually show the Axes in the new figure we have to make
    # a subplot with the positions etc for the Axes to go, so make a
    # subplot which will have a dummy Axes
    dummy_ax = fig.add_subplot(subplot_spec)

    # then copy the relevant data from the dummy to the ax
    ax.set_position(dummy_ax.get_position())

    # then remove the dummy
    dummy_ax.remove()

    # close the figure the original axis was bound to
    plt.close(old_fig)

def recenter_superimpose_traj(traj_fields, lig_id, rep_key):
    """Recenter the coordinates based on the ligand id and the
    representation key. Fields must only be 'positions' and
    'box_vectors'.

    """

    import numpy as np
    from wepy.util.util import traj_box_vectors_to_lengths_angles

    from geomm.superimpose import superimpose
    from geomm.grouping import group_pair
    from geomm.centering import center_around

    sel_idxs = lig_selection_idxs(lig_id)

    lig_idxs = sel_idxs['{}/ligand'.format(rep_key)]
    prot_idxs = sel_idxs['{}/protein'.format(rep_key)]
    bs_idxs = sel_idxs['{}/binding_site'.format(rep_key)]

    centered_traj_fields = get_centered_ref_state_traj_fields(lig_id, rep_key=rep_key)

    centered_ref_positions = centered_traj_fields['positions']

    box_lengths, _ = traj_box_vectors_to_lengths_angles(traj_fields['box_vectors'])

    ## regroup, center, and superimpose the frames

    # group the pair of ligand and binding site together in the same image
    grouped_positions = [group_pair(positions, box_lengths[idx],
                                        bs_idxs, lig_idxs)
                  for idx, positions in enumerate(traj_fields['positions'])]

    # center all the positions around the binding site
    centered_positions = [center_around(positions, bs_idxs)
                          for idx, positions in enumerate(grouped_positions)]

    # then superimpose the binding sites
    sup_positions = [superimpose(centered_ref_positions, pos, idxs=bs_idxs)[0]
                     for pos in centered_positions]

    return np.array(sup_positions), centered_ref_positions

def clean_gephi_gexf(gexf):

    from xmltodict import parse, unparse

    # parse the bad gexf file
    bad_gexf = parse(gexf)

    # the good tag we take things from
    good_gexf_tag = """
    <gexf version="1.2"
          xmlns:viz="http://www.gexf.net/1.2draft/viz"
          xmlns="http://www.gexf.net/1.2draft"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://www.w3.org/2001/XMLSchema-instance">
    </gexf>
    """

    # parse the good tag and get the values
    tag_attrs = {}
    for key, value in parse(good_gexf_tag)['gexf'].items():

        if not key.startswith('@'):
            continue
        else:
            tag_attrs[key] = value

    # read the bad xml and parse it

    # then replace the tag attributes for this with the good ones
    for key, value in tag_attrs.items():
        bad_gexf['gexf'][key] = value

    # then get the good xml string
    return unparse(bad_gexf)

def get_obs_shape(gexp):

    with get_gexp_wepy_h5(gexp) as wepy_h5:

        runs = []
        for run_idx in range(wepy_h5.num_runs):
            traj_lengths = []
            for traj_idx in range(wepy_h5.num_run_trajs(run_idx)):
                traj_lengths.append(wepy_h5.num_traj_frames(run_idx, traj_idx))

            runs.append(traj_lengths)

    return runs

def get_obs_traj_shape(gexp):

    with get_gexp_wepy_h5(gexp) as wepy_h5:

        traj_lengths = []
        for run_idx in range(wepy_h5.num_runs):
            for traj_idx in range(wepy_h5.num_run_trajs(run_idx)):
                traj_lengths.append(wepy_h5.num_traj_frames(run_idx, traj_idx))

    return traj_lengths

def are_features_scalar(features):

    return np.isscalar(features[0])

def scalarize_features(features):

    # if they are already scalars don't do anything
    if are_features_scalar(features):
        return features

    # first check that the features can be scalar
    feature_shape = len(features[0])
    assert features_shape == 1, f"Feature is length {features_shape} cannot be scalarized"

    n_features = len(features)

    return np.reshape(
        np.array(features),
        (n_features,)
    )

def unscalarize_features(features):

    # if its already not a scalar don't do anything
    if not are_features_scalar(features):
        return features

    n_features = len(features)

    return np.reshape(
        np.array(features),
        (n_features, 1)
    )

# TODO
def reshape_run_to_traj(run_obs, gexp):
    raise NotImplemented

def reshape_run_to_features(run_obs):

    # just a double concatenate

        return np.concatenate(
                 [np.concatenate(
                     [traj for traj in run])
                  for run in run_obs]
        )

# TODO
def reshape_traj_to_run(traj_obs, gexp):
    raise NotImplemented

def reshape_traj_to_features(traj_obs):

    # just a flattening
    return np.concatenate(traj_obs)

def reshape_features_to_run(features,
                            gexp,
                            run_idxs=None,
):

    # this is the "run" shape for all of the runs
    obs_shape = get_obs_shape(gexp)

    # if the run_idxs are given assume the features are for these runs
    # only.

    # otherwise we generate the run_idxs from the whole contigtree
    # structure
    if run_idxs is None:
        run_idxs = list(range(len(obs_shape)))

    start_idx = 0

    obs_runs = []
    for run_idx in run_idxs:

        # for the run_idx get the sub-shape structure for the run
        run = obs_shape[run_idx]

        obs_trajs = []
        for num_traj_frames in run:

            obs_trajs.append(features[start_idx : start_idx + num_traj_frames])

            start_idx += num_traj_frames

        obs_runs.append(obs_trajs)

    return obs_runs

def reshape_features_to_traj(features, gexp):

    traj_lengths = get_obs_traj_shape(gexp)

    start_idx = 0

    obs_trajs = []
    for num_frames in traj_lengths:

        obs_trajs.append(features[start_idx : start_idx + num_frames])

        start_idx += num_frames

    return obs_trajs

# TODO: move this to the section with other run-> feature conversions
def feature_idx_to_trace_idx(
        feature_idx,
        gexp,
):
    """For a feature major structure get the trace index (run, traj,
    cycle) from the feature index for the splitting/shuffling used in
    that classifier.

    """

    run_feature_idxs = make_feature_idx_observable(gexp)

    count = 0
    for run_idx, run in enumerate(run_feature_idxs):
        for traj_idx, traj in enumerate(run):
            if feature_idx in traj:
                frame_idx = np.searchsorted(traj, feature_idx)

                return (run_idx, traj_idx, frame_idx)

def make_feature_idx_observable(gexp):

    run_obs = []
    start_idx = 0
    for run in get_obs_shape(gexp):

        traj_obs = []
        for num_frames in run:

            traj_obs.append(list(range(start_idx, start_idx + num_frames)))

            start_idx += num_frames

        run_obs.append(traj_obs)

    return run_obs

#  translation of my names to the mdtraj names

MDJ_COL_NAMES_TRANSLATION = {
    'resname' : 'resName',
    'resid' : 'resSeq',
    'pdb_serial' : 'serial',
    'segid' : 'segmentID',
    'chain' : 'chainID',
}


BACKBONE_ATOM_NAMES = (
    'C',
    'CA',
    'N',
    'O',
)

### Names & Path Templates For Common Structures

# ROIs means Regions of Interest

# The reference structures are standard across all ligands and
# non-specific
REFERENCE_STRUCTURES = {

    # this is the raw structure from the RCSB crystal
    # structure PDB
    'raw_crystal' : {
        'pdb_path' : "data/seh/rcsb/4od0.pdb",
        'fasta_path' : "data/seh/rcsb/4od0.fasta.txt",
    },

    # this is the structure after some minor things were fixed with
    # PDBFixer, it includes the whole protein and the crystallographic
    # waters for the 'active_domain' ROI. This is used for assembling
    # the final complex for simulation. It is not used for docking.
    'fixed_crystallographic_waters' : {
        'pdb_path' : "data/seh/pdbfixer/4od0_prot_DOI_water.pdbfixer.pdb"
    },


    # this is just the protein structure after being fixed with
    # PDBFixer. This is the target for docking. It includes both
    # domains.
    'fixed_protein' : {
        'pdb_path' : "data/seh/pdbfixer/4od0_protein.pdbfixer.pdb"
    },

    ## This is the reference structure to use for alignments of other
    ## things.

    # this is the same as 'fixed_protein' but it also has the ligand
    'fixed_complex' : {
        'pdb_path' : "data/seh/pdbfixer/4od0_simplified.pdbfixer.pdb"
    },

}


# ways to select the different ROIs, if the field is there it can be
# used to get it uniquely, if its not it you can't. The access methods
# are ordered such that the first is the default.

REFERENCE_ROIS = {
    'raw_crystal' : {

        'ligand' : (
            ('resname', '2RV'),
            ('resid', (603,)),
            ('pdb_serial', (4336, 4360)),
        ),
        'protein' : (
            ('chain', 0),
            ('pdb_serial', (1, 4328)),
            ('idx', (0, 4327)),
        ),
        'active_domain' : (
            ('resid', (231, 547,)), # this is inclusive of the second index
        ),
        'crystallographic_waters' : (
            ('resname', 'HOH'),
            ('pdb_serial_range', (4361, 4376)),
        ),
        'phosphate' : (
            ('resname', 'PO4'),
        ),
        'magnesium' : (
            ('resname', 'MG'),
        ),
    },

    'fixed_crystallographic_waters' : {

        'ligand' : None,
        'protein' : (
            ('chain', 0),
            ('pdb_serial', (1, 4333)),
            ('idx', (0, 4327)),
        ),
        'active_domain' : (
            ('resid', (231, 547,)), # this is inclusive of the second index
        ),
        'crystallographic_waters' : (
            ('resname', 'HOH'),
            ('pdb_serial_range', (4334, 4343)),
        ),
        'phosphate' : None,
        'magnesium' : None,
    },

    'fixed_protein' : {

        'ligand' : None,
        'protein' : (
            ('chain', 0),
            ('pdb_serial', (1, 4333)),
            ('idx', (0, 4327)),
        ),
        'active_domain' : (
            ('resid', (231, 547,)), # this is inclusive of the second index
        ),
        'crystallographic_waters' : None,
        'phosphate' : None,
        'magnesium' : None,
    },

    'fixed_complex' : {

        'ligand' : (
            ('resname', '2RV',),
        ),
        'protein' : (
            ('chain', 0),
            ('pdb_serial', (1, 4333)),
            ('idx', (0, 4327)),
        ),
        'active_domain' : (
            ('resid', (231, 547,)), # this is inclusive of the second index
        ),
        'crystallographic_waters' : (
            ('resname', 'HOH'),
        ),
        'phosphate' : (
            ('resname', 'PO4'),
        ),
        'magnesium' : (
            ('resname', 'MG'),
        ),
    },

}


## raw loaders for reference structures
def load_ref_structure_pdb(ref_id):
    """Load the raw PDB for the given ref_id"""

    pdb_path = REFERENCE_STRUCTURES[ref_id]['pdb_path']
    traj = mdj.load_pdb(pdb_path)

    return traj


def get_ref_structure_mdj_top(ref_id):
    return load_ref_structure_pdb(ref_id).top


def get_ref_structure_top_df(ref_id):
    return get_ref_structure_mdj_top(ref_id).to_dataframe()[0]

def get_ref_structure_json_top(ref_id):

    from wepy.util.mdtraj import mdtraj_to_json_topology

    return mdtraj_to_json_topology(
        load_ref_structure_pdb(ref_id).top
    )



## methods to get subsets of a reference structure

def ref_selection_idxs(
        ref_id,
        roi,
):
    """Return the idxs of the ROI in the given reference structure"""


    top_df = get_ref_structure_top_df(ref_id)


    # get the key and value for selecting the ROI, use the first one
    # as a default

    key, val = REFERENCE_ROIS[ref_id][roi][0]

    # convert to the mdtraj key
    mdj_key = MDJ_COL_NAMES_TRANSLATION[key]

    # these are both single selections
    if key in ('resname', 'chain', 'segid',):

        idxs = top_df[top_df[mdj_key] == val].index.values

    # these are ranges
    elif key in ('resid', 'pdb_serial', 'idx',):

        start, stop = val

        idxs = top_df[
            (top_df[mdj_key] >= start) &
             (top_df[mdj_key] <= stop)
        ].index.values

    else:
        raise ValueError("Don't know how to handle the ROI access method")

    return idxs

def ref_selection_top(
        ref_id,
        roi
):

    from wepy.util.json_top import json_top_subset

    # get the indices to use for the subset
    sel_idxs = ref_selection_idxs(
        ref_id,
        roi,
    )


    # get the full JSON topology
    all_json_top = get_ref_structure_json_top(ref_id)

    sel_top = json_top_subset(
        all_json_top,
        sel_idxs,
    )

    return sel_top

def load_ref_structure_traj(
        ref_id,
        roi=Ellipsis,
):

    from wepy.util.mdtraj import json_to_mdtraj_topology

    # load the full raw structure
    traj = load_ref_structure_pdb(ref_id)

    # if it is Ellipsis just return the whole thing
    if roi is Ellipsis:
        return traj

    # get the subset of the topology
    sel_top = ref_selection_top(ref_id, roi)

    # get the selection idxs
    sel_idxs = ref_selection_idxs(ref_id, roi)

    # make an mtraj top out of this
    mdj_sel_top = json_to_mdtraj_topology(sel_top)


    # slice the coordinates
    sel_coords = traj.xyz[:, sel_idxs, :]

    # NOTE: we are ignoring the box here since this only becomes
    # important after we solvate which doesn't effect the reference
    # structures
    sel_traj = mdj.Trajectory(
        sel_coords,
        mdj_sel_top,
    )

    return sel_traj

def load_ref_structure_traj_fields(
        ref_id,
        roi=Ellipsis,
):

    traj = load_ref_structure_traj(
        ref_id,
        roi=roi,
    )

    return {
        'positions' : traj.xyz,
    }

# the complex structures are specific to each ligand. They names need
# string formatting to get a real path.

COMPLEX_STRUCTURES = {
    # this is the protein with the ligand and crystallographic waters
    # from the 'fixed_crystallographic_waters' model. It includs both domains.
    'docked_complex' : {
        'pdb_path' : "data/docking/assemblies/{lig_id}_4od0_docked.pdb"
    },

    # this is the result of solvating and truncating the
    # 'docked_complex' using the CHARMM-GUI servers. This includes the
    # 'active_domain', ligand, solvent ions, water solvent, as well as
    # being a fully hydrogenated model, unlike the other ones.
    #
    # This is the first structure that uses SEGIDs. It uses different
    # SEGIDs for the crystallographic waters over the solvation ones
    # which is nice.
    'chgui_solvated_complex' : {
        'pdb_path' : "data/charmm-gui_solvated_assemblies/{lig_id}/step2_solvator.pdb"
    },

    # This is a processed version of the 'chgui_solvated_complex'. It
    # doesn't change anything from CHARMM-GUI except some minor things
    # it got wrong for making the forcefields. It also adds complete
    # connectivity records.
    #
    # This is the model that is actually used for simulations and will
    # be minimized and equilibrated. The topology here is made into
    # the standard JSON topology that is used throughout. It also has
    # an detailed CSV file for the atoms that details their
    # relationship to the atom types in the force field etc.
    #
    # The primary changes are simply to normalize the structure so
    # that the indices etc. are more coherent for the actual
    # simulation model, rather than the entire protein itself. It also
    # renames the SOD ions as their own SEGID, rather than SOLV which
    # is reserved for the waters. It also reorders it so all the SOLV
    # waters are at the bottom.
    'fixed_complex' : {
        'pdb_path' : "data/ff_top/{lig_id}/sEH_lig-{lig_id}_system.FIXED.pdb"
    },

    # this is the 'fixed_complex' structure that has been minimized &
    # equilibrated.
    #
    # This also has positions in the state pickle which is used for
    # generating the initial simulation snapshots
    'equilibrated_complex' : {
        'pdb_path' : "data/equilibration/{lig_id}/sEH_lig-{lig_id}_equilibrated.pdb"
    },

}

COMPLEX_ROIS = {
    'docked_complex' : {

        'ligand' : (
            ('resname', 'UNL'),
        ),
        'protein' : (
            ('chain', 0),
        ),
        'active_domain' : (
            ('resid', (231, 547,)), # this is inclusive of the second index
        ),
        'crystallographic_waters' : (
            ('resname', 'HOH'),
        ),
    },

    'chgui_solvated_complex' : {

        'ligand' : (
            ('segid', 'HETA'),
        ),
        'protein' : (
            ('segid', 'PROA'),
        ),

        'crystallographic_waters' : (
            ('segid', 'WATA'),
        ),

        'sodium' : (
            ('segid', 'SOD'),
        ),
        'solvent_waters' : (
            ('segid', 'SOLV'),
        ),

    },

    'fixed_complex' : {

        'ligand' : (
            ('segid', 'HETA'),
        ),
        'protein' : (
            ('segid', 'PROA'),
        ),

        'crystallographic_waters' : (
            ('segid', 'WATA'),
        ),

        'sodium' : (
            ('segid', 'SOD'),
        ),
        'solvent_waters' : (
            ('segid', 'SOLV'),
        ),
    },

    'equilibrated_complex' : {

        'ligand' : (
            ('segid', 'HETA'),
        ),
        'protein' : (
            ('segid', 'PROA'),
        ),

        'crystallographic_waters' : (
            ('segid', 'WATA'),
        ),

        'sodium' : (
            ('segid', 'SOD'),
        ),
        'solvent_waters' : (
            ('segid', 'SOLV'),
        ),

    },

}

## raw loaders for reference structures
def load_complex_structure_pdb(
        complex_id,
        lig_id,
):
    """Load the raw PDB for the given complex_id"""

    pdb_path = COMPLEX_STRUCTURES[complex_id]['pdb_path'].format(
        lig_id=lig_id
    )
    traj = mdj.load_pdb(pdb_path)

    return traj


def get_complex_structure_mdj_top(
        complex_id,
        lig_id,
):
    return load_complex_structure_pdb(complex_id, lig_id).top


def get_complex_structure_top_df(
        complex_id,
        lig_id,
):
    return get_complex_structure_mdj_top(complex_id, lig_id).to_dataframe()[0]

def get_complex_structure_json_top(
        complex_id,
        lig_id,
):

    from wepy.util.mdtraj import mdtraj_to_json_topology

    return mdtraj_to_json_topology(
        load_complex_structure_pdb(complex_id, lig_id).top
    )

## methods to get subsets of a reference structure

def complex_selection_idxs(
        complex_id,
        lig_id,
        roi,
):
    """Return the idxs of the ROI in the given reference structure"""


    top_df = get_complex_structure_top_df(complex_id, lig_id)


    # get the key and value for selecting the ROI, use the first one
    # as a default

    key, val = COMPLEX_ROIS[complex_id][roi][0]

    # convert to the mdtraj key
    mdj_key = MDJ_COL_NAMES_TRANSLATION[key]

    # these are both single selections
    if key in ('resname', 'chain', 'segid',):

        idxs = top_df[top_df[mdj_key] == val].index.values

    # these are ranges
    elif key in ('resid', 'pdb_serial', 'idx',):

        start, stop = val

        idxs = top_df[
            (top_df[mdj_key] >= start) &
             (top_df[mdj_key] <= stop)
        ].index.values

    else:
        raise ValueError("Don't know how to handle the ROI access method")

    return idxs

def complex_selection_top(
        complex_id,
        lig_id,
        roi
):

    from wepy.util.json_top import json_top_subset

    # get the indices to use for the subset
    sel_idxs = complex_selection_idxs(
        complex_id,
        lig_id,
        roi,
    )


    # get the full JSON topology
    all_json_top = get_complex_structure_json_top(complex_id, lig_id)

    sel_top = json_top_subset(
        all_json_top,
        sel_idxs,
    )

    return sel_top

def load_complex_structure_traj(
        complex_id,
        lig_id,
        roi=Ellipsis,
):

    from wepy.util.mdtraj import json_to_mdtraj_topology

    # load the full raw structure
    traj = load_complex_structure_pdb(complex_id, lig_id)

    # if it is Ellipsis just return the whole thing
    if roi is Ellipsis:
        return traj

    # get the subset of the topology
    sel_top = complex_selection_top(complex_id, lig_id, roi)

    # get the selection idxs
    sel_idxs = complex_selection_idxs(complex_id, lig_id, roi)

    # make an mtraj top out of this
    mdj_sel_top = json_to_mdtraj_topology(sel_top)


    # slice the coordinates
    sel_coords = traj.xyz[:, sel_idxs, :]

    # NOTE: we are ignoring the box here since this only becomes
    # important after we solvate which doesn't effect the complexerence
    # structures
    sel_traj = mdj.Trajectory(
        sel_coords,
        mdj_sel_top,
        unitcell_lengths=traj.unitcell_lengths,
        unitcell_angles=traj.unitcell_angles,
    )

    return sel_traj

def load_complex_structure_traj_fields(
        complex_id,
        lig_id,
        roi=Ellipsis,
):

    traj = load_complex_structure_traj(
        complex_id,
        lig_id,
        roi=roi,
    )

    return {
        'positions' : traj.xyz,
    }

# just for reference all the keys we want to fill in
SELECTION_KEYS = ('all_atoms/ligand',
                  'all_atoms/ligand/homology',
                  'all_atoms/protein',
                  'all_atoms/binding_site',
                  'all_atoms/main_rep',
                  'all_atoms/main_rep/ligand',
                  'all_atoms/main_rep/ligand/homology',
                  'all_atoms/main_rep/protein',
                  'all_atoms/main_rep/binding_site',
                  'all_atoms/correct_rep',
                  'all_atoms/correct_rep/ligand',
                  'all_atoms/correct_rep/ligand/homology',

                  'all_atoms/image',

                  'main_rep/ligand',
                  'main_rep/ligand/homology',
                  'main_rep/protein',
                  'main_rep/binding_site',

                  'ligand/homology',

                  'image/ligand',
                  'image/binding_site',
                  'image/protein',

                  'correct_rep/ligand',
                  'correct_rep/ligand/homology',
                  'correct_rep/protein',
                  'correct_rep/binding_site',
                  'correct_rep/missing',
                  'correct_rep/main_rep',

                  # and for the legacy stuff
                  'real_rep',
                  'real_rep/ligand',
                  'real_rep/ligand/homology',
                  'real_rep/protein',
                  'real_rep/binding_site',

                  )

TOP_KEYS = ('all_atoms', 'main_rep', 'image', 'correct_rep', 'ligand')


#@jlmem.cache
def lig_selection_idxs_tops(gexp):
    """ compute both at once for convenience """

    from wepy.util.json_top import json_top_subset
    from wepy.util.mdtraj import (
        mdtraj_to_json_topology,
        json_to_mdtraj_topology,
        traj_fields_to_mdtraj,
    )

    from seh_prep.modules import (
        ligand_idxs,
        protein_idxs,
        binding_site_idxs,
        old_ligand_idxs,
        old_protein_idxs,
        old_binding_site_idxs,
    )
    from seh_prep.parameters import BINDING_SITE_CUTOFF

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    wepy_h5 = get_gexp_wepy_h5(gexp)

    with wepy_h5:
        main_rep_idxs = wepy_h5.main_rep_idxs


    sel_idxs = {}

    # the initial state
    init_state = get_ref_state(gexp)

    # get some things from this
    init_box_vectors = init_state['box_vectors'] # * init_state.box_vectors_unit

    # for the units if its legacy it won't have the units so we get
    # them from another gexp since they are the same
    if gexp in LEGACY_GEXPS:

        # get the init state for another known gexp
        tmp_state = get_ref_state('17')

        # use this for the units
        positions_unit = tmp_state.positions_unit
        box_vectors_unit = tmp_state.box_vectors_unit


    # otherwise just use the init state
    else:
        positions_unit = init_state.positions_unit
        box_vectors_unit = init_state.box_vectors_unit

    # the full topology
    all_atoms_top = get_topology(gexp)

    all_atoms_positions = init_state['positions'] # * positions_unit

    # get the indices of the ligand, protein and binding site
    all_atoms_ligand_idxs = ligand_idxs(all_atoms_top)
    all_atoms_protein_idxs = protein_idxs(all_atoms_top)
    all_atoms_bs_idxs = binding_site_idxs(all_atoms_top,
                                          all_atoms_positions * positions_unit,
                                          init_box_vectors * box_vectors_unit,
                                          BINDING_SITE_CUTOFF)


    # the all atoms subselections
    sel_idxs['all_atoms/ligand'] = all_atoms_ligand_idxs
    sel_idxs['all_atoms/protein'] = all_atoms_protein_idxs
    sel_idxs['all_atoms/binding_site'] = all_atoms_bs_idxs

    # get the idxs of the ligand homology atoms in terms of the ligand itself
    lig_hom_idxs = np.array(dict(LIGAND_HOMOLOGY_INDICES)[lig_id])
    sel_idxs['ligand/homology'] = lig_hom_idxs

    # set the selection for the ligand homology indices from the all atoms
    sel_idxs['all_atoms/ligand/homology'] = sel_idxs['all_atoms/ligand'][sel_idxs['ligand/homology']]


    # do the indices for the faulty main_rep

    # do for the indices within the all_atoms rep

    all_atoms_fields = {'positions' : np.array([all_atoms_positions]),
                       'box_vectors' : np.array([init_box_vectors])}

    # load the topology/reference structure PDB
    ref_traj = traj_fields_to_mdtraj(all_atoms_fields, all_atoms_top)

    # the ligand and protein indices for this ligand system
    all_atoms_main_rep_lig_idxs = old_ligand_idxs(ref_traj.top, "UNL")
    # we get the indices for the old way of getting protein indices,
    # this is missing atoms
    all_atoms_main_rep_prot_idxs = old_protein_idxs(ref_traj.top)

    all_atoms_main_rep_bs_idxs = old_binding_site_idxs(ref_traj.top, "UNL",
                                        all_atoms_positions * positions_unit,
                                        init_box_vectors * box_vectors_unit,
                                        BINDING_SITE_CUTOFF)


    # the indices of the main rep within the all_atoms rep
    sel_idxs['all_atoms/main_rep'] = np.concatenate((all_atoms_main_rep_lig_idxs,
                                                     all_atoms_main_rep_prot_idxs))
    sel_idxs['all_atoms/main_rep/ligand'] = all_atoms_main_rep_lig_idxs
    sel_idxs['all_atoms/main_rep/protein'] = all_atoms_main_rep_prot_idxs
    sel_idxs['all_atoms/main_rep/binding_site'] = all_atoms_main_rep_bs_idxs

    sel_idxs['all_atoms/main_rep/ligand/homology'] = sel_idxs['all_atoms/main_rep/ligand'][sel_idxs['ligand/homology']]



    # and the internal indices of the main rep


    main_rep_top = json_top_subset(all_atoms_top, main_rep_idxs)
    main_rep_positions = all_atoms_positions[main_rep_idxs]


    main_rep_fields = {'positions' : np.array([main_rep_positions]),
                       'box_vectors' : np.array([init_box_vectors])}

    # load the topology/reference structure PDB
    main_rep_traj = traj_fields_to_mdtraj(main_rep_fields, main_rep_top)

    # the ligand and protein indices for this ligand system
    main_rep_lig_idxs = old_ligand_idxs(main_rep_traj.top, "UNL")
    # we get the indices for the old way of getting protein indices,
    # this is missing atoms
    main_rep_prot_idxs = old_protein_idxs(main_rep_traj.top)

    main_rep_bs_idxs = old_binding_site_idxs(
        main_rep_traj.top,
        "UNL",
        main_rep_positions * positions_unit,
        init_box_vectors * box_vectors_unit,
        BINDING_SITE_CUTOFF,
    )


    # the indices of the main rep within the all_atoms rep
    sel_idxs['main_rep/ligand'] = main_rep_lig_idxs
    sel_idxs['main_rep/protein'] = main_rep_prot_idxs
    sel_idxs['main_rep/binding_site'] = main_rep_bs_idxs

    sel_idxs['main_rep/ligand/homology'] = sel_idxs['main_rep/ligand'][sel_idxs['ligand/homology']]

    # make a correct rep topology

    # the selections relative to the correct rep
    correct_rep_idxs = np.concatenate((all_atoms_ligand_idxs,
                                       all_atoms_protein_idxs))
    correct_rep_top = json_top_subset(all_atoms_top, correct_rep_idxs)
    correct_rep_positions = all_atoms_positions[correct_rep_idxs]

    # get the ligand, protein, and binding site indices for this
    correct_ligand_idxs = ligand_idxs(correct_rep_top)
    correct_protein_idxs = protein_idxs(correct_rep_top)
    correct_bs_idxs = binding_site_idxs(correct_rep_top,
                                        correct_rep_positions * positions_unit,
                                        init_box_vectors * box_vectors_unit,
                                        BINDING_SITE_CUTOFF)


    # determine the indices of the missing atoms from the main rep in
    # the all_atoms rep
    missing_atom_idxs = list(set(correct_rep_idxs).difference(main_rep_idxs))
    missing_atom_idxs.sort()

    # and the atoms they both have
    common_atom_idxs = list(set(correct_rep_idxs).intersection(main_rep_idxs))
    common_atom_idxs.sort()

    # then get them in terms of the correct rep
    correct_rep_missing_idxs = np.where(np.in1d(correct_rep_idxs, missing_atom_idxs))[0]
    correct_rep_common_idxs = np.where(np.in1d(correct_rep_idxs, common_atom_idxs))[0]

    # correct rep selections
    sel_idxs['all_atoms/correct_rep'] = correct_rep_idxs
    sel_idxs['all_atoms/correct_rep/missing'] = np.array(missing_atom_idxs)
    sel_idxs['all_atoms/correct_rep/main_rep'] = np.array(common_atom_idxs)

    # TODO
    # sel_idxs['all_atoms/correct_rep/ligand'] =
    # sel_idxs['all_atoms/correct_rep/ligand/homology'] = sel_idxs['all_atoms/correct_rep/ligand'][sel_idxs['ligand/homology']]

    sel_idxs['correct_rep/ligand'] = correct_ligand_idxs
    sel_idxs['correct_rep/protein'] = correct_protein_idxs
    sel_idxs['correct_rep/binding_site'] = correct_bs_idxs

    sel_idxs['correct_rep/missing'] = correct_rep_missing_idxs
    sel_idxs['correct_rep/main_rep'] = correct_rep_common_idxs
    sel_idxs['correct_rep/ligand/homology'] = sel_idxs['correct_rep/ligand'][sel_idxs['ligand/homology']]


    # image

    # combine the ligand and binding site to get the image idxs
    image_idxs = np.concatenate((all_atoms_ligand_idxs, all_atoms_bs_idxs))

    # get the topology for this
    image_top = json_top_subset(all_atoms_top, image_idxs)

    # get the ligand and protein from that
    image_ligand_idxs = ligand_idxs(image_top)
    image_protein_idxs = protein_idxs(image_top)

    sel_idxs['all_atoms/image'] = image_idxs
    sel_idxs['image/ligand'] = image_ligand_idxs
    sel_idxs['image/protein'] = image_protein_idxs
    # since the image is the binding site we also set this here which
    # is expected by other programs
    sel_idxs['image/binding_site'] = image_protein_idxs

    # ligand
    ligand_top = json_top_subset(all_atoms_top, all_atoms_ligand_idxs)

    sel_tops = {}

    sel_tops['all_atoms'] = all_atoms_top
    sel_tops['main_rep'] = main_rep_top
    sel_tops['ligand'] = ligand_top
    sel_tops['image'] = image_top
    sel_tops['correct_rep'] = correct_rep_top

    # add in stuff for legacy if necessary
    if gexp in LEGACY_GEXPS:

        leg_sel_idxs, leg_sel_tops = real_rep_selection_idxs_tops(
            sel_idxs,
            sel_tops,
        )

        sel_idxs.update(leg_sel_idxs)
        sel_tops.update(leg_sel_tops)

    return sel_idxs, sel_tops

def lig_selection_idxs(gexp):

    sels, tops = lig_selection_idxs_tops(gexp)

    return sels

def lig_selection_tops(gexp):

    sels, tops = lig_selection_idxs_tops(gexp)

    return tops

def ligs_selection_tops():
    d = {}
    for gexp in GEXPS:
        d[gexp] = lig_selection_tops(gexp)
    return d

def ligs_selection_idxs():
    d = {}
    for gexp in GEXPS:
        d[gexp] = lig_selection_idxs(gexp)

    return d


# make a dataframe summarizing the counts of different selections for each ligand
def ligs_selection_df():
    from collections import defaultdict

    selection_df = defaultdict(list)
    for gexp, selections in ligs_selection_idxs().items():

        selection_df['gexp'].append(gexp)
        for selection_key, idxs in selections.items():
            selection_df[selection_key].append(len(idxs))

    selection_df = pd.DataFrame(selection_df)

    return selection_df

REF_TOP_DIR = "data/ref_top"

# these are selections that can be made on 'ref_prot'. The name refers
# to the name from my previous paper.
BS_ROIS = {
    'Tyr236' : {
        'resid' : (235, 235,),
        'color' : 'green',
    },

    'Asp105' : {
        'resid' : (105, 105,),
        'color' : 'orange',
    },

    'Tyr153' : {
        'resid' : (152, 152,),
        'color' : 'blue',
    },
    'Met189' : {
        'resid' : (188, 188,),
        'color' : 'yellow',
    },
    'Ala134' : {
        'resid' : (134, 134,),
        'color' : 'pink',
    },

}

# TODO: double check and note here exactly which "rep" this is
# for. Hopefully shouldn't make a difference since it is for the
# ligand which was complete always

# the order of the atom indices here is very important as this is what
# maps them to each other. I.e. the index of the atom index is the key
# for the abstract homological atom we are comparing.
LIGAND_HOMOLOGY_INDICES = (
    ('17', (14, 5, 0, 19, 36, 37,)),
    ('3', (43, 40, 29, 0, 11, 13)),
    ('10', (10, 5, 0, 24, 35, 37)),
    ('18', (14, 11, 0, 19, 36, 37)),
    ('20', (46, 37, 32, 0, 11, 13)),
)

def get_lig_bs_bs_atom_idxs(gexp,
                            rep_key='main_rep'):
    """This will return the indices of the atoms that are considered for
    the feature vector.

    Basically this is a selection of specific points in the binding
    site to make a reduced feature vector size.

    These indices are relative to the BS_ANATOMY_TEMPLATE_NAME: ref_prot

    """

    # NOTE: gexp is ignored, they are all the same


    # always load the centered one
    # ref_traj = mdj.load_pdb(str(Path(REF_TOP_DIR) / f"{rep_key}_center_ref.pdb") )

    ref_traj = lig_centered_ref_selection_trajs(
        gexp,
    )[rep_key]

    ref_top_df = ref_traj.top.to_dataframe()[0]

    # we will use a single backbone C-alpha atom ('CA' type) from
    # specific amino acid residues of interest (ROIs)

    # these are the ROI names for the amino acids
    aa_rois = (
        'Tyr236',
        'Asp105',
        'Tyr153',
        'Met189',
        'Ala134',
    )

    # for each of them get the topology of them and get the index of
    # the CA
    aa_ca_idxs = []
    for aa_roi in aa_rois:

        # get the resid, this is encoded as a range, so discard the
        # end
        resid, _ = BS_ROIS[aa_roi]['resid']

        aa_ca_idx = ref_top_df[
            (ref_top_df['resSeq'] == resid) &
            (ref_top_df['name'] == 'CA')
        ].index.values[0]

        aa_ca_idxs.append(aa_ca_idx)

    # thats it, just return these
    return np.array(aa_ca_idxs)


def get_lig_bs_lig_atom_idxs(gexp,
                             rep_key='main_rep'):

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    # these are the indices chosen on the ligand. They are
    # "homological" because they are matched across the different
    # ligands to be the same
    hom_idxs = list(dict(LIGAND_HOMOLOGY_INDICES)[lig_id])

    # these are in terms of the ligand itself we need the relative to
    # the rep
    lig_sel_idxs = lig_selection_idxs(lig_id)[f'{rep_key}/ligand'][hom_idxs]

    return lig_sel_idxs


def lig_bs_atom_pairs(
        gexp,
        rep_key='main_rep',
):
    """Get the pairs of atoms between a selection of atoms on the ligand
    and the selected binding site atoms.

    This is always done in terms of the 'main_rep' since that is the
    only reliable metric.

    """

    from itertools import product

    # Get the atoms on the binding site that we are pairing with
    bs_sel_idxs = get_lig_bs_bs_atom_idxs(
        gexp,
        rep_key=rep_key
    )

    lig_sel_idxs = get_lig_bs_lig_atom_idxs(
        gexp,
        rep_key=rep_key,
    )

    # make the all to all pairs of them
    pairs = np.array(list(product(lig_sel_idxs, bs_sel_idxs)))

    return pairs


def vmd_query_lig_bs_atom_pairs(
        gexp,
        rep_key='main_rep',
):
    """Print the VMD query for selecting the selected atoms for the ligand
    and binding site pairing.
    """

    lig_sel_idxs = get_lig_bs_lig_atom_idxs(
        gexp,
        rep_key=rep_key
    )

    bs_sel_idxs = get_lig_bs_bs_atom_idxs(
        gexp,
        rep_key=rep_key
    )

    # verify visually so print out the VMD query for the selection of
    # atoms for both the ligand and the protein
    vmd_q_str = "index {}"

    print("Protein selection indices: ",
          vmd_q_str.format(' '.join([str(i) for i in bs_sel_idxs])))

    print("Ligand selection indices: ",
          vmd_q_str.format(' '.join([str(i) for i in lig_sel_idxs])))

def get_span_warps_df(
        gexp,
        span_id,
):

    contig = get_gexp_span_contig(gexp, span_id)

    with contig:
        return contig.warping_records_dataframe()


def get_gexp_warps_df(
        gexp,
):
    """Make a table of all warping events across a gexp.

    Assign an idx to each warp event and add a column for the span
    that it is in.

    """

    span_dfs = []
    for span_id in get_gexp_span_ids(gexp):

        contig = get_gexp_span_contig(gexp, span_id)

        # get the warps df
        with contig:
            span_warps_df = contig.warping_records_dataframe()

        if len(span_warps_df) > 0:
            # add a column for the span_id
            span_warps_df['span_id'] = span_id

            span_dfs.append(span_warps_df)

    # concatenate to the master table
    warps_df = pd.concat(span_dfs)

    # then reindex
    warps_df.index = range(len(warps_df))

    return warps_df

# @jlmem.cache
def get_warp_lineage_traces(lig_id, span_id):

    contig = get_gexp_span_contig(lig_id, span_id)

    with contig:
        warp_trace = contig.warp_contig_trace()
        lineage_traces = contig.lineages(warp_trace)

    return lineage_traces


def get_warp_lineage_traj_fields(gexp, span_id):
    """Generator for the warp lineages and do recentering if necessary.

    Trajectory fields only.

    Shouldn't have all in memory at once so we yield them.
    """

    lineage_traces = get_warp_lineage_traces(gexp, span_id)

    contig = get_gexp_span_contig(gexp, span_id)

    for warp_idx, lineage_trace in enumerate(lineage_traces):

        with contig:

            if gexp == '3':

                traj_fields = contig.wepy_h5.get_trace_fields(
                    lineage_trace,
                    ['positions', 'box_vectors']
                )

            else:
                traj_fields = contig.wepy_h5.get_trace_fields(
                    lineage_trace,
                    ['positions', 'box_vectors', 'alt_reps/missing']
                )


        # if we can convert the positions to the correct rep
        if gexp == '3':
            alt_rep = 'main_rep'

        else:
            alt_rep = 'correct_rep'
            # get the correct rep from this data
            traj_fields = traj_fields_to_correct_rep(traj_fields, gexp)
            gc.collect()

        ## do processing

        # 1. recentering
        traj_fields['positions'] = recenter_superimpose_traj(
            traj_fields,
            gexp,
            alt_rep,
        )[0]


        yield traj_fields

def get_warp_lineage_trajs(gexp, span_id):
    """Generator for the warp lineages and do recentering if necessary.

    mdtraj trajectories

    Shouldn't have all in memory at once so we yield them.
    """

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    if lig_id == '3':
        json_top = lig_selection_tops(lig_id)['main_rep']

    else:
        json_top = lig_selection_tops(lig_id)['correct_rep']

    for traj_fields in get_warp_lineage_traj_fields(gexp, span_id):
        yield traj_fields_to_mdtraj(
           traj_fields,
           json_top
        )

# @jlmem.cache
def get_final_lineage_traces(lig_id, span_id):

    contig = get_gexp_span_contig(lig_id, span_id)

    with contig:
        final_trace = contig.final_contig_trace()
        lineage_traces = contig.lineages(final_trace)

    return lineage_traces


def get_final_lineage_traj_fields(gexp, span_id):
    """Generator for the final walkers lineages and do recentering if necessary.

    Trajectory fields only.

    Shouldn't have all in memory at once so we yield them.
    """

    lineage_traces = get_final_lineage_traces(gexp, span_id)

    contig = get_gexp_span_contig(gexp, span_id)

    for lineage_trace in lineage_traces:


        if gexp == '3':
            alt_rep = 'main_rep'
            with contig:
                traj_fields = contig.wepy_h5.get_trace_fields(
                    lineage_trace,
                    ['positions', 'box_vectors',]
                )


        else:
            alt_rep = 'correct_rep'

            with contig:
                traj_fields = contig.wepy_h5.get_trace_fields(
                    lineage_trace,
                    ['positions', 'box_vectors', 'alt_reps/missing']
                )

            # get the correct rep from this data
            traj_fields = traj_fields_to_correct_rep(traj_fields, gexp)
            gc.collect()

        ## do processing

        # 1. recentering
        traj_fields['positions'] = recenter_superimpose_traj(
            traj_fields,
            gexp,
            alt_rep,
        )[0]


        yield traj_fields

def get_final_lineage_trajs(gexp, span_id):
    """Generator for the final lineages and do recentering if necessary.

    mdtraj trajectories

    Shouldn't have all in memory at once so we yield them.
    """

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    json_top = lig_selection_tops(lig_id)['correct_rep']

    for traj_fields in get_final_lineage_traj_fields(gexp, span_id):
        yield traj_fields_to_mdtraj(
                           traj_fields,
                           json_top
        )

def get_high_progress_traces(
        gexp,
        span_id,
        top_n,
        progress_key):

    contig_trace = get_high_progress_contig_traces(
        gexp,
        span_id,
        top_n,
        progress_key)

    contig = get_gexp_span_contig(gexp, span_id)
    with contig:
        # convert to a run trace for getting data out of HDF5
        run_trace = contig.walker_trace_to_run_trace(
            contig_trace
            )

    return run_trace

def get_high_progress_contig_traces(
        gexp,
        span_id,
        top_n,
        progress_key):

    contig = get_gexp_span_contig(gexp, span_id)

    with contig:

        # get the progress records for this contig as a big array
        prog_vals = np.array(
            [getattr(rec, progress_key)
             for rec in contig.progress_records()])

    # do a partial sort for top N values and slice the indices
    # of those off
    top_prog_val_flatidxs = np.argpartition(prog_vals.flat, -top_n)[-top_n:]

    # function to unflatten the idxs
    row_len = prog_vals.shape[1]
    unflatten_idx = (lambda flat_idx :
                     (flat_idx // row_len, flat_idx % row_len))


    with contig:

        # get the unflattenend contig walker trace indices of them
        # (walker, cycle)
        contig_trace = [(walker_idx, cycle_idx)
                        for cycle_idx, walker_idx in
                        [unflatten_idx(flatidx)
                         for flatidx in top_prog_val_flatidxs
                         ]
                        ]


    return contig_trace

def get_high_progress_traj_fields(
        gexp,
        span_id,
        top_n,
        progress_key
):

    run_trace = get_high_progress_traces(
                gexp,
                span_id,
                top_n,
                progress_key,
            )

    contig = get_gexp_span_contig(gexp, span_id)

    with contig:

        # now use the trace to get the fields and traj for this with
        # recentering
        traj_fields = contig.wepy_h5.get_trace_fields(
            top_prog_val_run_trace,
            ['positions', 'box_vectors', 'alt_reps/missing'],
        )


    # get the correct rep from this data
    traj_fields = traj_fields_to_correct_rep(traj_fields, gexp)
    gc.collect()


    traj_fields['positions'] = recenter_superimpose_traj(
        traj_fields,
        gexp,
        'correct_rep',
    )[0]

    return traj_fields


def get_high_progress_trajs(
        gexp,
        span_id,
        top_n,
        progress_key,
):

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    json_top = lig_selection_tops(lig_id)['correct_rep']


    traj_fields = get_high_progress_traj_fields(
                gexp,
                span_id,
                top_n,
                progress_key,
            )

    yield traj_fields_to_mdtraj(
            traj_fields,
            json_top
        )



### The lineages version


def get_high_progress_lineage_traces(
        gexp,
        span_id,
        top_n,
        progress_key,
):

    # the trace of the final walkers in a lineage
    child_trace = get_high_progress_contig_traces(
        gexp,
        span_id,
        top_n,
        progress_key)

    contig = get_gexp_span_contig(gexp, span_id)

    with contig:
        lineage_traces = contig.lineages(child_trace)

    return lineage_traces

def get_high_progress_lineages_traj_fields(
        gexp,
        span_id,
        top_n,
        progress_key,
):

    lineage_traces = get_high_progress_lineage_traces(
          gexp,
          span_id,
          top_n,
          progress_key,
      )

    contig = get_gexp_span_contig(gexp, span_id)

    for lineage_trace in lineage_traces:

        with contig:
            traj_fields = contig.wepy_h5.get_trace_fields(
                lineage_trace,
                ['positions', 'box_vectors', 'alt_reps/missing']
            )

        # get the correct rep from this data
        traj_fields = traj_fields_to_correct_rep(traj_fields, gexp)
        gc.collect()

        ## do processing

        # 1. recentering
        traj_fields['positions'] = recenter_superimpose_traj(
            traj_fields,
            gexp,
            'correct_rep'
        )[0]


        yield traj_fields

def get_high_progress_lineages_trajs(
        gexp,
        span_id,
        top_n,
        progress_key,
):

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    json_top = lig_selection_tops(lig_id)['correct_rep']

    for traj_fields in get_high_progress_lineages_traj_fields(
            gexp,
            span_id,
            top_n,
            progress_key,
    ):

        yield traj_fields_to_mdtraj(
            traj_fields,
            json_top
        )

SPANS_TABLE_COLS = (
    'gexp',
    'span_id',
    'n_runs',
    'n_warps',
    'n_cycles',
    'total_sampling_time_ps', # picoseconds
    'warps_per_ps', # picoseconds
    'rate_ps', # picoseconds
    'total_warped_weight',
    'warp_weight_mean',
    'warp_weight_median',
    'warp_weight_min',
    'warp_weight_max',
    'warp_wait-time_mean_ps',
    'warp_wait-time_median_ps',
    'warp_wait-time_min_ps',
    'warp_wait-time_max_ps',
)

@jlmem.cache
def get_span_stats_df(gexp):

    from wepy.analysis.rates import contig_warp_rates

    jobs_df = get_gexp_jobs_df(gexp)

    # extract the columns we need for a new dataframe
    spans_df = jobs_df[['span_id', 'segment_idx']]

    print(spans_df)
    # count the number of segments in the jobs
    spans_df = spans_df.groupby(['span_id']).count().reset_index().rename(
        columns={'segment_idx' : 'n_runs'})

    print(spans_df)

    # get the number of cycles and warps in the span by reading the
    # contig
    new_cols = {
        'n_warps' : [],
        'n_cycles' : [],

        'total_sampling_time_us' : [],
        'warps_per_ps' : [],
        'rate_s' : [],
        'rt_ms' : [],
        'total_warped_weight' : [],
        'warp_weight_mean' : [],
        'warp_weight_median' : [],
        'warp_weight_min' : [],
        'warp_weight_max' : [],
        'warp_wait-time_mean_us' : [],
        'warp_wait-time_median_us' : [],
        'warp_wait-time_min_us' : [],
        'warp_wait-time_max_us' : [],
    }

    for row_idx, row in spans_df.iterrows():
        print(row['span_id'])
        contig = get_gexp_span_contig(gexp, row['span_id'])

        with contig:

            n_cycles = contig.num_cycles

            # calculate rates and fluxes
            contig_rate = contig_warp_rates(
                contig,
                CYCLE_TIME
            )[0]

            try:
                total_weight, rate, total_sampling_time = contig_rate[0]

            except KeyError:

                print(f"No rates for {gexp}-{row['span_id']}")
                total_weight = 0.0
                rate = np.inf * (1 / CYCLE_TIME.unit)

                total_sampling_time = contig.num_cycles * \
                    CYCLE_TIME * \
                    contig.num_walkers(contig.num_cycles-1)

            rt = 1 / rate

            warps_df = contig.warping_records_dataframe()



        n_warps = len(warps_df)

        if n_warps > 0:

            # the unweighted rate of warp events
            warps_per_sampling = n_warps / total_sampling_time
            total_warped_weight = warps_df['weight'].sum()

            # descriptive statistics on the warps weights
            mean_warp_weight = warps_df['weight'].mean()
            median_warp_weight = warps_df['weight'].median()
            min_warp_weight = warps_df['weight'].min()
            max_warp_weight = warps_df['weight'].max()


            wait_times_us = warps_df['cycle_idx'] * CYCLE_TIME.in_units_of(tkunit.microsecond)

            # descriptive statistics on the warps waiting times
            mean_warp_wait_time = wait_times_us.mean()
            median_warp_wait_time = wait_times_us.median()
            min_warp_wait_time = wait_times_us.min()
            max_warp_wait_time = wait_times_us.max()

        else:
            # the unweighted rate of warp events
            warps_per_sampling = n_warps / total_sampling_time
            total_warped_weight = 0.0
            mean_warp_weight = 0.0
            median_warp_weight = 0.0
            min_warp_weight = 0.0
            max_warp_weight = 0.0


            # descriptive statistics on the warps waiting times
            mean_warp_wait_time = np.inf * tkunit.microsecond
            median_warp_wait_time = np.inf * tkunit.microsecond
            min_warp_wait_time = np.inf * tkunit.microsecond
            max_warp_wait_time = np.inf * tkunit.microsecond


        new_cols['n_cycles'].append(n_cycles)
        new_cols['rate_s'].append(rate.value_in_unit((1/tkunit.second).unit))
        new_cols['rt_ms'].append(rt.value_in_unit(tkunit.millisecond))
        new_cols['total_sampling_time_us'].append(total_sampling_time.value_in_unit(tkunit.microsecond))
        new_cols['n_warps'].append(n_warps)
        new_cols['total_warped_weight'].append(total_warped_weight)
        new_cols['warps_per_ps'].append(warps_per_sampling.value_in_unit((1/tkunit.picosecond).unit))

        # warp weights statistics
        new_cols['warp_weight_mean'].append(mean_warp_weight)
        new_cols['warp_weight_median'].append(median_warp_weight)
        new_cols['warp_weight_min'].append(min_warp_weight)
        new_cols['warp_weight_max'].append(max_warp_weight)

        # warp wait time statistics
        new_cols['warp_wait-time_mean_us'].append(mean_warp_wait_time)
        new_cols['warp_wait-time_median_us'].append(median_warp_wait_time)
        new_cols['warp_wait-time_min_us'].append(min_warp_wait_time)
        new_cols['warp_wait-time_max_us'].append(max_warp_wait_time)

    # add the columns
    for colname, col in new_cols.items():
        spans_df[colname] = col

    return spans_df

def span_stats_table_str(gexp):

    spans_df = get_span_stats_df(gexp)

    table_str = tabulate(spans_df,
                         headers=spans_df.columns,
                         tablefmt='orgtbl'
    )

    return table_str

def get_master_span_stats_df():

    dfs = []
    for gexp in GEXPS:
        dfs.append(get_span_stats_df(gexp))

    return pd.concat(dfs)

def get_master_span_stats_table_str():

    spans_df = get_master_span_stats_df(gexp)

    table_str = tabulate(spans_df,
                         headers=spans_df.columns,
                         tablefmt='orgtbl'
    )

    return table_str

def plot_agg_prob(
        run_target_weights_rates,
        cycle_time,
        num_walkers,
        title=None,
        target_idx=0,
        time_unit=tkunit.picosecond,
        ylim=None,
        logscale=False,
):

    import numpy as np
    import scipy.stats
    import scipy.stats.mstats

    RUN_LABEL_TEMPLATE = "Run: {}"

    # collate data
    runs_weights = {}
    runs_times = {}
    for run_idx, target_weights_rates in run_target_weights_rates.items():
        run_cum_weights = []
        run_times = []

        for cycle_idx, timepoint in enumerate(target_weights_rates):

            # if there were no records (and thus no weights and rates
            # i.e. an empty dictionary) we set nans here so we can
            # easily discard them later but still keep track of time
            if not timepoint:
                weight, _, total_sampling_time = 0.0, 0.0 * 1/time_unit, cycle_time * cycle_idx * num_walkers
            else:
                weight, _, total_sampling_time = timepoint[target_idx]

            run_cum_weights.append(weight)
            run_times.append(total_sampling_time)

        runs_weights[run_idx] = run_cum_weights
        runs_times[run_idx] = run_times

    # compute the average unbound probability and MFPT across all runs

    # aggregate all of the runs weights into a single array and fill
    # with nans where missing, these will be ignored later

    num_runs = len(runs_times)

    # get the longest runs number of time points
    max_num_timepoints = max([len(run_times) for run_times in runs_times.values()])
    min_num_timepoints = min([len(run_times) for run_times in runs_times.values()])


    # get an array of the actual times for the longest one
    longest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                             reverse=True)[0][1]
    longest_run_times = runs_times[longest_run_idx]

    # the shortest one
    shortest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                              reverse=False)[0][1]
    shortest_run_times = runs_times[shortest_run_idx]



    # make an array for all of them where the longest we have is the
    # shortest one
    all_weights = np.empty((num_runs, min_num_timepoints))
    all_weights.fill(np.nan)
    # add the weights vectors to the array
    for i, run_weights in enumerate(runs_weights.values()):
        all_weights[i, 0:] = run_weights[0:min_num_timepoints]

    # compute the arithmetic mean, ignoring nans
    mean_weights = np.mean(all_weights, axis=0)


    # compute standard error of means
    sem_weights = scipy.stats.sem(all_weights, axis=0, nan_policy='raise')

    # get the bounds of the SEM envelope around the mean so we can
    # transform them
    upper_sem_bounds = mean_weights + sem_weights
    lower_sem_bounds = mean_weights - sem_weights

    # TODO compute harmonic mean as well of the rates

    # TODO: Alex wants me to calculate this just as a transformation
    # of the computed mean and std errors of the probability. Not sure
    # of the validity of this.

    shortest_run_time_values = np.array([time.value_in_unit(time_unit)
                                         for time in shortest_run_times])


    ### Plot
    fig, axes = plt.subplots(1, 1, constrained_layout=True)

    prob_ax = axes

    for run_idx in run_target_weights_rates.keys():

        # do the weights plot

        # get the times and weights converting if necessary
        run_weights = runs_weights[run_idx]
        run_times = [time.value_in_unit(time_unit) for time in runs_times[run_idx]]

        # slice them to the smallest run length
        run_weights = run_weights[0:min_num_timepoints]
        run_times = run_times[0:min_num_timepoints]

        label = RUN_LABEL_TEMPLATE.format(run_idx)

        replicate_color = dict(dict(REPLICATE_COLORS)[run_idx])['base']

        prob_ax.plot(run_times, run_weights, label=label, linewidth='3',
                     color=replicate_color)

        # then do the residence time plots


    # choose colors
    mean_color = dict(dict(REPLICATE_COLORS)['agg'])['base']
    sem_color = dict(dict(REPLICATE_COLORS)['agg'])['light']

    prob_ax.plot(shortest_run_time_values,
                 mean_weights,
                 label='Mean',
                 linewidth='3',
                 color=mean_color)

    prob_ax.fill_between(shortest_run_time_values,
                         mean_weights-sem_weights,
                         mean_weights+sem_weights,
                         label='Std. Error',
                         color=sem_color)

    # get the units from the value of the rates
    # rate_unit = runs_rates[list(runs_rates.keys())[0]][0].get_symbol()

    # # convert to the given unit as a rate
    # if time_unit is not None:
    #     rate_unit = 1 / (1 / rate_unit).in_units_of(time_unit)

    ## labels and fonts
    label_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '12'
    }

    ticks_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '10'
    }

    title_font = {
        'family' : 'sans-serif',
        'style' : 'normal',
        'size' : '12'
    }

    if title is None:
        prob_ax.set_title("Aggregate Probability",
                          fontdict=title_font)
    else:
        prob_ax.set_title(title,
                          fontdict=title_font)

    prob_ax.set_yscale('log')

    prob_ax.set_ylabel('Log Agg. Probability',
                       fontdict=label_font)

    # examples
    # prob_ax.set_ylim([1e-14, 1e-8])
    # prob_ax.set_ylim([1e-12, 10])
    if ylim is not None:
        prob_ax.set_ylim(ylim)


    if logscale:
        prob_ax.set_xscale('log')

        # TODO: set this intelligently
        prob_ax.set_xlim([50e-4, 1])

    prob_ax.set_xlabel('Log Simulation time (${}$)'.format(time_unit.get_symbol()),
                       fontdict=label_font)



    # set the tick fonts
    # for label in (prob_ax.get_xticklabels() + prob_ax.get_yticklabels()):
    #     label.set_fontname("serif")
    #     label.set_fontsize(10)

    prob_ax.legend(loc='best')

    return mean_weights, sem_weights, shortest_run_times, (fig, prob_ax)

def gexp_plot_agg_prob(gexp):

    from wepy.analysis.rates import contig_warp_rates

    YLIMS = {
        '3' : [1e-12, 1e-4],
        '10' : None,
        '17' : None,
        '18' : None,
        '20' : None,
        'TPPU-legaxy' : None,
    }

    span_rates = {}
    for span_id in get_gexp_span_ids(gexp):
        contig = get_gexp_span_contig(gexp, span_id)

        num_walkers = get_gexp_span_num_walkers(gexp, span_id)

        with contig:

            # get the values as a series
            contig_rates = contig_warp_rates(contig, CYCLE_TIME, time_points=Ellipsis)

            span_rates[span_id] = contig_rates

    # make the plot for each span
    result = plot_agg_prob(span_rates,
                           CYCLE_TIME,
                           num_walkers,
                           target_idx=0,
                           title=gexp,
                           time_unit=tkunit.microsecond,
                           ylim=YLIMS[gexp],
                           logscale=False)

    return result

def gexp_show_plot_agg_prob(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, total_sampling_times, plots = \
                                                            gexp_plot_agg_prob(gexp)

    print("Gexp: {}".format(gexp))

    plt.show()

def save_gexp_plot_agg_prob(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, total_sampling_times, plots = \
                                                            gexp_plot_agg_prob(gexp)

    print("Gexp: {}".format(gexp))

    fig, axes = plots

    save_gexp_fig(gexp, "agg_weights", fig)

def plot_rts(gexp,
             run_target_weights_rates,
             cycle_time,
             num_walkers,
             experimental_halftime,
             target_idx=0,
             mfpt_unit=tkunit.second,
             time_unit=tkunit.picosecond,
             ylim=None,
             halftime=False,
             logscale=False):
    """Plot the residence times

    Parameters
    ----------

    halftime : report the y axis as the halftime or the residence time

    experimental_halftime : The halftime of the ligand ln(2)/rate

    mfpt_unit : Unit of the y-axis, will convert values to this

    time_unit : Unit of the x-axis

    logscale : make the time axis logarithmic if True

    """

    import numpy as np
    import scipy.stats
    import scipy.stats.mstats

    RUN_LABEL_TEMPLATE = "Run: {}"


    # convert the experimental halftime to what is necessary
    if halftime:
        # already in this unit
        experimental_mfpt = experimental_halftime
    else:
        # conver to 1/rate residence time
        experimental_mfpt = convert_halftime_to_rt(experimental_halftime)

    # collate data
    runs_weights = {}
    runs_rates = {}
    runs_rts = {}
    runs_times = {}
    for run_idx, target_weights_rates in run_target_weights_rates.items():
        run_cum_weights = []
        run_rates = []
        run_times = []

        for cycle_idx, timepoint in enumerate(target_weights_rates):

            # if there were no records (and thus no weights and rates
            # i.e. an empty dictionary) we set nans here so we can
            # easily discard them later but still keep track of time
            if not timepoint:
                #weight, rate, time = np.nan, np.nan, cycle_time * cycle_idx
                weight, rate, total_sampling_time = (0.0,
                                                     0.0 * 1/time_unit,
                                                     cycle_time * cycle_idx * num_walkers)
            else:
                weight, rate, total_sampling_time = timepoint[target_idx]

            run_cum_weights.append(weight)
            run_rates.append(rate)
            run_times.append(total_sampling_time)

        # calculate the residence times
        run_rts = []
        for rate in run_rates:

            # choose whether to plot the MFPT as halftime (ln(2)/rate)
            # or residence time (1/rate)
            if halftime:

                run_rt = convert_rate_to_halftime(rate, time_base=mfpt_unit)

            else:
                run_rt = convert_rate_to_rt(rate, time_base=mfpt_unit)

            # TODO remove. SHould be covered by our conversion functions

            # we can't divide by zero with simtk units so we get the
            # value and invert it with a numpy float64 then
            # redimensionalize it
            # run_rt = (np.float64(1.0) / rate.value_in_unit(rate.unit)) * (1/rate.unit)

            # if the rate is nan (which was 0.0 aggregated
            # probability) we catch this and explicitly set the
            # residence time to infinity
            # if np.isnan(rate):
            #     run_rt = np.inf * time_unit
            # # otherwise compute the mfpt/residence time
            # else:
            #     run_rt = 1/rate

            run_rts.append(run_rt)

        runs_weights[run_idx] = run_cum_weights
        runs_rates[run_idx] = run_rates
        runs_rts[run_idx] = run_rts
        runs_times[run_idx] = run_times


    # compute the average unbound probability and MFPT across all runs

    # aggregate all of the runs weights into a single array and fill
    # with nans where missing, these will be ignored later

    num_runs = len(runs_times)

    # get the longest runs number of time points
    max_num_timepoints = max([len(run_times) for run_times in runs_times.values()])
    min_num_timepoints = min([len(run_times) for run_times in runs_times.values()])


    # get an array of the actual times for the longest one
    longest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                             reverse=True)[0][1]
    longest_run_times = runs_times[longest_run_idx]

    # the shortest one
    shortest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                              reverse=False)[0][1]
    shortest_run_times = runs_times[shortest_run_idx]



    # make an array for all of them where the longest we have is the
    # shortest one
    all_weights = np.empty((num_runs, min_num_timepoints))
    all_weights.fill(np.nan)
    # add the weights vectors to the array
    for i, run_weights in enumerate(runs_weights.values()):
        all_weights[i, 0:] = run_weights[0:min_num_timepoints]

    # compute the arithmetic mean, ignoring nans
    mean_weights = np.mean(all_weights, axis=0)

    # compute standard error of means
    sem_weights = scipy.stats.sem(all_weights, axis=0, nan_policy='raise')

    # get the bounds of the SEM envelope around the mean so we can
    # transform them
    upper_sem_bounds = mean_weights + sem_weights
    lower_sem_bounds = mean_weights - sem_weights

    # TODO compute harmonic mean as well of the rates

    # TODO: Alex wants me to calculate this just as a transformation
    # of the computed mean and std errors of the probability. Not sure
    # of the validity of this.

    shortest_run_time_values = np.array([time.value_in_unit(time_unit)
                                         for time in shortest_run_times])


    # then compute MFPTs

    # mask the values for which no data was observed (nans)

    # calculate the rates for all the runs, masking the zero values
    # which need special tolerances due to the small values. These
    # small values are dependent upon the minimum probabilities
    # allowed by the simulation so this must be a parameter.
    #masked_mean_weights = np.ma.masked_invalid(mean_weights)

    # actually calculate the average mfpts from the mean weights over
    # time
    mfpts = shortest_run_time_values / mean_weights

    # compute the error in the rates using the relationship between
    # the MFPT and the weight: MFPT_error = weight_error * t * weight^2
    mfpt_sems = sem_weights * shortest_run_time_values * (mean_weights**-2)

    upper_sem_mfpts = mfpts + mfpt_sems
    lower_sem_mfpts = mfpts - mfpt_sems


    # clip the things that are close to zero or less than zero since
    # this is just a floating point subtraction error
    close_to_zero_idx = np.argwhere(np.logical_or(np.isclose(lower_sem_mfpts, 0.),
                                                  lower_sem_mfpts < 0.))
    lower_sem_mfpts[close_to_zero_idx] = 0.



    # then compute the bounds for the MFPTs
    # upper_sem_mfpts = shortest_run_time_values / upper_sem_bounds
    # lower_sem_mfpts = shortest_run_time_values / lower_sem_bounds

    # because the lower standard error bound is lower than the average
    # we want to make sure the call to isclose doesn't count them as 0
    # so we make a new min_prob for the masking that takes into
    # account the lowest value. So we take the minimum value that is
    # not identically 0
    # lower_min_prob = lower_sem_bounds[lower_sem_bounds > 0.].min() * 1e-2

    # actually we are going to consder anythin 2 orders of magnitude
    # lower than the minimum probability to be zero so:
    # lower_sem_bound_min_value = min_probability * 1e-2

    # zero_masked_upper_sem_weights = np.ma.masked_values(upper_sem_bounds, 0.,
    #                                                     atol=min_probability)

    # zero_masked_lower_sem_weights = np.ma.masked_values(lower_sem_bounds, 0.,
    #                                                     atol=lower_min_prob)

    # # we calculate the standard error of the mean weights as well
    # zero_masked_upper_sem_mfpts = (shortest_run_time_values /
    #                                zero_masked_upper_sem_weights)
    # zero_masked_lower_sem_mfpts = (shortest_run_time_values /
    #                                zero_masked_lower_sem_weights)

    # make the plots, one for aggregate probability the other for
    # RT
    fig, axes = plt.subplots(1, 1, constrained_layout=True)

    rt_ax = axes

    for run_idx in run_target_weights_rates.keys():

        # do the weights plot

        # get the times and weights converting if necessary
        run_weights = runs_weights[run_idx]
        run_times = [time.value_in_unit(time_unit) for time in runs_times[run_idx]]

        # slice them to the smallest run length
        run_weights = run_weights[0:min_num_timepoints]
        run_times = run_times[0:min_num_timepoints]

        label = RUN_LABEL_TEMPLATE.format(run_idx)

        # then do the residence time plots

        run_rts = runs_rts[run_idx][0:min_num_timepoints]

        # convert the residence times to the time unit specified
        run_rts = [rt.value_in_unit(mfpt_unit) for rt in run_rts]

        replicate_color = dict(dict(REPLICATE_COLORS)[run_idx])['base']

        rt_ax.plot(run_times,
                   run_rts,
                   label=label,
                   linewidth='3',
                   color=replicate_color)



    shortest_run_times_values = [time.value_in_unit(time_unit) for time in shortest_run_times]
    # prob_ax.plot(shortest_run_times_values,
    #              mean_weights,
    #              label='Mean', linewidth='3', color='black')
    # prob_ax.fill_between(shortest_run_times_values,
    #                      mean_weights-sem_weights,
    #                      mean_weights+sem_weights,
    #                      label='Std. Error')

    # get the units from the value of the rates
    # rate_unit = runs_rates[list(runs_rates.keys())[0]][0].get_symbol()

    # # convert to the given unit as a rate
    # if time_unit is not None:
    #     rate_unit = 1 / (1 / rate_unit).in_units_of(time_unit)




    # plot the average rate

    # get the times and convert to the proper units
    times = [time.value_in_unit(time_unit) for time in shortest_run_times]

    # first because they were dedimensionalized we redimensionalize
    # them to the specified time_unit (x-axis) then convert them to
    # the MFPT unit, then dedimensionalize them to raw values for
    # plotting
    mfpt_values = np.array([(mfpt * time_unit).value_in_unit(mfpt_unit)
                            for mfpt in mfpts])

    upper_sem_mfpt_values = np.array([(mfpt * time_unit).value_in_unit(mfpt_unit)
                                      for mfpt in upper_sem_mfpts])

    lower_sem_mfpt_values = np.array([(mfpt * time_unit).value_in_unit(mfpt_unit)
                                      for mfpt in lower_sem_mfpts])

    # clip the things that are close to zero or less than zero since
    # this is just a floating point error from multiplying...
    close_to_zero_idx = np.argwhere(np.logical_or(np.isclose(lower_sem_mfpt_values, 0.),
                                                  lower_sem_mfpts < 0.))
    lower_sem_mfpt_values[close_to_zero_idx] = 0.

    # choose colors
    mean_color = dict(dict(REPLICATE_COLORS)['agg'])['base']
    sem_color = dict(dict(REPLICATE_COLORS)['agg'])['light']

    rt_ax.plot(times,
               mfpt_values,
               label='Mean',
               linewidth='3',
               color=mean_color)

    rt_ax.fill_between(times,
                       upper_sem_mfpt_values,
                       lower_sem_mfpt_values,
                       label='Std. Error of Mean',
                       color=sem_color,
    )

    if logscale:
        rt_ax.set_yscale('log')


    ## labels and fonts
    label_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '12'
    }

    ticks_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '10'
    }

    title_font = {
        'family' : 'sans-serif',
        'style' : 'normal',
        'size' : '12'
    }

    # prob_ax.set_title("Aggregate Probability",
    #                   fontdict=title_font)
    # prob_ax.set_yscale('log')
    # prob_ax.set_ylabel('Aggregated unbound probability',
    #                    fontdict=label_font)

    # # TODO set this intelligently
    # prob_ax.set_ylim([1e-12, 10])

    # prob_ax.set_xlabel('Simulation time (${}$)'.format(time_unit.get_symbol()),
    #                    fontdict=label_font)

    # # set the tick fonts
    # # for label in (prob_ax.get_xticklabels() + prob_ax.get_yticklabels()):
    # #     label.set_fontname("serif")
    # #     label.set_fontsize(10)

    # prob_ax.legend(loc='best')

    rt_ax.set_title(f"Residence Time: {gexp}",
                    fontdict=title_font)

    # plt.xticks([400, 800, 1200, 1600], ['400', '800', '1200', '1600'])
    # plt.yticks([0.1, 10, 1000, 100000, 10000000, 1000000000])
    rt_ax.set_xlabel('Simulation time (${}$)'.format(time_unit.get_symbol()),
                     fontdict=label_font)

    if halftime:
        label_eq = "Half-life ($ln(2) / k_{off}$)"

    else:
        label_eq = "residence time ($1 / k_{off}$)"

    rt_ax.set_ylabel('Predicted {} ({})'.format(
        label_eq,
        mfpt_unit.get_symbol(),
    ),
                     fontdict=label_font)


    # suggestion
    # rt_ax.set_ylim([1e-14, 1e-8])
    if ylim is not None:
        rt_ax.set_ylim(ylim)


    # plot the experimental value
    experimental_points = [experimental_mfpt.value_in_unit(mfpt_unit) for _ in run_times]
    rt_ax.plot(times,
               experimental_points,
               color='red',
               label='Experimental')

    rt_ax.legend()

    # return the actual computed values for the average aggregated
    # weights and residence times
    mfpts = np.array([(mfpt * time_unit).value_in_unit(mfpt_unit)
                      for mfpt in mfpts]) * mfpt_unit

    mfpt_sems = np.array([(mfpt * time_unit).value_in_unit(mfpt_unit)
                          for mfpt in mfpt_sems]) * mfpt_unit



    return mean_weights, sem_weights, mfpts, mfpt_sems, shortest_run_times, (fig, rt_ax)

def plot_rates(gexp,
               run_target_weights_rates,
               cycle_time,
               num_walkers,
               experimental_rate,
               target_idx=0,
               rate_unit=(1/tkunit.second).unit,
               time_unit=tkunit.picosecond,
               ylim=None,
               show_runs=True,
               logscale=False):

    import numpy as np
    import scipy.stats
    import scipy.stats.mstats

    RUN_LABEL_TEMPLATE = "Run: {}"

    # collate data
    runs_weights = {}
    runs_rates = {}
    runs_rts = {}
    runs_times = {}
    for run_idx, target_weights_rates in run_target_weights_rates.items():
        run_cum_weights = []
        run_rates = []
        run_times = []

        for cycle_idx, timepoint in enumerate(target_weights_rates):

            # if there were no records (and thus no weights and rates
            # i.e. an empty dictionary) we set nans here so we can
            # easily discard them later but still keep track of time
            if not timepoint:
                #weight, rate, time = np.nan, np.nan, cycle_time * cycle_idx
                weight, rate, total_sampling_time = (0.0,
                                                     0.0 * (1/time_unit).unit,
                                                     cycle_time * cycle_idx * num_walkers)
            else:
                weight, rate, total_sampling_time = timepoint[target_idx]

            run_cum_weights.append(weight)
            run_rates.append(rate)
            run_times.append(total_sampling_time)

        # calculate the residence times
        run_rts = []
        for rate in run_rates:

            # we can't divide by zero with simtk units so we get the
            # value and invert it with a numpy float64 then
            # redimensionalize it
            run_rt = (np.float64(1.0) / rate.value_in_unit(rate.unit)) * (1/rate.unit)

            # if the rate is nan (which was 0.0 aggregated
            # probability) we catch this and explicitly set the
            # residence time to infinity
            # if np.isnan(rate):
            #     run_rt = np.inf * time_unit
            # # otherwise compute the mfpt/residence time
            # else:
            #     run_rt = 1/rate

            run_rts.append(run_rt)

        runs_weights[run_idx] = run_cum_weights
        runs_rates[run_idx] = run_rates
        runs_rts[run_idx] = run_rts
        runs_times[run_idx] = run_times


    # compute the average unbound probability and MFPT across all runs

    # aggregate all of the runs weights into a single array and fill
    # with nans where missing, these will be ignored later

    num_runs = len(runs_times)

    # get the longest runs number of time points
    max_num_timepoints = max([len(run_times) for run_times in runs_times.values()])
    min_num_timepoints = min([len(run_times) for run_times in runs_times.values()])


    # get an array of the actual times for the longest one
    longest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                             reverse=True)[0][1]
    longest_run_times = runs_times[longest_run_idx]

    # the shortest one
    shortest_run_idx = sorted([(len(run_times), run_idx)
                              for run_idx, run_times in runs_times.items()],
                              reverse=False)[0][1]
    shortest_run_times = runs_times[shortest_run_idx]



    # make an array for all of them where the longest we have is the
    # shortest one
    all_weights = np.empty((num_runs, min_num_timepoints))
    all_weights.fill(np.nan)
    # add the weights vectors to the array
    for i, run_weights in enumerate(runs_weights.values()):
        all_weights[i, 0:] = run_weights[0:min_num_timepoints]

    # compute the arithmetic mean, ignoring nans
    mean_weights = np.mean(all_weights, axis=0)

    # compute standard error of means
    sem_weights = scipy.stats.sem(all_weights, axis=0, nan_policy='raise')

    # get the bounds of the SEM envelope around the mean so we can
    # transform them
    upper_sem_bounds = mean_weights + sem_weights
    lower_sem_bounds = mean_weights - sem_weights

    # TODO compute harmonic mean as well of the rates

    # TODO: Alex wants me to calculate this just as a transformation
    # of the computed mean and std errors of the probability. Not sure
    # of the validity of this.

    shortest_run_time_values = np.array([time.value_in_unit(time_unit)
                                         for time in shortest_run_times])


    # then compute MFPTs

    # mask the values for which no data was observed (nans)

    # calculate the rates for all the runs, masking the zero values
    # which need special tolerances due to the small values. These
    # small values are dependent upon the minimum probabilities
    # allowed by the simulation so this must be a parameter.
    #masked_mean_weights = np.ma.masked_invalid(mean_weights)

    # actually calculate the average mfpts from the mean weights over
    # time
    rates = mean_weights / shortest_run_time_values

    # compute the error in the rates using the relationship between
    # the MFPT and the weight: MFPT_error = weight_error * t * weight^-2
    #rt_sems = sem_weights * shortest_run_time_values * (mean_weights**-2)
    rate_sems = sem_weights / shortest_run_time_values

    upper_sem_rates = rates + rate_sems
    lower_sem_rates = rates - rate_sems

    # clip the things that are close to zero or less than zero since
    # this is just a floating point subtraction error
    close_to_zero_idx = np.argwhere(np.logical_or(np.isclose(lower_sem_rates, 0.),
                                                  lower_sem_rates < 0.))
    lower_sem_rates[close_to_zero_idx] = 0.

    # then compute the bounds for the MFPTs
    # upper_sem_mfpts = shortest_run_time_values / upper_sem_bounds
    # lower_sem_mfpts = shortest_run_time_values / lower_sem_bounds

    # because the lower standard error bound is lower than the average
    # we want to make sure the call to isclose doesn't count them as 0
    # so we make a new min_prob for the masking that takes into
    # account the lowest value. So we take the minimum value that is
    # not identically 0
    # lower_min_prob = lower_sem_bounds[lower_sem_bounds > 0.].min() * 1e-2

    # actually we are going to consder anythin 2 orders of magnitude
    # lower than the minimum probability to be zero so:
    # lower_sem_bound_min_value = min_probability * 1e-2

    # zero_masked_upper_sem_weights = np.ma.masked_values(upper_sem_bounds, 0.,
    #                                                     atol=min_probability)

    # zero_masked_lower_sem_weights = np.ma.masked_values(lower_sem_bounds, 0.,
    #                                                     atol=lower_min_prob)

    # # we calculate the standard error of the mean weights as well
    # zero_masked_upper_sem_mfpts = (shortest_run_time_values /
    #                                zero_masked_upper_sem_weights)
    # zero_masked_lower_sem_mfpts = (shortest_run_time_values /
    #                                zero_masked_lower_sem_weights)

    # make the plots, one for aggregate probability the other for
    # RT
    fig, axes = plt.subplots(1, 1, constrained_layout=True)

    # prob_ax = axes[0]
    rate_ax = axes
    for run_idx in run_target_weights_rates.keys():

        # do the weights plot

        # get the times and weights converting if necessary
        run_weights = runs_weights[run_idx]
        run_times = [time.value_in_unit(time_unit) for time in runs_times[run_idx]]

        # slice them to the smallest run length
        run_weights = run_weights[0:min_num_timepoints]
        run_times = run_times[0:min_num_timepoints]

        label = RUN_LABEL_TEMPLATE.format(run_idx)


        # then do the residence time plots

        # get the rates values previously calculated and truncate
        run_rates = runs_rates[run_idx][0:min_num_timepoints]

        # convert the rates to the time unit specified
        run_rates = [rate.value_in_unit(rate_unit)
                     for rate in run_rates]

        replicate_color = dict(dict(REPLICATE_COLORS)[run_idx])['base']

        # plot the rate estimates for each run
        rate_ax.plot(run_times,
                     run_rates,
                     label=label,
                     linewidth='3',
                     color=replicate_color)


    shortest_run_times_values = [time.value_in_unit(time_unit) for time in shortest_run_times]
    # prob_ax.plot(shortest_run_times_values,
    #              mean_weights,
    #              label='Mean', linewidth='3', color='black')
    # prob_ax.fill_between(shortest_run_times_values,
    #                      mean_weights-sem_weights,
    #                      mean_weights+sem_weights,
    #                      label='Std. Error')

    # get the units from the value of the rates
    # rate_unit = runs_rates[list(runs_rates.keys())[0]][0].get_symbol()

    # # convert to the given unit as a rate
    # if time_unit is not None:
    #     rate_unit = 1 / (1 / rate_unit).in_units_of(time_unit)




    # plot the average rate

    # get the times and convert to the proper units
    times = [time.value_in_unit(time_unit) for time in shortest_run_times]

    # first because they were dedimensionalized we redimensionalize
    # them to the specified time_unit (x-axis) then convert them to
    # the rate unit, then dedimensionalize them to raw values for
    # plotting
    rate_values = np.array([(rate / time_unit).value_in_unit(rate_unit)
                            for rate in rates])

    upper_sem_rate_values = np.array([(rate / time_unit).value_in_unit(rate_unit)
                                      for rate in upper_sem_rates])

    lower_sem_rate_values = np.array([(rate / time_unit).value_in_unit(rate_unit)
                                      for rate in lower_sem_rates])

    # clip the things that are close to zero or less than zero since
    # this is just a floating point error from multiplying...
    close_to_zero_idx = np.argwhere(np.logical_or(np.isclose(lower_sem_rate_values, 0.),
                                                  lower_sem_rates < 0.))
    lower_sem_rate_values[close_to_zero_idx] = 0.


    # choose colors
    mean_color = dict(dict(REPLICATE_COLORS)['agg'])['base']
    sem_color = dict(dict(REPLICATE_COLORS)['agg'])['light']

    # plot means
    rate_ax.plot(times,
                 rate_values,
                 label='Mean',
                 linewidth='3',
                 color=mean_color)

    rate_ax.fill_between(times,
                         upper_sem_rate_values,
                         lower_sem_rate_values,
                         label='Std. Error of Mean',
                         color=sem_color,)

    if logscale:
        rate_ax.set_yscale('log')


    ## labels and fonts
    label_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '12'
    }

    ticks_font = {
        'family' : 'serif',
        'style' : 'normal',
        'size' : '10'
    }

    title_font = {
        'family' : 'sans-serif',
        'style' : 'normal',
        'size' : '12'
    }

    # prob_ax.set_title("Aggregate Probability",
    #                   fontdict=title_font)
    # prob_ax.set_yscale('log')
    # prob_ax.set_ylabel('Aggregated unbound probability',
    #                    fontdict=label_font)

    # # TODO set this intelligently
    # prob_ax.set_ylim([1e-12, 10])

    # prob_ax.set_xlabel('Simulation time (${}$)'.format(time_unit.get_symbol()),
    #                    fontdict=label_font)

    # # set the tick fonts
    # # for label in (prob_ax.get_xticklabels() + prob_ax.get_yticklabels()):
    # #     label.set_fontname("serif")
    # #     label.set_fontsize(10)

    # prob_ax.legend(loc='best')

    rate_ax.set_title(f"Probability Flux: {gexp}",
                    fontdict=title_font)

    # plt.xticks([400, 800, 1200, 1600], ['400', '800', '1200', '1600'])
    # plt.yticks([0.1, 10, 1000, 100000, 10000000, 1000000000])
    rate_ax.set_xlabel('Simulation time (${}$)'.format(time_unit.get_symbol()),
                     fontdict=label_font)

    rate_ax.set_ylabel('Probability Flux ($P\/{}$)'.format(rate_unit.get_symbol()),
                     fontdict=label_font)

    # was used as a default before [10e-4, 10e7]
    if ylim is not None:
        rate_ax.set_ylim(ylim)

    # plot the experimental value
    experimental_points = [experimental_rate.value_in_unit(rate_unit) for _ in run_times]
    rate_ax.plot(times,
                 experimental_points,
                 color='red',
                 label='Experimental')

    rate_ax.legend()

    # return the actual computed values for the average aggregated
    # weights and residence times
    rates = np.array([(rate / time_unit).value_in_unit(rate_unit)
                      for rate in rates]) * rate_unit

    rate_sems = np.array([(rate / time_unit).value_in_unit(rate_unit)
                          for rate in rate_sems]) * rate_unit



    return mean_weights, sem_weights, rates, rate_sems, shortest_run_times, (fig, rate_ax)

def gexp_plot_rates(gexp):

    from wepy.analysis.rates import contig_warp_rates

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    YLIMS = {
        '3' : None,
        '10' : None,
        '17' : None,
        '18' : None,
        '20' : None,
        'TPPU-legaxy' : None,
    }

    # get the experimental koff

    experimental_rate = dict(experimental_koffs())[lig_id]

    span_rates = {}
    for span_id in get_gexp_span_ids(gexp):
        contig = get_gexp_span_contig(gexp, span_id)

        num_walkers = get_gexp_span_num_walkers(gexp, span_id)

        with contig:

            # get the values as a series
            contig_rates = contig_warp_rates(contig, CYCLE_TIME, time_points=Ellipsis)

            span_rates[span_id] = contig_rates

    # make the plot for each span
    result = \
                              plot_rates(gexp,
                                         span_rates,
                                         CYCLE_TIME,
                                         num_walkers,
                                         experimental_rate,
                                         target_idx=0,
                                         rate_unit=(1/tkunit.second).unit,
                                         time_unit=tkunit.microsecond,
                                         ylim=YLIMS[gexp],
                                         logscale=False)

    return result

def gexp_plot_rts(gexp):

    from wepy.analysis.rates import contig_warp_rates

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    YLIMS = {
        '3' : None,
        '10' : None,
        '17' : None,
        '18' : None,
        '20' : None,
        'TPPU-legaxy' : None,
    }

    # get the experimental mfpt/rt for this ligand in the correct unit
    experimental_halftime = dict(experimental_halftimes())[lig_id]

    # these are the "rates" structs from the function ignore that this
    # is a function for rts
    span_rates = {}
    for span_id in get_gexp_span_ids(gexp):
        contig = get_gexp_span_contig(gexp, span_id)

        num_walkers = get_gexp_span_num_walkers(gexp, span_id)

        with contig:

            # get the values as a series
            contig_rates = contig_warp_rates(contig, CYCLE_TIME, time_points=Ellipsis)

            span_rates[span_id] = contig_rates

    # make the plot for each span
    result = \
                              plot_rts(gexp,
                                       span_rates,
                                       CYCLE_TIME,
                                       num_walkers,
                                       experimental_halftime,
                                       target_idx=0,
                                       mfpt_unit=tkunit.minute,
                                       time_unit=tkunit.microsecond,
                                       ylim=YLIMS[gexp],
                                       logscale=True)

    return result

def gexp_show_plot_rates(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, rates, rate_sems, shortest_run_times, plots = \
                                                            gexp_plot_rates(gexp)

    print("Gexp: {}".format(gexp))

    print("Final Rate: {}".format(rates[-1]))
    print("Final Rate SEM: {}".format(rate_sems[-1]))

    plt.show()

def save_gexp_plot_rates(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, rates, rate_sems, shortest_run_times, plots = \
                                                            gexp_plot_rates(gexp)

    print("Gexp: {}".format(gexp))

    fig, axes = plots

    save_gexp_fig(gexp, "rates", fig)


def gexp_show_plot_rts(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, rts, rt_sems, shortest_run_times, plots = \
                                                            gexp_plot_rts(gexp)

    # TODO: add option to get this as halftime

    print("Gexp: {}".format(gexp))

    print("Final Residence Time: {}".format(rts[-1]))
    print("Final Residence Time SEM: {}".format(rt_sems[-1]))

    plt.show()

def save_gexp_plot_rts(gexp):

    import matplotlib.pyplot as plt

    mean_weights, sem_weights, rts, rt_sems, shortest_run_times, plots = \
                                                            gexp_plot_rts(gexp)


    # TODO: add option to get this as halftime

    print("Gexp: {}".format(gexp))

    fig, axes = plots

    save_gexp_fig(gexp, "residence_times", fig)

def box_volume_observable(fields):

    import numpy as np

    from wepy.util.util import traj_box_vectors_to_lengths_angles

    all_lengths, _ = traj_box_vectors_to_lengths_angles(fields['box_vectors'])

    volumes = np.multiply.reduce(all_lengths, axis=1)

    print("Computing on a chunk")

    return volumes

def lig_bs_rmsd_observable(lig_id, ligand_idxs, fields):

    import numpy as np

    from geomm.rmsd import calc_rmsd

    from wepy.util.util import traj_box_vectors_to_lengths_angles

    print("starting chunk calculation")

    # recenter the whole protein-ligand complex into the center of
    # the periodic boundary conditions, this uses the whole set of
    # receptor atom indices and not just the binding site.

    # TODO: do we need main_rep here? doesn't matter too much and I
    # already computed this

    # then also superimpose everything to the reference state binding
    # site indices
    sup_positions, ref_positions = recenter_superimpose_traj(
        fields,
        lig_id,
        'main_rep'
    )

    # then calculate the RMSD of the ligand of all the frames to the

    # reference state
    lig_rmsds = np.array([
        calc_rmsd(
            ref_positions,
            frame,
            idxs=ligand_idxs,
        )
        for frame in sup_positions])

    print("ending chunk calculation")

    return lig_rmsds

def lig_sasa_observable(
        lig_id,
        ligand_idxs,
        json_top,
        n_sphere_points,
        traj_fields):

    import time

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    import mdtraj as mdj

    shape = traj_fields['positions'].shape

    print(f"Computing SASAs for a chunk of size: {shape}")

    start = time.time()

    # ALERT: this needs to be main rep here for the calculation
    sup_positions, ref_positions = recenter_superimpose_traj(
        traj_fields,
        lig_id,
        'main_rep'
    )

    traj_fields['positions'] = sup_positions

    traj = traj_fields_to_mdtraj(traj_fields, json_top)

    # compute sasas with all atoms present then slice the ligand idxs
    # off, then sum up the values for each atom in a frame
    sasas = mdj.shrake_rupley(
        traj,
        n_sphere_points=n_sphere_points,
    )[:,ligand_idxs].sum(axis=1)

    end = time.time()

    duration = end - start

    print(f"Finished Chunk, took: {duration}")

    return sasas

def lig_prot_atom_pair_observable(pair_idxs, top, fields):

    from msmbuilder.featurizer import AtomPairsFeaturizer

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    n_frames = fields['positions'].shape[0]

    print(f"Starting a chunk, with # of frames: {n_frames} ")

    featurizer = AtomPairsFeaturizer(pair_idxs, periodic=True)

    # convert the fields to an mdtraj trajectory
    traj = traj_fields_to_mdtraj(fields, top)

    features = featurizer.partial_transform(traj)

    print(f"finished a chunk, with # of frames: {n_frames}")

    return features

def pc_projections_observable(lig_id, hom_idxs, model, fields):

    from geomm.centroid import centroid

    # TODO: do we need the main reps here because this is about
    # homology and that is special
    hom_positions = recenter_superimpose_traj(
        fields,
        lig_id,
        'main_rep'
    )[0][:,hom_idxs,:]

    coms = np.array([centroid(frame) for frame in hom_positions])

    center_projections = model.transform(coms)

    return center_projections

def legacy_pc_projections_observable(lig_id, hom_idxs, model, fields):

    import numpy as np

    import mdtraj as mdj

    from wepy.util.util import traj_box_vectors_to_lengths_angles

    from geomm.superimpose import superimpose
    from geomm.grouping import group_pair
    from geomm.centering import center_around
    from geomm.centroid import centroid

    sel_idxs = lig_selection_idxs(lig_id)

    lig_idxs = sel_idxs['correct_rep/ligand']
    prot_idxs = sel_idxs['correct_rep/protein']
    bs_idxs = sel_idxs['correct_rep/binding_site']

    box_lengths, _ = traj_box_vectors_to_lengths_angles(fields['box_vectors'])

    ref_traj = mdj.load_pdb(osp.join(data_path(), 'top/{}/real_rep_center_ref.pdb'.format(lig_id)))

    centered_ref_positions = ref_traj.xyz[0]

    ## regroup, center, and superimpose the frames

    # group the pair of ligand and binding site together in the same image
    grouped_positions = [group_pair(positions, box_lengths[idx],
                                        bs_idxs, lig_idxs)
                  for idx, positions in enumerate(fields['positions'])]

    # center all the positions around the binding site
    centered_positions = [center_around(positions, bs_idxs)
                          for idx, positions in enumerate(grouped_positions)]

    # then superimpose the binding sites
    sup_positions = np.array([superimpose(centered_ref_positions, pos, idxs=bs_idxs)[0]
                     for pos in centered_positions])


    hom_positions = sup_positions[:,hom_idxs,:]

    coms = np.array([centroid(frame) for frame in hom_positions])

    center_projections = model.transform(coms)

    return center_projections

MAX_FE_BUFFER = 5
MIN_FE_AXIS_VALUE = 0.0


def plot_fe_profile(
        fe_profile,
        bin_edges,
        max_fe=None,
        title="Free Energy Profile",
        observable_label="Observable",
):
    """Plots a single, already calculated free energy profile and
    binning spec.

    Used only for configuring the actual plot detailes like axes and
    labels.

    See Also
    --------

    plot_fe_profiles :: for plotting multiple curves at once
    """

    import numpy as np
    import matplotlib.pyplot as plt



    # if no maximum free energy is given set it to the maximum valid
    # value
    if max_fe is None:
        max_fe = np.ma.masked_invalid(fe_profile).max()

    bin_centers = np.array([(bin_edges[i] + (bin_edges[i + 1] - bin_edges[i]))
                            for i in range(bin_edges.shape[0] - 1)])

    # construct the Figure and Axes objects. We want to be able to
    # collect the individual Axes objects into a single Figure later
    # or have a standalone figure
    fig, ax = plt.subplots()

    # add the figure Artists

    # add a super title to the figure
    title = fig.suptitle(title)

    # add the axes artists
    artists = ax.plot(bin_centers, fe_profile)

    ax.set_xlim(left=bin_edges[0], right=bin_edges[-1])
    ax.set_ylim(bottom=MIN_FE_AXIS_VALUE, top=max_fe+MAX_FE_BUFFER)

    ax.set_xlabel(observable_label)
    ax.set_ylabel("Free Energy ($-ln(p)$)")

    return fig, ax

def plot_fe_profiles(
        fe_profiles,
        bin_edges,
        max_fe=None,
        labels=None,
        title="Free Energy Profiles",
        observable_label="Observable",
):
    """Plot multiple already calculated FE curves.

    This is just for configuring details of presentation.

    """

    import numpy as np
    import matplotlib.pyplot as plt

    # if no maximum free energy is given set it to the maximum valid
    # value
    if max_fe is None:
        max_fe = np.ma.masked_invalid(np.concatenate(fe_profiles)).max()

    # get the bin centers for plotting a line graph over the bins
    bin_centers = np.array([(bin_edges[i] + (bin_edges[i + 1] - bin_edges[i]))
                            for i in range(bin_edges.shape[0] - 1)])


    # construct the Figure and Axes objects. We want to be able to
    # collect the individual Axes objects into a single Figure later
    # or have a standalone figure
    fig, ax = plt.subplots()
    title = fig.suptitle(title)

    # if this is the last curve make the label show the true number of
    # cycles as the upper limit
    for curve_idx, profile in enumerate(fe_profiles):

        if labels is None:
            label = str(curve_idx)
        else:
            label = labels[curve_idx]

        artists = ax.plot(bin_centers, profile,
                          label=label)


    ax.set_xlim(left=bin_edges[0], right=bin_edges[-1])
    ax.set_ylim(bottom=MIN_FE_AXIS_VALUE, top=max_fe + MAX_FE_BUFFER)

    ax.set_xlabel(observable_label)
    ax.set_ylabel("Free Energy ($-ln(p)$)")

    artists.append(ax.legend())

    return fig, ax

def plot_lig_obs_fe_spans(
        gexp,
        field_key,
        observable_label,
        bin_method='auto',
):
    """Calculate a FE profile for each replicate/span in a GEXP dataset."""

    from wepy.analysis.profiles import ContigTreeProfiler, contigtrees_bin_edges

    contigtree = get_contigtree(gexp)
    set_recursion_limit()
    with contigtree:

        profiler = ContigTreeProfiler(contigtree)

        # get the bin edges for the whole contig
        bin_edges = profiler.bin_edges(bin_method, field_key)
        print("Num bins {}".format(len(bin_edges) - 1))

        bin_centers = bin_centers_from_edges(bin_edges)

        # profile the finale FE profile spans for the ligand
        span_profiles = []
        for span_idx in profiler.contigtree.span_traces.keys():
            print("span {}".format(span_idx))
            # then make a fe profile with those bins for the first span
            fe_profile = profiler.fe_profile(span_idx, field_key, bins=bin_edges)

            span_profiles.append(fe_profile)

        fig, ax = plot_fe_profiles(span_profiles, bin_edges,
                                   title="GEXP {} span FE profiles".format(gexp),
                                   observable_label=observable_label)

    return fig, ax

def plot_contigtrees_spans_observable_convergence(
        contigtrees,
        field_key,
        obs_label,
        num_partitions=5,
        bin_method='sqrt',
):
    """Plot a timeseries of free energy profiles for every
    span/replicate.

    Specify how many timepoints you want with the 'num_partitions'.

    """

    import numpy as np
    import matplotlib.pyplot as plt

    from wepy.analysis.profiles import ContigTreeProfiler, contigtrees_bin_edges


    # get the bin edges for all of the contigtrees for all ligands
    bin_edges = contigtrees_bin_edges(contigtrees,
                                      bin_method,
                                      field_key)

    # then get the profiles for each span of each contig
    contigtree_span_convergences = []
    contigtree_span_cum_num_cycles = []
    for contigtree in contigtrees:

        with contigtree.wepy_h5:

            profiler = ContigTreeProfiler(contigtree)

            # save the convergences for each span
            spans_convergences = []
            spans_cum_num_cycles = []

            # iterate over the spans in sorted order
            span_idxs = list(contigtree.span_traces.keys())
            span_idxs.sort()

            for span_idx in span_idxs:

                fe_profiles, num_cycles = \
                    profiler.fe_cumulative_profiles(
                        span_idx,
                        field_key,
                        bins=bin_edges,
                        num_partitions=num_partitions,
                        ignore_truncate=True
                    )

                spans_convergences.append(fe_profiles)
                spans_cum_num_cycles.append(num_cycles)

            contigtree_span_convergences.append(spans_convergences)
            contigtree_span_cum_num_cycles.append(spans_cum_num_cycles)

    # get the maximum FE from all of the simulations
    max_fe = np.ma.masked_invalid(
        [np.concatenate(
            [np.concatenate(profiles) for profiles in spans])
         for spans in contigtree_span_convergences]
        ).max()

    figs = []
    for contigtree_idx, spans_convergences in enumerate(contigtree_span_convergences):

        for span_idx, fe_profiles in enumerate(spans_convergences):

            num_cycles = contigtree_span_cum_num_cycles[contigtree_idx][span_idx]

            labels = ["Cycles: {}-{}".format(0, num_cycles[curve_idx])
                      for curve_idx, fe_profile in enumerate(fe_profiles)]

            title = \
        f"Free Energy Convergence of {obs_label}: GEXP: {contigtree_idx} span {span_idx}"

            fig, ax = plot_fe_profiles(
                fe_profiles,
                bin_edges,
                max_fe=max_fe,
                title=title,
                labels=labels)

            figs.append((fig, ax))

    return figs

@jlmem.cache
def contigtrees_observable_fe(
        contigtrees,
        field_key,
        bin_method='auto',
):

    import numpy as np
    import matplotlib.pyplot as plt

    from wepy.analysis.profiles import (
        ContigTreeProfiler,
        contigtrees_bin_edges,
    )

    # DEBUG
    print("calculating bin edges across all the gexps")

    # get the bin edges for all of the contigtrees for all ligands
    bin_edges = contigtrees_bin_edges(
        contigtrees,
        bin_method,
        field_key,
    )

    gc.collect()

    # DEBUG
    print("finished calculating bin edges across all the gexps")


    # DEBUG
    print("Calculating FE profiles for each GEXP contigtree")

    # then get the profiles for each span of each contig
    fe_profiles = []
    for idx, contigtree in enumerate(contigtrees):

        # DEBUG
        print(f"Calculating for GEXP contigtree: {idx}")


        with contigtree.wepy_h5:

            profiler = ContigTreeProfiler(contigtree)

            fe_profile = profiler.fe_profile_all(field_key, bins=bin_edges)

            fe_profiles.append(fe_profile)

        # DEBUG
        print(f"Finished calculating for GEXP contigtree: {idx}")

    # DEBUG
    print(f"Finished calculating for all gexps")


    return fe_profiles, bin_edges

def plot_contigtrees_trace_observable_fe(
        contigtrees,
        traces,
        observable_key,
        bin_method='sqrt',
):
    """Generate a FE profile for each contigtree of only the provided
    trace.

    """

    import numpy as np
    import matplotlib.pyplot as plt

    from wepy.analysis.profiles import ContigTreeProfiler, contigtrees_bin_edges


    field_key = "observables/{}".format(observable_key)

    # get the bin edges for all of the contigtrees for all ligands
    bin_edges = contigtrees_bin_edges(contigtrees,
                                      bin_method, field_key)

    # then get the profiles for each span of each contig
    fe_profiles = []
    for tree_idx, contigtree in enumerate(contigtrees):

        with contigtree.wepy_h5:

            profiler = ContigTreeProfiler(contigtree)

            fe_profile = profiler.fe_profile_trace(traces[tree_idx], field_key,
                                                   bins=bin_edges)

            fe_profiles.append(fe_profile)

    return fe_profiles, bin_edges

# TODO: not sure this is even necessary as its own function

# def plot_contigtrees_trace_observable_fe(
#         fe_profiles,
#         bin_edges,
#         observable_key,
#         trace_title=None,
# ):
#     """Plot FE profiles for 
#     """

#     # get the maximum FE from all of the simulations
#     max_fe = np.ma.masked_invalid(np.concatenate(fe_profiles)).max()

#     for contigtree_idx, fe_profile in enumerate(fe_profiles):

#         labels = ["{}".format(contigtree_idx)]

#         title = "Free Energy {}".format(observable_key)

#         fig, ax = plot_fe_profiles(
#             fe_profiles,
#             bin_edges,
#             max_fe=max_fe,
#             title=title,
#             labels=labels,
#         )

OBS_SPAN_FE_OPTS = {
      'lig_rmsd' : {
          'label' : "Ligand RMSD $\AA$",
          'bin_method' : 'auto',
      },

      'lig_sasa' : {
          'label' : "Ligand SASA $\AA^{2}$",
          'bin_method' : 'auto',
      },
}

GEXP_SPAN_FE_OPTS = {
    '3' : {
        'convergence_n_partitions' : 5,
    },

    '10' : {
        'convergence_n_partitions' : 5,
    },

    '17' : {
        'convergence_n_partitions' : 5,
    },

    '18' : {
        'convergence_n_partitions' : 5,
    },

    '20' : {
        'convergence_n_partitions' : 5,
    },
}

def plot_gexp_obs_span_fe(
        gexp,
        observable,
):

    field_key = f'observables/{observable}'

    OBS_OPTS = OBS_SPAN_FE_OPTS[observable]

    figs, axs = plot_lig_obs_fe_spans(
        gexp,
        field_key,
        OBS_OPTS['label'],
        bin_method=OBS_OPTS['bin_method'],
    )

    return figs, axs

def plot_gexp_obs_span_convergence_fe(
        gexp,
        observable,
):

    field_key = f'observables/{observable}'

    OBS_OPTS = OBS_SPAN_FE_OPTS[observable]
    GEXP_OPTS = GEXP_SPAN_FE_OPTS[gexp]

    # get the contigtree
    contigtree = get_contigtree(gexp)

    figs = plot_contigtrees_spans_observable_convergence(
        [contigtree],
        field_key,
        OBS_OPTS['label'],
        num_partitions=GEXP_OPTS['convergence_n_partitions'],
        bin_method=OBS_OPTS['bin_method'],
    )

    # there is just one here so pop that out for return as: fig, ax
    return figs

FE_ALL_AGG_OPTS = {

    # which gexps to include
    'gexps' : ['3', '10', '17', '18', '20'],

    # options for each observable
    'obs_opts' : {

        'lig_rmsd' : {
            'bin_method' : 'auto',
        },
    },

}

def plot_all_agg_fe_obs(observable_key):

    field_key = "observables/{}".format(observable_key)

    # get the gexp contigtrees requested
    contigtrees = [get_contigtree(gexp) for gexp in FE_ALL_AGG_OPTS['gexps']]

    # get the options for this observable
    OBS_OPTS = FE_ALL_AGG_OPTS['obs_opts'][observable_key]

    # compute the profiles and the bin edges for all of GEXP
    # contigtrees
    fe_profiles, bin_edges = contigtrees_observable_fe(
        contigtrees,
        field_key,
        bin_method=OBS_OPTS['bin_method'],
    )

    # now plot the fe profiles

    title = "Free Energy {}".format(observable_key)
    labels = [f"{gexp}" for gexp in FE_ALL_AGG_OPTS['gexps']]

    # get the maximum FE from all of the simulations
    max_fe = np.ma.masked_invalid(np.concatenate(fe_profiles)).max()

    fig, ax = plot_fe_profiles(
        fe_profiles,
        bin_edges,
        max_fe=max_fe,
        title=title,
        labels=labels,
    )

    return fig, ax

def gexp_show_plot_spans_fe_obs(gexp, observable):
    """Render the Spans/Replicates chosen plot for the GEXP"""

    import matplotlib.pyplot as plt

    fig, axs = plot_gexp_obs_span_fe(gexp, observable)

    print(f"Obs: {observable}  Gexp: {gexp}")

    plt.show()

def gexp_save_plot_spans_fe_obs(gexp, observable):

    import matplotlib.pyplot as plt

    fig, axs = plot_gexp_obs_span_fe(gexp, observable)

    print(f"Obs: {observable}  Gexp: {gexp}")

    save_gexp_fig(gexp, f'fe_profiles/spans/{observable}', fig)

def gexp_show_plot_spans_convergence_fe_obs(gexp, observable):
    """Render the Spans/Replicates chosen plot for the GEXP"""

    import matplotlib.pyplot as plt

    figs = plot_gexp_obs_span_convergence_fe(gexp, observable)

    print(f"Obs: {observable}  Gexp: {gexp}")

    plt.show()

def gexp_save_plot_spans_convergence_fe_obs(gexp, observable):

    import matplotlib.pyplot as plt

    figs = plot_gexp_obs_span_convergence_fe(gexp, observable)

    print(f"Obs: {observable}  Gexp: {gexp}")

    for span_idx, figax in enumerate(figs):
        fig, ax = figax
        print(f"Span Idx: {span_idx}")

        save_gexp_fig(
            gexp,
            f'fe_profiles/spans_convergence/{observable}',
            fig,
            tags={'span' : f'{span_idx}'}
        )

def all_render_plot_agg_fe_obs(observable,
                               save=False,
                               show=True,
):
    """Render the Aggregate FEP for the observable of all GEXPs
    together.

    If save and show are both turned on will show first, then save
    that way you can Ctrl-C to interrupt saving if you don't like what
    it showed. (or just rerun).

    """


    fig, ax = plot_all_agg_fe_obs(observable)

    print(f"All Agg FE profile for Obs: {observable}")

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    if save:
        save_fig(
            f'fe_profiles/all_agg/{observable}',
            f"{observable}",
            fig,
        )

def plot_ts_pc0_span_fe(lig_id):

    import matplotlib.pyplot as plt

    field_key = 'observables/TS-warp+featdist-0.07_pc-0'
    bin_method = 'auto'
    observable_label = "TS PC-0"

    print("Plotting")
    figs, axs = plot_lig_obs_fe_spans(lig_id, field_key, observable_label, bin_method=bin_method)

    plt.show()

def plot_ts_pc1_span_fe(lig_id):

    import matplotlib.pyplot as plt

    field_key = 'observables/TS-warp+featdist-0.07_pc-1'
    bin_method = 'auto'
    observable_label = "TS PC-1"

    print("Plotting")
    figs, axs = plot_lig_obs_fe_spans(lig_id, field_key, observable_label, bin_method=bin_method)

    plt.show()

from sklearn.cluster import (
    MiniBatchKMeans,
    KMeans,
)
from msmbuilder.cluster.kcenters import _KCenters as KCenters


# CLUST_METHODS = {
#     'KCenters' : KCenters,
#     'MiniBatchKMeans' : MiniBatchKMeans,
# }

# CLUST_METRICS = {
#     'canberra' : 'canberra',
#     'euclidean' : 'euclidean',
# }

# this is the configuration used in the ACS poster etc.
OLD_CLUSTER_SPECS = {

    'subsampling' : {
        'test_proportion' : 0.5,
    },
    'clust_method' : 'KCenters',
    'params' : {
        'metric' : 'canberra',
        'random_state' : 1,
        'n_clusters' : 1000,
    },
}


# "Normal" Clustering Specs, the key is the classifier ID
CLUSTER_SPECS = {

    '0' : {
        'clf_class' : KCenters,
        'clf_kwargs' : {
            'n_clusters' : 1000,
            'metric' : 'canberra',
            'random_state' : 1,
        },
        'splitter_class' : None,
        'splitter_kwargs' : None,
    },
}

# this is used as a "null" splitting, so that we handle split and
# unsplit data the same
class NoSplit():
    """Splitter implementing API, but doesn't split at all. All features
    go to training."""

    def __init__(self, **kwargs):

        pass

    def split(self, X, **kwargs):

        train_idxs = list(range(len(X)))
        test_idxs = []

        # yield here to satisfy the splitter API, we only have one set
        # to do, unlike other splitters
        yield train_idxs, test_idxs

class StaticSplit():
    """Split according to a predetermined set of splits."""

    def __init__(self,
                 splits,
                 **kwargs):
        """

        Parameters
        ----------

        splits : list of (list of int, list of int)
            The splits that you want to have it yield.

        """

        self._splits = splits

    def split(self, X, **kwargs):

        # just yield the ones we stored
        for split in self._splits:
            yield split


def gexp_classify_observable(
        gexp,
        observable_key,
        clf,
        clf_id,
        splitter=None,
):
    """Classify an observable given a classifier.

    Optionally, provide a splitter to partition the data into train
    test. This isn't for testing, but just efficiency of not having to
    cluster everything. That means it will just grab the first splitting
    That the splitter generates.

    """

    # if no splitting was request we just use the NoSplit, which has
    # the same API.
    if splitter is None:
        splitter = NoSplit()

    # get the observables, reshape them to features, and vectorize them if needed
    all_features = unscalarize_features(
        reshape_run_to_features(
            get_observable(
                observable_key,
                gexp,
                source='h5',
            )))

    # we only use the first splitting.
    train_idxs, test_idxs = next(splitter.split(all_features))

    # we return it as a dict
    train_test_d = {
        'training_idxs' : train_idxs,
        'testing_idxs' : test_idxs,
    }

    # we need to clean up memory so remove all_features from memory
    del all_features
    gc.collect()


    # we get the train features into memory
    train_features = unscalarize_features(
        reshape_run_to_features(
            get_observable(
                observable_key,
                gexp,
                source='h5',
            )
        )
    )[train_idxs]

    ## do the classification
    print("Running the Fitting")
    start = time.time()

    clf.fit(train_features)

    end = time.time()

    duration_s = end - start
    duration_m = duration_s / 60
    duration_h = duration_m / 60

    print("Finished fitting")
    print(f"Took: {duration_s} s OR {duration_m} m OR {duration_h} h")

    # retrive and name training labels
    train_labels = clf.labels_

    del train_features
    gc.collect()

    # now we can get the test features and classify them
    test_features = unscalarize_features(
        reshape_run_to_features(
            get_observable(
                observable_key,
                gexp,
                source='h5',
            )
        )
    )[test_idxs]

    # Classify the test features

    print("Classifying Test Features")

    # only do this if there are any to do it on
    if len(test_features) > 0:
        test_labels = clf.predict(test_features)
    else:
        test_labels = []

    print("Finished Classifying Test Features")

    del test_features
    gc.collect()

    # we can unshuffle the observables using the indices and put them
    # back how we need them for them to be observables.
    labels = unshuffle_features(
        train_idxs,
        test_idxs,
        train_labels,
        test_labels,
    )

    return labels, clf, train_test_d

def do_gexp_classification(
        gexp,
        observable_key,
        clf_id,
):

    # get the spec for this classifier id
    clf_specs = CLUSTER_SPECS[clf_id]

    # create the classifier and splitter given the spec
    clf = clf_specs['clf_class'](**clf_specs['clf_kwargs'])

    if clf_specs['splitter_class'] is not None:
        splitter = clf_specs['splitter_class'](**clf_specs['splitter_kwargs'])
    else:
        splitter = None

    # perform the classification for an observable and a gexp
    assignments, clf, train_test_d = gexp_classify_observable(
        gexp,
        observable_key,
        clf,
        clf_id,
        splitter=splitter,
    )


    # make the assignment labels the proper shape to set as an
    # observable
    assignments_obs = reshape_features_to_run(
        assignments,
        gexp,
    )

    # save the assignments

    # as an observable file
    assg_obs_name = f"clust-assigs_clfid-{clf_id}"

    print(f"Saving assignments as observable {assg_obs_name}")
    save_observable(
        assg_obs_name,
        gexp,
        assignments_obs,
    )

    # then attach it to the HDF5 as an observable

    # first check if its in the HDF5, if it is delete it then add it
    wepy_h5 = get_gexp_wepy_h5(gexp)
    delete = False
    with wepy_h5:
        # see if its in the observable names
        if assg_obs_name in wepy_h5.observable_field_names:
            delete = True

    if delete:
        delete_wepy_h5_observable(
            gexp,
            assg_obs_name,
        )

    print("Saving to HDF5")
    attach_observable(
        assg_obs_name,
        gexp,
    )

    # serialize the splitting for this classifier so we can get the
    # cluster centers and stuff back out
    print("Saving train-test splits")
    save_clustering_traintest(
        clf_id,
        gexp,
        train_test_d['training_idxs'],
        train_test_d['testing_idxs'],
    )

    print("Saving clustering model")
    save_clustering_model(
        clf_id,
        gexp,
        clf,
    )

def unshuffle_features(
        train_idxs,
        test_idxs,
        train_features,
        test_features,
):
    """Unshuffle all of the features (which is both the full set of
    training and test features) features given the splitting indices.

    We need to shuffle our data for splitting training and testing
    data, this is to unshuffle them.

    """

    n_samples = len(train_idxs) + len(test_idxs)

    unshuffled_features = np.empty(
        (n_samples, *train_features.shape[1:]),
        dtype=train_features.dtype)

    for idx, shuf_idx in enumerate(train_idxs):
        unshuffled_features[shuf_idx] = train_features[idx]

    for idx, shuf_idx in enumerate(test_idxs):
        unshuffled_features[shuf_idx] = test_features[idx]

    return unshuffled_features

def unshuffle_feature_idxs(
        shuffle_mapping,
        shuffled_features_idxs,
):
    """For a list feature idxs find their unshuffled feature idxs."""

    # feature_idxs = [None for i in range(len(shuffled_features_idxs))]

    feature_idxs = []
    for shuffled_feature_idx in shuffled_features_idxs:

        # find the feature index from the shuffled feature
        feature_idx = shuffle_mapping.index(shuffled_feature_idx)

        feature_idxs.append(feature_idx)

    return feature_idxs

def get_unshuffled_feature_idxs(
        clf_id,
        gexp,
        clf_feature_idxs,
):

    # get the splitting and the training index, shuffling for this
    # classifier
    train_test_d = get_clustering_traintest(
        clf_id,
        gexp,
    )

    return unshuffle_feature_idxs(
        train_test_d['training_idxs'],
        clf_feature_idxs,
    )

def clf_features_to_run_trace_idxs(
        gexp,
        clf_id,
        clf_feature_idxs,
):
    """Convert a list of features relative to the classifier to a
    equivalent run idxs (i.e. (run, traj, cycle)).

    """

    # because the clf_feature_idxs are in terms of the shuffled
    # features we need to get them unshuffled
    feature_idxs = get_unshuffled_feature_idxs(
        clf_id,
        gexp,
        clf_feature_idxs,
    )

    trace = []
    for node_idx, feature_idx in enumerate(feature_idxs):

        # using the unshuffled feature idx, get the trace index:
        # i.e. (run, traj, cycle)
        trace_idx = feature_idx_to_trace_idx(feature_idx, gexp)

        trace.append(trace_idx)

    return trace

# MSM specs. Keys are 'msm_id'

MSM_SPECS = {

    'cutoff-10' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder-unweighted',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 10,
        },
        'transition_prob_method' : 'nonreversible-normalize',
        'transition_prob_kwargs' : {
            'reversible' : False,
        },
    },

    'cutoff-5' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder-unweighted',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 5,
        },
        'transition_prob_method' : 'nonreversible-normalize',
        'transition_prob_kwargs' : {
            'reversible' : False,
        },
    },

    'cutoff-2' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder-unweighted',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 2,
        },
        'transition_prob_method' : 'nonreversible-normalize',
        'transition_prob_kwargs' : {
            'reversible' : False,
        },
    },

    'cutoff-1' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder-unweighted',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 1,
        },
        'transition_prob_method' : 'nonreversible-normalize',
        'transition_prob_kwargs' : {
            'reversible' : False,
        },
    },


}

    ## weighted methods
OLD_MSM_SPECS= {

    '0' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder', #'csnanalysis',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : seh_params.PMIN, # 1e-16
        },
        'transition_prob_method' : 'mle-pyemma',
        'transition_prob_kwargs' : {
            'reversible' : True,
            'eps_mu' : 1e-20,
            },
    },

    'pmin_irreversible' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder', #'csnanalysis',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : seh_params.PMIN, # 1e-16
        },
        'transition_prob_method' : 'mle-pyemma',
        'transition_prob_kwargs' : {
            'reversible' : False,
            'eps_mu' : 1e-20,
            },
    },

    'pmin_transpose' : {
        'csn_id' : '0',
        'trim_method' : 'msmbuilder', #'csnanalysis',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : seh_params.PMIN, # 1e-16
        },
        'transition_prob_method' : 'transpose',
        'transition_prob_kwargs' : {},
    },

    'big_cutoff' : {
        'csn_id' : '0',
        'lag_time' : 2,
        'obs_name' : 'clust-assigs_clfid-0',
        'clf_id' : '0',
        'trim_method' : 'msmbuilder', #'csnanalysis',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 1.0,
        },
        'transition_prob_method' : 'mle-pyemma',
        'transition_prob_kwargs' : {
            'reversible' : True,
            'eps_mu' : 1e-6,
        },
    },

    'small-cutoff_0' : {
        'csn_id' : '0',
        'lag_time' : 2,
        'obs_name' : 'clust-assigs_clfid-0',
        'clf_id' : '0',
        'trim_method' : 'msmbuilder', #'csnanalysis',
        'trim_kwargs' : {
            'ergodic_trim_cutoff' : 0.5e-11,
        },
        'transition_prob_method' : 'mle-pyemma',
        'transition_prob_kwargs' : {
            'reversible' : True,
            'eps_mu' : 1e-15,
        },
    },

}

def make_macro_state_network_label(
        gexp,
        label_features,
        lag_time=2
):
    """Create a MacroStateNetwork for a gexp data given some features
    which are labels.

    This is a minimal function and tries not to compute anything that
    isn't necessary over the network.

    Parameters
    ----------

    gexp : str

    label_features : nested lists of shape n_runs x n_trajs x n_frames
                     of arraylikes of shape feature_shape.

    lag_time : int

    """

    from wepy.analysis.network import MacroStateNetwork

    contigtree = get_contigtree(gexp)

    label_features = reshape_features_to_run(label_features,
                                             gexp)

    run_labels = label_features
    del label_features

    msn = MacroStateNetwork(contigtree,
                            transition_lag_time=lag_time,
                            assignments=run_labels,
                            )

    return msn

def make_macro_state_network_obs(
        gexp,
        obs_name,
        lag_time=2
):
    """Create a MacroStateNetwork for a gexp data given a name of an
       observable that is a label.

    """

    from wepy.analysis.network import MacroStateNetwork

    assg_field_key = f"observables/{obs_name}"

    contigtree = get_contigtree(gexp)

    msn = MacroStateNetwork(contigtree,
                            transition_lag_time=lag_time,
                            assg_field_key=assg_field_key,
                            )

    return msn

CSN_SPECS = {
    '0' : {
        'lag_time' : 2,
        'obs_name' : 'clust-assigs_clfid-0',
        'clf_id' : '0',
    }
}

@jlmem.cache
def make_basic_csn(
        gexp,
        csn_id,
):
    """Get the basic MSN that is derived directly from the classifier and
    simulation data.

    This will not have any additional computed fields and analyses for
    the nodes that may have been saved, see the `get_network` function for
    getting this one.

    """

    observable = CSN_SPECS[csn_id]['obs_name']
    lag_time = CSN_SPECS[csn_id]['lag_time']
    clf_id = CSN_SPECS[csn_id]['clf_id']

    # using the parameters in the CSN_SPEC make the msn
    msn = make_macro_state_network_obs(
        gexp,
        observable,
        lag_time=lag_time,
    )

    # then we do some basic post-processing

    # set the macrostate weights
    with msn:
        msn.set_macrostate_weights()

    # set the index of the cluster centers for each macrostate
    #
    # both the feature index (from the raw vector of features) and the
    # trace-like index (run, traj, frame)
    clf = get_clustering_model(
        clf_id,
        gexp,
    )

    # TODO: should probably just functionalize this as a pass that is
    # dependent on the clustering algorithm

    # the feature vector index (the concatenated features)

    # note that the feature idx is for how the features were
    # clustered, which may have been split and shuffled
    msn.set_nodes_attribute(
        'center_feature_idx',
        {
            node_idx : feature_idx
            for node_idx, feature_idx
            in enumerate(clf.cluster_ids_)
        }
    )

    # get the run trace index of each cluster center.

    # this is an index over the actual data set, i.e. (run, traj, cycle)
    msn.set_nodes_attribute(
        'center_idx',
        {
            node_idx : trace_idx
            for node_idx, trace_idx
            in enumerate(
                clf_features_to_run_trace_idxs(
                    gexp,
                    clf_id,
                    clf.cluster_ids_
                )
            )
        }
    )

    # calculate the free energies
    node_fes = calc_msn_free_energies(msn)

    # set them as node attributes
    msn.set_nodes_attribute('free_energy', node_fes)


    # then we do some basic groups

    # make a group which is the node with the native state

    # find the native node
    native_node_id = classify_native_state_cluster(
        clf_id,
        gexp,
    )

    # make a group of it
    msn.set_node_group('native_state', [native_node_id])

    # TODO: make the trimmed network groups

    return msn.base_network

def calc_msn_free_energies(
        msn,
):

    from geomm.free_energy import free_energy

    node_weights = [
        (node_id, weight)
        for node_id, weight in
        msn.get_nodes_attribute('_observables/total_weight').items()
    ]

    node_ids = [node_id for node_id, weight in node_weights]
    weights = [weight for node_id, weight in node_weights]

    node_fes = {
        node_id : fe
        for node_id, fe
        in zip(
            node_ids,
            free_energy(weights)
        )
    }

    return node_fes

# @jlmem.cache
def classify_native_state_cluster(
        clf_id,
        gexp,
):
    """Determine which cluster the starting state is from.

    It does this by (re)computing the pair distances for the reference
    state and then classifying them with the trained model.

    """

    # REVD: It should be main rep, explicit above now

    # get the feature vector for the reference structure
    ref_fields = get_centered_ref_state_traj_fields(
        gexp,
        rep_key='main_rep',
    )

    lig_id = dict(GEXP_LIG_IDS)[gexp]
    top = lig_selection_tops(lig_id)['main_rep']

    # get the atom pairs for the reference state

    pair_idxs = lig_bs_atom_pairs(
        gexp,
        rep_key='main_rep',
    )
    ref_feature = lig_prot_atom_pair_observable(
        pair_idxs,
        top,
        ref_fields,
    )

    # predict which cluster it would be in
    clf = get_clustering_model(
        clf_id,
        gexp,
    )

    node_id = clf.predict(ref_feature)[0]

    return node_id

def get_cluster_micro_trace(model_name, lig_id, node_id):

    # use lag time 2 since it is cached and not relevant here
    net = get_network(model_name, lig_id, 2)

    # get the assignments from the nodes in the specified cluster
    return net.node_assignments(node_id)

def get_cluster_micro_traj(model_name, lig_id, node_id, downsample=None):

    import random

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    trace = get_cluster_micro_trace(model_name, lig_id, node_id)

    # downsample the trajectory if requested
    if downsample is not None or downsample < 1.0:
        assert downsample > 0

        n_samples = round(downsample * len(trace))

        trace = random.sample(trace, n_samples)

    with get_gexp_wepy_h5(lig_id) as wepy_h5:
        top = wepy_h5.get_topology()
        traj_fields = wepy_h5.get_trace_fields(
            trace,
            ['positions', 'box_vectors', 'alt_reps/missing'])

        traj_fields = traj_fields_to_correct_rep(traj_fields, lig_id)

    traj_fields['positions'] = recenter_superimpose_traj(
        traj_fields,
        lig_id,
        'correct_rep'
    )[0]

    traj = traj_fields_to_mdtraj(traj_fields, top)

    return traj

def state_weight_ranking(net):
    """Get the ordered ranking of the clusters according to their weights."""

    import numpy as np

    weights = [node_attrs['_observables/total_weight'] for node_id, node_attrs
               in net.graph.nodes.items()]

    ranking = np.array(net.graph.nodes)[np.argsort(weights)]

    return ranking

@jlmem.cache
def get_cluster_center_traj(
        gexp,
        csn_id,
        alt_rep=Ellipsis,
):
    """Get the cluster center trajs for a gexp. If 'alt_rep' is given as a
    valid value it will return the traj in that representation (if
    possible error otherwise). If Ellipsis is given it will choose the
    "best" one for visualization ('correct_rep' if available 'main_rep' if
    not).

    If None it will default to the common denominator for all of them (i.e. 'main_rep')

    """

    from wepy.util.mdtraj import traj_fields_to_mdtraj

    # dispatch behavior on alt_rep
    if alt_rep is None:
        alt_rep = 'main_rep'

    elif alt_rep is Ellipsis:
        if gexp in ('3',):
            alt_rep = 'main_rep'
            field_names = ['positions', 'box_vectors']

        else:
            alt_rep = 'correct_rep'
            field_names = ['positions', 'box_vectors', 'alt_reps/missing']

    else:

        if alt_rep == 'main_rep':
            field_names = ['positions', 'box_vectors']

        elif alt_rep == 'correct_rep':
            field_names = ['positions', 'box_vectors', 'alt_reps/missing']

    clf_id = CSN_SPECS[csn_id]['clf_id']

    clf = get_clustering_model(clf_id, gexp)

    # convert the feature idxs of the cluster_ids to a nice trace
    trace = clf_features_to_run_trace_idxs(
        gexp,
        clf_id,
        clf.cluster_ids_,
    )

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    sel_tops = lig_selection_tops(lig_id)
    top = sel_tops[alt_rep]

    # get the traj fields
    with get_gexp_wepy_h5(lig_id) as wepy_h5:

        traj_fields = wepy_h5.get_trace_fields(
            trace,
            field_names,
        )

    # if we are getting the correct rep we need to transform them
    # to combine the missing positions
    if alt_rep == 'correct_rep':

        traj_fields = traj_fields_to_correct_rep(
            traj_fields,
            lig_id,
        )

    # recenter them
    traj_fields['positions'] = recenter_superimpose_traj(
        traj_fields,
        lig_id,
        alt_rep,
    )[0]

    gc.collect()

    # write them
    traj = traj_fields_to_mdtraj(traj_fields, top)

    return traj

def mappable_obs_node_stats(
        msn,
        field_key,
        node_id,
        node_attrs,
        node_fields,
):
    """Function compatible with the 'node_fields_map' function for a
    MacroStateNetwork, after partially evaluated for the field key
    """

    assert field_key in node_fields

    field_data = node_fields[field_key]

    # idx of microstate as a trace for the center
    center_trace = [msn.get_node_attribute(
        0,
        'center_idx',
    )]

    # ALERT: notice the hdf5 must already be opened since this could
    # be done in parallel

    # the value of the field for it
    center_value = msn.wepy_h5.get_trace_fields(
        center_trace,
        [field_key,],
    )[field_key][0]

    data_min = np.min(field_data)
    data_max = np.max(field_data)

    node_stats = {
        'mean' : np.mean(field_data),
        'median' : np.median(field_data),
        'std-dev' : np.std(field_data),
        'min' : data_min,
        'max' : data_max,
        'range' : data_max - data_min,
        'cluster_center' : center_value,
    }

    return node_stats

def calc_msn_obs_stats(
        msn,
        observable_name,
        save=False,
):
    """Calculate the stats for an observable of a macrostate network.

    The MSN must be a full MacroStateNetwork so it has access to the
    microstates in the HDF5.

    If the 'save' parameter is True it will save the stats in a
    predetermined schema to the MacroStateNetwork.

    *This will not serialize to disk, it will just update the 'msn' object
    passed in.*

    If false this will be ignored.

    The schema is of the format (as argument to `set_nodes_observable`):

        'obs/{observable_name}/{stat_field_name}'

    Where 'stat_field_name' is the keys of the stats computed by the
    observable function, should be:

    - mean
    - median
    - std-dev
    - min
    - max
    - range
    - cluster_center

    """
    from functools import partial

    # get the full field key for the observable
    field_key = f'observables/{observable_name}'

    # partially evaluate the mapped function with the field_key
    stats_func = partial(
        mappable_obs_node_stats,
        msn,
        field_key,
    )

    # calculate the stats for each node

    # ALERT: must open for reading here and not inside the mapped
    # function since this could be in parallel
    with msn:
        node_stats = msn.node_fields_map(
            stats_func,
            [field_key,]
        )


    if save:

        # get the keys that were returned for each node
        stat_field_names = list(node_stats[list(node_stats.keys())[0]].keys())

        for stat_field_name in stat_field_names:

            # extract data for just one key
            stat_node_values = {
                node_id : stats_d[stat_field_name]
                for node_id, stats_d in node_stats.items()
            }

            # set that as an attribute
            attr_key = f'obs/{observable_name}/{stat_field_name}'

            msn.set_nodes_attribute(
                attr_key,
                stat_node_values,
            )

    return node_stats

def subset_trim(countsmat, subset_idxs):

    # K is the input states, N is the output states
    in_n_states = countsmat.shape[0]
    out_n_states = len(subset_idxs)

    # allocate trimmed matrix
    trim_countsmat = np.zeros(
        (out_n_states, out_n_states,),
        dtype=countsmat.dtype,
    )

    # copy values for the subset
    for out_idx, in_idx in enumerate(subset_idxs):

        trim_countsmat[out_idx, :] = countsmat[in_idx, subset_idxs]
        trim_countsmat[:, out_idx] = countsmat[subset_idxs, in_idx]

    # make the mapping
    mapping = {in_idx : out_idx
               for out_idx, in_idx in enumerate(subset_idxs)}

    return trim_countsmat, mapping

def ergodic_trim(
        msn,
        method='csnanalysis',
        method_kwargs=None,
):
    """Trim a counts matrix and get the largest connected component.

    Parameters
    ----------

    msn : MacroStateNetwork or BaseMacroStateNetwork

    Returns
    -------

    """

    from msmtools.estimation import \
        connected_sets as pyemma_connected_set

    from csnanalysis.csn import CSN

    from msmbuilder.msm.core import \
        _strongly_connected_subgraph as msmb_strongly_connected_subgraph

    ## methods using the unweighted counts
    if method == 'csnanalysis-unweighted':

        # get the counts matrix according to this method
        countsmat = msn.edge_attribute_to_matrix(
            'unweighted_counts',
            fill_value=0,
        )

        # for CSNAnalysis the matrix is the reverse orientation
        countsmat_T = countsmat.T

        # make a CSN
        csn = CSN(countsmat_T)

        # trim it to the main components
        csn.trim(
            min_count=method_kwargs['ergodic_trim_cutoff'],
        )

        # use the trimmed off idxs to get which states remain
        chosen_subset = csn.trim_indices

    elif method == 'msmbuilder-unweighted':

        # get the counts matrix according to this method
        countsmat = msn.edge_attribute_to_matrix(
            'unweighted_counts',
            fill_value=0,
        )

        # this counts matrix is in the correct orientation for
        # msmbuilder method

        trimmed_countsmat, mapping, _ = \
                msmb_strongly_connected_subgraph(
                    countsmat,
                    weight=method_kwargs['ergodic_trim_cutoff'],
                    verbose=False,
                    )

        # the mapping here is in terms of the node_idxs

        # chosen_subset is node_idx
        chosen_subset = list(mapping.keys())

    ### These are the weighted counts methods

    # The CSNAnalysis way
    elif method == "csnanalysis-weighted":

        # get the counts matrix according to this method
        countsmat = msn.edge_attribute_to_matrix(
            'weighted_counts',
            fill_value=0,
        )

        # for CSNAnalysis the matrix is the reverse orientation
        countsmat_T = countsmat.T

        # make a CSN
        csn = CSN(countsmat_T)

        # trim it to the main components
        csn.trim(
            min_count=method_kwargs['ergodic_trim_cutoff'],
        )

        # use the trimmed off idxs to get which states remain
        chosen_subset = csn.trim_indices

    # the MSMBuilder way
    elif method == "msmbuilder-weighted":

        # get the counts matrix according to this method
        countsmat = msn.edge_attribute_to_matrix(
            'weighted_counts',
            fill_value=0,
        )

        # this counts matrix is in the correct orientation for
        # msmbuilder method

        trimmed_countsmat, mapping, _ = \
                msmb_strongly_connected_subgraph(
                    countsmat,
                    weight=method_kwargs['ergodic_trim_cutoff'],
                    verbose=False,
                    )

        chosen_subset = list(mapping.keys())


    # ALERT: the pyemma way, this wasn't working for me
    elif method == "pyemma":

        raise NotImplementedError()

        subsets = pyemma_connected_sets(msn.countsmat)

        # get the index of largest one and the set itself
        subset_sizes = [len(subset) for subset in subsets]

        largest_subset_idx = np.argmax(subset_sizes)

        # choose the one to use
        chosen_subset = subsets[largest_subset_idx]

    # get the submatrix of this along with the mapping for the
    # countsmat trimmings
    trim_countsmat, trimming_mapping = subset_trim(countsmat, chosen_subset)

    # trimming_mapping is in node_idxs

    return trim_countsmat, trimming_mapping

# SNIPPET: the old way I was doing things to the CSNAnalysis
# def net_csn(msn, min_count=1):
#     """Convert a macrostate network to a CSN object."""

#     from csnanalysis.csn import CSN

#     # make a CSN
#     csn = CSN(msn.countsmat)

#     # trim it to the main components
#     csn.trim(min_count=min_count)

#     trimmed_group = [msn.node_idx_to_id(idx) for idx in
#                      set(range(csn.nnodes)).difference(csn.trim_indices)]

#     return csn, trimmed_group

## Transition Probabilities Method 1: Maximum Likelihood (MLE) PYEMMA
## version
def transprob_mle_pyemma(
        counts_mat,
        method='auto',
        reversible=True,
        eps_mu=1e-6,
):

    from msmtools.estimation \
        import transition_matrix \
        as pyemma_transmat_mle

    print("Running the pyemma transition matrix")
    # this uses the MLE method
    transprob_mat, stationary_distribution = \
        pyemma_transmat_mle(
            counts_mat,
            # this should be calculated by this method, so we set to None
            # here
            mu=None,

            # this is which algorithms to use based on the sparsity of the
            # matrix, which we just set to auto since we should get the
            # same answer either way and let it figure it out.
            method=method,

            # this has a bunch of options that go with it, but
            # this makes the matrix reversible, which is what we
            # want, rest of the options assume this is True
            reversible=reversible,

            # this is the epsilon value you want to have for finding
            # convergence
            eps_mu=eps_mu,

            # return the stationary distribution as well
            return_statdist=True,

        )

    return transprob_mat, stationary_distribution

## Transition Probabilities Method 2: Maximum Likelihood (MLE)
## MSMBuilder version
def transprob_mle_msmb(
        counts_mat,
        **kwargs,
):
    # MSMBuilder way

    from msmbuilder.msm._markovstatemodel \
        import _transmat_mle_prinz \
        as msmb_transmat_mle

    # no options here except the tolerance
    transprob_mat, stationary_distribution = \
        msmb_transmat_mle(counts_mat)


    return transprob_mat, stationary_distribution

## Transition Probabilities Method 3: Transpose method
def transprob_transpose(
        counts_mat,
        **kwargs,
):

    rev_counts = 0.5 * (counts_mat + counts_mat.T)

    populations = rev_counts.sum(axis=0)
    populations /= populations.sum(dtype=float)
    transmat = rev_counts.astype(float) / rev_counts.sum(axis=1)[:, None]

    return transmat, populations

## Transition Probabilities Method 4: Nonreversible Normalize
#
# This method simply normalizes the outgoing probabilities of each
# node so that they sum to 1 and are probabilities. Does not do
# symmetrization
def transprob_nonreversible_normalize(
        counts_mat,
        **kwargs,
):

    # the first axis is the sources, second axis is the targets

    # we want to normalize only the source->target transitions for
    # each source node

    # sum along the source axis
    out_total = counts_mat.sum(axis=0)

    # scale all the values down and renormalize such that the sum to 1
    # along the rows
    transmat = counts_mat / out_total

    # ALERT: not sure I'm doing this right, but I don't use it anyhow

    # compute the populations
    populations = out_total / out_total.sum()

    return transmat, populations

def compute_transprob(
        counts_mat,
        method='mle-pyemma',
        method_kwargs=None,
):

    if method_kwargs is None:
        method_kwargs = {}

    # dispatch on the method

    if method == 'nonreversible-normalize':

        transprob_mat, stationary_dist = transprob_nonreversible_normalize(
            counts_mat,
            **method_kwargs,
        )

    elif method == 'mle-pyemma':

        transprob_mat, stationary_dist = transprob_mle_pyemma(
            counts_mat,
            **method_kwargs,
        )

    elif method == 'mle-msmb':

        transprob_mat, stationary_dist = transprob_mle_msmb(
            counts_mat,
            **method_kwargs,
        )

    elif method == 'transpose':

        transprob_mat, stationary_dist = transprob_transpose(
            counts_mat,
            **method_kwargs,
        )


    else:
        raise ValueError(f"Unkown method {method}")

    return transprob_mat, stationary_dist

def make_pyemma_msm(
        msn,
        ergodic_trim_cutoff=1,
        trim_method='csnanalysis',
        trim_kwargs=None,
        transprob_method='transpose',
        transprob_kwargs=None,
):

    from pyemma.msm import MSM as MarkovStateModel

    # get the transition probability matrix of the strongly ergodic
    # component

    print(f"Doing Ergodic Trimming with method: {trim_method}")

    # Get the connected subset we are interested in.

    # trimming mapping is in terms of node_idx
    trim_countsmat, trimming_mapping = ergodic_trim(
        msn,
        method=trim_method,
        method_kwargs=trim_kwargs,
    )

    print("Making the Transition Probability Matrix")
    tprob_mat, populations = compute_transprob(
        trim_countsmat,
        method=transprob_method,
        method_kwargs=transprob_kwargs,
    )

    # for pyemma
    msm = MarkovStateModel(
        tprob_mat.T,
    )

    ## compose the state mapping for the MSN with the trimming mapping
    ## so we know that the idxs in the trimmed countsmat match the
    ## root state idxs
    dict_match = lambda dict1, dict2 : \
        {k: dict2.get(v) for k, v in dict1.items() if v in dict2}

    # the trimming_mapping is in terms of node_idx, this converts it to node_id
    mapping = dict_match(
        msn.node_id_to_idx_dict(),
        trimming_mapping,
    )

    # the mapping is: node_id -> trimmed_node_idx

    return msm, mapping

def add_msm_to_msn(
        msn,
        msm_id,
        pyemma_msm,
        trimming_mapping,
):
    """Adds basic information like the transition probability matrix
    values to the MSN model.

    Currently only the trimming information and the transition
    probabilities are stored.

    The trimming mapping is saved as two groups:

    - trimmed_nodes/msmid-{msm_id}/largest_component :: the nodes that are in the
      largest connected group, given by the 'trimming_mapping'
      dictionary.

    - trimmed_nodes/msmid-{msm_id}/not_largest_component :: the nodes that were left
      out of the largest group.

    The transition probability values are stored assymetrically for
    only those edges in the largest connected component. All edges not
    in this are set to np.nan. The edge attribute will be of the
    format: "transition_probability/msmid-{msm_id}"

    """

    # copy so that we don't screw up the original
    msn = deepcopy(msn)

    ## make the groups for the trimmings
    keep_node_ids = list(trimming_mapping.keys())
    trimmed_node_ids = set(msn.node_ids) - set(keep_node_ids)

    msn.set_node_group(
        f"trimmed_nodes/msmid-{msm_id}/largest_component",
        keep_node_ids,
    )

    msn.set_node_group(
        f"trimmed_nodes/msmid-{msm_id}/not_largest_component",
        trimmed_node_ids,
    )

    # ALERT: pay attention to this

    # get the transitions probability matrix, we have to transpose
    # since that is how they store it in PYEMMA
    transition_prob_mat = pyemma_msm.P.T

    prob_attr_key = f"transition_probability/msmid-{msm_id}"

    # iterate through every edge and see if there is a value for it
    for i, j in msn.graph.edges.keys():

        # if both nodes are in the trimming we will set the probability
        if ((i in trimming_mapping.keys()) and
            (j in trimming_mapping.keys())):

            # lookup the trimmed idxs
            trim_i = trimming_mapping[i]
            trim_j = trimming_mapping[j]

            prob = transition_prob_mat[trim_i, trim_j]

            # set it into the network
            msn.graph.edges[(i, j)][prob_attr_key] = prob

        # if it isn't assign nan
        else:
            msn.graph.edges[(i, j)][prob_attr_key] = np.nan

    return msn

## a little special error class
class MSMError(Exception):
    """Error for propagating bad GMRQ calculations."""

    pass


### Calculating GMRQ for training and test sets

# I've chopped it up into a number of different functions for clarity
# of understanding rather than performance. This really is fast
# compared to the clustering process so speed is not an issue here at
# all.

def map_eigenvectors(source_mapping,
                     source_eigvecs,
                     target_mapping,
):
    """Given a mapping of K states that has two sub-mappings N_s and N_t
    (N_s >= N_t) each of which has a set of eigenvectors we want to
    get the eigenvectors from N_s that are equivalent to the
    eigenvectors in N_t.

    Parameters
    ----------

    source_mapping : dict of int to int

    source_eigvecs : arraylike

    target_mapping : dict of int to int

    Returns
    -------

    mapped_eigvecs : arraylike

    """

    # def dict_compose
    dict_match = lambda dict1, dict2 : \
        {k: dict2.get(v) for k, v in dict1.items() if v in dict2}

    # make a mapping between the source states to the target states:
    # N_s : K match K : N_t --> N_s : N_t
    transform_mapping = dict_match(
        # reverse the source mapping so its N_s : K
        {v: k for k, v in source_mapping.items()},
        # this is K : N_t
        target_mapping
    )

    # get the source and value indices as two lists
    source_indices, dest_indices = zip(*transform_mapping.items())

    # make a new eigenvector array that is the size of the target
    mapped_eigvecs = np.zeros((
        len(target_mapping),
        source_eigvecs.shape[1]
    ))

    # if there are any eigenvectors in the source not in the target
    # than they will just be zero here

    # copy the source eigenvectors to the mapped eigenvectors
    mapped_eigvecs[dest_indices, :] = np.take(
        source_eigvecs,
        source_indices,
        axis=0,
    )

    return mapped_eigvecs


def gmrq_overlap_matrix(
        stationary_distribution,
):
    """Compute the overlap matrix (sometimes abbreviated as S) from the
    stationary probability distribution (typically abbreviated as pi)
    for a transition probability matrix.

    Parameters
    ----------

    stationary_distribution : numpy.array of shape (N_s,)
       The stationary probability distribution (pi) of the transition
       probability matrix (T) where N_s is the number of states in the
       state model.

    Returns
    -------

    overlap_matrix : np.array of shape (N_s, N_s)
       The overlap matrix (S) of the stationary probability distribution.

    """

    # S is the typical single letter abbreviation
    return np.diag(stationary_distribution)


def gmrq_correlation_matrix(
        overlap_matrix,
        transition_matrix,
):
    """Compute the correlation matrix for a transition matrix and the
    associated overlap matrix.

    N_s is the number of states in the MSM.

    Parameters
    ----------

    overlap_matrix : np.array of shape (N_s, N_s)
       The overlap matrix (S) of the stationary probability distribution.

    transition_matrix : np.array of shape (N_s, N_s)
       The transition probability matrix (T) of the MSM.

    Returns
    -------

    correlation_matrix : np.array of shape (N_s, N_s)
       The correlation matrix (C) of the stationary probability distribution.

    """

    # C is the common abbreviation here
    return overlap_matrix.dot(transition_matrix)


def gmrq_score_train_msm(train_msm):
    """Score the MSM based on the eigenvalues of its transition matrix."""

    return np.sum(train_msm.eigenvalues())


def gmrq_score_test_msm(
        train_msm,
        train_mapping,
        test_msm,
        test_mapping,
):
    """Score a test MSM given the training MSM eigenvectors.

    The test_msm must have the same state definitions as the train_msm.

    Parameters
    ----------

    train_msm : pyemma.MarkovStateModel

    train_mapping : dict of int to int

    test_msm : pyemma.MarkovStateModel

    test_mapping : dict of int to int

    Returns
    -------

    test_gmrq : float
        GMRQ score for the test set using the training model.

    """

    # calculate the eigenvectors for the training MSM which will be
    # used to score their performance on the testing MSM. This is A in
    # the McGibbon paper
    A = train_msm.eigenvectors_right()

    # we need to make sure we are comparing the correct eigenvectors
    # between the training and test MSMs because both can have
    # different sets of states which are included in their transition
    # probability matrices due to the need to do ergodic
    # trimming. That means that the root MSN will have K states as
    # specified from clustering but the training MSM has N_train and
    # the test MSM has N_test states which are both subsets of K but
    # are potentially disjoint sets between each other. To reconcile
    # this we have the mappings of indices from K to each N. We use
    # these to match states/eigenvectors from N_train to N_test. This
    # assumes that N_train >= N_test.

    # if the two mappings are identical we can just ignore this
    # mapping step
    if train_mapping != test_mapping:

        # otherwise go ahead and map them and replace the A
        # eigenvectors with the mapped one

        # ALERT: not sure this is correct here. Is it okay if N_train < N_test ?

        # this might fail if the training set is smaller than testing
        # set
        try:
            A = map_eigenvectors(
                train_mapping,
                A,
                test_mapping
            )

        # TODO: -inf or NaN here?

        # if it does fail we automatically score this as -inf, which
        # will result in an infinite score during optimization
        except MSMError:
            breakpoint()
            return np.nan

    # from the test set we compute the overlap and correlation
    # matrices which will be used to calculate the score on the
    # training eigenvectors
    S = gmrq_overlap_matrix(
        test_msm.pi
    )
    C = gmrq_correlation_matrix(
        S,
        test_msm.transition_matrix,
    )

    # compute the GMRQ score of the training eigenvectors on the
    # testing overlap and correlation matrices.

    # we handle linear algebra answers and return NaN if necessary.
    try:
        test_gmrq = np.trace(
            A.T.dot(
                C.dot(A)
            ).dot(
                np.linalg.inv(
                    A.T.dot(
                        S.dot(A)
                    )
                )
            )
        )
    except np.linalg.LinAlgError:
        test_gmrq = np.nan


    return test_gmrq

## Training stuff

# two layers to functions: one that applies a classfier to a
# particular split of samples, and one that runs all of the splits.
# another is for aggregating the split scores

def score_msm_gmrq_split_classifier(
        clf,
        train_rep_idxs,
        test_rep_idxs,
        gexp=None,
        lag_time=None,
        all_features=None,
        span_runs=None,
):
    """Score a classifier using GMRQ given a single splitting of training
    and sample replicates (i.e. a single cross-validation resampling).

    Parameters
    ----------

    clf : classifier following scikit learn API

    train_rep_idxs : list of int
       The simulation replicates in the training set.

    test_rep_idxs : list of int
       The simulation replicates in the testing set.


    gexp : str

    lag_time : int
        Lag time for computing transitions must be 2 or greater.

    all_features : dict of int to features
        A dictionary of all the clustering features by run. Its easier
        to read from memory once and then take from here since there
        is a lot of redundancy in the splits.

    span_runs : dict of int to list of int
        A mapping of a spanning contig (replicate) to the runs it
        contains.

    Returns
    -------

    test_gmrq : float
       The GMRQ score for the test dataset.

    """

    # make sure the "globals" are given
    assert gexp is not None
    assert lag_time is not None
    assert all_features is not None
    assert span_runs is not None

    # use the reps for each split to get the runs for each rep
    train_run_idxs = it.chain(*[span_runs[rep_idx] for rep_idx in train_rep_idxs])
    test_run_idxs = it.chain(*[span_runs[rep_idx] for rep_idx in test_rep_idxs])

    # ALERT sort them here. we need to sort and order the run idxs
    # in such a way that we can restructure them later after
    # desctructuring them here.

    # sort them to get a stable ordering
    train_run_idxs = sorted(train_run_idxs)
    test_run_idxs = sorted(test_run_idxs)

    # we grab the features for the runs here and inline
    # restructure them so we don't have to deal with the garbage
    # collector

    train_samples = unscalarize_features(
      reshape_run_to_features(
          [all_features[run_idx] for run_idx in train_run_idxs]
      )
    )

    test_samples = unscalarize_features(
      reshape_run_to_features(
          [all_features[run_idx] for run_idx in test_run_idxs]
      )
    )

    # perform the clustering on the training data
    clf.fit(train_samples)

    # get the labels from this and reshape into an assignments
    # shape. Restructure them to run_idx -> traj, cycles, 1
    train_labels = {
      run_idx : run_obs
      for run_idx, run_obs
      in zip(train_run_idxs,
             reshape_features_to_run(
                 clf.labels_,
                 gexp,
                 run_idxs=train_run_idxs,
             )
      )
    }

    # we can also get the labels for the training data
    test_labels = clf.predict(test_samples)

    # restructure them
    test_labels = {
      run_idx : run_obs
      for run_idx, run_obs
      in zip(test_run_idxs,
             reshape_features_to_run(
                 test_labels,
                 gexp,
                 run_idxs=test_run_idxs,
             )
      )
    }

    # then we use the clustering data to make a network

    # first get the contigtree which is a subset of runs we are using

    # Training

    # make a sub-contigtree of just the training runs/spans
    train_contigtree = get_contigtree(gexp,
                                    runs=train_run_idxs)

    # make a network with all the labels
    train_msn = MacroStateNetwork(train_contigtree,
                              transition_lag_time=lag_time,
                              assignments=train_labels
    )

    # Testing
    test_contigtree = get_contigtree(gexp,
                                  runs=test_run_idxs)

    test_msn = MacroStateNetwork(test_contigtree,
                              transition_lag_time=lag_time,
                              assignments=test_labels
    )

    # get the counts matrix for both the training data and the test
    # data
    train_countsmat = train_msn.countsmat

    test_countsmat = test_msn.countsmat

    ## Compute the transition probabilities and stationary
    ## distribution, there are different ways to do this, choose the
    ## appropriate one

    # the mappings are relative to their MSNs count matrix

    # make the MSM for the training data
    train_msm, train_mapping = make_msm(train_msn)

    # score the training msm for fun, but this isn't actually used
    train_gmrq = gmrq_score_train_msm(train_msm)

    test_msm, test_mapping = make_msm(test_msn)

    # then score the testing data from this model
    test_gmrq = gmrq_score_test_msm(
      train_msm,
      train_mapping,
      test_msm,
      test_mapping,
    )

#     # Report on results
#     msg = f"""
# Results for split
# -----------------

# Train   Test
# {train} {test}

# Total number of states: {train_msn.num_states}

# - Training MSM:
#   - num states :: {train_msm.nstates}
#   - GMRQ :: {train_gmrq}

# - Testing MSM:
#   - num states :: {test_msm.nstates}
#   - GMRQ :: {test_gmrq}
# """

    return test_gmrq

def score_msm_gmrq_crossval_classifier(
        gexp,
        observable_key,
        clf,
        splitter,
        lag_time=2,
):
    """Score a classifier using cross validation over the splitter using GMRQ."""

    # SNIPPET: these were some test parameters used

    ## Parameters
    # GEXP = '10'
    # LAG_TIME = 2
    # OBSERVABLE = 'lig_rmsd'
    # N_SPLITS = 3
    # splitter = KFold(n_splits=N_SPLITS)

    ## Code

    # split into training and testing based on the replicates (spans
    # in the contigtree)

    all_contigtree = get_contigtree(gexp)

    n_replicates = len(all_contigtree.span_traces)

    replicate_idxs = [i for i in range(n_replicates)]

    # split on the replicates using the splitter
    splits = list(splitter.split(replicate_idxs))

    # OPT: figure out where these should be loaded so as to reduce
    # loading from disk. Here it would would be once per
    # hyperparameter sample which is decent since it would get loaded
    # K times for each fold otherwise.

    # get all the labels for the runs and index them by the run in a
    # dict. Subsamples will draw from here
    all_features = {run_idx : run_obs
                  for run_idx, run_obs in
                  enumerate(get_observable(
                      observable_key,
                      gexp,
                      source='fs',))
    }

    # get a mapping of replicate/span to which runs it has
    span_runs = {rep_idx :
                 set((run_idx
                      for run_idx, cycle_idx
                      in span_trace))
                 for rep_idx, span_trace
                 in all_contigtree.span_traces.items()}

    split_scores = []
    # perform the process for each splitting
    for split_idx, split in enumerate(splits):

        train_rep_idxs = split[0]
        test_rep_idxs = split[1]


        test_gmrq = score_classifier_split(
            clf,
            train_rep_idxs,
            test_rep_idxs,
            gexp=gexp,
            lag_time=lag_time,
            all_features=all_features,
            span_runs=span_runs,
        )

        split_scores.append(test_gmrq)

    # aggregate the scores for each split in cross validation

    stats = {
        'mean' : np.mean(split_scores),
        'median' : np.median(split_scores),
        'min' : np.min(split_scores),
        'max' : np.max(split_scores),
        'std' : np.std(split_scores),
        'var' : np.var(split_scores),
        'scores' : split_scores,
    }

    # this is the final score for the classifier
    return stats

BASIN_SPECS = {}

## basin model for the topn bound and 0.7 unbound cutoff across all
## MSMs
_basin_top2_7 ={

    f'msm-{msm_id}_basin-bound-top2-unbound-0.7-cutoff' : {
        'msm_id' : f'{msm_id}',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.7 * tkunit.nanometer,
            },
        },
    }

    for msm_id in MSM_SPECS.keys()
}


#
_basin_top2_6 ={

    f'msm-{msm_id}_basin-bound-top2-unbound-0.6-cutoff' : {
        'msm_id' : f'{msm_id}',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.6 * tkunit.nanometer,
            },
        },
    }

    for msm_id in MSM_SPECS.keys()
}

#
_basin_top2_5 = {

    f'msm-{msm_id}_basin-bound-top2-unbound-0.5-cutoff' : {
        'msm_id' : f'{msm_id}',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.5 * tkunit.nanometer,
            },
        },
    }

    for msm_id in MSM_SPECS.keys()
}

bound_cutoffs = [
    0.1 * tkunit.nanometer,
    0.2 * tkunit.nanometer,
    0.3 * tkunit.nanometer,
]

unbound_cutoffs = [
    0.7 * tkunit.nanometer,
    0.6 * tkunit.nanometer,
    0.5 * tkunit.nanometer,
]

msm_ids = [
    'cutoff-1',
    'cutoff-2',
    'cutoff-5',
    'cutoff-10',
]

matrix_combos = [
    (cutoff, bound, unbound)
    for cutoff, bound, unbound
    in it.product(
        msm_ids,
        bound_cutoffs,
        unbound_cutoffs,
    )
]

#
_basin_matrix = {

    f'msm-{msm_id}_basin-bound-{bound.value_in_unit(bound.unit)}-unbound-{unbound.value_in_unit(unbound.unit)}-cutoff' : {
        'msm_id' : f"{msm_id}",
        'bound' : {
            'method' : 'rmsd',
            'method_kwargs' : {
                'cutoff' : bound,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : unbound,
            },
        },
    }

    for msm_id, bound, unbound in matrix_combos

}

# add all the specs into the main one
for basin_specs in [
        _basin_top2_7,
        _basin_top2_6,
        _basin_top2_5,
        _basin_matrix,
]:
    BASIN_SPECS.update(basin_specs)


def print_basin_spec_shell():

    print(' '.join([f"'{basin_id}'" for basin_id in BASIN_SPECS.keys()]))

OLD_BASIN_SPECS = {
    '0' : {
        'msm_id' : '0',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.7 * tkunit.nanometer,
            },
        },
    },

    '0_big_cutoff' : {
        'msm_id' : 'big_cutoff',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.7 * tkunit.nanometer,
            },
        },
    },


    '1' : {
        'msm_id' : '0',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'lig_dist',
            'method_kwargs' : {
                'cutoff' : 0.24 * tkunit.nanometer,
            },
        },
    },


    'choose-unbound-0' : {
        'msm_id' : '0',
        'bound' : {
            'method' : 'weights',
            'method_kwargs' : {
                'top_n' : 2,
            },
        },
        'unbound' : {
            'method' : 'chosen',
            'method_kwargs' : {
                'node_ids' : [691,],
            },
        },
    },

}

def idxs_in_msm(
        msn,
        msm_id,
        query_node_ids,
):
    """Function checks that all the node idxs for an MSM and determines
    whether they are in the largest connected component or not.

    Uses the node groups from the MacroStateNetwork of the format:

        'trimmed_nodes/msmid-{msm_id}/largest_component'

    Returns
    -------

    in_msm_idxs
        The indices in the MSM

    not_in_msm_idxs
        The indices not in the MSM largest connected component

    Raises
    ------

    MSMError : if no largest component is available to use.

    """

    node_group_key = f'trimmed_nodes/msmid-{msm_id}/largest_component'

    # the nodes that are in the largest component
    try:
        good_nodes = msn.node_groups[node_group_key]
    except KeyError:
        raise MSMError(f"No node group for largest component for the MSM model: {msm_id}.")

    in_nodes = [node_id
                for node_id in query_node_ids
                if node_id in good_nodes]

    out_nodes = list(set(query_node_ids) - set(in_nodes))

    return in_nodes, out_nodes

def compute_csn_bound_basin(
        basin_id,
        gexp,
        csn_id,
):

    # parse the basin specs
    method = BASIN_SPECS[basin_id]['bound']['method']
    basin_kwargs = BASIN_SPECS[basin_id]['bound']['method_kwargs']

    msm_id = BASIN_SPECS[basin_id]['msm_id']

    msn = get_msn(
        csn_id,
        gexp,
    )

    if method == 'weights':

        bound_basin_idxs = compute_csn_bound_basin_weights(
            msn,
            **basin_kwargs
        )

    elif method == 'rmsd':

        bound_basin_idxs = compute_csn_bound_basin_rmsd(
            msn,
            csn_id,
            gexp,
            **basin_kwargs,
        )


    # elif method == 'natives':

    #     bound_basin_idxs = compute_csn_bound_basin_natives(
    #         clustering_model,
    #         lig_id,
    #     )

    # elif method == 'feature_dist':

    #     bound_basin_idxs = compute_csn_bound_basin_feature_dist(
    #         clustering_model,
    #         lig_id,
    #         cutoff,
    #     )


    else:

        raise ValueError(f"method not recognized for basin_id: {basin_id}")

    # remove the idxs that are not in the largest connected component
    # for the MSM
    bound_basin_idxs, _ = idxs_in_msm(
        msn,
        msm_id,
        bound_basin_idxs,
    )

    return bound_basin_idxs

def compute_csn_unbound_basin(
        basin_id,
        gexp,
        csn_id,
):

    # parse the basin specs
    method = BASIN_SPECS[basin_id]['unbound']['method']
    basin_kwargs = BASIN_SPECS[basin_id]['unbound']['method_kwargs']

    msm_id = BASIN_SPECS[basin_id]['msm_id']

    msn = get_msn(
        csn_id,
        gexp,
    )

    if method == 'warp':

        unbound_basin_idxs = compute_csn_unbound_basin_warp(
            msn,
            gexp,
            **basin_kwargs,
        )

    elif method == 'lig_dist':

        unbound_basin_idxs = compute_csn_unbound_basin_lig_dist(
            gexp,
            csn_id,
            **basin_kwargs,
        )

    elif method == 'chosen':

        unbound_basin_idxs = BASIN_SPECS[basin_id]['unbound']['method_kwargs']['node_ids']

    # elif method == 'lig_dist':

    #     unbound_basin_idxs = compute_csn_unbound_basin_lig_dist(
    #         clustering_model,
    #         lig_id,
    #         cutoff,
    #     )

    else:

        raise ValueError(f"method not recognized for basin_id: {basin_id}")

    # remove the idxs that are not in the largest connected component
    # for the MSM
    unbound_basin_idxs, _ = idxs_in_msm(
        msn,
        msm_id,
        unbound_basin_idxs,
    )

    return unbound_basin_idxs


def basins_to_msn(
        basin_id,
        gexp,
        msn,
        bound_basin_idxs,
        unbound_basin_idxs,
):

    new_msn = deepcopy(msn)

    # remove old groups if they are in there
    if f'basin-{basin_id}/bound_basin' in new_msn._node_groups:
        del new_msn._node_groups[f'basin-{basin_id}/bound_basin']

    if f'basin-{basin_id}/unbound_basin' in new_msn._node_groups:
        del new_msn._node_groups[f'basin-{basin_id}/unbound_basin']

    # save the things in the groups
    new_msn.set_node_group(f'basin-{basin_id}/bound_basin', bound_basin_idxs)
    new_msn.set_node_group(f'basin-{basin_id}/unbound_basin', unbound_basin_idxs)

    return new_msn

def compute_csn_bound_basin_weights(
        msn,
        top_n=1,
):
    """This function ranks the weights of all the nodes in the network and
    gives the top N of them as the bound basin.
    """

    # just give the biggest weighted nodes according to the state
    # weight ranking
    return state_weight_ranking(msn)[-top_n:]


def compute_csn_bound_basin_natives(
        msn,
        clf_id=None,
        gexp=None,
):
    """This function will simply get the native state as the bound
    basin.

    This may not be appropriate for docked or otherwise modelled
    structures.

    """

    assert clf_id is not None
    assert gexp is not None

    native_state = classify_native_state_cluster(
        clf_id,
        gexp,
    )

    return [native_state,]




@jlmem.cache
def compute_csn_bound_basin_rmsd(
        msn,
        csn_id,
        gexp,
        cutoff=None,
        ):
    """This returns all cluster centers that are within a certain cutoff
    to a certain state as the bound basin.
    """

    from geomm.rmsd import calc_rmsd

    assert cutoff is not None

    lig_id = dict(GEXP_LIG_IDS)[gexp]
    clf_id = CSN_SPECS[csn_id]['clf_id']

    sel_idxs = lig_selection_idxs(lig_id)

    lig_idxs = sel_idxs['main_rep/ligand']

    # coordinates are in nanometers, convert cutoff to this
    cutoff_unitless = cutoff.value_in_unit(tkunit.nanometer)

    # get the features from the model
    clf = get_clustering_model(
        clf_id,
        gexp,
    )

    # compute RMSDs of the ligand to the native state as an ad hoc
    # selection criterion

    # get the cluster center structures. Use the common denomitor
    # representation (which is 'main_rep') here for compatibility,
    # this isn't used for visualization.
    cluster_center_traj = get_cluster_center_traj(
        gexp,
        csn_id,
        alt_rep='main_rep',
    )

    # use the highest weight state for the bound state
    bound_state_id = state_weight_ranking(msn)[-1]

    dists_to_bound = []
    for center_idx in range(len(clf.cluster_centers_)):

        rmsd = calc_rmsd(cluster_center_traj[bound_state_id].xyz[0][lig_idxs],
                         cluster_center_traj[center_idx].xyz[0][lig_idxs])

        dists_to_bound.append(rmsd)

    matched_centers = [l[0] for l in
                       np.argwhere(
                           np.array(dists_to_bound) < cutoff_unitless
                       )
                       ]

    # add the centers that are within the cutoff along with the bound state
    bound_basin_ids = np.array(matched_centers)

    return bound_basin_ids

# TODO: I don't want to fix this right now, I won't be using it.
@jlmem.cache
def compute_csn_bound_basin_feature_dist(clustering_model, lig_id, cutoff):

    raise NotImplementedError

    import numpy as np
    import scipy.spatial.distance as spdist

    # use lag time of 2 since it should be cached and isn't actually used.
    net = lig_state_network(lig_id, clustering_model, 2)

    # get the features from the model
    clf = get_clustering_model(clustering_model, lig_id)

    # use the highest weight state for the bound state
    bound_state_id = state_weight_ranking(net)[-1]
    other_idxs = list(set(range(len(clf.cluster_centers_))) -
                      set([bound_state_id]))


    bound_state_center = clf.cluster_centers_[bound_state_id]
    other_centers = clf.cluster_centers_[other_idxs]

    dists_to_bound = (spdist.cdist([bound_state_center], other_centers, 'canberra') /\
                      bound_state_center.shape[0])[0]

    matched_centers = [l[0] for l in np.argwhere(dists_to_bound < cutoff)]

    # add the centers that are within the cutoff along with the bound state
    bound_basin_ids = np.array([bound_state_id] + matched_centers)

    return bound_basin_ids

@jlmem.cache
def compute_csn_unbound_basin_warp(
        msn,
        gexp,
):
    """Report the nodes which had any of the warping events for the gexp
    simulations."""

    contigtree = get_contigtree(gexp)

    # get the trace of the frames which were warped
    trace = contigtree.warp_trace()

    warped_nodes = []
    # then find which cluster these were in
    for node_id in msn.graph.nodes:

        assignments = msn.get_node_attribute(node_id, 'assignments')

        # check if any of the trace is in the node assignments and if
        # they are record this node index and remove them from the
        # list since they won't be there in the future
        matches = [rec for rec in trace if rec in assignments]

        if len(matches) > 0:
            warped_nodes.append(node_id)
            for match in matches:
                trace.remove(match)

    return warped_nodes

@jlmem.cache
def compute_csn_unbound_basin_lig_dist(
        gexp,
        csn_id,
        cutoff=None,
):
    """This gets the walkers which are beyond a certain distance much like
    the boundary conditions themselves but only takes into account the
    cluster centers.

    This is good if you didn't get any walkers unbound.

    Cutoff will be converted to nanometer

    """

    # compute the distances of the ligand to the ligand to the protein
    # for the cluster centers and then take the minimum and find if it
    # is above the cutoff

    assert cutoff is not None
    assert tkunit.is_quantity(cutoff)

    # convert cutoff to nanometer
    cutoff_nm = cutoff.value_in_unit(tkunit.nanometer)

    # get the cluster center structures
    cluster_center_traj = get_cluster_center_traj(
        gexp,
        csn_id,
        alt_rep='main_rep',
    )

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    # get the ligand and protein idxs
    sel_idxs = lig_selection_idxs(lig_id)

    lig_idxs = sel_idxs['main_rep/ligand']
    prot_idxs = sel_idxs['main_rep/protein']

    atom_pairs = np.array(list(it.product(lig_idxs, prot_idxs)))

    cluster_pair_distances = mdj.compute_distances(cluster_center_traj, atom_pairs)

    unbound_basin_ids = [node_id for node_id, node_dists
                         in enumerate(cluster_pair_distances)
                         if node_dists.min() > cutoff_nm]

    return unbound_basin_ids

committor_methods = [
    'msmb',
    # 'pyemma',
    # 'csnanalysis',
]

TS_SPECS = {}
for method in committor_methods:

    for basin_id in BASIN_SPECS.keys():

        TS_SPECS[f'basin-{basin_id}_committor-{method}'] = {
            'committor_method' : method,
            'basin_id' : f'{basin_id}',
            'upper' : 0.6,
            'lower' : 0.4,
        }

def print_ts_spec_shell():

    print(' '.join([f"'{ts_id}'" for ts_id in TS_SPECS.keys()]))

def forward_committor_probabilities(
        transition_probability_matrix,
        source_basin,
        sink_basin,
        tolerance=1e-6,
        maxstep=20,
):

    raise NotImplementedError

    # make the sinks on the matrix by setting these columns to 0 & 1
    # for the source & sink respectively
    sink_prob_mat = copy(transition_probability_matrix)

    # for each idx in the sink basin set the column to identity vector
    for sink_node_idx in sink_basin:

        # set column to 0.0
        sink_prob_mat[:,sink_prob_mat] = 0.0

        # set the diagonal to 1.0
        sink_prob_mat[sink_node_idx, sink_node_idx] = 1.0


        pass


    return forward_committor_probs

def committors_csnanalysis(
        countsmat,
        bound_basin,
        unbound_basin,
):

    from csnanalysis.csn import CSN

    csn = CSN(countsmat)

    committors = csn.calc_committors(
        [bound_basin, unbound_basin],
        labels=['bound', 'unbound']
    )

    forward_probs = committors[:,0]
    backward_probs = committors[:,1]

    return forward_probs, backward_probs

def committors_pyemma(
        msm,
        bound_basin,
        unbound_basin,
):

    # this is a pyemma MSM which has a method
    forward_probs = msm.committor_forward(
        bound_basin,
        unbound_basin,
    )

    backward_probs = msm.committor_backward(
        bound_basin,
        unbound_basin,
    )

    return forward_probs, backward_probs

def committors_msmb(
        transition_probability_mat,
        bound_basin,
        unbound_basin,
):
    from msmbuilder.tpt.committor import \
        _committors as msmb_forward_committors

    print("running msmbuilder function")

    forward_probs = msmb_forward_committors(
        bound_basin,
        unbound_basin,
        transition_probability_mat,
    )

    return forward_probs

def calc_msm_committors(
        msm_id,
        gexp,
        bound_basin,
        unbound_basin,
        method='pyemma',
):
    """

    Parameters
    ----------

    bound_basin : list of node_id
        The node_id

"""

    if method == 'pyemma':

        pyemma_msm, trimming_mapping = load_msm(
            msm_id,
            gexp,
        )

        unbound_committors, bound_committors = committors_pyemma(
            pyemma_msm,
            bound_basin,
            unbound_basin,
        )

    elif method == 'msmb':

        print("Using MSMBuilder method")

        # this trimming mapping is a mapping of node_id -> trimmed_node_idx
        pyemma_msm, trimming_mapping = load_msm(
            msm_id,
            gexp,
        )

        # this matrix is trimmed
        transition_probability_mat = pyemma_msm.P

        # the bound and unbound basin are in terms of node_ids, so we
        # need to convert these to the trimmed node_idxs of the
        # trimmed matrix

        trim_bound_basin = [trimming_mapping[node_id]
                           for node_id in bound_basin]

        trim_unbound_basin = [trimming_mapping[node_id]
                             for node_id in unbound_basin]

        print("Running wrapper function")
        # here the basins are in node_idx and are in terms of the
        # transition probability matrix which is trimmed.
        unbound_committors = committors_msmb(
            transition_probability_mat,
            trim_bound_basin,
            trim_unbound_basin,
        )

        # just invert the unbound committors for the bound ones
        bound_committors = 1 - unbound_committors

        # NOTE: committors are in terms of node_idx

    elif method == 'csnanalysis':

        msn = get_msn(
            msm_id,
            gexp,
        )

        pyemma_msm, trimming_mapping = load_msm(
            msm_id,
            gexp,
        )

        subset_idxs = list(trimming_mapping.keys())

        trim_countsmat = np.zeros(
            (len(subset_idxs), len(subset_idxs),),
            dtype=msn.countsmat.dtype,
        )

        for full_idx, trimmed_idx in trimming_mapping.items():

            trim_countsmat[trimmed_idx, :] = msn.countsmat[full_idx, subset_idxs]
            trim_countsmat[:, trimmed_idx] = msn.countsmat[subset_idxs, full_idx]

        # calculate the committor probabilities
        unbound_committors, bound_committors = committors_csnanalysis(
            trim_countsmat,
            bound_basin,
            unbound_basin,
        )

    else:
        raise ValueError(f"Unkown method {method}")

    # in terms of node_idx
    return unbound_committors, bound_committors

def predict_ts(
        msn,
        gexp,
        ts_id,
        committor_probs,
        trimming_mapping,
):

    ts_upper_bound = TS_SPECS[ts_id]['upper']
    ts_lower_bound = TS_SPECS[ts_id]['lower']

    # get nodes which are in the transition state, and convert to ids

    # boolean selection
    ts_sel_idxs = np.argwhere(
        (committor_probs >= ts_lower_bound) &
        (committor_probs <= ts_upper_bound)
    ).flatten()

    # ALERT: the committor probs are of the trimmed indices, so the
    # trimming mapping must be used to unconvert
    trimming_mapping_rev = {trim_idx : node_id
                            for node_id, trim_idx
                            in trimming_mapping.items()
                            }

    # convert these idxs to the node ids
    ts_node_ids = [
        trimming_mapping_rev[idx]
        for idx in ts_sel_idxs
    ]

    return ts_node_ids

def msn_committors_ts(
        msn,
        gexp,
        ts_id,
):

    basin_id = TS_SPECS[ts_id]['basin_id']
    committor_method = TS_SPECS[ts_id]['committor_method']

    msm_id = BASIN_SPECS[basin_id]['msm_id']

    pyemma_msm, trimming_mapping = load_msm(
        msm_id,
        gexp,
    )

    # get the trimmed nodes and largest component nodes
    trimmed_node_ids = msn.node_groups[f'trimmed_nodes/msmid-{msm_id}/not_largest_component']
    untrimmed_node_ids = msn.node_groups[f'trimmed_nodes/msmid-{msm_id}/largest_component']

    # get the basins from the MSN; these are in node_id
    bound_basin = msn.node_groups[f'basin-{basin_id}/bound_basin']
    unbound_basin = msn.node_groups[f'basin-{basin_id}/unbound_basin']

    # calculate the committors

    # the unbound and bound committors are in terms of the node_idx
    unbound_committors, bound_committors = calc_msm_committors(
        msm_id,
        gexp,
        bound_basin,
        unbound_basin,
        method=committor_method,
    )

    # get the values for all the nodes as a dictionary
    # first initialize them with the trimmed nodes as nans
    bound_node_values = {
        node_id : np.nan
        for node_id
        in trimmed_node_ids
    }

    unbound_node_values = {
        node_id : np.nan
        for node_id
        in trimmed_node_ids
    }

    # then get the actual committors

    # convert the node_idx (given by the natural indexing)

    for node_id in untrimmed_node_ids:

        # then get the trimmed node id of them
        trim_node_idx = trimming_mapping[node_id]

        bound_node_values[node_id] = bound_committors[trim_node_idx]
        unbound_node_values[node_id] = unbound_committors[trim_node_idx]


    # unbound committors should be in terms of the node_idx.
    # the returned ts_node_ids are in terms of node_id
    ts_node_ids = predict_ts(
        msn,
        gexp,
        ts_id,
        unbound_committors,
        trimming_mapping,
    )

    # return as forward, backward, ts
    return unbound_node_values, bound_node_values, ts_node_ids



def add_ts_to_msn(
        msn,
        ts_id,
        forward_node_probs,
        backward_node_probs,
        ts_node_ids,
):
    """Add data to the MSN for a transition state model.

    This will store a node group for the transition state model in a
    group of the format:

        'tsid-{ts_id}/TS'

    The forward/backwards committor probabilities will be saved as node
    attributes of the form:

        'tsid-{ts_id}/forward_committors'
        'tsid-{ts_id}/backwards_committors'

    We also make a convenience node attribute for coloring the network
    called the 'UBT'. This field is string valued and is one of:

        U : Unbound Basin
        B : Bound Basin
        T : Transisition State Ensemble Prediction
        N : Not in any of the other ones, but still in the MSM
        X : Not ergodically connected to the MSM

    This will be saved in a node attribute of the form:

        'tsid-{ts_id}/UBT'

    """

    new_msn = deepcopy(msn)

    # things we will be setting
    ts_node_group_key = f"tsid-{ts_id}/TS"
    forward_attr_key = f'tsid-{ts_id}/forward_committors'
    backward_attr_key = f'tsid-{ts_id}/backwards_committors'
    UBT_attr_key = f'tsid-{ts_id}/UBT'

    # get the basin id
    basin_id = TS_SPECS[ts_id]['basin_id']

    msm_id = BASIN_SPECS[basin_id]['msm_id']

    # set the group for the transition state nodes
    new_msn.set_node_group(
        ts_node_group_key,
        ts_node_ids,
    )

    # set the node attributes for committor probabilities

    # TODO: remove this. Should be done at the site of creation not
    # after passing it in here
    # # we have to pad it for the ones that don't have committors with
    # # nans.
    # missing_probs_node_ids = list(set(list(new_msn.graph.nodes.keys())) -
    #                           set(forward_node_probs.keys()))
    # for node_id in missing_probs_node_ids:
    #     forward_node_probs[node_id] = np.nan
    #     backward_node_probs[node_id] = np.nan

    new_msn.set_nodes_observable(
        forward_attr_key,
        forward_node_probs,
    )
    new_msn.set_nodes_observable(
        backward_attr_key,
        backward_node_probs,
    )

    # then we construct and set the UBT node attributes
    UBT_node_values = {}
    for node_id in new_msn.graph.nodes:

        if node_id in new_msn.node_groups[f'basin-{basin_id}/unbound_basin']:
            UBT_node_values[node_id] = 'U'

        elif node_id in new_msn.node_groups[f'basin-{basin_id}/bound_basin']:
            UBT_node_values[node_id] = 'B'

        # HACK: we can use the TS node group since we set it above
        elif node_id in new_msn.node_groups[ts_node_group_key]:
            UBT_node_values[node_id] = 'T'

        elif node_id in new_msn.node_groups[f'trimmed_nodes/msmid-{msm_id}/not_largest_component']:
            UBT_node_values[node_id] = 'X'

        else:
            UBT_node_values[node_id] = 'N'


    new_msn.set_nodes_observable(
        UBT_attr_key,
        UBT_node_values,
    )


    return new_msn

def save_gexp_basin_eval_plots(
        gexp,
        case_0_key,
        case_1_key,
        **figs,
):

    for fig_key, fig in figs.items():

        save_gexp_fig(
            gexp,
            f"basin-eval/{fig_key}/case-0-{case_0_key}_case-1-{case_1_key}",
            fig,
        )

@jlmem.cache
def compute_legacy_all_lig_coms(
        node_ids,
        gexp,
):

    from geomm.centroid import centroid

    clust_assgs = get_legacy_cluster_canonical_traj_assignments()
    ts_microstates_trace = list(it.chain(*[
        clust_assgs[clust_id]
        for clust_id in node_ids
    ]))

    # get the HDF5
    wepy_h5 = get_legacy_h5()

    # get the microstates in batches by chunk and slice
    # recenter, slice homology, and get coms
    CHUNK_SIZE = 500

    chunks = [
        ts_microstates_trace[i:i + CHUNK_SIZE]
        for i in range(
                0,
                len(ts_microstates_trace),
                CHUNK_SIZE,
        )
        ]

    com_chunks = []
    for i, chunk_trace in enumerate(chunks):

        print(f"Doing chunk {i}")

        with wepy_h5:

            traj_fields = wepy_h5.get_trace_fields(
                chunk_trace,
                [
                    'positions',
                    'box_vectors',
                ],
            )

        chunk_coms = compute_traj_fields_hom_coms(
            gexp,
            traj_fields,
        )

        com_chunks.append(chunk_coms)


    coms = np.concatenate(com_chunks)


    return coms

@jlmem.cache
def compute_msn_all_lig_coms(
        gexp,
        csn_id,
        node_ids,
        homology_idxs,
):
    """Compute the centers of mass (coms) of the homologous ligand indices for
    a MacroStateNetwork.

    This computes from scratch rather than loading them.


    Returns
    -------

    coms : array of shape (N, 3)
       Where N is the number of all microstates across the node_ids.

    """

    from geomm.centroid import centroid

    # load the MSN with the H5
    print("Getting MSN H5")

    msn_h5 = get_h5_msn(
        csn_id,
        gexp,
    )

    # get the microstates from the TSE states, recenter and
    # superimpose them, then slice out only the homology indices
    # of the ligand from this
    with msn_h5:

        # we do this one node at a time since this will blow memory if we don't
        all_chunks = []
        for node_id in node_ids:

            print(f"Getting HOM COMS for node {node_id}")

            traj_fields = msn_h5.states_to_traj_fields([node_id])

            chunk_coms = compute_traj_fields_hom_coms(
                gexp,
                traj_fields,
            )

            all_chunks.append(chunk_coms)

    # combine them
    coms = np.concatenate(all_chunks)

    return coms

def compute_tspca_model(
        coms,
        test_size=0.25,
):
    """Compute the PC model, centers of mass, and the score for the PCA
    model"""

    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split

    # split into training and test set
    train, test = train_test_split(coms, test_size=0.25)

    model = PCA(n_components=3)

    print(f"Fitting PCA model on {len(train)} samples")

    # train and get the mode projections
    train_projections = model.fit_transform(train)

    # keep a record of the individual mode scores
    mode_scores = []


    # loop over the modes and how much variance they explain
    for c_idx, c_var in enumerate(model.explained_variance_):

        # get the total explained variance of the mode
        var_ratio = model.explained_variance_ratio_[c_idx]

        mode_score = {
            'percent_explained_variance' : 100 * var_ratio,
            'percent_total_explained_variance' : 100 * c_var,
        }

        mode_scores.append(mode_score)

    print("Scoring the training set")

    # explain the model on the test set
    score = model.score(test)


    return model, score, mode_scores

def compute_gexps_tspca_model(
        gexps,
        csn_id,
        ts_id,
        test_size=0.25,
):

    from geomm.centroid import centroid

    print("Getting the COMs for each gexp")

    # for each gexp get the COMs for the TSs
    all_ts_coms = []
    for gexp in gexps:

        print(gexp)

        sel_idxs = lig_selection_idxs(gexp)

        # UGLY: treat the legacy one differently
        if gexp in LEGACY_GEXPS:

            # get the microstate trace for the TSE nodes
            ts_cluster_idxs = get_legacy_ts_cluster_idxs()

            gexp_ts_coms = compute_legacy_all_lig_coms(
                ts_cluster_idxs,
                gexp,
            )

        # the modern gexps get handled this way
        else:

            # get the node_ids for the ts model
            msn_h5 = get_msn(
                csn_id,
                gexp,
            )

            tse_node_ids = msn_h5.node_groups[f"tsid-{ts_id}/TS"]

            print("Computing the COMs for the TSE")

            # calculate the COMs of the homology atoms in the ligand for all
            # of the microstates in the TS nodes
            gexp_ts_coms = compute_msn_all_lig_coms(
                gexp,
                csn_id,
                tse_node_ids,
            )

        all_ts_coms.append(gexp_ts_coms)

    # concatenate the arrays of coms for each gexp
    all_ts_coms = np.concatenate(all_ts_coms)

    print("Computing the TSPCA model")

    pca_model, model_score, mode_scores = compute_tspca_model(
        all_ts_coms,
        test_size=test_size,
    )

    return pca_model, model_score, mode_scores

def ts_pca_model_selection(
        gexps,
        csn_id,
        ts_id,
        test_size=0.25,
):

    # then generate the test matrix of subsamples
    num_samples = len(gexps)
    subsamples = []
    for subsample_size in range(1, num_samples + 1):

        combs = it.combinations(gexps, subsample_size)

        subsamples.extend(combs)

    # then we will make the tspca_id for each of them
    tspca_id_tmpl = "tsid-{ts_id}_gexps-{gexp_list}_testsize-{test_size}"

    tspca_ids = []
    for subsample in subsamples:

        gexp_list = '-'.join([str(gexp) for gexp in subsample])

        tspca_id = tspca_id_tmpl.format(
            ts_id=ts_id,
            gexp_list=gexp_list,
            test_size=test_size,
        )

        tspca_ids.append(tspca_id)

    # then we go through and make/retrieve the models for the tspca_ids
    model_results = {}
    for tspca_id, gexps in zip(tspca_ids, subsamples):

        # ALERT: don't do manual "caching" just let joblib take care of it

        print("Getting the PCA model for tspca_id:")
        print(tspca_id)

        # the 'coms' are the centers of mass
        model, model_score, mode_scores = compute_gexps_tspca_model(
            gexps,
            csn_id,
            ts_id,
            test_size=test_size,
        )

        model_results[tspca_id] = {
            'model' : model,
            'gexps' : gexps,
            'test_size' : test_size,
            'model_score' : model_score,
            'mode_scores' : mode_scores,
        }


    return model_results

def plot_tspca_model_score(
        model_results,
):

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1,1)

    scores = [res['model_score'] for res in model_results.values()]

    axes.boxplot(scores)

    return fig, axes

def make_tspca_model_score_table(
        model_results,
):

    from collections import defaultdict

    df_d = defaultdict(list)
    for tspca_id, result_d in model_results.items():

        df_d['tspca_id'].append(
            tspca_id,
        )

        df_d['model_score'].append(
            result_d['model_score'],
        )

        df_d['gexps'].append(
            ','.join([gexp for gexp in result_d['gexps']]),
        )

        df_d['test_size'].append(
            result_d['test_size'],
        )

        for i, mode_score_d in enumerate(result_d['mode_scores']):

            df_d[f'mode_{i}_per-explained-variance'].append(
                mode_score_d['percent_explained_variance'],
            )

            df_d[f'mode_{i}_per-total-explained-variance'].append(
                mode_score_d['percent_total_explained_variance'],
            )

    model_score_df = pd.DataFrame(df_d)

    return model_score_df

def observable_traj_fields_hom_com_projections(
        model,
        gexp,
        traj_fields
):
    """Function for computing the projections of a TSPCA model over the
    coms of the ligand homology atoms.

    """

    coms = compute_traj_fields_hom_coms(
        gexp,
        traj_fields,
    )

    pc_projections = compute_traj_hom_com_projections(
        model,
        coms,
    )

    return pc_projections

def compute_traj_fields_hom_coms(
        gexp,
        traj_fields,
):
    """Returns the centers of mass (COMS) for trajectory fields for a
    given gexp."""

    from geomm.centroid import centroid

    if gexp in LEGACY_GEXPS:
        rep = 'real_rep'
    else:
        rep = 'main_rep'

    # get the ligand homology atoms for the rep chosen
    sel_idxs = lig_selection_idxs(gexp)
    homology_idxs = sel_idxs[f'{rep}/ligand/homology']

    hom_positions = recenter_superimpose_traj(
        traj_fields,
        gexp,
        rep,
    )[0][:,homology_idxs,:]

    gc.collect()

    coms = np.array([centroid(frame) for frame in hom_positions])

    return coms

def compute_traj_com_projections(
        model,
        coms,
):
    """For a model and some centers of mass make a trajectory of the
    centers of mass using a phony topology and then calculate the
    projection values of these points for each PC.

    Parameters
    ----------

    model : scikit learn PCA model

    coms : array of shape (N, 3)
        Center of mass positions.


    Returns
    -------

    traj : mdtraj.Trajectory
        A trajectory of the COMs for the ligands.

    pc_projections : list of array of float
        Each element of the list is for each of the modes (in order of
        the model). The values of each array correspond to the frames
        of the trajectory and are scalar values which are the
        projected value of that frame onto the PC.

    """

    # first we get their projections which is just easier to re-transform
    # all of them
    com_projections = model.transform(coms)

    pc_projections = []
    # make a pdb traj for each component
    for pc_idx in range(model.n_components_):

        # reshape the projection value for this PC, which can be used
        # to color things
        pc_projection = com_projections[:, pc_idx].reshape((1, com_projections.shape[0]))

        pc_projections.append(pc_projections)


    return pc_projections


def make_coms_traj(coms):
    """Make an mdtraj trajectory object from single particle center of mass positions."""

    # we get a topology for the coms as atoms in a single topology
    top = n_atom_mdj_top(coms.shape[0])

    # then make a trajectory using the coms xyzs, after we reshape to trajectory style
    coms_frames = np.reshape(
        coms,
        (1, coms.shape[0], 3)
    )

    traj = mdj.Trajectory(coms_frames, top)

    return traj



def compute_gexp_tspca_projections_csn(
        msn,
        tse_node_ids,
        contigtree,
        homology_idxs,
        lig_id,
        pca_model,
):

    ## Part 2. Getting projections onto TS cluster centers
    #---------------------------------------------------------------------------
    # get the actual cluster centers only so we can project them

    # get the trace of all of the node centers
    node_center_trace_d = msn.get_nodes_attribute('center_idx')


    # First get only the TS centers

    # reshaping this into two lists ordered the same

    # zip the recs and their nodes
    ts_center_recs = [
        (node_id, rec)
        for node_id, rec
        in node_center_trace_d.items()
        if node_id in tse_node_ids
    ]

    # then unzip so we can reference them separately
    ts_center_trace_node_ids = [node_id for node_id, rec in ts_center_recs]
    ts_center_trace = [rec for node_id, rec in ts_center_recs]

    # get the trajectory fields, this is 'main_rep'
    with contigtree.wepy_h5 as wepy_h5:

        ts_center_traj_fields = wepy_h5.get_trace_fields(
            ts_center_trace,
            ['positions', 'box_vectors'],
        )

    ts_center_hom_positions = recenter_superimpose_traj(
        ts_center_traj_fields,
        lig_id,
        'main_rep'
    )[0][:,homology_idxs,:]

    # get their ligand COMs
    ts_center_coms = np.array([
        centroid(frame)
        for frame
        in ts_center_hom_positions
    ])

    # project them onto the modes
    ts_center_projections = model.transform(ts_center_coms)


    # then write out a trajectory for all the COMs as single atoms

    # first we get their projections which is just easier to re-transform
    # all of them
    ts_center_com_projections = model.transform(ts_center_coms)

    # we get a topology for the coms as atoms in a single topology
    ts_center_top = n_atom_mdj_top(ts_center_com_projections.shape[0])

    # then make a trajectory using the coms xyzs, after we reshape to trajectory style
    ts_center_coms_frames = np.reshape(ts_center_coms, (1, ts_center_coms.shape[0], 3))
    traj = mdj.Trajectory(ts_center_coms_frames, ts_center_top)

    # make a pdb traj for each component
    for pc_idx in range(model.n_components_):

        # add in the extra field for each mode so we can visualize, in pdb
        # mode so we can have the bfactors
        traj.save_pdb('data/ts_pca/ts_center_coms_pc_{}_lig-{}.pdb'.format(pc_idx, LIG_ID),
                      bfactors=ts_center_projections[:, pc_idx].reshape((1, ts_center_projections.shape[0])))



    # Option B: All of the centers

    # zip the recs and their nodes
    center_recs = [(node_id, rec) for node_id, rec in node_center_recs.items()]

    # then unzip so we can reference them separately
    center_trace_node_ids = [node_id for node_id, rec in center_recs]
    center_trace = [rec for node_id, rec in center_recs]

    # get the trajectory fields
    with contigtree.wepy_h5 as wepy_h5:
        center_traj_fields = wepy_h5.get_trace_fields(center_trace,
                                                      ['positions', 'box_vectors'])

    # TODO: audit main_rep here
    center_hom_positions = recenter_superimpose_traj(center_traj_fields,
                                                     LIG_ID, 'main_rep')[0][:,hom_idxs,:]
    # get their ligand COMs
    center_coms = np.array([centroid(frame) for frame in center_hom_positions])

    # project them onto the modes
    center_projections = model.transform(center_coms)


    return

def junk():
    ## Params

    ## Part 3a: Render onto network
    #---------------------------------------------------------------------------

    # now set these as observables in the network
    # pcs_obs = {}
    for pc_idx in range(model.n_components_):
        pc_obs = {}
        for idx, node_id in enumerate(center_trace_node_ids):

            # # if it is not a ts node then just set it as a nan
            # if node_id not in ts_center_trace_node_ids:
            #     pc_obs[node_id] = 0. #np.nan

            # # get the projections of the cluster centers from the transition state
            # else:
            #     ts_idx = ts_center_trace_node_ids.index(node_id)
            #     pc_obs[node_id] = ts_center_projections[ts_idx, pc_idx]

            pc_obs[node_id] = center_projections[idx, pc_idx]


        # pcs_obs[pc_idx] = pc_obs

        # TODO use the key for the TS model, i'm not doing it here because
        # it is too busy for this data exploration task

        # set as an observable in the network
        net.set_nodes_observable('ts_PC-{}'.format(pc_idx), pc_obs)


    # save the network with this visualization
    gephi_graph = load_gephi_graph(MODEL, LIG_ID, LAG_TIME, name_pattern='main')
    net = update_network_from_gephi(net, gephi_graph, layout_name='main')

    save_gexf(MODEL, LIG_ID, LAG_TIME, net)
