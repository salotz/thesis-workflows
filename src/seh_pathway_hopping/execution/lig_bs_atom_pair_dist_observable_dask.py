def observable_lib_bs_atom_pair_dist(client, **kwargs):

    from functools import partial
    import os.path as osp

    import joblib
    import numpy as np

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    from seh_pathway_hopping._tasks import (
        data_path,
        gexp_wepy_h5_path,
        lig_selection_idxs,
        lig_selection_tops,
        lig_prot_atom_pair_observable,
        lig_prot_atom_pairs,
    )


    FIELD = 'bs_lig_pair_dists'
    FIELDS = ['positions', 'box_vectors']
    RUN_IDXS = Ellipsis
    CHUNK_SIZE = 100

    lig_id = kwargs['gexp']

    wepy_h5_path = gexp_wepy_h5_path(kwargs['gexp'])
    # TODO: audit main_rep here
    ligand_idxs = lig_selection_idxs(lig_id)['main_rep/ligand']

    # do a partial evaluation of the function so we can use
    # the multiprocessing map which doesn't accept multiple
    # iterables for arguments
    func = partial(lig_prot_atom_pair_observable,
                   lig_prot_atom_pairs(lig_id),
                   lig_selection_tops(lig_id)['main_rep'])

    print("starting calculation")
    observable = compute_observable(func,
                                    wepy_h5_path,
                                    client,
                                    ['positions', 'box_vectors'],
                                    chunk_size=100,
                                    run_idxs=Ellipsis)

    print("finished calculation")

    with WepyHDF5(wepy_h5_path, mode='r+') as wepy_h5:

        wepy_h5.add_observable(FIELD, observable)


if __name__ == '__main__':
    from seh_pathway_hopping.execution import execute

    execute(observable_lib_bs_atom_pair_dist)
