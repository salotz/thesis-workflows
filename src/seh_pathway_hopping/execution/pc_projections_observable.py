def observable_pc_projections(client, lig_id, **kwargs):

    import functools

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    from seh_pathway_hopping._tasks import (
        gexp_wepy_h5_path,
        lig_selection_idxs,
        get_model,
    )

    from seh_pathway_hopping._tasks import pc_projections_observable

    TS_MODEL = 'TS/warp+featdist-0.07'

    wepy_h5_path = gexp_wepy_h5_path(lig_id)

    # TODO: audit main_rep here
    sel_idxs = lig_selection_idxs(lig_id)
    hom_idxs = sel_idxs['main_rep/ligand/homology']

    model = get_model(TS_MODEL, lig_id)

    obs_func = functools.partial(pc_projections_observable, lig_id, hom_idxs, model)


    with client:

        observable = compute_observable(obs_func, wepy_h5_path, client,
                                        ['positions', 'box_vectors'])

    # split the observable vector into 3 different observables
    pc_observables = []
    for pc_idx in range(model.n_components_):
        pc_observables.append([])

        for run_idx, run in enumerate(observable):
            pc_observables[pc_idx][run_idx].append([])

            for traj in run:

                pc_observables[pc_idx][run_idx].append(traj[:,pc_idx])





    # save the observable in the hdf5
    with WepyHDF5(wepy_h5_path, mode='r+') as wepy_h5:
        for pc_idx, pc_obs in enumerate(pc_observables):

            model_obs_name = TS_MODEL.replace('/', '-')
            obs_name = "{}_pc-{}".format(model_obs_name, pc_idx)

            # save the object for the observable
            save_observable(obs_name, lig_id, pc_obs)

            # and save to the HDF5
            wepy_h5.add_observable(obs_name, pc_obs)


if __name__ == "__main__":

    from seh_pathway_hopping.execution import execute

    execute(observable_pc_projections)
