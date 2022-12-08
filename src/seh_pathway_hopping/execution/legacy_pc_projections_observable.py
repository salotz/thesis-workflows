def observable_pc_projections(client, lig_id, **kwargs):

    import functools
    import os.path as osp

    import joblib

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    from seh_pathway_hopping._tasks import (
        gexp_wepy_h5_path,
        lig_selection_idxs,
        get_model,
        legacy_linker_file_path,
        save_observable,
        data_path,
    )

    from seh_pathway_hopping._tasks import legacy_pc_projections_observable

    TS_MODEL = 'TS/warp+featdist-0.07'
    LEGACY_LIG_ID = 17

    wepy_h5_path = legacy_linker_file_path()

    sel_idxs = lig_selection_idxs(LEGACY_LIG_ID)
    hom_idxs = sel_idxs['real_rep/ligand/homology']

    model = get_model(TS_MODEL, lig_id)

    obs_func = functools.partial(legacy_pc_projections_observable, LEGACY_LIG_ID,
                                 hom_idxs, model)


    # DEBUG: for debugging
    with WepyHDF5(wepy_h5_path, mode='r') as wepy_h5:
        observable = wepy_h5.compute_observable(obs_func, ['positions', 'box_vectors'],
                                                (),
        )

    # DEBUG: failsafe so I don't lose work
    failsafe_path = osp.join(data_path(), 'tmp/legacy_projections_observable_failsafe.jl.pkl')
    joblib.dump(observable, failsafe_path)
    print("failsafe worked")


    # with client:

    #     observable = compute_observable(obs_func, wepy_h5_path, client,
    #                                     ['positions', 'box_vectors'])

    # split the observable vector into 3 different observables
    pc_observables = []
    for pc_idx in range(model.n_components_):
        pc_observables.append([])

        for traj_idx, traj in enumerate(observable):

            pc_observables[pc_idx].append(traj[:,pc_idx])

    for pc_idx, pc_obs in enumerate(pc_observables):

        model_obs_name = TS_MODEL.replace('/', '-')
        obs_name = "lig-{}_{}_pc-{}".format(lig_id, model_obs_name, pc_idx)

        # save the object for the observable
        save_observable(obs_name, LEGACY_LIG_ID, pc_obs)


if __name__ == "__main__":

    from seh_pathway_hopping.execution import execute

    execute(observable_pc_projections)
