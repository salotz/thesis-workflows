def observable_box_volume(client, **kwargs):

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    from seh_pathway_hopping._tasks import gexp_wepy_h5_path

    from seh_pathway_hopping._tasks import box_volume_observable

    wepy_h5_path = gexp_wepy_h5_path(kwargs['gexp'])


    with client:

        observable = compute_observable(
            box_volume_observable,
            str(wepy_h5_path),
            client,
            ['box_vectors']
        )

    # save the observable in the hdf5
    with WepyHDF5(wepy_h5_path, mode='r+') as wepy_h5:
        wepy_h5.add_observable('box_volume', observable)


if __name__ == "__main__":

    from seh_pathway_hopping.execution import execute

    execute(observable_box_volume)
