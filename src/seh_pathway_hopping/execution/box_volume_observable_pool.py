def observable_box_volume(map_func, **kwargs):

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    from seh_pathway_hopping._tasks import get_gexp_wepy_h5

    from seh_pathway_hopping._tasks import box_volume_observable

    wepy_h5 = get_gexp_wepy_h5(kwargs['gexp'])

    with wepy_h5:

        observable = wepy_h5.compute_observable(
            box_volume_observable,
            ['box_vectors'],
            (),
            map_func=map_func,
        )

    # save the observable in the hdf5
    # with WepyHDF5(wepy_h5_path, mode='r+') as wepy_h5:
    #     wepy_h5.add_observable('box_volume', observable)

#----------------------------------------
import click

from seh_pathway_hopping.execution import parse_kwargs

from multiprocessing.pool import Pool as MPPool
from ray.util.multiprocessing import Pool as RayPool

func = observable_box_volume

@click.command(context_settings=dict(
ignore_unknown_options=True,))
@click.option('--n-cores', '-n', type=int, default=1)
@click.option('--mapper', '-m',
              type=click.Choice(['mp', 'ray']),
              default='mp')
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(n_cores, mapper, specs):

    kwargs = parse_kwargs(specs)

    if mapper is 'mp':
        pool = MPPool(n_cores)
    elif mapper is 'ray':
        pool = RayPool(n_cores)

    func(pool.map, **kwargs)



if __name__ == "__main__":

    cli()
