N_SPHERE_POINTS = 960

def observable_lig_sasa(map_func, **kwargs):

    from functools import partial

    from wepy.hdf5 import WepyHDF5
    from wepy.analysis.distributed import compute_observable

    # the task or flow
    from seh_pathway_hopping._tasks import (
        GEXP_LIG_IDS,
        get_gexp_wepy_h5,
        lig_sasa_observable,
        lig_selection_idxs,
        lig_selection_tops,
    )

    gexp = kwargs['gexp']

    # params for the sasa func
    if 'n_sphere_points'in kwargs:
        n_sphere_points = int(kwargs['n_sphere_points'])

    else:
        n_sphere_points = N_SPHERE_POINTS

    FIELD = f'lig-sasa_npoints-{n_sphere_points}'
    FIELDS = ['positions', 'box_vectors']

    # RUN_IDXS = [0]
    RUN_IDXS = Ellipsis


    wepy_h5 = get_gexp_wepy_h5(gexp,
                                mode='r+')

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    # generate the function we need to make this run
    # TODO: audit main_rep here
    # REVD: main rep is good here since it has to do data for all of the runs.
    ligand_idxs = lig_selection_idxs(lig_id)['main_rep/ligand']
    top = lig_selection_tops(lig_id)['main_rep']

    # do a partial evaluation of the function so we can use
    # the multiprocessing map which doesn't accept multiple
    # iterables for arguments
    func = partial(
        lig_sasa_observable,
        lig_id,
        ligand_idxs,
        top,
        n_sphere_points,
    )

    print(f"starting calculation for gexp: {gexp}")

    with wepy_h5:
        observable = wepy_h5.compute_observable(
            func,
            FIELDS,
            (),
            map_func=map_func,
            save_to_hdf5=FIELD,
         )

    print("finished calculation")

# ----------------------------------------
import time

import click

from seh_pathway_hopping.execution import parse_kwargs

from multiprocessing.pool import Pool as MPPool
from ray.util.multiprocessing import Pool as RayPool

func = observable_lig_sasa

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

        CHUNKSIZE = 1
        def map_func(func, *args):
            return pool.imap(func,
                             *args,
                             CHUNKSIZE,
            )

    elif mapper is 'ray':
        pool = RayPool(n_cores)
        map_func = pool.map

    start = time.time()

    func(map_func, **kwargs)

    end = time.time()

    duration = end - start

    print(f"Took a total of: {duration} s")



if __name__ == "__main__":

    cli()
