def observable_lig_bs_rmsd(map_func, **kwargs):

    from functools import partial

    from wepy.hdf5 import WepyHDF5

    from seh_pathway_hopping._tasks import get_gexp_wepy_h5

    # the task or flow
    from seh_pathway_hopping._tasks import (
        lig_bs_rmsd_observable,
        lig_selection_idxs,
    )

    FIELD = 'lig_rmsd'
    FIELDS = ['positions', 'box_vectors']
    RUN_IDXS = Ellipsis

    wepy_h5 = get_gexp_wepy_h5(kwargs['gexp'])

    # generate the function we need to make this run

    # TODO: audit main_rep here
    ligand_idxs = lig_selection_idxs(kwargs['gexp'])['main_rep/ligand']

    # do a partial evaluation of the function so we can use
    # the multiprocessing map which doesn't accept multiple
    # iterables for arguments
    func = partial(lig_bs_rmsd_observable, kwargs['gexp'], ligand_idxs)

    print("starting calculation")
    wepy_h5.set_mode('r+')
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
from functools import partial

import click

from seh_pathway_hopping.execution import parse_kwargs

from multiprocessing.pool import Pool as MPPool
from ray.util.multiprocessing import Pool as RayPool

func = observable_lig_bs_rmsd

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

        CHUNKSIZE = 10
        def map_func(func, *args):
            return pool.imap(func,
                             *args,
                             CHUNKSIZE,
            )

    elif mapper is 'ray':
        pool = RayPool(n_cores)
        map_func = pool.map

    func(map_func, **kwargs)



if __name__ == "__main__":

    cli()
