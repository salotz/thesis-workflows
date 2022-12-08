def observable_lig_bs_atom_pair_dist(map_func, **kwargs):

    from functools import partial
    import os.path as osp

    from seh_pathway_hopping._tasks import (
        get_gexp_wepy_h5,
        lig_selection_idxs,
        lig_selection_tops,
        lig_prot_atom_pair_observable,
        lig_bs_atom_pairs,
        GEXP_LIG_IDS,
    )


    FIELD = 'bs_lig_pair_dists'
    FIELDS = ['positions', 'box_vectors']
    RUN_IDXS = Ellipsis

    gexp = kwargs['gexp']

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    wepy_h5 = get_gexp_wepy_h5(gexp)

    atom_pairs = lig_bs_atom_pairs(
        gexp,
        rep_key='main_rep',
    )

    # do a partial evaluation of the function so we can use
    # the multiprocessing map which doesn't accept multiple
    # iterables for arguments
    func = partial(
        lig_prot_atom_pair_observable, # the function
        atom_pairs,
        lig_selection_tops(lig_id)['main_rep'], # topology
    )

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

import click

from seh_pathway_hopping.execution import parse_kwargs

from multiprocessing.pool import Pool as MPPool
from ray.util.multiprocessing import Pool as RayPool

func = observable_lig_bs_atom_pair_dist

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
