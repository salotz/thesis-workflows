from functools import partial

import click
from seh_pathway_hopping.execution import parse_kwargs

from multiprocessing.pool import Pool as MPPool

from wepy.hdf5 import WepyHDF5

from seh_pathway_hopping._tasks import get_gexp_wepy_h5

from seh_pathway_hopping._tasks import (
    get_gexp_wepy_h5,
    load_tspca_model,
    observable_traj_fields_hom_com_projections,
)


def obs_func(map_func, **kwargs):

    gexp = kwargs['gexp']
    tspca_id = kwargs['tspca_id']

    FIELD = f"tspca_projection/{tspca_id}"
    FIELDS = ['positions', 'box_vectors']
    RUN_IDXS = Ellipsis


    wepy_h5 = get_gexp_wepy_h5(
        gexp,
        mode='r+',
    )

    # load the tspca model
    tspca_model, _, _, _, _ = load_tspca_model(tspca_id)


    chunk_func = partial(
        compute_traj_fields_hom_com_projections,
        tspca_model,
        gexp,
    )

    with wepy_h5:

        observable = wepy_h5.compute_observable(
            func,
            FIELDS,
            (),
            map_func=map_func,
            save_to_hdf5=FIELD,
        )

    print("finished calculation")


@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.option('--n-cores', '-n', type=int, default=1)
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        load_tspca_model,
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    tspca_id = kwargs['tspca_id']
    gexp = kwargs['gexp']


    # load the contigtree for the gexp
    contigtree = get_contigtree(gexp)

    pool = MPPool(n_cores)

    CHUNKSIZE = 500
    def map_func(chunk_func, *args):

        return pool.imap(chunk_func,
                         *args,
                         CHUNKSIZE,
                         )

    obs_func(map_func, **kwargs)


if __name__ == "__main__":

    cli()
