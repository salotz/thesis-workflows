import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        get_msn,
        save_all_csn_stuff,
        MSM_SPECS,
        BASIN_SPECS,
        compute_csn_bound_basin,
        compute_csn_unbound_basin,
        basins_to_msn,
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    # msm_id = kwargs['msm_id']
    gexp = kwargs['gexp']
    csn_id = kwargs['csn_id']
    layout_id = kwargs['layout_id']

    # if there was a basin_id only compute for that otherwise do all
    # of the basin specs
    if 'basin_id' in kwargs:
        basin_ids = [kwargs['basin_id']]

        print("Getting the basins for basin_id:")
        print(basin_ids[0])

    for basin_id in basin_ids:

        msm_id = BASIN_SPECS[basin_id]['msm_id']

        ## Load
        msn = get_msn(
            csn_id,
            gexp,
        )

        ## Transformations

        # get the bound basin according to the 'basin_id' specs
        bound_basin_idxs = compute_csn_bound_basin(
            basin_id,
            gexp,
            csn_id,
        )

        unbound_basin_idxs = compute_csn_unbound_basin(
            basin_id,
            gexp,
            csn_id,
        )

        print(f"Basins for basin_id: {basin_id}")
        print(f"Found {len(bound_basin_idxs)} nodes for bound basin")
        print(" ".join([str(i) for i in bound_basin_idxs]))
        print(f"Found {len(unbound_basin_idxs)} nodes for unbound basin")
        print(" ".join([str(i) for i in unbound_basin_idxs]))

        msn = basins_to_msn(
            basin_id,
            gexp,
            msn,
            bound_basin_idxs,
            unbound_basin_idxs,
        )


        ## Save it
        save_all_csn_stuff(
            csn_id,
            gexp,
            msn,
            layout_id=layout_id,
            overwrite=True,
        )

if __name__ == "__main__":

    cli()
