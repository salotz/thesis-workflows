import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        get_msn,
        load_msm,
        save_all_csn_stuff,
        CSN_SPECS,
        MSM_SPECS,
        BASIN_SPECS,
        TS_SPECS,
        msn_committors_ts,
        add_ts_to_msn,
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    # msm_id = kwargs['msm_id']
    gexp = kwargs['gexp']
    csn_id = kwargs['csn_id']
    layout_id = kwargs['layout_id']

    # if we have a single basin id requested only do that
    if 'ts_id' in kwargs:
        ts_ids = [kwargs['ts_id']]
    else:
        ts_ids = list(TS_SPECS.keys())

    clf_id = CSN_SPECS[csn_id]['clf_id']

    # load the msn
    msn = get_msn(
        csn_id,
        gexp,
    )

    # do it for each ts spec
    for ts_id in ts_ids:

        basin_id = TS_SPECS[ts_id]['basin_id']

        msm_id = BASIN_SPECS[basin_id]['msm_id']

        ## Load

        msm = load_msm(
            msm_id,
            gexp,
        )

        ## Transformations

        # calculate committors and make TS prediction
        forward_probs, backwards_probs, ts_node_ids = msn_committors_ts(
            msn,
            gexp,
            ts_id,
        )

        print("Adding data to the MSN")
        # save data to the msn
        msn = add_ts_to_msn(
            msn,
            ts_id,
            forward_probs,
            backwards_probs,
            ts_node_ids,
        )

        print("Saving All CSN stuff")
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
