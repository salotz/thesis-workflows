import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        save_all_csn_stuff,
        get_h5_msn,
        calc_msn_obs_stats,
    )

    kwargs = parse_kwargs(specs)

    csn_id = kwargs['csn_id']
    gexp = kwargs['gexp']
    layout_id = kwargs['layout_id']

    observable_name = kwargs['obs_name']

    if 'tag' in kwargs:
        tag = kwargs['tag']
    else:
        tag = None

    ### load the MSN with the ContigTree backing
    msn = get_h5_msn(
        csn_id,
        gexp,
        tag=tag,
    )

    # then compute the stats for the observable adding it to the
    # network
    _ = calc_msn_obs_stats(
        msn,
        observable_name,
        save=True,
    )

    # we want to only save the base MSN so make sure to only do that part

    # save it as the current set of network files that are used for
    # analysis and visualization etc.
    save_all_csn_stuff(
        csn_id,
        gexp,
        msn.base_network,
        tag=tag,
        layout_id=layout_id,
        overwrite=True,
    )

if __name__ == "__main__":

    cli()
