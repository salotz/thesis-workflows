import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        load_gephi_graph,
        update_network_layout_from_gephi,
        save_all_csn_stuff,
        get_msn,
    )

    kwargs = parse_kwargs(specs)

    csn_id = kwargs['csn_id']
    gexp = kwargs['gexp']
    layout_id = kwargs['layout_id']

    if 'tag' in kwargs:
        tag = kwargs['tag']
    else:
        tag = None

    ## Basically we load the gexf file and then dump the results to
    ## the rest of the un-updated files

    # get the gephi gexf networkx graph
    gephi_nx = load_gephi_graph(
        csn_id,
        gexp,
        tag=tag,
        layout_id=layout_id,
    )

    # load the msn to update
    msn = get_msn(
        csn_id,
        gexp,
        tag=tag,
    )

    # update the msn with the gexf data
    update_msn = update_network_layout_from_gephi(
        msn,
        gephi_nx,
        layout_name=layout_id,
    )

    # save it as the current set of network files that are used for
    # analysis and visualization etc.
    save_all_csn_stuff(
        csn_id,
        gexp,
        update_msn,
        tag=tag,
        layout_id=layout_id,
        overwrite=True,
    )

if __name__ == "__main__":

    cli()
