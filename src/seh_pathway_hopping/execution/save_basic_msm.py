import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        make_basic_csn,
        save_all_csn_stuff,
    )

    kwargs = parse_kwargs(specs)

    gexp = kwargs['gexp']
    csn_id = kwargs['csn_id']

    if 'tag' in kwargs:
        tag = kwargs['tag']
    else:
        tag = None

    if 'layout_id' in kwargs:
        layout_id = kwargs['layout_id']
    else:
        layout_id = None


    # make a basic CSN from the specs
    msn = make_basic_csn(
        gexp,
        csn_id,
    )

    # save it as the current set of network files that are used for
    # analysis and visualization etc.
    save_all_csn_stuff(
        csn_id,
        gexp,
        msn,
        tag=tag,
        layout_id=layout_id,
        overwrite=True,
    )

if __name__ == "__main__":

    cli()
