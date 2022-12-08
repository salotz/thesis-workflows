import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        do_gexp_classification,
    )

    kwargs = parse_kwargs(specs)

    do_gexp_classification(
        kwargs['gexp'],
        kwargs['observable'],
        kwargs['clf_id'],
    )

if __name__ == "__main__":

    cli()
