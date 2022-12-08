import click
from seh_pathway_hopping.execution import parse_kwargs

import matplotlib.pyplot as plt

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    # msm_id = kwargs['msm_id']
    csn_id = kwargs['csn_id']
    tspca_id = kwargs['tspca_id']
    gexp = kwargs['gexp']
    layout_id = kwargs['layout_id']

    ## Analysis

if __name__ == "__main__":

    cli()
