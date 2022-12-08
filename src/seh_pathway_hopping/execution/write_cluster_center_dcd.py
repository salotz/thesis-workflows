import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        save_cluster_centers,
        get_cluster_center_traj,
    )


    kwargs = parse_kwargs(specs)

    msm_id = kwargs['msm_id']
    gexp = kwargs['gexp']

    centers_traj = get_cluster_center_traj(
        gexp,
        csn_id,
    )

    save_cluster_centers(
        msm_id,
        gexp,
        centers_traj,
    )

if __name__ == "__main__":

    cli()
