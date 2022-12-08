import click

@click.command()
@click.option('--top-n', '-n', type=int, default=5)
@click.argument('gexp')
def gen_high_progress_walker_lineages(
        top_n,
        gexp,
):

    from seh_pathway_hopping._tasks import (
        save_gexp_high_progress_walkers_lineages_dcds,
    )
    save_gexp_high_progress_walkers_lineages_dcds(
        gexp,
        top_n,
        progress_key='min_distances',
    )

if __name__ == "__main__":

    gen_high_progress_walker_lineages()
