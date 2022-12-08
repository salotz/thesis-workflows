import click

@click.command()
@click.argument('gexp')
def gen_warp_lineages_multi(gexp):

    from seh_pathway_hopping._tasks import (
        save_gexp_warp_lineages_dcds,
    )

    save_gexp_warp_lineages_dcds(gexp)

if __name__ == "__main__":

    gen_warp_lineages_multi()
