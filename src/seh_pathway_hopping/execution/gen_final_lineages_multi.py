import click

@click.command()
@click.argument('gexp')
def gen_final_lineages_multi(gexp):

    from seh_pathway_hopping._tasks import (
        save_gexp_final_lineages_dcds,
    )

    save_gexp_final_lineages_dcds(gexp)

if __name__ == "__main__":

    gen_final_lineages_multi()
