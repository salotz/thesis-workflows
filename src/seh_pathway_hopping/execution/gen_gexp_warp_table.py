import click

@click.command()
@click.argument('gexp')
def gen_gexp_warp_table(gexp):


    from seh_pathway_hopping._tasks import (
        save_gexp_warp_table,
    )

    save_gexp_warp_table(gexp, overwrite=True)


if __name__ == "__main__":

    gen_gexp_warp_table()
