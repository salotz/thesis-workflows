import click

@click.command()
@click.argument('gexp')
def gen_span_warp_table(gexp):


    from seh_pathway_hopping._tasks import (
        save_span_warp_table,
        get_gexp_span_ids,
    )

    for span_id in get_gexp_span_ids(gexp):

        save_span_warp_table(gexp, span_id, overwrite=True)


if __name__ == "__main__":

    gen_span_warp_table()
