import click

@click.command()
@click.argument('gexp')
def gen_contig_stats_table(gexp):


    from seh_pathway_hopping._tasks import (
        save_span_stats_table,
        span_stats_table_str,
    )

    save_span_stats_table(gexp, overwrite=True)

    click.echo(span_stats_table_str(gexp))


if __name__ == "__main__":

    gen_contig_stats_table()
