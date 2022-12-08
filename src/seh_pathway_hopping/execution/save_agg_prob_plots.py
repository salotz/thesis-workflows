import click

@click.command()
@click.argument('gexp')
def save_agg_plots(gexp):

    from seh_pathway_hopping._tasks import (
        save_gexp_plot_agg_prob,
    )

    save_gexp_plot_agg_prob(gexp)

if __name__ == "__main__":

    save_agg_plots()
