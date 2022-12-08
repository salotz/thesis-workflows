import click

@click.command()
@click.argument('gexp')
def show_agg_prob_plots(gexp):

    from seh_pathway_hopping._tasks import (
        gexp_show_plot_agg_prob,
    )

    gexp_show_plot_agg_prob(gexp)

if __name__ == "__main__":

    show_agg_prob_plots()
