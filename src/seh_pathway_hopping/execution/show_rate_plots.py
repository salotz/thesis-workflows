import click

@click.command()
@click.argument('gexp')
def show_rate_plots(gexp):

    from seh_pathway_hopping._tasks import (
        gexp_show_plot_rates,
    )

    gexp_show_plot_rates(gexp)

if __name__ == "__main__":

    show_rate_plots()
