import click

@click.command()
@click.argument('gexp')
def save_rate_plots(gexp):

    from seh_pathway_hopping._tasks import (
        save_gexp_plot_rates,
    )

    save_gexp_plot_rates(gexp)

if __name__ == "__main__":

    save_rate_plots()
