import click

@click.command()
@click.argument('gexp')
def save_rt_plots(gexp):

    from seh_pathway_hopping._tasks import (
        save_gexp_plot_rts,
    )

    save_gexp_plot_rts(gexp)

if __name__ == "__main__":

    save_rt_plots()
