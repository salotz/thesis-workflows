import click

@click.command()
@click.argument('gexp')
def show_rt_plots(gexp):

    from seh_pathway_hopping._tasks import (
        gexp_show_plot_rts,
    )

    gexp_show_plot_rts(gexp)

if __name__ == "__main__":

    show_rt_plots()
