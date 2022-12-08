import click

@click.command()
@click.argument('gexp')
@click.argument('observable')
def main(gexp, observable):

    from seh_pathway_hopping._tasks import (
        gexp_save_plot_spans_fe_obs,
    )

    gexp_save_plot_spans_fe_obs(gexp, observable)

if __name__ == '__main__':
    main()
