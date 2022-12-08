import click

@click.command()
@click.option('--save', is_flag=True)
@click.option('--show', is_flag=True)
@click.argument('observable')
def main(save, show, observable):

    from seh_pathway_hopping._tasks import (
        all_render_plot_agg_fe_obs,
    )

    all_render_plot_agg_fe_obs(
        observable,
        save=save,
        show=show,
    )

if __name__ == '__main__':
    main()
