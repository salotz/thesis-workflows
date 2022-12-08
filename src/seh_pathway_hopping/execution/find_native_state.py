import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        get_msn,
        save_all_csn_stuff,
        MSM_SPECS,
        classify_native_state_cluster,
    )

    kwargs = parse_kwargs(specs)

    ## Load
    msn = get_msn(
        kwargs['msm_id'],
        kwargs['gexp'],
        csn_id=kwargs['csn_id'],
    )

    clf_id = MSM_SPECS[kwargs['msm_id']]['clf_id']

    ## Transformations
    native_node_id = classify_native_state_cluster(
        clf_id,
        kwargs['gexp'],
    )

    print(f"native node", native_node_id)

    msn.set_node_group('native_state_id', [native_node_id])

    ## Save it
    save_all_csn_stuff(
        kwargs['msm_id'],
        kwargs['gexp'],
        msn,
        csn_id=kwargs['csn_id'],
        layout_id=kwargs['layout_id'],
        overwrite=True,
    )

if __name__ == "__main__":

    cli()
