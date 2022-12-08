import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        MSM_SPECS,
        get_msn,
        save_all_csn_stuff,
        save_msm,
        make_pyemma_msm,
        add_msm_to_msn,
    )


    kwargs = parse_kwargs(specs)

    msm_id = kwargs['msm_id']
    gexp = kwargs['gexp']
    layout_id = kwargs['layout_id']

    if 'tag' in kwargs:
        tag = kwargs['tag']
    else:
        tag = None

    msm_spec = MSM_SPECS[msm_id]
    csn_id = msm_spec['csn_id']

    print("MSM Spec:", msm_spec)

    trim_method = msm_spec['trim_method']
    trim_kwargs = msm_spec['trim_kwargs']

    transprob_method = msm_spec['transition_prob_method']
    transprob_kwargs = msm_spec['transition_prob_kwargs']

    print("trim method: ", trim_method)
    print("trim kwargs: ", trim_kwargs)

    print("transprob method: ", transprob_method)
    print("transprob kwargs: ", transprob_kwargs)

    print("Loading the MSN")
    ## Load
    msn = get_msn(
        csn_id,
        gexp,
        tag=tag,
    )

    print("Making the MSM")
    ## Transformation
    pyemma_msm, trimming_mapping = make_pyemma_msm(
        msn,
        trim_method=trim_method,
        trim_kwargs=trim_kwargs,
        transprob_method=transprob_method,
        transprob_kwargs=transprob_kwargs,
    )

    print(f"Trimmed network has {len(trimming_mapping)} nodes")

    msm_msn = add_msm_to_msn(
        msn,
        msm_id,
        pyemma_msm,
        trimming_mapping,
    )

    ## Save the network and MSM

    print("Saving the MSM")
    save_msm(
        msm_id,
        gexp,
        pyemma_msm,
        trimming_mapping,
    )

    print("Saving the CSN stuff")
    save_all_csn_stuff(
        csn_id,
        gexp,
        msm_msn,
        tag=tag,
        layout_id=layout_id,
        overwrite=True,
    )

if __name__ == "__main__":

    cli()
