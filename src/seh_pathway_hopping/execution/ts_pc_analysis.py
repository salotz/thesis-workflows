# size of test set for scoring the PCA model
TEST_SIZE = 0.25

import click
from seh_pathway_hopping.execution import parse_kwargs

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        get_msn,
        save_all_csn_stuff,
        CSN_SPECS,
        MSM_SPECS,
        BASIN_SPECS,
        TS_SPECS,
        compute_gexp_ts_pca,
        compute_traj_com_projections,
        save_tspca_model,
        save_com_trajs,
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    # msm_id = kwargs['msm_id']
    gexp = kwargs['gexp']
    csn_id = kwargs['csn_id']
    layout_id = kwargs['layout_id']
    ts_id = kwargs['ts_id']

    # load the msn
    msn = get_msn(
        csn_id,
        gexp,
    )

    ## Analysis

    # Make the PCA model and get scores

    pca_model, ts_coms, mode_scores, model_score = compute_gexp_ts_pca(
        gexp,
        csn_id,
        ts_id,
        test_size=TEST_SIZE,
    )

    # make data for 3D visualization of the COMs

    # then we want to visualize the Transition State coms, make the trajectory and PC
    # projection values for each mode
    ts_com_traj, ts_com_pc_projections = compute_traj_com_projections(
        pca_model,
        ts_coms,
    )

    # TODO
    # apply the projections values for the whole network
    # compute_gexp_tspca_projections_csn(
    #     msn_h5,
    #     pca_model,
    # )


    ## Do all the saving at once so its more likely to be consistent

    # Save the PCA model and associated data
    # _ = save_tspca_model(
    #     ts_id,
    #     gexp,
    #     pca_model,
    #     ts_coms,
    #     mode_scores,
    #     model_scores,
    # )

    # # save the TS COMs as a trajectory with the projection values as the color
    # _ = save_com_trajs(
    #     ts_id,
    #     gexp,
    #     ts_com_traj,
    #     pc_projections,
    # )

    # # TODO
    # # add the data of the PC mode projections to the network




    # # Plot the 1D free energy


    # # Plot the 2D free energy


    # print("Saving All CSN stuff")
    # ## Save it
    # save_all_csn_stuff(
    #     csn_id,
    #     gexp,
    #     msn,
    #     layout_id=layout_id,
    #     overwrite=True,
    # )

if __name__ == "__main__":

    cli()
