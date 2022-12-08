import click


def match_run_idx(wepy_h5, start_hash, end_hash):
    """Get the run idx of the run in the HDF5 using the start and end hash"""

    with wepy_h5:

        for run_idx in wepy_h5.run_idxs:

            if (start_hash == wepy_h5.run_start_snapshot_hash(run_idx) and
                end_hash == wepy_h5.run_end_snapshot_hash(run_idx)):

                # a match, short circuit
                return run_idx

    return None

@click.argument('gexp')
@click.command()
def main(gexp):

    from seh_pathway_hopping._tasks import (
        get_gexp_jobs_df,
        get_gexp_wepy_h5,
        get_gexp_master_orch,
    )

    orch = get_gexp_master_orch(gexp)
    jobs_df = get_gexp_jobs_df(gexp)
    wepy_h5 = get_gexp_wepy_h5(gexp)

    # we need to add continuations so set this here
    wepy_h5.set_mode('r+')

    # clear the continuations
    with wepy_h5:

        # delete the dset
        del wepy_h5.h5['_settings/continuations']

        # re-initialize
        wepy_h5._init_continuations()

    for row_idx, row in jobs_df.iterrows():

        # use the real start hash to get the continuations
        continued_run = orch.run_continues(
             row['start_hash'],
             row['end_hash'])


        if continued_run is None:
            continue

        cont_patched_start_hash, cont_end_hash = continued_run

        # get the run idxs in the HDF5 file so we can add the continuation

        # the idx of this run
        curr_run_idx = match_run_idx(wepy_h5,
                                     row['patched_start_hash'],
                                     row['end_hash'])

        # the run that was being continued
        cont_run_idx = match_run_idx(wepy_h5,
                                     cont_patched_start_hash,
                                     cont_end_hash)

        print("----------------------------------------")
        print("Continuing Run")
        print("start_hash:", row['start_hash'])
        print("patched_start_hash:", row['patched_start_hash'])
        print("end_hash:", row['end_hash'])
        print("HDF5 run idx", curr_run_idx)
        print("--------------------")
        print("Continued Run")
        # print("start_hash:", row['start_hash'])
        print("patched_start_hash:", cont_patched_start_hash)
        print("end_hash:", cont_end_hash)
        print("HDF5 run idx", cont_run_idx)
        print("--------------------")

        print("Continuation is:", curr_run_idx, cont_run_idx)
        print("Setting to HDF5")


        with wepy_h5:
            wepy_h5.add_continuation(curr_run_idx, cont_run_idx)


        print("----------------------------------------")

        print("Continuations after")
        with wepy_h5:
            print(wepy_h5.continuations)




if __name__ == "__main__":

    main()
