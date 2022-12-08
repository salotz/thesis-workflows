if __name__ == "__main__":

    from seh_pathway_hopping._tasks import (
        save_real_rep_top,
        make_legacy_results_linker_h5,
    )

    # create the reference topology from the source materials in the
    # current project
    save_real_rep_top()

    # NOTE: Don't do this
    # and the centered one (must come after the first one)
    # save_centered_real_rep_top()

    # using that topology and some other stuff make a linker HDF5 with
    # the proper headers and such for the new dataset
    make_legacy_results_linker_h5()
