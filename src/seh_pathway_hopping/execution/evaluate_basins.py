import click
from seh_pathway_hopping.execution import parse_kwargs

import matplotlib.pyplot as plt
import itertools as it

## specs

AXES = [
    'unbound',
    'bound',
    'trim',
]

CASES_SPECS = {
    'unbound' : {
        'label' : "Unbound Cutoff (nm)",
        'cases' : [
            '0.7',
            '0.6',
            '0.5',
        ],
    },

    'bound' : {
        'label' : "Bound Cutoff (nm)",
        'cases' : [
            '0.1',
            '0.2',
            '0.3',
        ],
    },

    'trim' : {
        'label' : "Trim Cutoff (num. transitions)",
        'cases' : [
            '10',
            '5',
            '2',
            '1'
        ],
    },

}

CHOOSE_SPECS = {
    'unbound' : '0.6',
    'bound' : '0.3',
    'trim' : '1',
}

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        get_msn,
        load_msm,
        MSM_SPECS,
        BASIN_SPECS,
        TS_SPECS,
        calc_msm_committors,
        save_gexp_basin_eval_plots,
    )

    kwargs = parse_kwargs(specs)

    # get the show and save args
    if 'save' in kwargs:
        save = bool(kwargs['save'])
    else:
        save = False

    if 'show' in kwargs:
        show = bool(kwargs['show'])
    else:
        show = False


    # This is specified in the basin id since it is dependent on it for trimming
    gexp = kwargs['gexp']
    csn_id = kwargs['csn_id']

    case_0_key = kwargs['case_0']
    case_1_key = kwargs['case_1']

    case_0_label = CASES_SPECS[case_0_key]['label']
    case_1_label = CASES_SPECS[case_1_key]['label']

    # Load the msn
    msn = get_msn(
        csn_id,
        gexp,
    )

    cases_0 = CASES_SPECS[case_0_key]['cases']
    cases_1 = CASES_SPECS[case_1_key]['cases']

    # find which axis is not being used
    no_case_key = [axis for axis in AXES
                    if axis not in (case_0_key, case_1_key)]

    assert len(no_case_key) == 1
    no_case_key = no_case_key[0]

    # get the pinned value for the unused parameter
    no_case_val = CHOOSE_SPECS[no_case_key]


    # number of nodes in committors
    fig, axes = plt.subplots(
        len(cases_0),
        len(cases_1),
        constrained_layout=True,
        figsize=(10,10),
    )

    # sasas in committors
    sasa_fig, sasa_axes = plt.subplots(
        len(cases_0),
        len(cases_1),
        constrained_layout=True,
        figsize=(10,10),
    )

    def find_no_case(
            c0_case,
            c1_case,
    ):

        vals = {}

        for axis in AXES:

            # find the case that is the trim
            panel_case = list(it.filterfalse(
                    lambda x: x != axis,
                    (case_0_key, case_1_key,)
            ))

            if len(panel_case) == 0:

                vals[axis] = no_case_val

            else:

                if axis == case_0_key:
                    vals[axis] = c0_case

                elif axis == case_1_key:
                    vals[axis] = c1_case

        return vals['trim'], vals['bound'], vals['unbound']

    # get the ts ids for these cases
    ts_tups = []
    for c0_idx, c0_case in enumerate(cases_0):
        for c1_idx, c1_case in enumerate(cases_1):

            trim_val, bound_val, unbound_val = find_no_case(
                c0_case,
                c1_case,
            )


            # use these to build query strings for matching the TS Spec key
            msm_id = f"cutoff-{trim_val}"

            basin_id = f"msm-{msm_id}_basin-bound-{bound_val}-unbound-{unbound_val}-cutoff"

            ts_id = f"basin-{basin_id}_committor-msmb"

            plot_idx = (c0_idx, c1_idx,)
            ts_tups.append((
                ts_id,
                trim_val,
                unbound_val,
                bound_val,
                plot_idx,
            ))

    for plot_idx, ts_tup in enumerate(ts_tups):

        ts_id, trim_val, unbound_val, bound_val, plot_id = ts_tup

        ts_spec = TS_SPECS[ts_id]

        basin_id = ts_spec['basin_id']
        committor_method = ts_spec['committor_method']

        print("----------------------------------------")
        print(f"Doing Analysis for plot: {plot_id}")
        print(f"Trimming Value: {trim_val}")
        print(f"Unbound Value: {unbound_val}")
        print(f"Bound Value: {bound_val}")

        msm_id = BASIN_SPECS[basin_id]['msm_id']

        bound_basin = msn.node_groups[f'basin-{basin_id}/bound_basin']
        unbound_basin = msn.node_groups[f'basin-{basin_id}/unbound_basin']

        print("Bound Basin", bound_basin)
        print("Unbound Basin", unbound_basin)

        print("Calculating committors")
        unbound_committors, bound_committors = calc_msm_committors(
            msm_id,
            gexp,
            bound_basin,
            unbound_basin,
            method=committor_method,
        )
        print("Finished Calculating committors")



        unbound_method = BASIN_SPECS[basin_id]['unbound']['method']
        unbound_params = BASIN_SPECS[basin_id]['unbound']['method_kwargs']

        bound_method = BASIN_SPECS[basin_id]['bound']['method']
        bound_params = BASIN_SPECS[basin_id]['bound']['method_kwargs']

        title_str = f"Trim: {trim_val}; Unbound: {unbound_val}; Bound {bound_val}"


        ## Plot histogram of num nodes over the committor values

        axes[plot_id].hist(unbound_committors)
        axes[plot_id].set_xlabel("Unbound Committor (P)")
        axes[plot_id].set_ylabel("Node Counts")
        axes[plot_id].set_title(title_str)


        ## histogram of total weights over committors

        nodes_total_weights = msn.get_nodes_attribute('_observables/total_weight')

        # get them associated with the committors we have for plotting
        total_weights = []
        for node_idx, committor in enumerate(unbound_committors):

            node_id = msn.node_idx_to_id(node_idx)

            total_weights.append(nodes_total_weights[node_id])

        # plot the weights and set the range to around the TS
        sasa_axes[plot_id].hist(
            unbound_committors,
            weights=total_weights,
            range=(0.3, 1.0),
        )

        # also draw a vertical line at 0.5 committor
        sasa_axes[plot_id].axvline(x=0.5, color='k')

        sasa_axes[plot_id].set_xlabel("Unbound Committor (P)")
        sasa_axes[plot_id].set_ylabel("Total Weight")
        sasa_axes[plot_id].set_title(title_str)


        # OLD: don't do sasas now
        # ## Plot chart of SASAs over the committor probabilities

        # nodes_mean_sasa = msn.get_nodes_attribute('obs/lig-sasa_npoints-100/mean')

        # # get them associated with the committors we have for plotting
        # mean_sasas = []
        # for node_idx, committor in enumerate(unbound_committors):

        #     node_id = msn.node_idx_to_id(node_idx)

        #     mean_sasas.append(nodes_mean_sasa[node_id])

        # sasa_axes[plot_id].scatter(
        #     mean_sasas,
        #     unbound_committors,
        # )
        # sasa_axes[plot_id].set_xlabel("Ligand SASA")
        # sasa_axes[plot_id].set_ylabel("Committor Probability")
        # sasa_axes[plot_id].set_title(title_str)

    if save:

        save_gexp_basin_eval_plots(
            gexp,
            case_0_key,
            case_1_key,
            num_nodes_hist=fig,
            weights_hist=sasa_fig,
        )

    if show:
        print("SHOWING")
        plt.show()

    # save the plots


if __name__ == "__main__":

    cli()
