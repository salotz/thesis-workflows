# size of test set for scoring the PCA model
TEST_SIZE = 0.25

import click
from seh_pathway_hopping.execution import parse_kwargs

import matplotlib.pyplot as plt
from tabulate import tabulate

@click.command(context_settings=dict(
    ignore_unknown_options=True,))
@click.argument('specs', nargs=-1, type=click.UNPROCESSED)
def cli(specs):

    from seh_pathway_hopping._tasks import (
        ts_pca_model_selection,
        save_tspca_model,
        save_tspca_score_table,
        make_tspca_model_score_table,
        plot_tspca_model_score,
    )

    kwargs = parse_kwargs(specs)

    # This is specified in the basin id since it is dependent on it for trimming
    # msm_id = kwargs['msm_id']
    csn_id = kwargs['csn_id']
    ts_id = kwargs['ts_id']

    ## Analysis

    # choose which gexps have a TS model that can actually be used for
    # the model selection, i.e. if some simulations couldn't have a TS
    # predicted exclude them
    TS_PCA_MODEL_GEXPS = (
        'TPPU-legacy',
        '3',
        '10',
        '17',
    )

    # calculate the various subsamplings etc for the PCA model and
    # return them all with their scores
    model_results = ts_pca_model_selection(
        TS_PCA_MODEL_GEXPS,
        csn_id,
        ts_id,
        test_size=TEST_SIZE,
    )


    # Save the PCA model and associated data
    for tspca_id, result_d in model_results.items():

        _ = save_tspca_model(
            tspca_id,
            result_d['model'],
            result_d['model_score'],
            result_d['mode_scores'],
            result_d['test_size'],
            result_d['gexps'],
        )

    # make the table of the scores:
    model_score_df = make_tspca_model_score_table(
        model_results,
    )

    # print the table nicely
    table_str = tabulate(
        model_score_df,
        headers=model_score_df.columns,
        tablefmt='orgtbl',
    )

    print("Results of the model selection")
    print(table_str)

    # save it
    save_tspca_score_table(model_score_df)

    # plot scores of the model
    fig, axes = plot_tspca_model_score(model_results)

    # show the plots
    plt.show()

if __name__ == "__main__":

    cli()
