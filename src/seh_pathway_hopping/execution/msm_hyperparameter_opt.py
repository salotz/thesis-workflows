import click
@click.command()
@click.argument('gexp')
@click.argument('observable')
def main(gexp, observable):

      from seh_pathway_hopping._tasks import (
          score_msm_gmrq_crossval_classifier,
      )

    space = [
        sp.Integer(100, 110, name="n_clusters"),
    ]



    @use_named_args(space)
    def objective(**params):

        # the classifier to use
        clf = MiniBatchKMeans(
            n_clusters=params['n_clusters'],
            batch_size=100,
        )

        # the cross-validation scheme as a splitter object
        splitter = KFold(n_splits=n_splits)

        # do the cross-validation score for the classifier
        stats = score_msm_gmrq_crossval_classifier(
            gexp,
            observable,
            clf,
            splitter,
            lag_time=lag_time,
        )


        # the final score is the negative of the mean
        return -stats['mean']


    ## There are different ways to control optimization with skopt

    method = 'highlevel'

    # the high level function and the sklearn compatible wrapper both
    # use callbacks to get info. Otherwise you can control the loop
    # with the optimizer.

    # callbacks for reporting
    def print_result_cb(res):

        msg = f"""
Current Results
---------------

- last objective function score :: {res['fun']}
- solution :: {res['x']}
"""
        print(msg)

    verbose_cb = skopt.callbacks.VerboseCallback(1)

    DEADLINE = 108000 # seconds
    deadline_cb = skopt.callbacks.DeadlineStopper(DEADLINE)

    # use the high level function
    if method == 'highlevel':

        res_gp = skopt.gp_minimize(
            objective,
            space,
            n_calls=10,
            callback=[
                print_result_cb,
                verbose_cb,
                deadline_cb,
            ],
        )

    elif method == 'sklearn':

        pass

    elif method == 'asktell':

        optimizer = skopt.Optimizer(
            space,
        )

        for step in range(3):

            # get the next sample of hyperparameters from the optimizer
            next_sample = optimizer.ask()

            print(next_sample)

            # compute the objective function
            obj_val = objective(next_sample)

            print(f"Objective Function Score: {obj_val}")

            # Report the score for the sample to the optimizer

            optimizer.tell(
                next_sample,
                obj_val,
            )



if __name__ == "__main__":

    main()

    import itertools as it

    import numpy as np

    import skopt
    import skopt.space as sp
    from skopt.utils import use_named_args

    from sklearn.cluster import MiniBatchKMeans
    from sklearn.model_selection import (
        ShuffleSplit,
        KFold,
    )

    from wepy.analysis.network import MacroStateNetwork

    from seh_pathway_hopping._tasks import *

    ## perform the clustering on the observable using only the
    ## training runs

    # space = [
    #     sp.Integer(2, 1000, name="n_clusters"),
    #     sp.Integer(100, 10000, name="batch_size"),
    # ]
