import sys

import click

from dask.distributed import Client, LocalCluster

def parse_kwargs(kwarg_strings):
    """Parse kwargs on command line of the form key=value.

    Parameters
    ----------

    kwarg_strings : list of str
        The kwarg strings each is like 'key=value'

    Returns
    -------
    kwargs : dict of str : str

    """

    kwargs = {}
    for kwarg_str in kwarg_strings:
        key, value = kwarg_str.split('=')
        kwargs[key] = value

    return kwargs


def execute(ctx, func):
    """Function to execute the work function. This is what accepts command
    line arguments to connect a function to an executor.

    You may run immediately with this in the __main__ block:

    execute(func)

    Or partially evaluate as a thunk and invoke later in a __main__ block:

    functools.partial(execute, func)

    """

    if sys.argv[-1] == '-h' or sys.argv[-1] == '--help':
        print("Usage: execute <method> [key=value, ...]")

    cluster_address = sys.argv[1]
    kwargs = parse_kwargs(sys.argv[2:])

    DASHBOARD_PORT = 9998
    N_WORKERS = 2
    PROCESSES = False

    worker_kwargs = {
        'memory_limit' : '8GB',
        'memory_spill_fraction' : 1.0,
    }

    # if the address is None just start a local cluster, with the default options
    if cluster_address == ':local':

        # start a local cluster
        cluster = LocalCluster(processes=PROCESSES,
                               n_workers=N_WORKERS,
                               dashboard_address=f":{DASHBOARD_PORT}",
                               **worker_kwargs)
        print(f"Ad hoc cluster online. Dashboard on port {DASHBOARD_PORT}")

        client = Client(cluster)

    # otherwise just connect
    else:
        client = Client(cluster_address)


    func(client, **kwargs)
