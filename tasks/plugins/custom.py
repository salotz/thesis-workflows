"""Put user defined tasks in the plugins folder. You can start with
some customizations in this file which is included by default."""

from invoke import task

import platform
from invoke.vendor.six.moves import input

def confirm(question, assume_yes=True):
    """
    Ask user a yes/no question and return their response as a boolean.
    ``question`` should be a simple, grammatically complete question such as
    "Do you wish to continue?", and will have a string similar to ``" [Y/n] "``
    appended automatically. This function will *not* append a question mark for
    you.
    By default, when the user presses Enter without typing anything, "yes" is
    assumed. This can be changed by specifying ``assume_yes=False``.
    .. note::
        If the user does not supplies input that is (case-insensitively) equal
        to "y", "yes", "n" or "no", they will be re-prompted until they do.
    :param str question: The question part of the prompt.
    :param bool assume_yes:
        Whether to assume the affirmative answer by default. Default value:
        ``True``.
    :returns: A `bool`.
    """
    # Set up suffix
    if assume_yes:
        suffix = "Y/n"
    else:
        suffix = "y/N"
    # Loop till we get something we like
    # TODO: maybe don't do this? It can be annoying. Turn into 'q'-for-quit?
    while True:
        # TODO: ensure that this is Ctrl-C friendly, ISTR issues with
        # raw_input/input on some Python versions blocking KeyboardInterrupt.
        response = input("{0} [{1}] ".format(question, suffix))
        response = response.lower().strip()  # Normalize
        # Default
        if not response:
            return assume_yes
        # Yes
        if response in ["y", "yes"]:
            return True
        # No
        if response in ["n", "no"]:
            return False
        # Didn't get empty, yes or no, so complain and loop
        err = "I didn't understand you. Please specify '(y)es' or '(n)o'."
        print(err, file=sys.stderr)



## linking and referencing project resources

## project resources go in
## $TREE/$DOMAIN/resources/project-resources/$project_name

# these are for each project in $TREE/$DOMAIN/projects/$project_name

# they must have exactly the same structure as the project directory
# except with the resources folders filled in. Here we will simply
# create a bunch of symlinks in the locations of the resource
# files/folders in the project resource dir

PROJECT_NAME = "seh.pathway_hopping"

# STUB: choose the root paths

HPCC_SCRATCH_DIR = "/mnt/gs18/scratch/users/lotzsamu"

SCRATCH_PROJECT_DIR = f"{HPCC_SCRATCH_DIR}/tree/lab/projects/{PROJECT_NAME}"
SCRATCH_RESOURCE_DIR = f"{HPCC_SCRATCH_DIR}/tree/lab/resources/project-resources/{PROJECT_NAME}"

HOME_PROJECT_DIR = "$HOME/tree/lab/projects/{}".format(PROJECT_NAME)
HOME_RESOURCE_DIR = "$HOME/tree/lab/resources/project-resources/{}".format(PROJECT_NAME)


if platform.node() in ['ostrich', 'superior',]:
    PROJECT_DIR = HOME_PROJECT_DIR
    RESOURCE_DIR = HOME_RESOURCE_DIR
elif platform.node() in [
        'dev-intel16-k80',
        'dev-intel18',
        'dev-intel16',
        'dev-intel14',
        'gateway-03',
        'gateway-02',
        'gateway-01',
        'gateway-00',
        'globus-01', # rsync.hpcc.msu.edu
        'globus-00',
]:
    # then use scratch
    PROJECT_DIR = SCRATCH_PROJECT_DIR
    RESOURCE_DIR = SCRATCH_PROJECT_DIR


# relative paths to $project_dir of things considered resources
RESOURCES = [
    'cache',
    'data',
    'db',
    'results',
]

@task
def link_resources(ctx):
    """Make links to the project resource folders in this project"""

    for resource in RESOURCES:

        command = "ln -s -r -f -T {res}/{resource} {proj}/{resource}".format(
            res=RESOURCE_DIR,
            proj=PROJECT_DIR,
            resource=resource)

        print(command)
        if confirm("Is this command okay?"):
            ctx.run(command)

@task
def exe_scripts(cx):

    for p in [f"{PROJECT_DIR}/prep", f"{PROJECT_DIR}/hpcc/scripts"]:
        cx.run(f'chmod ug+x {p}/*')

@task(post=[exe_scripts])
def org_tangle(ctx):

    ctx.run("emacs -Q --batch -l org project.org -f org-babel-tangle")



@task
def org_clean(ctx):

    # specify extensions so we don't remove the .keep file (TODO
    # probably a better way to do this)

    # prep
    ctx.run("rm -f prep/*.sh")
    ctx.run("rm -f prep/*.bash")
    ctx.run("rm -f prep/*.xsh")
    ctx.run("rm -f prep/*.py")

    # hpcc scripts
    ctx.run("rm -f hpcc/scripts/*.sh")
    ctx.run("rm -f hpcc/scripts/*.bash")
    ctx.run("rm -f hpcc/scripts/*.xsh")
    ctx.run("rm -f hpcc/scripts/*.py")

    # templates
    ctx.run("rm -f hpcc/templates/*.j2")


# HPCC job generation stuff

LIG_IDS = [3, 10, 17, 18, 20]
RES_ID = "UNL"

@task
def init_lig_sim_dir(cx, lig=None):

    if lig is None:
        lig_ids = LIG_IDS
    else:
        lig_ids = [lig]

    for lig_id in lig_ids:
        sim_dir = f"{PROJECT_DIR}/hpcc/simulations/{lig_id}_simulations"
        cx.run(f'mkdir -p "{sim_dir}"')

        cx.run(f'mkdir -p "{sim_dir}/logs"')
        cx.run(f'touch "{sim_dir}/logs/.keep"')

        cx.run(f'mkdir -p "{sim_dir}/submissions"')
        cx.run(f'touch "{sim_dir}/submissions/.keep"')

        cx.run(f'mkdir -p "{sim_dir}/tasks"')
        cx.run(f'touch "{sim_dir}/tasks/.keep"')

        # create the resource directories
        res_dir = f"{RESOURCE_DIR}/hpcc/simulations/{lig_id}_simulations"
        cx.run(f'mkdir -p "{res_dir}/configurations"')
        cx.run(f'mkdir -p "{res_dir}/input"')
        cx.run(f'mkdir -p "{res_dir}/orchs"')
        cx.run(f'mkdir -p "{res_dir}/results"')

        # link to the resource directories
        cx.run(f'ln -r -s -f -T "{res_dir}/configurations" "{sim_dir}/configurations"')
        cx.run(f'ln -r -s -f -T "{res_dir}/input" "{sim_dir}/input"')
        cx.run(f'ln -r -s -f -T "{res_dir}/orchs" "{sim_dir}/orchs"')
        cx.run(f'ln -r -s -f -T "{res_dir}/results" "{sim_dir}/results"')


@task(pre=[org_tangle])
def gen_configurations(cx):

    worker_types = ["Worker", "OpenMMGPUWorker", "OpenMMCPUWorker"]

    # this is not really useful to have be specific since it is always set at runtime
    num_workers = 1

    for worker_type in worker_types:
        cx.run(f"bash -i {PROJECT_DIR}/prep/gen_configuration.sh {worker_type} 1")

@task
def clean_configurations(cx):

    for lig_id in LIG_IDS:
        cx.run(f'rm -rf "{PROJECT_DIR}/data/configurations/{lig_id}"/*')
        cx.run(f'rm -rf "{PROJECT_DIR}/hpcc/simulations/{lig_id}_simulations/input"/*')
        cx.run(f'rm -rf "{PROJECT_DIR}/hpcc/simulations/{lig_id}_simulations/configurations"/*')


@task(pre=[org_tangle])
def hpcc_gen_tasks(c):

    for lig_id in LIG_IDS:
        print("Generating tasks for Ligand {}".format(lig_id))

        sim_dir = './hpcc/simulations/{}_simulations'.format(lig_id)

        # remove the exisiting tasks
        c.run('rm -f "{}/tasks/*"'.format(sim_dir))

        # regenerate the scripts
        c.run('python ./hpcc/scripts/gen_tasks.py '
              './hpcc/templates {sim_dir}/tasks {sim_dir}/run_spec.toml '.format(
            sim_dir=sim_dir))


@task()
def hpcc_clean_tasks(cx):

    for lig_id in LIG_IDS:

        print("removing tasks for: {}".format(lig_id))

        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(PROJECT_DIR, lig_id)

        # remove the exisiting tasks
        cx.run('rm -f "{}"/tasks/*'.format(sim_dir))

@task(pre=[org_tangle, hpcc_gen_tasks])
def hpcc_gen_submissions(c, tag=None):

    for lig_id in LIG_IDS:
        print("Generating submissions for Ligand {}".format(lig_id))

        if tag is not None:
            label = "lig_{}_sims_{}".format(lig_id, tag)
        else:
            label = "lig_{}_sims".format(lig_id)

        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(SCRATCH_PROJECT_DIR, lig_id)

        command = \
"""slurmify \
        --config ./hpcc/run_settings.toml \
        --context {sim_dir}/context_settings.toml \
        --batch-in {sim_dir}/tasks \
        --batch-out {sim_dir}/submissions \
        {label}
""".format(
    sim_dir=sim_dir, label=label)

        # regenerate the scripts
        c.run(command)

@task()
def hpcc_clean_submissions(cx):

    for lig_id in LIG_IDS:


        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(PROJECT_DIR, lig_id)
        print("removing submissions: {}".format(sim_dir))

        # remove the exisiting tasks
        cx.run('rm -f "{}"/submissions/*'.format(sim_dir))

@task
def hpcc_test_submit(cx):

    for lig_id in LIG_IDS:

        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(PROJECT_DIR, lig_id)
        print("Submitting test run for: {}".format(sim_dir))

        # the sbatch command must be run from the simulation directory
        # since it is context sensitive
        command =\
'''(cd {sim_dir}; sbatch submissions/lig-{lig_id}_test.sh.slurm)
'''.format(
    lig_dir=lig_id,
    sim_dir=sim_dir
)
        # remove the exisiting tasks
        cx.run(command)


@task(pre=[hpcc_gen_tasks])
def hpcc_test_local(cx, lig=None):

    if lig is None:
        lig_ids = LIG_IDS
    else:
        lig_ids = [lig]

    for lig_id in lig_ids:

        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(PROJECT_DIR, lig_id)
        print("Running local test run for: {}".format(sim_dir))

        task_script = "{}/tasks/lig-{}_test-local.sh".format(sim_dir, lig_id)

        test_job_dir = '{}/test_jobs/test'.format(sim_dir)
        inputs_dir = '{}/input'.format(sim_dir)

        # make the job directory, cleaning it if it had stuff
        cx.run('mkdir -p "{}"'.format(test_job_dir))
        cx.run('rm -rf "{}"/*'.format(test_job_dir))

        # copy the inputs to the directory
        cx.run('cp "{}"/* "{}"/'.format(inputs_dir, test_job_dir))

        command =\
'''(cd "{job_dir}"; {task})
'''.format(
    task=task_script,
    job_dir=test_job_dir,
)

        print(command)
        if confirm("Run this command?"):
            cx.run(command)



@task
def hpcc_clean_tests(cx):

    for lig_id in LIG_IDS:

        sim_dir = '{}/hpcc/simulations/{}_simulations'.format(PROJECT_DIR, lig_id)
        print("Cleaning tests for: {}".format(sim_dir))

        cx.run('rm -rf "{sim_dir}"/test_jobs')

@task(hpcc_gen_submissions)
def hpcc_submit(cx):
    # TODO
    pass

@task(org_tangle)
def hpcc_push_scripts(cx):

    # the actual tangled scripts
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/hpcc/scripts/ \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/hpcc/scripts
""")

    # the tempalates
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/hpcc/templates/ \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/hpcc/templates
""")

    # the tasks
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/tasks/ \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/tasks
""")

    # the analysis file
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/project.org \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/project.org
""")

    # env stuff
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/envs/ \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/envs
""")


    # the run settings
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/hpcc/run_settings.toml \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/hpcc/
""")

    ## analysis jobs

    # dask_server
    cx.run(
"""rsync -ravvhhiz $HOME/tree/lab/projects/seh.pathway_hopping/hpcc/analysis/dask_server/submissions/ \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/hpcc/analysis/dask_server/submissions/
""")

    ## src/seh_pathway_hopping
    cx.run(
"""rsync -ravvhhiz \
$HOME/tree/lab/projects/seh.pathway_hopping/src \
lotzsamu@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/lotzsamu/tree/lab/projects/seh.pathway_hopping/
""")



@task
def hpcc_clean_scripts(cx):

    # prep scripts
    cx.run(f'rm -rf "{PROJECT_DIR}"/prep/*')
    cx.run(f'touch "{PROJECT_DIR}"/prep/.keep')

    # hpcc/script scripts
    cx.run(f'rm -rf "{PROJECT_DIR}"/hpcc/scripts/*')
    cx.run(f'touch "{PROJECT_DIR}"/hpcc/scripts/.keep')


@task(pre=[hpcc_clean_tasks,
           hpcc_clean_submissions,
           hpcc_clean_tests,
           hpcc_clean_scripts])
def hpcc_clean(cx):
    pass

@task(pre=[org_clean, hpcc_clean])
def clean(cx):

    # also clean the tmp dir
    cx.run('rm -rf "{}"/tmp'.format(PROJECT_DIR))
    cx.run('touch "{}"/tmp/.keep'.format(PROJECT_DIR))

gephi_stuff = (
    'gexf',
    'edges',
    'nodes',
    'gephi',
)

@task
def push_to_superior(cx):

    for thing in gephi_stuff:
        cx.run(f"rsync -ravhiz ~/tree/lab/projects/seh.pathway_hopping/data/csn/{thing} salotz@superior.bch.msu.edu:~/local/seh.pathway_hopping/data/csn/")

@task
def pull_from_superior(cx):

    for thing in gephi_stuff:
        cx.run(f"rsync -ravhiz salotz@superior.bch.msu.edu:~/local/seh.pathway_hopping/data/csn/{thing} ~/tree/lab/projects/seh.pathway_hopping/data/csn/")

