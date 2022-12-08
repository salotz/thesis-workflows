from invoke import task

from ..config import (
    EXE_DIRS,
    TANGLE_DIRS,
)

@task(default=True)
def tangle(cx):
    """Tangle all the code blocks for dirs in TANGLE_DIRS."""

    # tangle them
    cx.run("emacs -Q --batch -l org project.org -f org-babel-tangle")

    for tangle_dir in EXE_DIRS:
        # make them executable
        cx.run(f'chmod ug+x ./{tangle_dir}/*')

@task
def clean(cx):
    """Clean all the tangled documents in TANGLE_DIRS."""

    for tangle_dir in TANGLE_DIRS:
        cx.run(f"rm -f {tangle_dir}/*")
        cx.run(f"touch {tangle_dir}/.keep")
