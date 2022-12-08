"""User settings for a project."""

# load the system configuration. You can override them in this module,
# but beware it might break stuff
from .sysconfig import *

### Envs

# which virtual environment tool to use: venv or conda
ENV_METHOD = 'conda'

# which env spec to use by default
DEFAULT_ENV = 'common_dev'


PROJECT_SLUG = "seh_pathway_hopping"

### Project

PROJECT_DIRS = None
RESOURCES = [
    'cache',
    'data',
    'db',
]
RESOURCE_DIR = "$HOME/tree/lab/resources/"
PROJECT_DIR = "$HOME/tree/lab/projects/seh.pathway_hopping"
