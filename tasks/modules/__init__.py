
# should be copied in by the installation process
from . import core
from . import clean
from . import env
from . import git
from . import docs
from . import py
from . import lxd

from . import cache
from . import project
from . import analytics

MODULES = [
    core,
    clean,
    env,
    git,
    docs,
    py,
    lxd,
    cache,
    project,
    analytics,
]
