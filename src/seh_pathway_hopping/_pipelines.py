"""Generated file from the analysis.org file. Do not edit directly."""
from copy import copy
import inspect

from prefect import Flow
import prefect

import seh_pathway_hopping._tasks as tasks_module

# these helper functions are for automatically listing all of the
# functions defined in the tasks module
def is_mod_function(mod, func):
    return inspect.isfunction(func) and inspect.getmodule(func) == mod

def get_functions(mod):

    # get only the functions that aren't module functions and that
    # aren't private
    return {func.__name__ : func for func in mod.__dict__.values()
            if (is_mod_function(mod, func) and
                not func.__name__.startswith('_')) }

# get the task functions and wrap them as prefect tasks
tasks = {name : prefect.task(func)
         for name, func in get_functions(tasks_module).items()}

@prefect.task
def say_hello(person: str) -> None:
    print("hello world, {}".format(person))

@prefect.task
def add(x, y=1):
    return x + y

test_flow = prefect.Flow("flow for testing")

with test_flow as flow:
    a = add(1, 2)
    b = add(100, 100)
    c = say_hello(str(a))
    d = say_hello(str(b))
