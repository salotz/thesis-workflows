from invoke import task

import joblib

jlmem = joblib.Memory(f"./cache/joblib")

@task
def clear(cx):

    jlmem.clear()
