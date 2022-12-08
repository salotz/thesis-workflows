# Thesis Workflows

This repo is a port of the source code only from the code I used in my thesis
workflows.

It is written such that most of the tasks and entrypoints are written in a
literate Org-mode document and then tangled to a source code structure under
`src`.

You needn't worry about these details however when looking at the code except
that the output files are very large and not how I would normally write modules.

The `jigs` folder contains some other helper setups I used for running
Prometheus monitoring that got piped over port-forwarded via SSH to locally
hosted Grafana dashboards on our HPC systems.
