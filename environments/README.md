# Separate environments

`cpu.txt` pins the CPU reproduction dependency. Python 3.13 with NumPy 2.2.6 has
been tested from a fresh virtual environment; the package supports Python 3.10+.

Do not install all predictors into this environment. Their CUDA, PyTorch,
JAX, ESM, and Transformers dependencies differ. Each model is invoked through its
own Python interpreter or command-line executable, supplied in `resources.json`.

The three predictor constraint files record core versions from the installations
referenced by the launchers. Source revisions are recorded in
`reference/predictor_installations.json`, including the ESM and Transformers forks.
Use constraints with the provider's source installation procedure; these files
are not standalone installers or complete dependency locks. Small inference tests
passed in these existing GPU environments; clean installation tests
and immutable historical checkpoint verification remain pending. Installing
current unpinned packages would not establish reproduction of the paper's settings.
See [prediction instructions](../docs/prediction.md) for supported parameters and
required assets.

`opendde.txt` installs the exact OpenDDE source revision with GPU extras.
Use an independent environment and verify the checkpoint hash from
`configs/predictors/opendde.json`. The general-purpose checkpoint, common runtime
assets, and local MSA files are required; weights are not distributed here.
