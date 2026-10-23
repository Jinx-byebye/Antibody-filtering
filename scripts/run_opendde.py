#!/usr/bin/env python3
"""Run OpenDDE with lossless confidence output for downstream scoring."""
from abfilter.predictors.opendde_output import install_lossless_dumper

install_lossless_dumper()

from runner.inference import run

if __name__ == '__main__':
    run()
