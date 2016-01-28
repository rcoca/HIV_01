#!/usr/bin/env python
"""
This script illustrates batching topgether multiple simulation to run one after the other.
Due to multiprocessing streching machine to the limit, there's no point running thema at the same time.
"""
import importlib

if __name__=='__main__':
    simulationScripts=["runsim_orig","runsim_sw_art","runsim_miner_prep"]
    for module in simulationScripts:
        importlib.import_module(module).run()
