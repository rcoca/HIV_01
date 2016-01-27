"""
Provides code to run multiple simulations.
Uses sim01, sim, concurrency, concurrency01
"""
from sim_miner_prep import run_simulations,  run_simulations_mp
#, onesim

# Run model with medium sex worker population
from params import ehg_staged05_csw as model_params, phi_ehg

### model globals
model_params = model_params.copy()
model_params.update(
        model_phi=phi_ehg,
        outfile_pattern='eps??sim??.out',
        n_sim=10,  # the number of replications for each parameter set
        pop_size=1 * 10 ** 4,  # equal number of males & females
)

if __name__ == '__main__':
    outfolder = 'temp'
    try:
        import os
        os.mkdir(outfolder)
    except OSError:
        pass
    model_name = 'test'
    model_params['outfolder'] = outfolder
    model_params['model_name'] = model_name

    run_simulations_mp(model_params)
