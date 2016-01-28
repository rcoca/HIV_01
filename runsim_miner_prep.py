"""
Provides code to run multiple simulations.
Uses sim_miner_prep, sim, concurrency, person_miner_prep
Also exposes a run method for enqueueing from outside in a queue of simulations.
"""
from sim_miner_prep import run_simulations,  run_simulations_mp
#, onesim

# Run model with medium sex worker population
from params import phi_ehg
from params_miner_prep import ehg_staged05_csw as model_params
import time

### model globals
model_params = model_params.copy()
model_params.update(
        model_phi=phi_ehg,
        outfile_pattern='eps??sim??.out',
        n_sim=10,  # the number of replications for each parameter set
        pop_size=1 * 10 ** 3,  # equal number of males & females
        sim_days=365*3
)
def run():
    outfolder = 'temp02'
    try:
        import os
        os.mkdir(outfolder)
    except OSError:
        pass
    model_name = 'test'
    model_params['outfolder'] = outfolder
    model_params['model_name'] = model_name
    T0=time.time()
    run_simulations_mp(model_params)
    print "simulation took {minutes:.2f} minutes".format(minutes=(time.time()-T0)/60.0)
    
if __name__ == '__main__':
    run()
