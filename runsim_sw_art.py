"""
Clone of runsim01. exposes a run method which is the handle by which to enqueue it from outside.
"""
from sim01 import run_simulations_mp, onesim


# Run model with medium sex worker population 
from params import ehg_staged03_csw as model_params, phi_ehg

### model globals
model_params = model_params.copy()
model_params.update(
    model_phi = phi_ehg,
    outfile_pattern = 'eps??sim??.out',
    n_sim = 10, #the number of replications for each parameter set
    pop_size = 1 * 10**3, #equal number of males & females
    sim_days = 3*365
    )

def run():
    outfolder = 'temp01'
    try:
        import os
        os.mkdir(outfolder)
    except OSError:
        pass
    model_name = 'test'
    model_params['outfolder'] = outfolder
    model_params['model_name'] = model_name
    run_simulations_mp(model_params)


    
if __name__ == '__main__':
    run()

