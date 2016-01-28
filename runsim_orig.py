"""
Clone of runsim. wraps the body of code in a run method to be called later.
"""


from sim01 import run_simulations_mp, onesim

## IMPORTANT: GET PARAMS FOR THIS **MODEL**:  SHARED BY SIMULATIONS!
from params import ehg_staged01 as model_params, phi_ehg

# Base scenario
### model globals
model_params = model_params.copy()
model_params.update(
        model_phi=phi_ehg,
        outfile_pattern='eps??sim??.out',
        n_sim=10,  # the number of replications for each parameter set
        pop_size=1 * 10 ** 3,  # equal number of males & females
)

def run():
    outfolder = 'temp00'
    try:
        import os
        os.mkdir(outfolder)
    except OSError:
        pass
    model_name = 'test00'
    model_params['outfolder'] = outfolder
    model_params['model_name'] = model_name
    run_simulations_mp(model_params)
    

if __name__ == '__main__':
    run()
