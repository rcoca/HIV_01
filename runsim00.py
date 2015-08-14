"""
Provides code to run multiple simulations.
Uses sim01, sim, concurrency, concurrency01
"""
from sim01 import run_simulations, onesim

## IMPORTANT: GET PARAMS FOR THIS **MODEL**:  SHARED BY SIMULATIONS!
from params import ehg_staged01 as model_params, phi_ehg


# Base scenario
### model globals
model_params = model_params.copy()
model_params.update(
	model_phi = phi_ehg,
	outfile_pattern = 'eps??sim??.out',
	n_sim = 100, #the number of replications for each parameter set
	pop_size = 2 * 10**4, #equal number of males & females 
	)

if __name__ == '__main__':
	outfolder = 'temp00'
	model_name = 'test00'
	model_params['outfolder'] = outfolder 
	model_params['model_name'] = model_name
	run_simulations(model_params)


