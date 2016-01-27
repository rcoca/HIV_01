"""
Provides constants and parameter set for simulation.
It does so by importing relevant parts from the pool
of parameters params.py and overriding/updating them
as needed by the current simulation.
This simulation aims to include two extra types of individuals
in the population: miners and PREP-treated persons.
"""
from params import shared_params, ehg_staged03_csw

# miners
p_miners  = 0.08

#under prep treatment
p_PREP    = 0.01
beta_PREP = (0.0,0.0,0.0,0.0)

#Base case scenario; should be the same as running sim.py and concurrency.py 
#get the shared params
ehg_staged01 = shared_params.copy()


#update the sim adding miners inside
ehg_staged05_csw = ehg_staged03_csw.copy()
ehg_staged05_csw.update(
        p_miners = p_miners,
        p_PREP   = p_PREP, 
        beta_PREP= beta_PREP
	)
