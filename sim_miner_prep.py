""" 
This code uses sim.py to add commercial sex workers to the model

"""
#import functools, random
#import numpy as np
# from concurrency01 import CommercialPartnership, ARTPartnership
import multiprocessing as mp
import  glob, math, shutil, tempfile, time, pprint
import os.path as path
from itertools import izip as zip  # Python 2.7
from collections import defaultdict
import datetime, logging
try:
    import os
    os.mkdir("out")
except:
    pass
logging.basicConfig(level=logging.INFO, filename='out/temp.log', filemode='w')
logging.info(datetime.datetime.today().isoformat())

from numpy.random import RandomState
import logging
import sim
from concurrency import Person, Partnership, stagedHIVfactory
import gc
from concurrency01 import Person01,Person02
from person_miner_prep import  Person03,Person04

#  11 different values of epsilon (the concurrency parameter)
#  EHG's R version: epsilon = 1-1.02*(exp(3.9*(seq(0,1,.1)-1))-.02)  
#  :note: lowest index is biggest (approx 1) epsilon
ehg_epsilons = tuple(1 - 1.02 * (math.exp(3.9 * (xi / 10. - 1)) - .02) for xi in range(11))
nOutGroups = 12  # based on EHG's Rsim.cpp; otherwise arbitrary


class Scheduler(sim.Scheduler):
    def coresim(self, males, females, day, params):
        """Return None; modify `params`.
        Run the core stages of the simulation (one iteration):
        form partnerships, do HIV transmissions,
        dissolve partnerships, do HIV deaths.
        These are run every iteration (including burn days).
        """
        self.form_partnerships(males, females, day)
        self.hiv_transmissions(day)
        self.dissolve_partnerships(day)
        self.hiv_deaths(day)

    def form_partnerships(self, males, females, day):  # chkchkchk schedule partnerships
        params = self.params
        ##### form partnerships
        logging.debug('coresim: BEGIN form partnerships')
        old_partnership_count = self.n_partnerships  # chkchkchk
        assert self.count_partnerships() == old_partnership_count
        # population sizes
        yearday = day%365 # the day in the year
        print day,yearday
        if yearday < Person03.activeRange[0] or  yearday>= Person03.activeRange[-1]:            
            filtered_males = filter(lambda x:not isinstance(x,Person03),males)
            msg="day: {day}, year day: {yday} Filtered out {nminers} from partnership formation".\
              format(nminers=len(males)-len(filtered_males),day=day,yday=yearday)
            logging.info(msg)
            print msg
        else:
            filtered_males = males

        nM = len(filtered_males);
        nF = len(females);
        # we follow EHG and impose a hard ceiling on the number of partnerships!
        # (there a lots of things not to like about this)
        max_new_partnerships = (nM + nF) // 2 - self.n_partnerships
        assert type(max_new_partnerships) is int
        # determine how many new partnerships to form (binomial draw with probability rho)
        prng = params['prng']
        rho = params['rho']  # aggregate partnership formation paramter
        assert type(rho) is float
        nFormPartnerships = 0  # needed to match EHG behavior (R's rbinom(0, rho)=0)
        if max_new_partnerships > 0:
            nFormPartnerships = prng.binomial(n=max_new_partnerships, p=rho)
        # form exactly nFormPartnerships new partnerships
        npairings = 0  # counter for testing that `random_pairings` returns the right number
        #:note: random pairings uses phi, which is in `params`
        for male, female in random_pairings(nFormPartnerships, filtered_males, females, params):
            newpartnership = Partnership(male, female, day, registry=self, params=params)  # form partnership
            npairings += 1
        assert npairings == nFormPartnerships
        logging.debug('coresim: END form partnerships')
        logging.info('\t{0} partnerships formed'.format(nFormPartnerships))
        assert self.count_partnerships() == old_partnership_count + nFormPartnerships


def random_pairings(n_pairs, males, females, params):
    """Return n_pairs random male-female pairs (**with** replacement).

    n_pairs : int
        the number of pairs to form
    males : list
        the potential male partners
    females : list
        the potential female partners
    params : dict
        params['prng'] is the random number generator;
        params['phi'] is a function that returns partnership formation probability
    """
    _phi = params['sim_phi']
    _prng = params['prng']
    _randint = _prng.randint
    _random_sample = _prng.random_sample
    if n_pairs < 0:
        raise ValueError('cannot form a negative number of pairs')
    nM = len(males)
    nF = len(females)
    assert n_pairs < nM and n_pairs < nF
    while n_pairs:
        # `n_pairs` random males and females (*with* replacement)
        males2pair = (males[idxM] for idxM in _randint(0, nM, n_pairs))
        females2pair = (females[idxF] for idxF in _randint(0, nF, n_pairs))
        draws = _random_sample(n_pairs)
        # must form pairs sequentially due to concurrency resistence
        for (male, female, draw) in zip(males2pair, females2pair, draws):
            # form partnership with probability phi, which varies with concurrency
            if test_random_pair(male, female, draw, _phi):
                n_pairs -= 1
                yield male, female


def test_random_pair(male, female, draw, phi):
    pairup = False
    if male.comsex and female.comsex:
        pairup = True
    elif (draw < phi(male, female)):  # ai: EHG use binomial; why?
        # ai: check that partnership doesn't already exist
        if female not in male.partners:
            assert male not in female.partners
            pairup = True
    return pairup

# Identical copy with the one in sim.py. onesim creates namespace collision being present in both
# maybe this needs refactoring/splitting in generic/specific
def run_simulations_mp(params, workerpool=None):
    """
    Doesn't return.
    Runs sumulations, same as run_simulation.
    Expected input is the dictionary of input parameters for simulation run.
    It uses multiprocessing to distribute work to more processors.
    """
    import sys, os
    if 'linux' in sys.platform:
        logging.info("Installing interrupt handler")
        import signal
        def signal_handler(signal, frame):
            logging.warn('Killing worker %d' % os.getpid())
            sys.exit(0)

        # only for *nix-based
        signal.signal(signal.SIGINT, signal_handler)

    if workerpool is None:
        n_processors = mp.cpu_count()
        print '{0} cores available'.format(n_processors)
        workerpool = mp.Pool(n_processors)

    if 'linux' in sys.platform:
        signal.signal(signal.SIGINT, signal.SIG_DFL)
    try:
        workerpool.map(onesim, [p for p in get_param_set(params)])
    except KeyboardInterrupt:
        logging.exception("INTERRUPTED")
        pass
    workerpool.close()
    workerpool.join()

def run_simulations(params):
    """Return None. Run all pending simulation for this model.
    Input should be a dict of model-specific parameters.
    """
    for simparams in get_param_set(params):
        onesim(simparams)


def get_param_set(params):
    """Yield dict, a replicate specific set of parameters.
    There are n_sim (e.g., 100) replicates for each epsilon,
    and each one has its own random seed, sim_name, and outfilename.
    This generator function should yield dicts with **pickleable** elements only,
    so that it can be used with multiprocessing.
    """
    outfolder = params['outfolder']
    outfile_pattern = 'eps??sim??.out'
    search_pattern = path.normpath(path.join(outfolder, outfile_pattern))
    previously_done = glob.glob(search_pattern)
    for i, eps in enumerate(ehg_epsilons):  # ehg_epsilons is a global constant
        for n in range(params['n_sim']):            
            assert type(eps) is float
            sim_name = 'eps{0:02d}sim{1:02d}'.format(i, n)
            outfilename = path.normpath(path.join(outfolder, sim_name + '.out'))
            if path.normpath(outfilename) in previously_done:
                print 'Skip simulation ' + sim_name + ' (already completed).'
                continue
            else:
                seed = params['rndseed'], n, int(eps * 10 ** 5)  # replicate(n)-specific seed

                simparams = params.copy()  # fresh copy for each replicate (IMPORTANT!)
                simparams.update(
                        prng=RandomState(seed=seed),
                        # get phi for this simulation (to be passed to `random_pairings` via `coresim`)
                        epsilon=eps,  # `phi` needs this!
                        sim_name=sim_name,
                        outfilename=outfilename,
                )
                yield simparams


def onesim(params):
    """Return None; setup and run one replicate
    (e.g., all daily iterations for 1 out of 100 replicates, for 1 value of epsilon)
    """
    try:
        sim_name = '{0}: {1}'.format(params['model_name'], params['sim_name'])
        print '\n\nBegin ' + sim_name + '\n' + pprint.pformat(params)
        # create a temp file to hold the results of one simulation
        outfolder = params['outfolder']
        tempfh = tempfile.NamedTemporaryFile(mode='w', suffix='.out', dir=outfolder, delete=False)
        assert params.get('fout', None) is None
        params['fout'] = tempfh  # used by `record_output`
        tempfname = tempfh.name  # we'll rename this if the simulation runs to completion
        # time the simulation
        t0 = time.clock()
        schedule = Scheduler(params=params)  # event scheduler (transmissions, deaths, dissolutions)
        pop_size = params['pop_size']
        nM = nF = pop_size // 2
        sim_days = params['sim_days']
        burn_days = params['burn_days']
        out_interval = params['out_interval']  # usu. 365 (i.e., once a year)
        p_nclients = params['p_nclients']
        p_nsexworkers = params['p_nsexworkers']
        p_nM_ART = params['p_nM_ART']
        p_nF_ART = params['p_nF_ART']
        p_miners = params['p_miners']
        p_PREP   = params['p_PREP']
        beta_PREP      = params["beta_PREP"]
        
        # type check parameters
        assert type(pop_size) is int
        assert type(sim_days) is int
        assert type(burn_days) is int
        assert type(out_interval) is int

        # create the disease for this model
        # note: don't add factory to model_params as it won't pickle -> MP problems
        params['Disease'] = stagedHIVfactory(durations=params['dur'],
                                             transM2F=params['beta_M2F'],
                                             transF2M=params['beta_F2M'],
                                             transM2F_ART=params['beta_M2F_ART'],
                                             transF2M_ART=params['beta_F2M_ART']
                                             )
        # create the phi for this model
        # note: don't add to model_params as it won't pickle -> MP problems
        # params['sim_phi'] = functools.partial(params['model_phi'], epsilon=params['epsilon'])
        params['sim_phi'] = lambda male, female: params['model_phi'](male, female, params['epsilon'])

        # counters used to tally incidence transmission by stage of infection
        # note: don't make this global (MP!)
        params['counters'] = counters = defaultdict(int)

        prep_params=params.copy()
        prep_params["Disease"] = stagedHIVfactory(durations=params['dur'],
                                             transM2F=prep_params['beta_PREP'],
                                             transF2M=prep_params['beta_PREP'],
                                             transM2F_ART=prep_params['beta_PREP'],
                                             transF2M_ART=prep_params['beta_PREP']
                                             )
        
        # prepare outfile by writing header
        header = 'nMinfect,nFinfect,'
        header += 'MPrimaryTrans,FPrimaryTrans,MAsymptomaticTrans,FAsymptomaticTrans,MSymptomaticTrans,FSymptomaticTrans,'
        header += 'MPships,,,,,,,,,,,,FPships,,,,,,,,,,,,iMPships,,,,,,,,,,,,iFPships,,,,,,,,,,,,'
        header += 'MPrimary,FPrimary'
        tempfh.write(header)

        # Prepare variables (note that `males` and `females` not reused across simulations)
        nminers     = int(p_miners * nM)
        nmpreps     = int(p_PREP * nM)
        nfpreps     = int(p_PREP  * nF)
        males       = list(Person(sex='M', registry=schedule, params=params) for i in range(nM-nminers-nmpreps))
        males      += list(Person03(sex='M', registry=schedule, params=params) for i in range(nminers))
        males      += list(Person04(sex='M', registry=schedule, params=prep_params) for i in range(nmpreps))
        females     = list(Person(sex='F', registry=schedule, params=params) for i in range(nF-nfpreps))
        females    += list(Person04(sex='F', registry=schedule, params=prep_params) for i in range(nfpreps))
        prng        = params['prng']
        # ai: look for data!
        # jk: vandepitte (2006); carael (2006) data
        nclients    = int(p_nclients * len(males))
        nsexworkers = int(p_nsexworkers * len(females))
        
        clients = prng.choice(males, nclients)
        sexworker = prng.choice(females, nsexworkers)
        # upgrade selected sexworkers from Person to Person01 (sexworkers)        
        for sw in sexworker:
            Person01.CastFromPerson(sw, True)

        # jk: look for more data
        nF_ART = int(p_nF_ART * len(females))
        nM_ART = int(p_nM_ART * len(males))
        f_ARTs = prng.choice(females, nF_ART)
        m_ARTs = prng.choice(males, nM_ART)

        # upgarde selected fameles_ARTs and males_ARTs to Person02
        for fa in f_ARTs:
            Person02.CastFromPerson(fa, True)
        for ma in m_ARTs:
            Person02.CastFromPerson(ma, True)

        # begin simulation loop /* Do the simulations */
        for day in range(sim_days + burn_days):
            logging.info('\nBEGIN ITERATION for day {0}.'.format(day))
            logging.debug(schedule.show_one_day(day))

            # Seed infections after burn_days have passed (ONCE)
            if (day == burn_days):
                assert schedule.count_scheduled_deaths() == 0
                diseases = sim.seed_infections(males, females, day, schedule=schedule, params=params)
                assert schedule.count_scheduled_deaths() == len(diseases)  # for now only have fatal disease
                assert all(deathday >= day for deathday in schedule.deaths)  # equal only possible on day==burn_days

            # run the core of the simulation (runs even during burn days)
            # params holds counters and fout
            schedule.coresim(
                    males=males,
                    females=females,
                    day=day,
                    params=params
            )

            # Record the output once a "period" (i.e., every out_interval days)
            # :note: this won't record after last year is run (it needs one more day to pass the test).
            #        We keep it this way just to match EHG.
            if (day >= burn_days and (day - burn_days) % out_interval == 0):
                print '.',
                outIndex = (day - burn_days) / out_interval
                # ai: recording strategy
                #    params holds counters and fout
                sim.record_output(males, females, params)
                # reset counters
                counters.clear()

            gc.collect()
        # END of simulation; just need to clean up: reset static elements
        # ai: mostly don't need to clean up in Python unless we intend to reuse objects
        # (and EHG don't even reuse the persons, so the rest of object creation is a trivial cost)
        # prepare classes for reuse
        schedule.clear_partnerships()  # clears the `partnerships` and `transmissions` multimaps
        schedule.deaths.clear()

        tempfh.close()
        dt = time.clock() - t0
        # since the simulation ran to completion, we can use the output file
        outfilename = params['outfilename']
        shutil.move(tempfname, outfilename)
        msg = """
	{0} completed successfully in {1:d} minutes.
	{0} output written to {2}.
	""".format(sim_name, int(dt) // 60, outfilename)
        logging.info(msg)
    except KeyboardInterrupt:
        logging.exception("Interrupted")
        raise
