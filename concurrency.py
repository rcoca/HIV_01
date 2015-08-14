"""Provides key classes for concurrency simulations.

This code is replication of Sawers, Isaac, Stillwaggon (2011)
"HIV and concurrent sexual partnerships: modelling the role of coital dilution"

"""

import logging
from collections import defaultdict
schedule = dict(
	deaths = defaultdict(list), # multimap for scheduling deaths: map date -> list of Person
	transmissions = defaultdict(list),  # map **upcoming** transmission date to list of Partnership
	dissolutions = defaultdict(list)  # map each **end** date list of Partnership
	)
	


def stagedHIVfactory(durations, transM2F, transF2M):
	"""Return class for staged HIV infection.
	"""
	class HIV(Infection):
		_durations = tuple(durations)
		_duration = sum(_durations)
		_beta_M2F = tuple(transM2F)
		_beta_F2M = tuple(transF2M)
		assert len(_durations)==len(_beta_F2M)==len(_beta_M2F)
	return HIV

class Infection(object):
	"""
	Provides an **ABSTRACT** infection class.
	DO NOT USE THIS CLASS.
	Instead, use the class factory `stagedHIVfactory`.
	This class is missing the following class variables:
	_durations : duration of each infection stage
	_duration : total duration of the infection
	_beta_M2F  : staged transmission rates for infected males
	_beta_F2M  : staged transmission rates for infected females
	"""
	_durations = NotImplemented
	_beta_M2F = NotImplemented
	_beta_F2M = NotImplemented
	def __init__(self, host, day):
		"""Return None.
		Record start date and host.
		Schedule host death.
		"""
		self._host = host
		self.start_date = day
		self.end_date = day + self.duration
		self.schedule_death()
	def schedule_death(self):
		dod = self.end_date
		host = self._host
		schedule['deaths'][dod].append(host)
	def stage(self, day):
		"""Return str, the name of the current infection stage.
		(Utilized to tally transmissions.)
		Called by `Person.get_HIV_stage`,
		which used to hold this code.
		"""
		assert self._host.is_infected
		dur_p, dur_a, dur_s, dur_0 = self._durations
		days_infected = day - self.start_date
		if(days_infected <= dur_p): # primary infection
			stage = 'primary'
		elif( days_infected <= (dur_p + dur_a) ): # asymptomatic infection
			stage = 'asymptomatic'
		elif( days_infected <= (dur_p + dur_a + dur_s) ): # sympotmatic infection
			stage = 'sympotmatic'
		else: # pre-death, no transmission
			stage = 'dying'
		return stage
	def expose(self, partnership, day): #called by Partnership::expose_transmission
		"""Return None; expose a *partnership* and possibly schedule transmission event.
		Expose partnership to infection by sampling a candidate transmission date during
		each stage based on the transmission rate during that stage, checking to see if
		the candidate date is before the end of that stage and before the end of the partnership.
		if so, transmission is scheduled for that date.  if not, proceed to the next stage.
		If transmission is not scheduled during any stage, no transmission occurs during the
		partnership.

		Called by `Partnership.expose_transmission`, which is called upon partnership
		formation or upon infection.

		:comment:
			We only schedule a transmission if exactly one of the two is infected,
			but with concurrency another partner might still infect the uninfected partner
			before this scheduled transmission happens.
		"""
		male, female = partnership.partners
		assert male.is_infected != female.is_infected
		sexfreq = partnership.sexfreq
		if male.is_infected:
			beta = self._beta_M2F
		if female.is_infected:
			beta = self._beta_F2M
		beta_p, beta_a, beta_s, beta_0 = ( sexfreq * beta_i for beta_i in beta )
		dur_p, dur_a, dur_s, dur_0 = self._durations
		candidate = eligDate = day  #candidate transmission day
		if(eligDate >= partnership.end_date):
			#ai: partnership ends beforehand; I think this only happens if eligDate==partnership.end_date
			# (because o/w partnership shd have dissolved) chk
			assert eligDate==partnership.end_date
			return

		stageEndDate = self.start_date #initialize to date of infection
		#ai: small change from EHG: loop over common code
		#ai: primary, asymptomatic, symptomatic infection stage
		for (dur_i,beta_i) in [(dur_p,beta_p),(dur_a,beta_a),(dur_s,beta_s),]:
			stageEndDate += dur_i  # infection period (increment then compare)
			if(eligDate < stageEndDate):
				if(beta_i > 0):  #daily transmission probability during this stage
					candidate = eligDate + partnership._prng.geometric(beta_i)  #chk
					if(candidate <= stageEndDate and candidate <= partnership.end_date ):
						#schedule a transmission
						partnership.transmission_scheduled = True
						partnership.transmission_date = candidate
						schedule['transmissions'][partnership.transmission_date].append(partnership)
						logging.info('Infection.expose: Transmission scheduled')
						return
				eligDate = stageEndDate  #if beta_i=0 or candidate too far out
			if(eligDate >= partnership.end_date):
				return

		#ai: final days (could reuse the above code)
		stageEndDate += dur_0
		if(eligDate < stageEndDate and beta_0 > 0):
			candidate = eligDate + partnership._prng.geometric(beta_0)
			if(candidate <= stageEndDate and candidate <= partnership.end_date):
				partnership.transmission_scheduled = True
				partnership.transmission_date = candidate
				schedule['transmissions'][partnership.transmission_date].append(partnership)
				return
	@property  #ai: really a class property
	def duration(self):
		return self._duration

class Person(object):
	"""Provides a basic person class.
	A person has a sex and a set of partnerships.
	A person may have a disease and an associated date of death.
	"""
	def __init__(self, sex, params):
		"""Return None;  initialize a Person.
		"""
		self._sex = sex  #'M' or 'F' (EHG use integer: 1=male, 2=female)
		self._prng = params['prng']
		self.reset()

	def reset(self):
		"""Return None.
		Prepare a Person instance for reuse after death.
		"""
		self.n_partners = 0
		self.partnerships = list() #use list (not set) to ensure order of iteration
		self.is_infected = False  #ai: rendundant? (Just use doi)
		self.disease = None
		self.iDOD = None   #ai: impose max (adult) lifespan? chk (need day!)

	def get_HIV_stage(self, day):
		"""Return str, the current stage of infection
		"""
		return self.disease.stage(day)

	def add_partnership(self, partnership):
		"""Return None; add a new partnership to self.partnerships list.
		"""
		self.partnerships.append(partnership)
		self.n_partners += 1
		assert (self.n_partners == len(self.partnerships) == len(set(self.partners)))
		return

	def remove_partnership(self, partnership):
		"""Return None; remove a partnership from self.partnerships list.
		"""
		logging.debug('ENTER: Person.remove_partnership')
		self.partnerships.remove(partnership)
		self.n_partners -= 1
		assert (self.n_partners == len(self.partnerships) == len(set(self.partners)))
		logging.debug('EXIT: Person.remove_partnership')

	def infect(self, Disease, day):  #compare EHG's HIVInfect
		"""Return None,
		infect an individual with HIV, schedule death, and expose other partners to transmission 
		"""
		logging.debug('ENTER: Person.infect')
		#change to infected state
		self.is_infected = True
		self.disease = Disease(self, day) #fatal disease will schedule death
		# expose current partners to disease transmission
		for pship in self.partnerships:
			pship.expose_transmission(day)
		logging.debug('EXIT: Person.infect')

	def HIVSeedInfect(self, Disease, day):  #compare EHG's HIVSeedInfect in person.cpp
		"""Return None.  seed infection
		:comment: There is a low probably of an offset of 0,
			which would mean the day of death is the day seeded.
			This is to match EHG and is only mildly annoying.
		"""
		logging.debug('ENTER: Person.HIVSeedInfect')
		#change to infected state
		self.is_infected = True
		duration = Disease._duration   #chk works but ugly
		# schedule current date uniformly during the duration of infection
		offset = self._prng.randint(0, duration)  #same as random.randrange
		doi = day - (duration - offset)
		self.disease = Disease(self, doi) #fatal disease will schedule death
		for pship in self.partnerships:
			pship.expose_transmission(day)  #odd to use day instead of DOI, but simpler, and matches EHG
		logging.debug('EXIT: Person.HIVSeedInfect')
	
	def has_primary(self):
		return any(p.tag=='primary' for p in self.partnerships)

	def die(self):  #compare EHG's `Kill`
		"""Return None; "kills" the individual (i.e., reset, so population is constant).
		Terminate partnerships and reset individual to unpartnered susceptible.
		"""
		logging.info('ENTER: Person.die')
		#terminate all partnerships
		for pship in tuple(self.partnerships):
			pship.DeletePartnership()
		assert len(self.partnerships)==0
		self.reset()  #reborn anew, uninfected, same sex
		logging.info('EXIT: Person.die')


	###################### ai: convenience additions ##############################
	@property
	def sex(self):
		return self._sex
	@property
	def date_of_infection(self):
		return getattr(self.disease, 'start_date', None)

	@property
	def partners(self):
		return (indiv for pship in self.partnerships for indiv in pship.partners if indiv is not self)
	@staticmethod
	def count_scheduled_deaths():
		return sum(len(lst) for lst in schedule['deaths'].values())

  
class Partnership(object):
	"""Provides a paired heterosexual "partnership" (i.e., sexual relationship).
	A Partnership has a predetermined random duration and schedules HIV transmission events.
	"""
	# initialize partnership counter
	n_partnerships = 0
	def __init__(self, male, female, start_date, params):
		"""Return None. Initialize partnership,
		schedule partnership dissolution,
		schedule disease transmission.
		"""
		logging.debug('ENTER: Partnership.__init__')
		assert (male.sex=='M' and female.sex=='F')
		self._mf = male, female
		#allow only one primary partner (check union of partnerships sets)
		has_primary = any(p.tag=='primary' for p in (male.partnerships + female.partnerships))
		if has_primary:
			self.tag = 'secondary'
			#set the relative frequency for non 'primary' partners 
			self.sexfreq = params.get('secondary_sexfreq', 1)
		else:
			self.tag = 'primary'
			self.sexfreq = 1  #primary partnership normal frequency
		self.start_date = start_date
		self._prng = prng = params['prng']
		self.end_date = start_date + prng.geometric(params['sigma'])  #add partnership duration (min of 1)
		self.transmission_scheduled = False	# T/F whether HIV transmission is scheduled for the partnership
		self.transmission_date = None      #date for scheduled transmission
		# add new partnership to each partner's partnership list
		male.add_partnership(self)
		female.add_partnership(self)
		# add partnership to partnership multimap with end date
		schedule['dissolutions'][self.end_date].append(self)
		self.__class__.n_partnerships += 1
		# expose partnership to infection transmission (if exactly one is infected)
		self.expose_transmission(start_date)
		logging.debug('EXIT: Partnership.__init__')
		
	def expose_transmission(self, day):
		"""Return None; expose and possibly schedule partnership for transmission event.
		Called by `Partnership.__init__`, so upon formation a partnership is
		checked for HIV exposure.

		:comment:
			We only schedule a transmission if exactly one of the two is infected,
			but with concurrency another partner might still infect the uninfected partner
			before this scheduled transmission happens.
		"""
		logging.debug('ENTER: Partnership.expose_transmission')
		male, female = self._mf

		if male.is_infected and not female.is_infected:
			infected = male
		elif not male.is_infected and female.is_infected:
			infected = female
		else:                #partners concordant -> exit the function
			return

		infected.disease.expose(self, day)  #expose the partnership to infection

		logging.debug('EXIT: Partnership.expose_transmission')

	def transmit(self, Disease, day): #compare EHG's HIVTransmission in partnership.cpp
		"""Return Person; transmit the HIV infection as **previously** scheduled
		from the infected partner to the uninfected partner.
		(If both are infected, then the originally uninfected partner got
		infected by someone else since this partnership formed.)
		"""
		logging.debug('ENTER: Partnership.transmit_HIV')
		assert self.transmission_scheduled, 'Only scheduled transmissions shd call this'
		assert self.transmission_date==day, 'Only scheduled transmissions shd call this'
		male, female = self._mf
		transmitter = None
		if(male.is_infected and not female.is_infected): # M --> F transmission
			transmitter = male
			female.infect(Disease, day)   # infect the female
		elif(not male.is_infected and female.is_infected): # F --> M transmission
			transmitter = female
			male.infect(Disease, day)   # infect the male
		else: #partner who was uninfected when partnership formed subsequently became infected
			assert (male.is_infected and female.is_infected)
		self.transmission_scheduled = False
		logging.debug('EXIT: Partnership.transmit_HIV')
		return transmitter

	def DeletePartnership(self):
		"""Return None; delete a partnership from ``partnerships`` multimap
		*before* scheduled end date (e.g., when partner dies from HIV).
		Also remove any scheduled transmissions.
		"""
		schedule['dissolutions'][self.end_date].remove(self)
		self.__class__.n_partnerships -= 1
		#destroy the partnership (match EHG's destructor for the partnership class)
		male, female = self._mf
		male.remove_partnership(self) #removes from set *and* decrements this person's count
		female.remove_partnership(self)
		# check if a future transmission is scheduled, and if so remove it
		if(self.transmission_scheduled):
			schedule['transmissions'][self.transmission_date].remove(self)

	#properties
	@property
	def partners(self):
		return self._mf
		
	#chk this may work very differently than the EHG code, but that shd't matter, it's just cleanup
	@classmethod
	def clear(cls): # delete all partnerships (prepare for next simulation); compare EHG's DumpPartnerships
		for lst in schedule['dissolutions'].values():
			for pship in lst[:]:
				pship.DeletePartnership()
		assert cls.n_partnerships == 0
		schedule['dissolutions'].clear()
		schedule['transmissions'].clear()

	######################  ai: new convenience methods #######################################
	@classmethod
	def count_partnerships(cls):
		return sum(len(lst) for lst in schedule['dissolutions'].values())

	@classmethod
	def count_transmissions(cls):
		return sum(len(lst) for lst in schedule['transmissions'].values())

