"""
This code uses concurrency.py to add commercial sex workers, and people on ART
to the model. 
"""

import logging
import random 
from concurrency import * 

class Person01(Person): 
	"""
    Person01 extends concurrency.py person class 
    New Attribute: sexworker
    New Methods: add_commercial_parnership, remove_commercial_partnership 
    Overridden Methods: 
    """
	def __init__(self, sex, registry, params): 
		Person.__init__(self, "F", registry, params)
		self.sexworker = False


	def add_commercial_partnership(self, commercial_partnership): 
		self.commercial_partnerships.append(commercial_partnership)
		self.n_commercial_partners += 1


	def remove_commercial_partnership(self, commercial_partnership):
		"""Return None; remove a commercial partnership from self.partnerships list.
		"""
		logging.debug('ENTER: Person.remove_partnership')
		self.commercial_partnerships.remove(commercial_partnership)
		self.n_partners -= 1

class CommercialPartnership(Partnership):
	"""
	CommercialPartnership subclasses concurrency.py Partnership.
	New Attribute:
	New Methods:
	Overriden Methods:
	"""
	def set_tag(self,params):
		#allow only one primary partner (check union of partnerships sets)
		self.tag = 'commercial'
		self.sexfreq = 5  #Five times more than normal primary partnership frequency
	def pshipduration(self, params):
		return 1


# Do we need different proportion of male&female who are immune/receive ART.
# Arguably they may have different epsilon values but I would keep the same. 
# Set all transmission rates to zero or very low. Look for data on this. 
# Look for data on number of immune/who received ART. 

class Person02(Person):
	"""
	Person02 extends concurrency.py Person class
	new attribute: ART (meaning receiving ART treatment)
	new methods: add_ART_partnership, remove_ART_partnership
	overridden methods:
	"""
	def __init__(self, sex, registry, params, ART):
		Person.__init__(self, sex, registry, params)
		self.ART = False 

	def add_ART_partnership(self, ART_partnership):
		self.ART_partnership.append(ART_partnership)
		self.n_ART_partners += 1 

	def remove_ART_partnership(self, ART_partnership):
		logging.debug('Enter: Person.remove_partnership')
		self.ART_partnership.remove(ART_partnership)
		self.n_partners -= 1


class ARTPartnership(Partnership):
	"""
	ImmunePartnership subclasses concurrency.py Partnership.
	New Attribute:
	New Methods:
	Overriden Methods:
	"""
	def set_tag(self,params):
		#allow only one primary partner (check union of partnerships sets)
		self.tag = 'ART'
		self.sexfreq = 1  #primary partnership normal frequency

	def pshipduration(self, params):
		return 1
'''
#not needed
	def infect_ART(self, Disease, day, registry):  
		logging.debug('ENTER: Person02.infect')
		#change to infected state
		self.is_infected = True
		#fatal disease will schedule death (chkchk change to list?)
		disease = self.disease = Disease(self, day, registry=registry)
		# expose current partners to disease transmission
		for pship in self.partnerships:
			pship.expose_ART_transmission(day)  #chk!
		logging.debug('EXIT: Person.infect')
		return disease

class Infection01(Infection):
	"""
	Infection01 extends concurrency.py Infection class. 
	New Attribute: 
		New Methods: 
	Overriden Mehods: 
	"""
	_beta_M2F_ART = NotImplemented 
	_beta_F2M_ART = NotImplemented 


	def expose_ART(self, ART_partnership, day): #called by ARTPartnership::expose_ART_transmission
		male, female = partnership.partners
		assert male.is_infected != female.is_infected
		sexfreq = partnership.sexfreq
		if male.is_infected:
			beta = self._beta_M2F_ART
		if female.is_infected:
			beta = self._beta_F2M_ART
		beta_p_ART, beta_a_ART, beta_s_ART, beta_0 = ( sexfreq * beta_i for beta_i in beta )
		dur_p, dur_a, dur_s, dur_0 = self._durations
		candidate = eligDate = day  #candidate transmission day
		if(eligDate >= partnership.end_date):
			assert eligDate==partnership.end_date
			return

		stageEndDate = self.start_date 
		schedule = self.registry  
		for (dur_i,beta_i) in [(dur_p,beta_p_ART),(dur_a,beta_a_ART),(dur_s,beta_s_ART),]:
			stageEndDate += dur_i  
			if(eligDate < stageEndDate):
				if(beta_i > 0):  
					candidate = eligDate + partnership._prng.geometric(beta_i)  
					if(candidate <= stageEndDate and candidate <= partnership.end_date ):
						#schedule a transmission
						partnership.transmission_scheduled = True
						partnership.transmission_date = candidate
						schedule.register_transmission(partnership)
						logging.info('Infection.expose: Transmission scheduled')
						return
				eligDate = stageEndDate  #if beta_i=0 or candidate too far out
			if(eligDate >= partnership.end_date):
				return

		stageEndDate += dur_0
		if(eligDate < stageEndDate and beta_0 > 0):
			candidate = eligDate + partnership._prng.geometric(beta_0)
			if(candidate <= stageEndDate and candidate <= partnership.end_date):
				partnership.transmission_scheduled = True
				partnership.transmission_date = candidate
				schedule.register_transmission(partnership)
				return
		
'''
