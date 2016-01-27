"""
This code uses concurrency.py to add commercial sex workers, and people on ART
to the model. 
"""

import logging
#import random
from concurrency import Person,Partnership
from concurrency01 import Person01,Person02

class Person03(Person):
    """
    This class represents the miner class, with higher sexfreq (3) only active ~ 200 days out of 365/year
    static data members:
    activeRange  = (0,200) : the range of days throughout the year when they can form partnerships
    sexfreq      = 3 : sex frequency parameter, native to the individual type, for all instances
    The type Person03 is used in scheduler to discriminate behavior in forming parnerships (isinstance)
    """
    activeRange  = (0,200)
    sexfreq      = 3
    def __init__(self, sex, registry, params):
        Person.__init__(self, sex, registry, params)
    
class Person04(Person):
    """
    This class models individuals taking PREP drugs.
    Its only influence on the simulation is the beta parameters (set to 0,0,0,0).
    """
    def __init__(self, sex, registry, params):
        Person.__init__(self, sex, registry, params)

