#######################################################
# Author: Timothy Tickle
# Description: Class to manage all tests in the regression suite
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libraries
import unittest

#Import test libraries
import AbundanceTableTest
import CladogramTest
import DiversityTest
import MLPYDistanceAdaptorTest
#import PCoATest
#import MicroPITATest
#import SVMTest
##import TimerTest
import Utility_MathTest

suite = unittest.TestSuite()
suite.addTest(AbundanceTableTest.suite())
suite.addTest(CladogramTest.suite())
suite.addTest(DiversityTest.suite())
#suite.addTest(MLPYDistanceAdaptorTest.suite())
#suite.addTest(PCoATest.suite())
#suite.addTest(MicroPITATest.suite())
#suite.addTest(SVMTest.suite())
##suite.addTest(TimerTest.suite())
suite.addTest(Utility_MathTest.suite())

runner = unittest.TextTestRunner()
runner.run(suite)
