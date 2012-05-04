#!/usr/bin/env python

"""
Author: Timothy Tickle
Description: Class to Generate InSilico Datasets for micropita
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libaries
import argparse
from Constants import Constants
from Constants_Arguments import Constants_Arguments
from CommandLine import CommandLine
import logging
import numpy as np
import os
import random
from Utility_Data import Utility_Data

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperConstructDataSets.py", description = """Creates In-silico data sets.""" )
#Arguments
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
#Outputfile
argp.add_argument( "strOutputDataFile", metavar = "output.txt", help = Constants_Arguments.c_genericOutputDataFileHelp )
#Data key to determine which data set to generate
argp.add_argument( "strDataKey", metavar = "Data Set key", help = Constants_Arguments.c_dataSetKeyHelp )
#Actual sample selection
argp.add_argument( "strActualFile", action="store", metavar = "Actual_file", help = Constants_Arguments.c_strMicropitaActualFileHelp, default = None, nargs = '?')

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    print("Main")

    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutputDataFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperConstructDataSets")
    print("Start MicropitaPaperConstructDataSets")

    if args.strDataKey =="Diversity":
        print("Diversity")
        Utility_Data.generateDiversityAbundanceTable(args.strOutputDataFile)
    elif args.strDataKey =="Unbalanced":
        print("Unbalanced")
        Utility_Data.generateAbundanceTable(args.strOutputDataFile,args.strActualFile)
    else:
        logging.error("".join(["The following key is not valid. Did not generate data set. Key=",args.strDataKey]))

    print("Stop MicropitaPaperConstructDataSets")
    logging.info("Stop MicropitaPaperConstructDataSets")

if __name__ == "__main__":
    _main( )
