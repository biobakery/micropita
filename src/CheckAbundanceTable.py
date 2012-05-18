#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Script checks abundance table files befroe reading in 
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

from AbundanceTable import AbundanceTable
import argparse
from Constants import Constants
from Constants_Arguments import Constants_Arguments
import logging
import os

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "CheckAbundanceTable.py", 
    description = """Checks the file that will generate an abundance table.""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataName, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)

#Data file
argp.add_argument( "strAbundanceFile", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument( "strOutFile", metavar = "OutputFile", help = Constants_Arguments.c_strOutputFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start CheckAbundanceTable")
    logging.info("CheckAbundanceTable. The following arguments were passed.")
    logging.info(str(args))

    if not args.sLastMetadataName:
      logging.error("CheckAbundanceTabe::Did not receive a value for sLastMetadataName. ")
      logging.info("Stop CheckAbundanceTable")
      #TODO does this stop further analysis, not if there is an already checked file....
      return False

    #Read abundance file
    #Abundance table object to read in and manage data
    totalData = AbundanceTable.funcCheckRawDataFile(strReadDataFileName=args.strAbundanceFile,
                                                    sLastMetadataName=args.sLastMetadataName, strOutputFileName=args.strOutFile)
    if not totalData:
      logging.error("CheckAbundanceTabe::Received error when AbundanceTable.funcCheckRawDataFile was called")
    logging.info("Stop CheckAbundanceTable")

if __name__ == "__main__":
    _main( )

