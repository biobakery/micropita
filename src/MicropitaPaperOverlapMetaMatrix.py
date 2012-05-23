#!/usr/bin/env python

"""
Author: Timothy Tickle
Description: Class to generate an overlap matrix based on several projects of micropita selection (Metaversion)
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
import logging
from MicroPITA import MicroPITA
import numpy as np
import os
from PlotMatrix import PlotMatrix
import re

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperOverlapMetaMatrix.py", description = """Creates an overlap matrix from multiple microPITA output.""" )
#Arguments
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
#Outputfile
argp.add_argument( "strOutputFigure", metavar = "output.txt", help = Constants_Arguments.c_genericOutputFigureFileHelp )
#Selection parameter
argp.add_argument("strSelectionMethods", metavar = "Selection_Methods", help = Constants_Arguments.c_strSelectionMethodsCommaHelp)
#Micropita selection file (predicted samples)
argp.add_argument( "strProjects", action="store", metavar = "Select_files", help = Constants_Arguments.c_strMicropitaProjectsHelp, nargs="+")

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutputFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    c_fInvert = (args.fInvert == "True")
    c_fFlipYLabels = False

    logging.info("Start MicropitaPaperOverlapMetaMatrix")

    #Split the selection methods into a set
    setSelectionMethods = set([sMethod for sMethod in args.strSelectionMethods.split(Constants.COMMA)])

    #Holds the universal y axis
    lsYLabels = None

    #Holds the universal x axis
    lsXLabels = None

    #Universal overlap counts
    npUniversalOverlap = None

    for strSelectionFile in args.strProjects:

        #Read and parse predicted data into a dictionary
        dictPredicted = MicroPITA.funcReadSelectionFileToDictionary(strSelectionFile)

        #Get labels
        lsLabels = dictPredicted.keys()

        #Subset them to just what is interested in looking at.
        lsLabels = list(set(lsLabels) & setSelectionMethods)
        iLabelLength = len(lsLabels)

        #Make sure there are labels to look at
        if iLabelLength == 0:
            logging.error("".join(["MicropitaPaperOverlapMetaMatrix. No labels were found so an overlap matrix could not be generated. Predicted Labels=",
                               str(dictPredicted.keys())]))
        elif iLabelLength < 2:
            logging.error("".join(["MicropitaPaperOverlapMetaMatrix. Only one label was found so an overlap matrix could not be generated. Predicted Labels=",
                               str(dictPredicted.keys())]))

        #Either set up universal axes or check to make sure the current one has the same composition as the incoming file
        if not lsYLabels:
            lsXLabels = lsLabels
            lsYLabels = list(lsXLabels)
            npUniversalOverlap = np.zeros(len(lsXLabels)*len(lsXLabels))
            if c_fFlipYLabels:
                lsYLabels.reverse()
        else:
            if not (len(set(lsLabels) ^ set(lsYLabels)) == 0):
                logging.error("".join(["MicropitaPaperOverlapMetaMatrix. Labels from one study to another ar not consistent. Labels before=",str(lsXLabels)," Labels after=",str(lsLabels)]))

        #Setup data for overlap matrix
        liOverlap = []
        for xLabel in lsXLabels:
            setPredicted = set(dictPredicted[xLabel])
            iPredictedPosition = None
            iArrayLengthTraveled = 0
            setPotentialSharedSamples = set()
            for iindex, yLabel in enumerate(lsYLabels):
                if yLabel == xLabel:
                    iPredictedPosition = iindex + iArrayLengthTraveled
                    liOverlap.append(None)
                else:
                    setOtherSamples = set(dictPredicted[yLabel])
                    liOverlap.append(len(setPredicted & setOtherSamples))
                    setPotentialSharedSamples = setOtherSamples | setPotentialSharedSamples
                iArrayLengthTraveled = iArrayLengthTraveled + iLabelLength
            liOverlap[iPredictedPosition] = len(setPredicted - setPotentialSharedSamples)

        #Cumulate results
        npUniversalOverlap = npUniversalOverlap + np.array(liOverlap)

    #Plot matrix
    npUniversalOverlap = npUniversalOverlap / float(len(args.strProjects))
    npUniversalOverlap.shape=iLabelLength,iLabelLength
    PlotMatrix.funcPlotMatrix(npMatrix=npUniversalOverlap, lsXLabels=lsYLabels, strOutputFigurePath=args.strOutputFigure, strXTitle="Predicted", strYTitle="Actual", fFlipYLabels=c_fFlipYLabels, fInvert=c_fInvert)

    logging.info("Stop MicropitaPaperOverlapMetaMatrix")

if __name__ == "__main__":
    _main( )
