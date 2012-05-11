#!/usr/bin/env python

"""
Author: Timothy Tickle
Description: Class to generate an overlap matrix based on micropita selection
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
argp = argparse.ArgumentParser( prog = "MicropitaPaperOverlapMatrix.py", description = """Creates an overlap matrix from microPITA output.""" )
#Arguments
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
#Micropita selection file (predicted samples)
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
#Outputfile
argp.add_argument( "strOutputFigure", metavar = "output.txt", help = Constants_Arguments.c_genericOutputFigureFileHelp )
#Selection parameter
argp.add_argument("strSelectionMethods", metavar = "Selection_Methods", help = Constants_Arguments.c_strSelectionMethodsHelp, nargs="*")

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

    logging.info("Start MicropitaPaperOverlapMatrix")

    #Read and parse predicted data into a dictionary
    dictPredicted = MicroPITA.funcReadSelectionFileToDictionary(args.strSelectionFile)

    #Get labels
    lsLabels = dictPredicted.keys()

    #Subset them to just what is interested in looking at.
    lsLabels = list(set(lsLabels) & set(args.strSelectionMethods))
    iLabelLength = len(lsLabels)

    if iLabelLength == 0:
        logging.error("".join(["MicropitaPaperOverlapMatrix. No labels were found so an overlap matrix could not be generated. Predicted Labels=",
                               str(dictPredicted.keys())]))
    elif iLabelLength < 2:
        logging.error("".join(["MicropitaPaperOverlapMatrix. Only one label was found so an overlap matrix could not be generated. Predicted Labels=",
                               str(dictPredicted.keys())]))

    #make a y axis
    yLabels = lsLabels
    if c_fFlipYLabels:
        yLabels.reverse()

    #Setup data for overlap matrix
    liOverlap = []
    for xLabel in lsLabels:
        setPredicted = set(dictPredicted[xLabel])
        iPredictedPosition = None
        iArrayLengthTraveled = 0
        setPotentialSharedSamples = set()
        for iindex, yLabel in enumerate(yLabels):
            if yLabel == xLabel:
                iPredictedPosition = iindex + iArrayLengthTraveled
                liOverlap.append(None)
            else:
                setOtherSamples = set(dictPredicted[yLabel])
                liOverlap.append(len(setPredicted & setOtherSamples))
                setPotentialSharedSamples = setOtherSamples | setPotentialSharedSamples
            iArrayLengthTraveled = iArrayLengthTraveled + iLabelLength
        liOverlap[iPredictedPosition] = len(setPredicted - setPotentialSharedSamples)

    #Plot confusion matrix
    liOverlap = np.array(liOverlap)
    liOverlap.shape=iLabelLength,iLabelLength
    PlotMatrix.funcPlotMatrix(npMatrix=liOverlap, lsLabels=lsLabels, strOutputFigurePath=args.strOutputFigure, strXTitle="Predicted", strYTitle="Actual", fFlipYLabels=c_fFlipYLabels, fInvert=c_fInvert)

    logging.info("Stop MicropitaPaperOverlapMatrix")

if __name__ == "__main__":
    _main( )
