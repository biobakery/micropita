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
import numpy as np
import os

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

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutputDataFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    c_fInvert = (args.fInvert == "True")

    logging.info("Start MicropitaPaperOverlapMatrix")
    print("Start MicropitaPaperOverlapMatrix")

    #Read and parse predicted data
    #Get predicted
    lsPredictedFileContents = []
    with open(args.strSelectionFile, 'r') as f:
        lsPredictedFileContents = f.read()
        lsPredictedFileContents = filter(None,re.split("\n",lsPredictedFileContents))
        f.close()

    #Parse selection Predicted
    dictPredicted = dict()
    for strSelectionMethod in lsPredictedFileContents:
        #Get method name
        astrSelectionMethod = strSelectionMethod.split(Constants.COLON)
        sCurSelectionMethodName = astrSelectionMethod[0]

        #Parse samples selected by the method
        astrSelectedSamples = astrSelectionMethod[1].split(Constants.COMMA)
        astrSelectedSamples = [strSelectedSample.strip() for strSelectedSample in astrSelectedSamples]

        #Add to dict
        dictPredicted[sCurSelectionMethodName]=astrSelectedSamples
    
    #Get labels
    lsLabels = intersect(set(dictPredicted.keys()))
    if len(lsLabels) == 0:
        logging.error("".join(["MicropitaPaperOverlapMatrix. No actual and predicted labels were in common so an overlap matrix could not be generated. Actual Labels=",
                               str(dictPredicted.keys())," Predicted Labels=",str(dictActual.keys())]))
    elif len(lsLabels) < 2:
        logging.error("".join(["MicropitaPaperOverlapMatrix. Only one actual and predicted labels were in common so an overlap matrix could not be generated. Actual Labels=",
                               str(dictPredicted.keys())," Predicted Labels=",str(dictActual.keys())]))

    #make a y axis
    yLabels = lsLabels
    if fFlipYLabels:
        yLabels.reverse()

    #Setup data for overlap matrix
    liOverlap = []
    for xLabel in lsLabels:
        setPredicted = set(dictPredicted[xLabel])
        iPredictedPosition = None
        setPotentialSharedSamples = {}
        for iindex, yLabel in enumerate(yLabels):
            if yLabel == xLabel:
                iPredictedPosition = iindex
                liOverlap.append(None)
            else:
                setOtherSamples = set(dictPredict[yLabel])
                liOverlap.append(len(setPredicted & setOtherSamples))
                setPotentialSharedSamples = setOtherSamples | setPotentialSharedSamples
        liOverlap[iPredictPosition] = setPredicted - setPotentialSharedSamples

    #Plot confusion matrix
    funcPlotMatrix(npMatrix=np.array(liConfusion), lsLabels=lsLabels, strOutputFigurePath=args.strOutputFigure, strXTitle="Predicted", strYTitle="Actual", fFlipYLabels=True, fInvert=c_fInvert):

    print("Stop MicropitaPaperOverlapMatrix")
    logging.info("Stop MicropitaPaperOverlapMatrix")

if __name__ == "__main__":
    _main( )
