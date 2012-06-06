#!/usr/bin/env python

"""
Author: Timothy Tickle
Description: Class to generate a confusion matrix based on micropita selection in multiple studies/samples
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

def funcReadActualClassFile(strInputFile):

    #Read in file
    strSelection = ""
    with open(strInputFile,'r') as fHndlInput:
        strSelection = fHndlInput.read()
    fHndlInput.close()

    #Break down content into lines
    lsSelection = filter(None,strSelection.split(Constants.ENDLINE))
    strSelection = None

    #Store information by classes
    #{"ClassName":["sample"...]
    #Will eventually change the list to sets
    dictClasses = {}
    iClassPrefixLength = len(Constants.c_strClassPrefix)
    for sSelection in lsSelection:
        if sSelection[0:iClassPrefixLength] == Constants.c_strClassPrefix:

            lsSelectionPieces = sSelection.split(Constants.COLON)
            if not len(lsSelectionPieces) == 2:
                logging.error("A line in the actual file was not formated corrected. Expected two pieces divided by a :. (Parsing class data)")
                return False
            sClass = lsSelectionPieces[0][iClassPrefixLength:].strip()
            dictClasses[sClass]= [s.strip() for s in lsSelectionPieces[1].split(Constants.COMMA)]

    #Get class count and make sure class information was parsed
    iClassCount = len(dictClasses.keys())
    if iClassCount < 1:
        logging.error("Did not read a class from the actual file.")
        return False

    #Given we have class information change the lists to sets
    for sClass in dictClasses:
        dictClasses[sClass] = set(dictClasses[sClass])

    #Will hold the amount of overlapp in sample groups (how many classes they come from)
    dictOverlap = dict()

    #Get correct sample selection
    dictActualClassSamples = dict()
    for sSelection in lsSelection:
        if not sSelection[0:iClassPrefixLength] == Constants.c_strClassPrefix:

            lsSelectionPieces = sSelection.split(Constants.COLON)
            if not len(lsSelectionPieces) == 2:
                logging.error("A line in the actual file was not formated corrected. Expected two pieces divided by a :. (Parsing sample data)")
                return False

            sSelectionMethodology = lsSelectionPieces[0].strip()
            dictActualClassSamples[sSelectionMethodology] = []
            for sSampleGrouping in lsSelectionPieces[1].split(Constants.COMMA):
                sSampleGrouping = sSampleGrouping.strip()
                if len(sSampleGrouping) > iClassPrefixLength:
                    if sSampleGrouping[0:iClassPrefixLength] == Constants.c_strClassPrefix:
                        dictActualClassSamples[sSelectionMethodology].extend(dictClasses[sSampleGrouping[iClassPrefixLength:]])
                        lcurSelection = dictOverlap.get(sSelectionMethodology,[])
                        lcurSelection.append(sSampleGrouping[iClassPrefixLength:])
                        dictOverlap[sSelectionMethodology] = lcurSelection
                    else:
                        dictActualClassSamples[sSelectionMethodology].append(sSampleGrouping)
                else:
                    dictActualClassSamples[sSelectionMethodology].append(sSampleGrouping)

    #Get sample grouping count and make sure information was parsed
    iGroupCount = len(dictActualClassSamples.keys())
    if iGroupCount < 1:
        logging.error("Did not read sample groupings from the actual file.")
        return False

    #Given we have sample grouping information change the lists to sets
    for sGroup in dictActualClassSamples:
        dictActualClassSamples[sGroup] = set(dictActualClassSamples[sGroup])

    #return the count of classes and the sample groupings per methodlogy
    #[int, {"sClass":[(sSample,sSample,sSample...)]}, {"sMethodology":[(sSample,sSample,sSample...)]}]
    return [iClassCount, dictActualClassSamples, dictClasses, dictOverlap]

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperConfusionMetaMatrix.py", description = """Creates a confusion matrix from multiple microPITA output.""" )
#Arguments
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
#Actual sample selection
argp.add_argument( "strActualFile", action="store", metavar = "Actual_file", help = Constants_Arguments.c_strMicropitaActualFileHelp)
#Outputfile
argp.add_argument( "strOutputFigure", metavar = "output.txt", help = Constants_Arguments.c_genericOutputFigureFileHelp )
#Selection parameter
argp.add_argument("strSelectionMethods", metavar = "Selection_Methods", help = Constants_Arguments.c_strSelectionMethodsCommaHelp)
#Micropita selection file (predicted samples)
argp.add_argument( "lsProjects", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp, nargs="+")


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

    logging.info("Start MicropitaPaperConfusionMetaMatrix")

    #Read and parse actual data into a dictionary
    lsActualData = funcReadActualClassFile(args.strActualFile)
    if not lsActualData:
        logging.error("Did not receive data from parsing the actual data file. Returning without generating a plot.")
        return False
    iClassCount,dictActual,dictClasses,dictOverlapMeasure = lsActualData

    #Split the selection methods into a set
    setSelectionMethods = set([sMethod for sMethod in args.strSelectionMethods.split(Constants.COMMA)])

    #Holds the universal y axis
    lsYLabels = None

    #Holds the universal x axis
    lsXLabels = None

    #Universal confusion counts
    npUniversalConfusion = None

    for strSelectionFile in args.lsProjects:

        #Read and parse predicted data into a dictionary
        dictPredicted = MicroPITA.funcReadSelectionFileToDictionary(strSelectionFile)

        logging.debug("dictPredicted")
        logging.debug(dictPredicted)
        logging.debug("dictActual")
        logging.debug(dictActual)

        #Get labels, make sure the labels are in common between actual predicted and those indicated to be plotted (intersection of all three classes)
        lsLabels = list(set(dictPredicted.keys()) & set(dictActual.keys()) & setSelectionMethods)
        iLabelLength = len(lsLabels)
        if iLabelLength == 0:
            logging.error("".join(["MicropitaPaperConfusionMetaMatrix. No actual and predicted labels were in common so a confusion matrix could not be generated. Actual Labels=",
                               str(dictPredicted.keys())," Predicted Labels=",str(dictActual.keys())]))
        elif iLabelLength < 2:
            logging.error("".join(["MicropitaPaperConfusionMetaMatrix. Only one actual and predicted labels were in common so a confusion matrix could not be generated. Actual Labels=",
                               str(dictPredicted.keys())," Predicted Labels=",str(dictActual.keys())]))

        #Either set up universal axes or check to make sure the current one has the same composition as the incoming file
        if not lsYLabels:
            lsXLabels = lsLabels
            lsYLabels = list(lsXLabels)
            npUniversalConfusion = np.zeros(len(lsXLabels)*len(lsXLabels))
            if c_fFlipYLabels:
                lsYLabels.reverse()

        dictMisclassifiedClasses = {}
        dictMisclassfiedSamples = {}

        #Setup data for confusion matrix
        liConfusion = []
        #Go through each label
        #Get the predicted and actual samples for the label
        #For all classes but representative:
        #Indicate which are correct by intersecting the predicted and actual results
        #In missclassification
        for yLabel in lsYLabels:
            setYActual = set(dictActual[yLabel])
            for xLabel in lsXLabels:
                setXPredicted = set(dictPredicted[xLabel])
                setXPredictedPossible = set(dictActual[xLabel])
                iCountPredictedSamples = len(setXPredicted)
                iMaxSamplesPerCategory = int(round(iCountPredictedSamples/float(iClassCount)))

                #Representative techniques are different because they do not have specific samples to
                #to select in a priori classes. Correct selection is defined by the shape of the data and how many data are sampled.
                #Given this it is impossible to know before hand what samples should be selected.
                #It is possible to know that the samples should be equally sampled between sample classes (as long as the sample
                #classes are truly different and have some defining characteristic that the algorithm will use to recognize.
                #So here we evaluate representative by how even it samples in the known sample groups
                if yLabel == xLabel:
                    #If this is an evenly sampling method, check counts per class
                    if Constants.c_strEvenSelection in setYActual:
                        iMisclassificationCount = 0
                        for sClass in dictClasses:
                            iClassSelectionCount = len(dictClasses[sClass]&setXPredicted)
                            if iClassSelectionCount > iMaxSamplesPerCategory:
                                iMisclassificationCount = iClassSelectionCount - iMaxSamplesPerCategory
                        liConfusion.append(iCountPredictedSamples - iMisclassificationCount)
                    #Treat as a normal sample with a defined a prior classification
                    #If comparing the amount correct just look at actual and predicted classes
                    else:
                        liConfusion.append(len(setXPredicted & setYActual))

                else:
                    #If this is an evenly sampling method,
                    #Return 0, these do not have specific samples associated with them
                    if Constants.c_strEvenSelection in setYActual:
                        liConfusion.append(0)
                    #If the predicted label is the evenly sampled group
                    #Return if it samples more than it is allowed with the actual group
                    #This means there could be a column or row total of more than the total samples
                    #sampled if the actual groups overlapp.
                    elif Constants.c_strEvenSelection in setXPredictedPossible:
                        lsAssociatedClasses = dictOverlapMeasure.get(yLabel, [])
                        iMissclassifiedCount = 0
                        for strClass in lsAssociatedClasses:
                            lAreadyCountedMistakes = dictMisclassifiedClasses.get(xLabel, [])
                            if strClass not in lAreadyCountedMistakes:
                                iCommonWithClass = len(dictClasses[strClass]&setXPredicted)
                                if iCommonWithClass > iMaxSamplesPerCategory:
                                    iMissclassifiedCount += (iCommonWithClass-iMaxSamplesPerCategory)
                                    lAreadyCountedMistakes.append(strClass)
                                    dictMisclassifiedClasses[xLabel] = lAreadyCountedMistakes 
                        liConfusion.append(iMissclassifiedCount)

                    #Treat as a normal sample with a defined a prior classification
                    #If looking at missclassification remember that the classes here are overlapping.
                    #It is expected that A diversity sample is selected by representative methodology for example
                    #So make sure to first remove samples that may be in other classes but are ok to be selected
                    #Before determining error.
                    else:
                        setMisclassified = setXPredicted & (setYActual - setXPredictedPossible)
                        iMisclassified = len(setMisclassified)
                        if iMisclassified > 0:
                            lsAlreadyMistaken = dictMisclassfiedSamples.get(xLabel, [])
                            setNewMisclassified = setMisclassified - set(lsAlreadyMistaken)
                            iMisclassified = len(setNewMisclassified)
                            lsAlreadyMistaken.extend(list(setNewMisclassified))
                            dictMisclassfiedSamples[xLabel] = lsAlreadyMistaken

                        liConfusion.append(iMisclassified)

        #Cumulate results
        npUniversalConfusion = npUniversalConfusion + np.array(liConfusion)

    #Plot confusion matrix
    npUniversalConfusion = npUniversalConfusion / float(len(args.lsProjects))
    npUniversalConfusion.shape=iLabelLength,iLabelLength
    PlotMatrix.funcPlotMatrix(npMatrix=npUniversalConfusion, lsXLabels=lsXLabels, strOutputFigurePath=args.strOutputFigure, strXTitle="Predicted", strYTitle="Actual", fFlipYLabels=c_fFlipYLabels, fInvert=c_fInvert)

    logging.info("Stop MicropitaPaperConfusionMetaMatrix")

if __name__ == "__main__":
    _main( )
