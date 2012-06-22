#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Validates selection by generating PCoA plots of selection from one data set into another.
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
from Constants_Figures import Constants_Figures
from Diversity import Diversity
import itertools
import logging
import matplotlib.pyplot as plt
from MicroPITA import MicroPITA
import os
from PCoA import PCoA

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidatePCoA.py", 
    description = """Validates selection by generating PCoA plots of selection from one data set into another.""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strValidationIDNameArgument, dest="sValidationIDName", metavar= "SampleRowNameValidation", default=None, help= Constants_Arguments.c_strValidationIDNameHelp)
argp.add_argument(Constants_Arguments.c_strValidationLastMetadataNameArgument, dest="sValidationLastMetadataName", metavar= "FirstDataRowValidation", default=None, help= Constants_Arguments.c_strValidationLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strValidationIsSummedArgument, dest="fValidationIsSummed", action = "store", metavar= "flagIndicatingSummationForValidationFile", help= Constants_Arguments.c_strValidationIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strValidationIsNormalizedArgument, dest="fValidationIsNormalized", action = "store", metavar= "flagIndicatingNormalizationForValidationFile", 
                  help= Constants_Arguments.c_strValidationIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strPairingMetadataArgument, dest="sPairedMetadata", action = "store", metavar= "sMetadataUsedInPairing", help= Constants_Arguments.c_strPairingMetadataHelp)
argp.add_argument(Constants_Arguments.c_strMetricArgument, dest="sMetric", action = "store", metavar= "sMetric_For_Validation", help= Constants_Arguments.c_strMetricHelp)
argp.add_argument(Constants_Arguments.c_strSupervisedLabelArgument, dest="sStratifyMetadata", action = "store", metavar= "sMetadata_For_Stratification", help= Constants_Arguments.c_strSupervisedLabelHelp)

#Data file
argp.add_argument( "strValidationAbundanceFile", metavar = "Validation_Abundance_file", help = Constants_Arguments.c_strValidationAbundanceFileHelp)
argp.add_argument( "strSelectionAbundanceFile", metavar = "Selection_Abundance_file", help = Constants_Arguments.c_strSelectionAbundanceFileHelp)
argp.add_argument( "strSelectionFile", metavar = "Selection_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)

#Outputfile
argp.add_argument( "strOutFigure", metavar = "OutputFile", help = Constants_Arguments.c_genericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperValidatePCoA")
    logging.info("MicropitaPaperValidatePCoA. The following arguments were passed.")
    logging.info(str(args))

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fValidationIsSummed = (args.fValidationIsSummed.lower() == "true")
    fValidationIsNormalized = (args.fValidationIsNormalized.lower() == "true")

    #Number of dimensions to explore
    iDimensionCount = 4

    #Analysis object
    analysis = PCoA()

    #Standardize figures
    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape
    #Alpha
    dAlpha = Constants_Figures.c_dAlpha
    #Size
    iShapeSize = Constants_Figures.iMarkerSize

    #Read abundance file
    #Abundance table object to read in and manage data
    #Validation table
    abndValidationData = AbundanceTable.makeFromFile(strInputFile=args.strValidationAbundanceFile, fIsNormalized=fValidationIsNormalized,
                                            fIsSummed=fValidationIsSummed, sMetadataID=args.sValidationIDName, sLastMetadata=args.sValidationLastMetadataName)
    if not fValidationIsSummed:
        abndValidationData.funcSumClades()
    if not fValidationIsNormalized:
        abndValidationData.funcNormalize()

    #Selection table
    abndSelectionTable = AbundanceTable.makeFromFile(strInputFile=args.strSelectionAbundanceFile, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)
    if not fIsSummed:
        abndSelectionTable.funcSumClades()
    if not fIsNormalized:
        abndSelectionTable.funcNormalize()

    #Get sample names as a set
    setsSampleNames = set(abndValidationData.funcGetSampleNames())

    #Read in selection file
    dictAllSelectionStudies = MicroPITA.funcReadSelectionFileToDictionary(args.strSelectionFile)

    #Get a set of sample names from the validation set
    setValidationSampleIds = set(abndValidationData.funcGetSampleNames())

    #For each diversity methodology in the selection
    for sMethod in dictAllSelectionStudies.keys():
        if sMethod == args.sMetric:

            #Get the selection
            setsMetricSelection = set(dictAllSelectionStudies[sMethod])

            #Make sure the selection is in this data and there is other data to compare to
            if len(setsMetricSelection) < len(setsSampleNames):

                #Check to make sure the paired key is primary in both tables
                if (not abndSelectionTable.funcIsPrimaryIdMetadata(args.sPairedMetadata)) or (not abndValidationData.funcIsPrimaryIdMetadata(args.sPairedMetadata)):
                    logging.error("".join(["MicropitaPaperValidatePCoA:: tried to validate on a none unique key:",args.sPairedMetadata]))
                    return False

                #In the selection file go from the sampleID to the paired value
                lsPairedSelected = abndSelectionTable.funcTranslateIntoMetadata(lsValues=setsMetricSelection, sMetadataFrom=abndSelectionTable.funcGetIDMetadataName(),
                                                             sMetadataTo=args.sPairedMetadata, fFromPrimaryIds=True)
                if not lsPairedSelected:
                    logging.error("MicropitaPaperValidatePCoA:: Did not recieve lsPairedSelected.")
                    return False

                #In the validation file go from the paired value to the sampleID
                lsSelectedInValidation = abndValidationData.funcTranslateIntoMetadata(lsValues=lsPairedSelected, sMetadataFrom=args.sPairedMetadata,
                                                             sMetadataTo=abndValidationData.funcGetIDMetadataName(), fFromPrimaryIds=True)
                if not lsSelectedInValidation:
                    logging.error("MicropitaPaperValidatePCoA:: Did not recieve lsSelectedInValidation.")
                    return False

                if len(lsSelectedInValidation) == len(setsMetricSelection):

                    #Generate PCoA
                    #LoadData
                    analysis.loadData(xData=abndValidationData, fIsRawData=True)
                    #Make distance matrix
                    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS, iDims = iDimensionCount)

                    #Colors
                    acharColors = []
                    #Labels
                    acharSelection = []

                    #Color for selected samples in PCoA based on selection method
                    charSelectedColor = ""
                    objColors = Constants_Figures()
                    objColors.invertColors(fInvert=fInvert)
                    charSelectedColor = objColors.dictConvertMethodToHEXColor[sMethod]

                    #Generate colors
                    for strSample in abndValidationData.funcGetMetadata(abndValidationData.funcGetIDMetadataName()):
                        if(strSample in lsSelectedInValidation):
                            acharColors.append(charSelectedColor)
                            acharSelection.append(sMethod)
                        else:
                            acharColors.append(objColors.c_charNoSelect)
                            acharSelection.append(objColors.c_strPCOANotSelected)

                    #If given a valid label, stratify otherwise normal PCoA
                    for iDim1, iDim2 in itertools.combinations(range(1,iDimensionCount+1),2):
                        lsFilePieces = os.path.splitext(args.strOutFigure)
                        strOutputFigure = "".join([lsFilePieces[0],str(iDim1+1),"-",str(iDim2+1),lsFilePieces[1]])
                        if args.sStratifyMetadata:
                            lsStratify = abndValidationData.funcGetMetadata(args.sStratifyMetadata)
                            if lsStratify:
                                #Stratified PCoA
                                analysis.plotList(lsLabelList=lsStratify,
                                    strOutputFileName=strOutputFigure, iSize=iShapeSize, dAlpha=dAlpha,
                                    charForceColor=[acharColors,acharSelection], fInvert=fInvert, iDim1=iDim1, iDim2=iDim2)
                            else:
                                #Normal PCoA
                                analysis.plot(tempPlotName=strOutputFigure, tempColorGrouping=acharColors, tempShape=acharShape,
                                  tempLabels=acharSelection, tempShapeSize=iShapeSize, tempAlpha=dAlpha, tempLegendLocation="lower left", tempInvert = fInvert, iDim1=iDim1, iDim2=iDim2)
                        else:
                            #Normal PCoA
                            analysis.plot(tempPlotName=strOutputFigure, tempColorGrouping=acharColors, tempShape=acharShape,
                              tempLabels=acharSelection, tempShapeSize=iShapeSize, tempAlpha=dAlpha, tempLegendLocation="lower left", tempInvert = fInvert, iDim1=iDim1, iDim2=iDim2)

                else:
                    logging.error("MicropitaPaperValidatePCoA::Not all samples selected for diversity are in the given sample.")
            else:
                logging.error("MicropitaPaperValidatePCoA::Not all samples selected for diversity are in the given sample.")

    logging.info("Stop MicropitaPaperValidatePCoA")

if __name__ == "__main__":
    _main( )

