#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Validate diversity selection in a sanmple set
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
import logging
import matplotlib.pyplot as plt
from MicroPITA import MicroPITA
import random
import os

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidateDiversity.py", 
    description = """Generates boxplots showing the distribution of diversity of samples in different selected groupings.""" )

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

#Data file
argp.add_argument( "strValidationAbundanceFile", metavar = "Validation_Abundance_file", help = Constants_Arguments.c_strValidationAbundanceFileHelp)
argp.add_argument( "strSelectionAbundanceFile", metavar = "Selection_Abundance_file", help = Constants_Arguments.c_strSelectionAbundanceFileHelp)
argp.add_argument( "strSelectionFile", metavar = "Selection_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)

#Outputfile
argp.add_argument( "strOutFigure", metavar = "BoxPlotOutputFile", help = Constants_Arguments.c_strGenericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperValidateDiversity")
    logging.info("MicropitaPaperValidateDiversity. The following arguments were passed.")
    logging.info(str(args))

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fValidationIsSummed = (args.fValidationIsSummed.lower() == "true")
    fValidationIsNormalized = (args.fValidationIsNormalized.lower() == "true")

    #Read abundance file
    #Abundance table object to read in and manage data
    #Validation table
    abndValidationData = AbundanceTable.funcMakeFromFile(strInputFile=args.strValidationAbundanceFile, fIsNormalized=fValidationIsNormalized,
                                            fIsSummed=fValidationIsSummed, sMetadataID=args.sValidationIDName, sLastMetadata=args.sValidationLastMetadataName)
    if not fValidationIsSummed:
        abndValidationData.funcSumClades()
    if not fValidationIsNormalized:
        abndValidationData.funcNormalize()

    #Selection table
    abndSelectionTable = AbundanceTable.funcMakeFromFile(strInputFile=args.strSelectionAbundanceFile, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)
    if not fIsSummed:
        abndSelectionTable.funcSumClades()
    if not fIsNormalized:
        abndSelectionTable.funcNormalize()

    #Get sample names as a set
    setsValidationMetadata = set(abndValidationData.funcGetMetadata(abndValidationData.funcGetIDMetadataName()))

    #Read in selection file
    dictAllSelectionStudies = MicroPITA.funcReadSelectionFileToDictionary(args.strSelectionFile)

    #For each diversity methodology in the selection
    for sMethod in dictAllSelectionStudies.keys():
        if sMethod == MicroPITA.c_strDiversity1:

            #Get the alpha metric being used and normalize
            sMetric = MicroPITA.dictConvertMicroPITAToAMetric[sMethod]
            
            #Get the selection
            setsDiversitySelection = set(dictAllSelectionStudies[sMethod])

            #Make sure the selection is in this data and there is other data to compare to
            if len(setsDiversitySelection) < len(setsValidationMetadata):

                #Check to make sure the paired key is primary in both tables
                if (not abndSelectionTable.funcIsPrimaryIdMetadata(args.sPairedMetadata)) or (not abndValidationData.funcIsPrimaryIdMetadata(args.sPairedMetadata)):
                    logging.error("".join(["MicropitaPaperValidateDiversity:: tried to validate on a none unique key:",args.sPairedMetadata]))
                    return False

                #In the selection file go from the sampleID to the paired value
                lsPairedSelected = abndSelectionTable.funcTranslateIntoMetadata(lsValues=setsDiversitySelection, sMetadataFrom=abndSelectionTable.funcGetIDMetadataName(),
                                                             sMetadataTo=args.sPairedMetadata, fFromPrimaryIds=True)
                if not lsPairedSelected:
                    logging.error("MicropitaPaperValidateDiversity:: Did not recieve lsPairedSelected.")
                    return False

                #In the validation file go from the paired value to the sampleID
                lsSelectedInValidation = abndValidationData.funcTranslateIntoMetadata(lsValues=lsPairedSelected, sMetadataFrom=args.sPairedMetadata,
                                                             sMetadataTo=abndValidationData.funcGetIDMetadataName(), fFromPrimaryIds=True)

                if not lsSelectedInValidation:
                    logging.error("MicropitaPaperValidateDiversity:: Did not recieve lsSelectedInValidation.")
                    return False

                if len(lsSelectedInValidation) == len(setsDiversitySelection):

                    #Start plot
                    #Get plot object
                    imgFigure = plt.figure()

                    #Get plot colors
                    objFigureControl = Constants_Figures()
                    objFigureControl.invertColors(fInvert=fInvert)

                    #Color/Invert figure
                    imgFigure.set_facecolor(objFigureControl.c_strBackgroundColorWord)
                    imgSubplot = imgFigure.add_subplot(111,axisbg=objFigureControl.c_strBackgroundColorLetter)
                    imgSubplot.set_xlabel("Sample Population")
                    imgSubplot.set_ylabel("Diversity (Inverse Simpson)")
                    imgSubplot.spines['top'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['bottom'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['left'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['right'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.xaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.yaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
                    #Adds light grid for numbers and puts them in the background
                    imgSubplot.yaxis.grid(True, linestyle='-', which='major', color=objFigureControl.c_strGridLineColor, alpha=objFigureControl.c_dAlpha)
                    imgSubplot.set_axisbelow(True)
                    imgSubplot.tick_params(axis='x', colors=objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.tick_params(axis='y', colors=objFigureControl.c_strDetailsColorLetter)
                    charMarkerEdgeColor = objFigureControl.c_strDetailsColorLetter

                    #Create selected and not selected groupings
                    #Measure diversity of selected and not selected groupings
                    ldSelectedDiversity = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abndValidationData.funcGetAbundanceCopy(), lsSampleNames = lsSelectedInValidation, lsDiversityMetricAlpha = sMetric)
                    ldNotSelectedDiversity = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abndValidationData.funcGetAbundanceCopy(), lsSampleNames = (setsValidationMetadata-set(lsSelectedInValidation)), lsDiversityMetricAlpha = sMetric)

                    #Make box plot
                    bp = plt.boxplot(x=[ldSelectedDiversity[0],ldNotSelectedDiversity[0]], notch=1, patch_artist=True)
                    ldJitteredX1 = [1.0+random.uniform(-.05,.05) for x in xrange(len(ldSelectedDiversity[0]))]
                    ldJitteredX2 = [2.0+random.uniform(-.05,.05) for x in xrange(len(ldNotSelectedDiversity[0]))]
                    plt.scatter(x=ldJitteredX1,y=ldSelectedDiversity[0],c=objFigureControl.dictConvertMethodToHEXColor[sMethod],marker="o",alpha=objFigureControl.c_dAlpha)
                    plt.scatter(x=ldJitteredX2,y=ldNotSelectedDiversity[0],c=objFigureControl.dictConvertMethodToHEXColor[sMethod],marker="o",alpha=objFigureControl.c_dAlpha)

                    #Color boxes
                    plt.setp(bp['boxes'], color=objFigureControl.c_strDetailsColorLetter, facecolor=objFigureControl.dictConvertMethodToHEXColor[sMethod], alpha=objFigureControl.c_dAlpha)
                    plt.setp(bp['whiskers'], color=objFigureControl.c_strDetailsColorLetter)

                    #Set ticks and title
                    xtickNames = plt.setp(imgSubplot, xticklabels=["".join(["Selected (",str(len(ldSelectedDiversity[0])),")"]),
                                                                   "".join(["Not Selected (",str(len(ldNotSelectedDiversity[0])),")"])])
                    imgSubplot.set_title("Maximum diversity shown in validation data set.")

                    #End plot
                    #Save to a file
                    imgFigure.savefig(args.strOutFigure, facecolor=imgFigure.get_facecolor())

                else:
                    logging.error("MicropitaPaperValidateDiversity::Not all samples selected for diversity are in the given sample.")
            else:
                logging.error("MicropitaPaperValidateDiversity::Not all samples selected for diversity are in the given sample.")

    logging.info("Stop MicropitaPaperValidateDiversity")

if __name__ == "__main__":
    _main( )

