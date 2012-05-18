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
import os

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidateDiversity.py", 
    description = """Generates boxplots showing the distribution of diversity of samples in different selected groupings.""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDName, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDName)
argp.add_argument(Constants_Arguments.c_strLastMetadataName, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strValidationIDName, dest="sValidationIDName", metavar= "SampleRowNameValidation", default=None, help= Constants_Arguments.c_strValidationIDName)
argp.add_argument(Constants_Arguments.c_strValidationLastMetadataName, dest="sValidationLastMetadataName", metavar= "FirstDataRowValidation", default=None, help= Constants_Arguments.c_strValidationLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strValidationIsSummedArgument, dest="fValidationIsSummed", action = "store", metavar= "flagIndicatingSummationForValidationFile", help= Constants_Arguments.c_strValidationIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strValidationIsNormalizedArgument, dest="fValidationIsNormalized", action = "store", metavar= "flagIndicatingNormalizationForValidationFile", 
                  help= Constants_Arguments.c_strValidationIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strPairingMetadata, dest="sPairedMetadata", action = "store", metavar= "sMetadataUsedInPairing", help= Constants_Arguments.c_strPairingMetadataHelp)

#Data file
argp.add_argument( "strValidationAbundanceFile", metavar = "Validation_Abundance_file", help = Constants_Arguments.c_strValidationAbundanceFileHelp)
argp.add_argument( "strSelectionAbundanceFile", metavar = "Selection_Abundance_file", help = Constants_Arguments.c_strSelectionAbundanceFileHelp)
argp.add_argument( "strSelectionFile", metavar = "Selection_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)

#Outputfile
argp.add_argument( "strOutFigure", metavar = "BoxPlotOutputFile", help = Constants_Arguments.c_genericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    print("args")
    print(args)

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
    totalData = AbundanceTable.makeFromFile(strInputFile=args.strValidationAbundanceFile, fIsNormalized=fValidationIsNormalized,
                                            fIsSummed=fValidationIsSummed, sMetadataID=args.sValidationIDName, sLastMetadata=args.sValidationLastMetadataName)
    totalData.funcNormalize()

    #Selection table
    abndSelectionTable = AbundanceTable.makeFromFile(strInputFile=args.strSelectionAbundanceFile, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)
    abndSelectionTable.funcNormalize()

    #Get sample names as a set
    setsSampleNames = set(totalData.funcGetSampleNames())

    #Read in selection file
    dictAllSelectionStudies = MicroPITA.funcReadSelectionFileToDictionary(args.strSelectionFile)

    #For each diversity methodology in the selection
    for sMethod in dictAllSelectionStudies.keys():
        if sMethod == MicroPITA.c_DIVERSITY_1:

            #Get the alpha metric being used and normalize
            sMetric = MicroPITA.convertMicroPITAToAMetric[sMethod]
            
            #Get the selection
            setsDiversitySelection = set(dictAllSelectionStudies[sMethod])

            #Make sure the selection is in this data and there is other data to compare to
            if len(setsDiversitySelection) < len(setsSampleNames):

                #Check to make sure the paired key is primary in both tables
                if (not abndSelectionTable.funcIsPrimaryIdMetadata(args.sPairedMetadata)) or (not totalData.funcIsPrimaryIdMetadata(args.sPairedMetadata)):
                    logging.error("".join(["MicropitaPaperValidateDiversity:: tried to validate on a none unique key:",args.sPairedMetadata]))
                    return False

                #In the selection file go from the sampleID to the paired value
                lsPairedSelected = abndSelectionTable.funcTranslateIntoMetadata(lsValues=setsDiversitySelection, sMetadataFrom=abndSelectionTable.funcGetIDMetadataName(),
                                                             sMetadataTo=args.sPairedMetadata, fFromPrimaryIds=True)
                if not lsPairedSelected:
                    logging.error("MicropitaPaperValidateDiversity:: Did not recieve lsPairedSelected.")
                    return False

                #In the validation file go from the paired value to the sampleID
                lsSelectedInValidation = totalData.funcTranslateIntoMetadata(lsValues=lsPairedSelected, sMetadataFrom=args.sPairedMetadata,
                                                             sMetadataTo=totalData.funcGetIDMetadataName(), fFromPrimaryIds=True)

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
                    imgSubplot.set_ylabel("Relative Abundance")
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
                    ldSelectedDiversity = Diversity.buildAlphaMetricsMatrix(tempSampleAbundance = totalData.funcGetAbundanceCopy(), tempSampleNames = lsSelectedInValidation, tempDiversityMetricAlpha = sMetric)
                    ldNotSelectedDiversity = Diversity.buildAlphaMetricsMatrix(tempSampleAbundance = totalData.funcGetAbundanceCopy(), tempSampleNames = setsSampleNames-set(lsSelectedInValidation), tempDiversityMetricAlpha = sMetric)

                    #Make box plot
                    bp = plt.boxplot(x=[ldSelectedDiversity[0],ldNotSelectedDiversity[0]], notch=1, patch_artist=True)

                    #Color boxes
                    plt.setp(bp['boxes'], color=objFigureControl.c_strDetailsColorLetter, facecolor=objFigureControl.dictConvertMethodToHEXColor[sMethod], alpha=objFigureControl.c_dAlpha)
                    plt.setp(bp['whiskers'], color=objFigureControl.c_strDetailsColorLetter)

                    #Set ticks and title
                    xtickNames = plt.setp(imgSubplot, xticklabels=["Selected", "Not Selected"])
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

