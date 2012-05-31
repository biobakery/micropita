#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Validate feature selection in a sanmple set
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
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidateFeatureAbundance.py", 
    description = """Generates boxplots showing the distribution of feature abundance in samples in different selected groupings.""" )

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
argp.add_argument(Constants_Arguments.c_strTargetedFeatureMethodArgument, dest="sFeatureSelection", metavar= "Feature Selection Method", default=Constants_Arguments.lsTargetedFeatureMethodValues[0], 
                  choices=Constants_Arguments.lsTargetedFeatureMethodValues, help= Constants_Arguments.c_strTargetedFeatureMethodHelp)

#Data file
argp.add_argument( "strValidationAbundanceFile", metavar = "Validation_Abundance_file", help = Constants_Arguments.c_strValidationAbundanceFileHelp)
argp.add_argument( "strSelectionAbundanceFile", metavar = "Selection_Abundance_file", help = Constants_Arguments.c_strSelectionAbundanceFileHelp)
argp.add_argument( "strSelectionFile", metavar = "Selection_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
argp.add_argument( "strValidateFeatureFile", metavar = "Feature_file", help = Constants_Arguments.c_strTargetedSelectionFileHelp)

#Outputfile
argp.add_argument( "strOutFigure", metavar = "BoxPlotOutputFile", help = Constants_Arguments.c_genericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperValidateFeatureAbundance")
    logging.info("MicropitaPaperValidateFeatureAbundance. The following arguments were passed.")
    logging.info(str(args))


    if not args.sFeatureSelection:
        logging.error("MicropitaPaperValidateFeatureAbundance::Did not receive a method to measure the feature so the plot was not generated.")
        return False

    c_PlotAbundance = None
    if args.sFeatureSelection.lower() == Constants_Arguments.c_TARGETED_METHOD_RANKED.lower():
        c_PlotAbundance = False
    elif args.sFeatureSelection.lower() == Constants_Arguments.c_TARGETED_METHOD_ABUNDANCE.lower():
        c_PlotAbundance = True
    else:
        logging.error("MicropitaPaperValidateFeatureAbundance::Did not receive a valid value for method to measure the feature so the plot was not generated. Received:"+str( args.sFeatureSelection))
        return False

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

    #Get features from file
    sContents = ""
    with open(args.strValidateFeatureFile, 'r') as f:
        sContents = f.read()
    f.close()
    lsFeatures =  filter(None, sContents.split(Constants.ENDLINE))
    if not lsFeatures:
        logging.error("".join(["MicropitaPaperValidateFeatureAbundance:: Could not read features form the file:",args.strValidateFeatureFile]))
        return False

    #For each diversity methodology in the selection
    for sMethod in dictAllSelectionStudies.keys():
        if sMethod == MicroPITA.c_strTaxa:

            #Get the selection
            setsSelection = set(dictAllSelectionStudies[sMethod])

            #Make sure the selection is in this data and there is other data to compare to
            if len(setsSelection) < len(setsSampleNames):

                #Check to make sure the paired key is primary in both tables
                if (not abndSelectionTable.funcIsPrimaryIdMetadata(args.sPairedMetadata)) or (not abndValidationData.funcIsPrimaryIdMetadata(args.sPairedMetadata)):
                    logging.error("".join(["MicropitaPaperValidateFeatureAbundance:: tried to validate on a none unique key:",args.sPairedMetadata]))
                    return False

                #In the selection file go from the sampleID to the paired value
                lsPairedSelected = abndSelectionTable.funcTranslateIntoMetadata(lsValues=setsSelection, sMetadataFrom=abndSelectionTable.funcGetIDMetadataName(),
                                                             sMetadataTo=args.sPairedMetadata, fFromPrimaryIds=True)
                if not lsPairedSelected:
                    logging.error("MicropitaPaperValidateFeatureAbundance:: Did not recieve lsPairedSelected.")
                    return False

                #In the validation file go from the paired value to the sampleID
                lsSelectedInValidation = abndValidationData.funcTranslateIntoMetadata(lsValues=lsPairedSelected, sMetadataFrom=args.sPairedMetadata,
                                                             sMetadataTo=abndValidationData.funcGetIDMetadataName(), fFromPrimaryIds=True)
                if not lsSelectedInValidation:
                    logging.error("MicropitaPaperValidateFeatureAbundance:: Did not recieve lsSelectedInValidation.")
                    return False

                if len(lsSelectedInValidation) == len(setsSelection):

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
                    if not c_PlotAbundance:
                        imgSubplot.set_ylabel("Ranked Abundance")
                    else:
                        imgSubplot.set_ylabel("Relative Abundance")
                    imgSubplot.spines['top'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['bottom'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['left'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.spines['right'].set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.xaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
                    #Adds light grid for numbers and puts them in the background
                    imgSubplot.yaxis.grid(True, linestyle='-', which='major', color=objFigureControl.c_strGridLineColor, alpha=objFigureControl.c_dAlpha)
                    imgSubplot.set_axisbelow(True)
                    imgSubplot.yaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.tick_params(axis='x', colors=objFigureControl.c_strDetailsColorLetter)
                    imgSubplot.tick_params(axis='y', colors=objFigureControl.c_strDetailsColorLetter)
                    charMarkerEdgeColor = objFigureControl.c_strDetailsColorLetter

                    #Create selected and not selected groupings
                    #Measure feature abundance of selected and not selected populations
                    #Reduce the table to just the features
                    #Rank the data if need be
                    if not c_PlotAbundance:
                        abndValidationData = abndValidationData.funcRankAbundance()
                    abndFeatureTable = abndValidationData.funcGetFeatureAbundanceTable(lsFeatures)
                    if not abndFeatureTable:
                        logging.error("MicropitaPaperValidateFeatureAbundance:: did not receive a reduced feature abundance table.")
                        return False

                    #Get a copy of the abundance table with everything but the features of interest.
                    abndOppositeFeatureTable = abndValidationData.funcGetFeatureAbundanceTable(set(abndValidationData.funcGetFeatureNames())-set(lsFeatures))
                    if not abndOppositeFeatureTable:
                        logging.error("MicropitaPaperValidateFeatureAbundance:: did not receive a abundance table with all but the feature.")
                        return False

                    #Get average abundance per sample selection group and plot
                    setsNotSelected = setsSampleNames-set(lsSelectedInValidation)
                    dFeatureSampleLength = float(abndFeatureTable.funcGetFeatureCount())
                    ldAverageSelectedAbundance = [sum(abndFeatureTable.funcGetSample(sSample))/dFeatureSampleLength for sSample in lsSelectedInValidation]
                    ldAverageNotSelectedAbundance = [sum(abndFeatureTable.funcGetSample(sSample))/dFeatureSampleLength for sSample in setsNotSelected]

                    #Make box plot
                    print "ValidateFeature::","Selected",ldAverageSelectedAbundance,"\nNot Selected",ldAverageNotSelectedAbundance,"\nSamples",setsSelection
                    bp = plt.boxplot(x=[ldAverageSelectedAbundance,ldAverageNotSelectedAbundance], notch=1, patch_artist=True)                       

                    #Color boxes
                    plt.setp(bp['boxes'], color=objFigureControl.c_strDetailsColorLetter, facecolor=objFigureControl.dictConvertMethodToHEXColor[sMethod], alpha=objFigureControl.c_dAlpha)
                    plt.setp(bp['whiskers'], color=objFigureControl.c_strDetailsColorLetter)

                    #Set ticks and title
                    #xtickNames = plt.setp(imgSubplot, xticklabels=["Selected", "Not Selected", "Diff Selected", "Diff Not Selected","x1","x2","x3","x4"])
                    xtickNames = plt.setp(imgSubplot, xticklabels=["".join(["Selected (",str(len(ldAverageSelectedAbundance)),")"]),
                                                                   "".join(["Not Selected (",str(len(ldAverageNotSelectedAbundance)),")"])])
                    plt.scatter(x=[1]*len(ldAverageSelectedAbundance),y=ldAverageSelectedAbundance,c=objFigureControl.dictConvertMethodToHEXColor[sMethod],marker="o",alpha=objFigureControl.c_dAlpha)
                    plt.scatter(x=[2]*len(ldAverageNotSelectedAbundance),y=ldAverageNotSelectedAbundance,c=objFigureControl.dictConvertMethodToHEXColor[sMethod],marker="o",alpha=objFigureControl.c_dAlpha)
                    imgSubplot.set_title("Targeted feature selection shown in validation data.")

                    if not c_PlotAbundance:
                        ax = plt.gca()
                        ax.set_ylim(ax.get_ylim()[::-1])

                    #End plot
                    #Save to a file
                    imgFigure.savefig(args.strOutFigure, facecolor=imgFigure.get_facecolor())

                else:
                    logging.error("MicropitaPaperValidateFeatureAbundance::Not all samples selected for diversity are in the given sample.")
            else:
                logging.error("MicropitaPaperValidateFeatureAbundance::Not all samples selected for diversity are in the given sample.")

    logging.info("Stop MicropitaPaperValidateFeatureAbundance")

if __name__ == "__main__":
    _main( )

