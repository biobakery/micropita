#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Validate Supervised selection in secondary dataset using
Distance from the other label group.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

import argparse
from BoxPlot import BoxPlot
from Constants import Constants
from Constants_Arguments import Constants_Arguments
from Diversity import Diversity
import logging
import math
from MicroPITA import MicroPITA
import numpy as np
from SVM import SVM
import os
from ValidationFigure import ValidationFigure


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidateHyperplaneByOrthogonalDistance.py", 
    description = """Generates boxplots showing the distribution of distances of selected samples from the average of the other label.""" )

c_PredictFile = "PredictFile"
c_sLabel = "Labels"

def funcGetAveragePopulation(abndTable, lfCompress):
  #Get the average populations
  lAverage = []

  for sFeature in abndTable.funcGetAbundanceCopy():
    sFeature = list(sFeature)[1:]
    sFeature=np.compress(lfCompress,sFeature,axis=0)
    #If there are no samples then return empty list.
    if len(sFeature) < 2:
      lAverage.append[0]
    #Return average of feature
    lAverage.append(sum(sFeature)/float(len(sFeature)))
  return lAverage

def funcGetDistanceFromAverage(abndTable,ldAverage,lsSamples,lfSelected,lfNotSelected):
  #Get the distance from label 1 of all samples in label0 splitting into selected and not selected lists
  ldSelectedDistances = []
  ldNotSelectedDistances = []

  for sSampleName in [sSample for iindex, sSample in enumerate(lsSamples) if lfSelected[iindex]]:
    #Get the sample measurements
    ldSelectedDistances.append(Diversity.funcGetBrayCurtisDissimilarity(np.array([abndTable.funcGetSample(sSampleName),ldAverage])))
  for sSampleName in [sSample for iindex, sSample in enumerate(lsSamples) if lfNotSelected[iindex]]:
    #Get the sample measurements
    ldNotSelectedDistances.append(Diversity.funcGetBrayCurtisDissimilarity(np.array([abndTable.funcGetSample(sSampleName),ldAverage])))
  return [ldSelectedDistances,ldNotSelectedDistances]

def funcPlot(abndRet,lsSelectedRet,sMetric,sOutputFigureName,objFigureControl,dictMisc):

  #Hold data for combined boxplots
  llBoxplotData = []
  lsBoxplotLabels = []

  #Run supervised micropita on the validation data set
  strFileNameBase = os.path.split(sOutputFigureName)[0]
  strInputPredictFile = "".join([strFileNameBase, "ValInputPrediction-select-", str(len(lsSelectedRet)), "-", str(dictMisc[c_sLabel]),".txt"])
  strPredictFile = "".join([strFileNameBase, "ValPrediction-select-", str(len(lsSelectedRet)), "-", str(dictMisc[c_sLabel]),".txt"])
  MicroPITA().funcRunMLPYSVM(abndAbundanceTable=abndRet,
                             sMetadataForLabel=dictMisc[c_sLabel],
                             strInputSVMFile=strInputPredictFile,
                             strPredictionFile=strPredictFile)

  #Get prediction file data
  strSVMSelection = ""
  with open(strPredictFile,'r') as fHndlInput:
    strSVMSelection = fHndlInput.read()

  #Get label and prediction data
  lsSVMFileLines = filter(None,strSVMSelection.split(Constants.ENDLINE))
  lsLabels = lsSVMFileLines[0].split(Constants.WHITE_SPACE)[1:]
  lsPredictions = [sLine.split(Constants.WHITE_SPACE)[0] for sLine in lsSVMFileLines[1:]]

  #Get label data
  if not len(lsLabels) == 2:
    logging.info("".join(["MicropitaPaperValidateHyperplaneByOrthogonalDistance. This only handles 2 label cases currently, was given a different about of labels (",str(len(lsLabels)),")."]))
    return

  #Create original labels
  liOriginalLabels = SVM.funcMakeLabels(abndRet.funcGetMetadata(dictMisc[c_sLabel]))

  #Indicate which labels were properly classified
  lfNotMisclassified = [lsPredictions[iindex] == ilabel for iindex, ilabel in enumerate(liOriginalLabels)]

  #Selected Validation samples (as list or strings of float)
  lsValidationSamples = abndRet.funcGetSampleNames()
  lfValidationSamples = [sSample in lsSelectedRet for sSample in lsValidationSamples]

  #Get the correctly classified sample populations by label
  #Label[0] Selected and not misclassified
  curLabel = lsLabels[0]
  lfLabel0Selected = [(lfNotMisclassified[iindex] and (liOriginalLabels[iindex] == curLabel) and lfValidationSamples[iindex]) for iindex, sSample in enumerate(lsValidationSamples)]
  #Label[0] Not Selected and not misclassified
  lfLabel0NotSelected = [(lfNotMisclassified[iindex] and (liOriginalLabels[iindex] == curLabel) and (not lfValidationSamples[iindex])) for iindex, sSample in enumerate(lsValidationSamples)]
  curLabel = lsLabels[1]
  #Label[1] Selected and not misclassified
  lfLabel1Selected = [(lfNotMisclassified[iindex] and (liOriginalLabels[iindex] == curLabel) and lfValidationSamples[iindex]) for iindex, sSample in enumerate(lsValidationSamples)]
  #Label[1] Not Selected and not misclassified
  lfLabel1NotSelected = [(lfNotMisclassified[iindex] and (liOriginalLabels[iindex] == curLabel) and (not lfValidationSamples[iindex])) for iindex, sSample in enumerate(lsValidationSamples)]

  #Get the average sample for label 0 and label 1 full populations
  lfFullPopulation0 = [lfLabel0Selected[iindex] or lfLabel0NotSelected[iindex] for iindex in xrange(len(lfLabel0Selected))]
  lfFullPopulation1 = [lfLabel1Selected[iindex] or lfLabel1NotSelected[iindex] for iindex in xrange(len(lfLabel1Selected))]

  #Get average populations
  lAverage0 = funcGetAveragePopulation(abndTable=abndRet, lfCompress=lfFullPopulation0)
  lAverage1 = funcGetAveragePopulation(abndTable=abndRet, lfCompress=lfFullPopulation1)

  ldSelectedDistances,ldNotSelectedDistances = funcGetDistanceFromAverage(abndTable = abndRet,ldAverage=lAverage1,
                                                    lsSamples=lsValidationSamples,lfSelected=lfLabel0Selected,
                                                    lfNotSelected=lfLabel0NotSelected)
  llBoxplotData.append(ldSelectedDistances)
  llBoxplotData.append(ldNotSelectedDistances)
  lsBoxplotLabels.append("".join(["Selected ",str(lsLabels[0]),"(",str(sum(lfLabel0Selected)),")"]))
  lsBoxplotLabels.append("".join(["Not Selected ",str(lsLabels[0]),"(",str(sum(lfLabel0NotSelected)),")"]))
  ldSelectedDistances,ldNotSelectedDistances = funcGetDistanceFromAverage(abndTable = abndRet,ldAverage=lAverage0,
                                                    lsSamples=lsValidationSamples,lfSelected=lfLabel1Selected,
                                                    lfNotSelected=lfLabel1NotSelected)
  llBoxplotData.append(ldSelectedDistances)
  llBoxplotData.append(ldNotSelectedDistances)
  lsBoxplotLabels.append("".join(["Selected ",str(lsLabels[1]),"(",str(sum(lfLabel0Selected)),")"]))
  lsBoxplotLabels.append("".join(["Not Selected ",str(lsLabels[1]),"(",str(sum(lfLabel0NotSelected)),")"]))

  #Start plot
  asFilePathPieces = os.path.splitext(sOutputFigureName)
  BoxPlot.funcPlot(ly=llBoxplotData, strOutputFigurePath="".join([asFilePathPieces[0],"-",sMetric,"-Boxplot-",asFilePathPieces[1]]),
                   lsLabels = lsBoxplotLabels, strTitle = "Distance from the other label ("+sMetric+").",
                   strXTitle= "Sample Population", strYTitle="Distance (Brays-Curtis)",
                   strColor = objFigureControl.dictConvertMethodToHEXColor[sMetric], fJitter=True, fInvert=objFigureControl.c_fInverted)

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
argp.add_argument(Constants_Arguments.c_strMetricArgument, dest="sMetric", action = "store", metavar= "sMetric_For_Validation", help= Constants_Arguments.c_strMetricHelp)
argp.add_argument(Constants_Arguments.c_strPairingMetadataArgument, dest="sPairedMetadata", action = "store", metavar= "sMetadataUsedInPairing", help= Constants_Arguments.c_strPairingMetadataHelp)
argp.add_argument(Constants_Arguments.c_strSupervisedLabelArgument, dest="sLabel", action = "store", metavar= "sMetadata_For_Label", help= Constants_Arguments.c_strSupervisedLabelHelp)
argp.add_argument(Constants_Arguments.c_strPredictFilePathArgument, dest = "sSVMPrediction", action = "store", default=None,
	help = Constants_Arguments.c_strPredictFilePathHelp)

#Data file
argp.add_argument( "strValidationAbundanceFile", metavar = "Validation_Abundance_file", help = Constants_Arguments.c_strValidationAbundanceFileHelp)
argp.add_argument( "strSelectionAbundanceFile", metavar = "Selection_Abundance_file", help = Constants_Arguments.c_strSelectionAbundanceFileHelp)
argp.add_argument( "strSelectionFile", metavar = "Selection_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)

#Outputfile
argp.add_argument( "strOutFigure", metavar = "BoxPlotOutputFile", help = Constants_Arguments.c_strGenericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Files
    strValidationSVMInputFile = os.path.splitext(args.strOutFigure)[0]+"-Val-Input-SVM.txt"
    strValidationSVMPredictFile = os.path.splitext(args.strOutFigure)[0]+"-Val-Predict-SVM.txt"

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperValidateHyperplaneByOrthogonalDistance")
    logging.info("MicropitaPaperValidateHyperplaneByOrthogonalDistance. The following arguments were passed.")
    logging.info(str(args))

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fValidationIsSummed = (args.fValidationIsSummed.lower() == "true")
    fValidationIsNormalized = (args.fValidationIsNormalized.lower() == "true")

    ValidationFigure().funcPlotFigure(funcPlot = funcPlot,
                         strValidationFile = args.strValidationAbundanceFile,
                         fNormalizeValidation = fValidationIsNormalized,
                         fSumValidation = fValidationIsSummed,
                         sSampleIDValidation = args.sValidationIDName,
                         sLastMetadataValidation = args.sValidationLastMetadataName,
                         strSelectionFile = args.strSelectionAbundanceFile,
                         fNormalizeSelection = fIsNormalized,
                         fSumSelection = fIsSummed,
                         sSampleIDSelection = args.sIDName,
                         sLastMetadatSelection = args.sLastMetadataName,
                         sPairingMetadata = args.sPairedMetadata,
                         sLabel = args.sLabel,
                         strMicropitaSelectionFile = args.strSelectionFile,
                         fInvert = fInvert,
                         lsMethods = [MicroPITA.c_strSVMClose, MicroPITA.c_strSVMFar],
                         sOutputFigureName = args.strOutFigure,
                         dictMisc = {c_PredictFile:args.sSVMPrediction,c_sLabel:args.sLabel})

    logging.info("Stop MicropitaPaperValidateHyperplaneByOrthogonalDistance")

if __name__ == "__main__":
    _main( )

