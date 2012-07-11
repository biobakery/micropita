#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Validate Supervised selection in secondary dataset
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
import logging
import math
from MicroPITA import MicroPITA
from SVM import SVM
import os
from ValidationFigure import ValidationFigure


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperValidateHyperplaneByProbability.py", 
    description = """Generates boxplots showing the distribution of probabilities of selected and not selected per class.""" )

c_PredictFile = "PredictFile"
c_sLabel = "Labels"

def funcPlot(abndRet,lsSelectedRet,sMetric,sOutputFigureName,objFigureControl,dictMisc):

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

  #Get label data
  iLabelLength = len(lsLabels)
  dCentralProbability = 1.0/iLabelLength
  if not (iLabelLength == 2):
    logging.error("".join(["MicropitaPaperValidationPCoA::Expected only 2 labels. Received ",str(len(lsLabels))," labels. Figure not generated."]))
    return False

  #Get distance from hyperplane, predictions, and selection
  lsPredictions = [sPredLine.split(Constants.WHITE_SPACE)[0] for sPredLine in lsSVMFileLines[1:]]

  lsDistances = []
  for sPredLine in lsSVMFileLines[1:]:
    lsProbabilities = sPredLine.split(Constants.WHITE_SPACE)[1:]
    lsDistances.append(sum([(dCentralProbability-float(sProbability))*(dCentralProbability-float(sProbability)) for sProbability in lsProbabilities]))
  lfSelection = [sSample in lsSelectedRet for sSample in abndRet.funcGetSampleNames()]

  #Create original labels
  liOriginalLabels = SVM.funcMakeLabels(abndRet.funcGetMetadata(dictMisc[c_sLabel]))

  #Indicate which labels were properly classified
  lfNotMisclassified = [lsPredictions[iindex] == ilabel for iindex, ilabel in enumerate(liOriginalLabels)]

  #Hold data for combined boxplots
  llBoxplotData = []
  lsBoxplotLabels = []

  #Stratify by label
  for strLabel in lsLabels:
    lsCurProbabilities = [float(lsDistances[iindex]) for iindex in xrange(len(lsDistances)) if ((lsPredictions[iindex] == strLabel) and (lfNotMisclassified[iindex]))]
    lsCurPredictions = [int(lsPredictions[iindex]) for iindex in xrange(len(lsDistances)) if ((lsPredictions[iindex] == strLabel) and (lfNotMisclassified[iindex]))]
    lsCurSelection = [lfSelection[iindex] for iindex in xrange(len(lsDistances)) if ((lsPredictions[iindex] == strLabel) and (lfNotMisclassified[iindex]))]

    #Create selected and not selected groupings
    #Stratify to selected and not selected groupings
    ldSelected = []
    ldNotSelected = []
    [ldNotSelected.append(lsCurProbabilities[iindex]) if lsCurSelection[iindex] == False else ldSelected.append(lsCurProbabilities[iindex])
       for iindex in xrange(len(lsCurProbabilities))]

    if len(ldSelected) or len(ldNotSelected):
      llBoxplotData.extend([ldSelected,ldNotSelected])
      lsBoxplotLabels.extend(["".join(["Selected Label "+str(strLabel)+" (",str(len(ldSelected)),")"]),
                                    "".join(["Not Selected Label "+str(strLabel)+" (",str(len(ldNotSelected)),")"])])
  #Start plot
  asFilePathPieces = os.path.splitext(sOutputFigureName)
  BoxPlot.funcPlot(ly=llBoxplotData, strOutputFigurePath="".join([asFilePathPieces[0],"-",sMetric,"-Boxplot-",asFilePathPieces[1]]),
                   lsLabels = lsBoxplotLabels, strTitle = "Probabilistic Distance from the Hyperplane ("+sMetric+").",
                   strXTitle= "Sample Population", strYTitle="Variance from the hyperplane.",
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

    logging.info("Start MicropitaPaperValidateHyperplaneByProbability")
    logging.info("MicropitaPaperValidateHyperplaneByProbability. The following arguments were passed.")
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

    logging.info("Stop MicropitaPaperValidateHyperplaneByProbability")

if __name__ == "__main__":
    _main( )

