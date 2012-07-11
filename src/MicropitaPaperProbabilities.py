#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
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
import os
from SVM import SVM
from ValidationFigure import ValidationFigure


c_PredictFile = "PredictFile"
c_sLabel = "Labels"

def funcPlot(abndRet,lsSelectedRet,sMetric,sOutputFigureName,objFigureControl,dictMisc):

  #Get prediction file data
  strSVMSelection = ""
  with open(dictMisc[c_PredictFile],'r') as fHndlInput:
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

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperSupervisedPCoA.py", description = """Creates PCoA plots for MicroPITA results.""" )
#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, 
                  help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False",
	help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False",
	help = Constants_Arguments.c_strInvertHelp)
argp.add_argument(Constants_Arguments.c_strPredictFilePathArgument, dest = "sSVMPrediction", action = "store", default=None,
	help = Constants_Arguments.c_strPredictFilePathHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strSumDataArgument, dest="fSumData", action = "store", metavar= "WouldlikeDataSummed", help= Constants_Arguments.c_strSumDataHelp)
argp.add_argument(Constants_Arguments.c_strSupervisedLabelArgument, dest="sLabel", action = "store", metavar= "sMetadata_For_Label", help= Constants_Arguments.c_strSupervisedLabelHelp)

#Abundance file
argp.add_argument("strFileAbund", action="store", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Select file
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.pdf", nargs = "?", help = Constants_Arguments.c_strOptionalOutputDataFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)
    logging.info("Start MicropitaPaperSupervisedPCoA")

    c_fInvert = (args.fInvert.lower() == "true")
    c_Normalize = (args.fNormalize.lower() == "true")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fSumData = (args.fSumData.lower() == "true")

    ValidationFigure().funcPlotFigure(funcPlot = funcPlot,
                         strValidationFile = None,
                         fNormalizeValidation = None,
                         fSumValidation = None,
                         sSampleIDValidation = None,
                         sLastMetadataValidation = None,
                         strSelectionFile = args.strFileAbund,
                         fNormalizeSelection = fIsNormalized,
                         fSumSelection = fIsSummed,
                         sSampleIDSelection = args.sIDName,
                         sLastMetadatSelection = args.sLastMetadataName,
                         sPairingMetadata = None,
                         sLabel = args.sLabel,
                         strMicropitaSelectionFile = args.strSelectionFile,
                         fInvert = c_fInvert,
                         lsMethods = [MicroPITA.c_strSVMClose, MicroPITA.c_strSVMFar],
                         sOutputFigureName = args.strOutFile,
                         dictMisc = {c_PredictFile:args.sSVMPrediction,c_sLabel:args.sLabel})

    logging.info("Stop MicropitaPaperSupervisedPCoA")

if __name__ == "__main__":
    _main( )

