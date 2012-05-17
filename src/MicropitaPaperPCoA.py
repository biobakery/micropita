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

from AbundanceTable import AbundanceTable
import argparse
from Constants import Constants
from Constants_Arguments import Constants_Arguments
from Constants_Figures import Constants_Figures
import logging
import matplotlib.cm as cm
from MicroPITA import MicroPITA
import os
from PCoA import PCoA
from ValidateData import ValidateData

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperPCoA.py", description = """Creates PCoA plots for MicroPITA results.""" )
#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, 
                  help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDName, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDName)
argp.add_argument(Constants_Arguments.c_strLastMetadataName, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False",
	help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False",
	help = Constants_Arguments.c_strInvertHelp)
argp.add_argument(Constants_Arguments.c_strPredictFilePath, dest = "sSVMPrediction", action = "store", default=None,
	help = Constants_Arguments.c_strPredictFilePathHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strSumDataArgument, dest="fSumData", action = "store", metavar= "WouldlikeDataSummed", help= Constants_Arguments.c_strSumDataHelp)

#Abundance file
argp.add_argument("strFileAbund", action="store", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Select file
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = Constants_Arguments.c_strOptionalOutputDataFile)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperPCoA")
    print("Start MicropitaPaperPCoA")
    print(args)

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    c_fCheckFile = False
    c_NotSelected = Constants_Figures.c_strPCOANotSelected

    c_fInvert = (args.fInvert == "True")
    c_Normalize = (args.fNormalize == "True")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fSumData = (args.fSumData.lower() == "true")

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable.makeFromFile(strInputFile=args.strFileAbund, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)

    #Normalize if needed and sum clades
    if fSumData:
        rawData.funcSumClades()
    if c_Normalize:
        rawData.funcNormalize()

    metadata = rawData.funcGetMetadataCopy()
    sampleNames = rawData.funcGetSampleNames()

    #Standardize figures
    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape
    #Alpha
    dAlpha = Constants_Figures.c_dAlpha
    #Size
    c_shapeSize = Constants_Figures.iMarkerSize
    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(xData=rawData, fIsRawData=True)
    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Draw known truths
    #Draw labeling from metadata
    for asMetadata in metadata:
      analysis.plotList(lsLabelList=metadata[asMetadata],strOutputFileName="".join([asFilePathPieces[0],"-Truth-",str(asMetadata),"-",asFilePathPieces[1]]),iSize=c_shapeSize, dAlpha =dAlpha, fInvert=c_fInvert)

    #Read in prediction file is supplied
    lsPredictions = list()
    if(not args.sSVMPrediction in [None, "None", "none"]):
        fHndlInput = open(args.sSVMPrediction,'r')
        strSVMSelection = fHndlInput.read()
        fHndlInput.close()

        #Make a list of predictions
        for strSVMSelectionLine in filter(None,strSVMSelection.split(Constants.ENDLINE)):
            lsPredictElements = strSVMSelectionLine.split(Constants.WHITE_SPACE)
            lsPredictions.append(lsPredictElements[0])
        analysis.plotList(lsLabelList=lsPredictions[1:],strOutputFileName="".join([asFilePathPieces[0],"-SVMPredictions",asFilePathPieces[1]]),iSize=c_shapeSize,dAlpha=dAlpha,fInvert=c_fInvert)

    #Draw selections
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))
    for strSelectionMethod in lstrSelection:
        #Colors
        acharColors = []
        #Labels
        acharSelection = []

        #Color for selected samples in PCoA based on selection method
        charSelectedColor = ""
        objColors = Constants_Figures()
        objColors.invertColors(fInvert=c_fInvert)
        astrSelectionMethod = strSelectionMethod.split(Constants.COLON)
        charSelectedColor = objColors.dictConvertMethodToHEXColor[astrSelectionMethod[0]]

        #Parse samples
        astrSelectedSamples = astrSelectionMethod[1].split(Constants.COMMA)
        astrSelectedSamples = [strSelectedSample.strip() for strSelectedSample in astrSelectedSamples]

        #Generate colors
        for strSample in sampleNames:
            if(strSample in astrSelectedSamples):
                acharColors.append(charSelectedColor)
                acharSelection.append(astrSelectionMethod[0])
            else:
                acharColors.append(objColors.c_charNoSelect)
                acharSelection.append(c_NotSelected)

        #Draw PCoA
        if astrSelectionMethod[0] in [MicroPITA.c_SVM_CLOSE, MicroPITA.c_SVM_FAR]:
          analysis.plotList(lsLabelList=lsPredictions[1:],strOutputFileName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]),
              iSize=c_shapeSize, dAlpha=dAlpha, charForceColor=[acharColors,acharSelection], fInvert=c_fInvert)
        else:
          analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]), tempColorGrouping=acharColors,
              tempShape=acharShape, tempLabels=acharSelection, tempShapeSize=c_shapeSize, tempAlpha=dAlpha, tempLegendLocation="lower left", tempInvert = c_fInvert)

    logging.info("Stop MicropitaPaperPCoA")

if __name__ == "__main__":
    _main( )

