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
argp.add_argument(Constants_Arguments.c_strSampleNameRowArgument, dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= Constants_Arguments.c_strSampleNameRowHelp)
argp.add_argument(Constants_Arguments.c_strFirstDataRow, dest="iFirstDataRow", metavar= "FirstDataRow", default=1, 
                  help= Constants_Arguments.c_strFirstDataRowHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False",
	help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False",
	help = Constants_Arguments.c_strInvertHelp)
argp.add_argument(Constants_Arguments.c_strPredictFilePath, dest = "sSVMPrediction", action = "store", default=None,
	help = Constants_Arguments.c_strPredictFilePathHelp)
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
    print(args)

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    c_fCheckFile = False
    c_NotSelected = "Not_Selected"
    c_shapeSize = 40

    c_fInvert = (args.fInvert == "True")
    c_Normalize = (args.fNormalize == "True")

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable()
    abundance,metadata = rawData.textToStructuredArray(tempInputFile=args.strFileAbund, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=c_Normalize)
    sampleNames = abundance.dtype.names[1:]

    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(tempReadData=args.strFileAbund, tempIsRawData=True, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=c_Normalize, tempCheckFile=c_fCheckFile)
    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Draw known truths
    #Draw labeling from metadata
    iMetadataIndex = 0

    for asMetadata in metadata:
      analysis.plotList(lsLabelList=asMetadata,strOutputFileName=str(iMetadataIndex),iSize=c_shapeSize,charForceColor='k',fInvert=c_fInvert)
      iMetadataIndex = iMetadataIndex + 1

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
        analysis.plotList(lsLabelList=lsPredictions[1:],strOutputFileName="".join([asFilePathPieces[0],"-SVMPredictions",asFilePathPieces[1]]),iSize=c_shapeSize,fInvert=c_fInvert)

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
              iSize=c_shapeSize, charForceColor=[acharColors,acharSelection], fInvert=c_fInvert)
        else:
          analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]), tempColorGrouping=acharColors,
              tempShape=acharShape, tempLabels=acharSelection, tempShapeSize=c_shapeSize, tempLegendLocation="lower left", tempInvert = c_fInvert)

    logging.info("Stop MicropitaPaperPCoA")

if __name__ == "__main__":
    _main( )

