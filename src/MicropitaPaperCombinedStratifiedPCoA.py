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
argp = argparse.ArgumentParser( prog = "MicropitaPaperCombinedStratifiedPCoA.py", description = """Creates PCoA plots for stratified MicroPITA results (Both supervised and unsupervised).""" )
#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, 
                  help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strUnsupervisedStratifyMetadataArgument, dest="strUnsupervisedStratify", metavar= "UnsupervisedStratify", default=None, 
                  help= Constants_Arguments.c_strUnsupervisedStratifyMetadataHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False", help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strSumDataArgument, dest="fSumData", action = "store", metavar= "WouldlikeDataSummed", help= Constants_Arguments.c_strSumDataHelp)

#Abundance file
argp.add_argument( "strFileAbund", action="store", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Select file
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", help = Constants_Arguments.c_strOptionalOutputDataFileHelp)
#Selection parameter
argp.add_argument("strSelectionMethods", metavar = "Selection_Methods", help = Constants_Arguments.c_strSelectionMethodsHelp, nargs="*")

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperCombinedStratifiedPCoA")

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    c_fCheckFile = False


    c_fInvert = (args.fInvert == "True")
    c_Normalize = (args.fNormalize == "True")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fSumData = (args.fSumData.lower() == "true")

    #Figure standardization
    objFigureControl = Constants_Figures()
    objFigureControl.invertColors(fInvert=c_fInvert)
    c_NotSelected = objFigureControl.c_strPCOANotSelected
    dAlpha = objFigureControl.c_dAlpha

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable.makeFromFile(strInputFile=args.strFileAbund, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)

    #Normalize if needed and sum clades
    if fSumData:
        rawData.funcSumClades()
    if c_Normalize:
        rawData.funcNormalize()

    sampleNames = rawData.funcGetSampleNames()

    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(xData=rawData, fIsRawData=True)

    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Colors
    acharColors = [objFigureControl.c_charNoSelect for iindex in xrange(len(sampleNames))]
    #Labels
    acharSelection = [c_NotSelected for iindex in xrange(len(sampleNames))]
    #Labels for shape
    astrSelectionLabels = rawData.funcGetMetadata(args.strUnsupervisedStratify)
    #Sizes
    aiSizes = [objFigureControl.iMarkerSize for iindex in xrange(len(sampleNames))]

    #Draw selections
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))
    for strSelectionMethod in lstrSelection:

        #Color for selected samples in PCoA based on selection method
        charSelectedColor = ""
        astrSelectionMethod = strSelectionMethod.split(Constants.COLON)
        sCurSelectionMethodName = astrSelectionMethod[0]

        #If the method parsed from the selection file is a method that is passed in as an argument and indicated as a method to plot
        if sCurSelectionMethodName in args.strSelectionMethods:
            #Get color to draw
            charSelectedColor = objFigureControl.dictConvertMethodToHEXColor[sCurSelectionMethodName]

            #Parse samples
            astrSelectedSamples = astrSelectionMethod[1].split(Constants.COMMA)
            astrSelectedSamples = [strSelectedSample.strip() for strSelectedSample in astrSelectedSamples]
            #Indicate colors
            for iindex in xrange(len(sampleNames)):

                #Go through samples and see if the sample was selected by the method
                if(sampleNames[iindex] in astrSelectedSamples):
                    # If it was selected by the method check to see if the color has already changed from
                    # the default, if it has then the sample was selected twice or more and should be
                    # indicated by a pie chart marker which is a list of colors with the first 
                    # element a char indicating the polygon type to be plotted
                    if not acharColors[iindex] == objFigureControl.c_charNoSelect:
                        #Update colors (this is where the pie cut marker is indicated to the PCoA software)
                        curColor = acharColors[iindex]
                        if not ValidateData.isValidList(curColor):
                          acharColors[iindex] = [curColor]
                        acharColors[iindex].append(charSelectedColor)
                        #Update selection which will be used to name selection
                        acharSelection[iindex]=sCurSelectionMethodName
                        #Update size
                        aiSizes[iindex] = objFigureControl.iPieMarkerSize
                    else:
                        acharColors[iindex]=charSelectedColor
                        acharSelection[iindex]=sCurSelectionMethodName

    #Draw Stratified PCoA
    analysis.plotList(lsLabelList=astrSelectionLabels,strOutputFileName=args.strOutFile,
        iSize=aiSizes, dAlpha=dAlpha, charForceColor=[acharColors,acharSelection], fInvert=c_fInvert)
    logging.info("Stop MicropitaPaperCombinedStratifiedPCoA")

if __name__ == "__main__":
    _main( )

