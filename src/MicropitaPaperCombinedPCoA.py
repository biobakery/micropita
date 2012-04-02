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
argp = argparse.ArgumentParser( prog = "MicropitaPaperCombinedPCoA.py", description = """Creates Combined PCoA plots for MicroPITA results.""" )
#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strSampleNameRowArgument, dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= Constants_Arguments.c_strSampleNameRowHelp)
argp.add_argument(Constants_Arguments.c_strFirstDataRow, dest="iFirstDataRow", metavar= "FirstDataRow", default=1, help= Constants_Arguments.c_strFirstDataRowHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False", help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
#Selection Methods to plot, max 4
#Abundance file
argp.add_argument( "strFileAbund", action="store", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Select file
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.png", help = Constants_Arguments.c_genericOutputFigureFileHelp)
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

    logging.info("Start MicropitaPaperCombinedPCoA")

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    #Parse possible selection methods
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))

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

    print("sampleNames")
    print(sampleNames)
    print("args.strSelectionMethods")
    print(args.strSelectionMethods)

    #Figure colors
    objColors = Constants_Figures()
    objColors.invertColors(fInvert=c_fInvert)

    #Shape for a sample that is selected twice
    acharMultSelectShape = objColors.c_charPCOAMultSelectionShape

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(tempReadData=args.strFileAbund, tempIsRawData=True, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=c_Normalize, tempCheckFile=c_fCheckFile)
    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Colors
    acharColors = [objColors.c_charNoSelect for iindex in xrange(len(sampleNames))]
    #Labels
    acharSelection = [c_NotSelected for iindex in xrange(len(sampleNames))]
    #Shapes
    acharShapes = [objColors.c_charPCOAShape for iindex in xrange(len(sampleNames))]

    #Selection methods
    lsSelectionMethods = []

    #Draw selections
    for strSelectionMethod in lstrSelection:

        #Color for selected samples in PCoA based on selection method
        astrSelectionMethod = strSelectionMethod.split(Constants.COLON)
        sCurSelectionMethodName = astrSelectionMethod[0]
        print("sCurSelectionMethodName")
        print(sCurSelectionMethodName)

        #If the method parsed fromthe selection file is a method that is passed in as an argument and indicated as a method to plot
        if sCurSelectionMethodName in args.strSelectionMethods:

            #Get the correct color for the method
            charSelectedColor = objColors.dictConvertMethodToHEXColor[sCurSelectionMethodName]

            #Parse samples selected by the method
            astrSelectedSamples = astrSelectionMethod[1].split(Constants.COMMA)
            astrSelectedSamples = [strSelectedSample.strip() for strSelectedSample in astrSelectedSamples]
            print("astrSelectedSamples")
            print(astrSelectedSamples)
            #Indicate colors
            for iindex in xrange(len(sampleNames)):
                #Go through samples and see if the sample was selected by the method
                print("sampleNames[iindex]")
                print(sampleNames[iindex])
                if(sampleNames[iindex] in astrSelectedSamples):
                    print("IS IN")
                    # If it was selected by the method check to see if the color has already changed from
                    # the default, if it has then the sample was selected twice or more and should be
                    # indicated by shape for caution. So change shape.
                    if not acharColors[iindex] == objColors.c_charNoSelect:
                        acharShapes[iindex] = acharMultSelectShape
                    acharColors[iindex]=charSelectedColor
                    acharSelection[iindex]=sCurSelectionMethodName
                print("acharColors[iindex]")
                print(acharColors[iindex])
                print("acharSelection[iindex]")
                print(acharSelection[iindex])

    #Draw PCoA
    print("acharSelection")
    print(acharSelection)
    print(len(acharSelection))
    print("acharShapes")
    print(acharShapes)
    print(len(acharShapes))
    analysis.plotList(lsLabelList=acharSelection, strName="-".join(sCurSelectionMethodName), asFilePathPieces=asFilePathPieces, iSize=c_shapeSize, charForceColor=[acharColors,acharSelection], charForceShape=acharShapes, fInvert=c_fInvert)
#        else:
#          analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]), tempColorGrouping=acharColors, tempShape=acharShape, tempLabels=acharSelection, tempShapeSize=c_shapeSize, tempLegendLocation="lower left", tempInvert = c_fInvert)

    logging.info("Stop MicropitaPaperCombinedPCoA")

if __name__ == "__main__":
    _main( )

