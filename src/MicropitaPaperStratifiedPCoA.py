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
from Constants_Figures import Constants_Figures
import logging
import matplotlib.cm as cm
from MicroPITA import MicroPITA
import os
from PCoA import PCoA
from ValidateData import ValidateData

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperStratifiedPCoA.py", description = """Creates PCoA plots for stratified MicroPITA results.""" )
#Arguments
#Logging
argp.add_argument("-l", dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], 
                  help= "Logging level which will be logged to a .log file with the same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
argp.add_argument("-n", dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering.")
argp.add_argument("-d", dest="iFirstDataRow", metavar= "FirstDataRow", default=1, 
                  help= "The row in the abundance file that is the first row to contain abundance data. This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata.")
argp.add_argument("-u", dest="strUnsupervisedStratify", metavar= "UnsupervisedStratify", default=None, 
                  help= "The metatdata to stratify unsupervised analysis.")
argp.add_argument( "-r", dest = "fNormalize", action = "store", default="False",
	help = "Normalize the abundance data before working with it (default=False)." )
argp.add_argument( "-i", dest = "fInvert", action = "store", default="False",
	help = "Invert the image to a black background (default=False)." )
#Abundance file
argp.add_argument( "strFileAbund", action="store", metavar = "Abundance_file", help = "An abundance table." )
#Select file
argp.add_argument( "strSelectionFile", action="store", metavar = "Select_file",
    help = "A file containing the samples selected which will be visualized." )

#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = "An optional output file" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperStratifiedPCoA")

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

        #Draw Stratified PCoA
        plotList(objPCOA=analysis,lsLabelList=metadata[args.strUnsupervisedStratify],strName="-".join([args.strUnsupervisedStratify,astrSelectionMethod[0]]),asFilePathPieces=asFilePathPieces,iSize=c_shapeSize, charForceColor=[acharColors,acharSelection], fInvert=c_fInvert)

    logging.info("Stop MicropitaPaperStratifiedPCoA")

#charForceColor if set, automatic coloring will not occur but will occur based on the charForceColor
#CharForceColor should be a list of selection methods or not selected which will be automatically broken into colors
#charForceColor must contain the same elements as the lsLabeList
#Currently can be a list (1 color per marker in the order of the data), or 1 char to force all markers to
#charForceShapes if set, automatic shapes will not occur
#Currently can only be a char (forcing effects all markers equally)
def plotList(objPCOA,lsLabelList,strName,asFilePathPieces,iSize,charForceColor=None,charForceShape=None, fInvert=False):
    #Get uniqueValues for labels
    acharUniqueValues = list(set(lsLabelList))
    iCountUniqueValues = len(acharUniqueValues)

    #Set colors
    atupldLabelColors = None

    #Set shapes
    alLabelShapes = None
    if charForceShape == None:
      #Get shapes
      acharShapes = PCoA.getShapes(iCountUniqueValues)
      if acharShapes == None:
        return False
      #Make label shapes
      alLabelShapes = [ acharShapes[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
    else:
      alLabelShapes = acharShapes

    #If the coloring is not forced, color so it is based on the labels
    if charForceColor == None:
      #Get colors based on labels
      atupldColors = [PCoA.RGBToHex(cm.jet(float(iUniqueValueIndex)/float(iCountUniqueValues))) for iUniqueValueIndex in xrange(0,iCountUniqueValues)]
      #Make label coloring
      atupldLabelColors = [ atupldColors[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
    #If the coloring is forced, color so it is based on the charForcedColor list
    elif(ValidateData.isValidList(charForceColor)):
      atupldLabelColors = charForceColor[0]
      if not len(lsLabelList) == len(charForceColor[1]):
        logging.error("MicropitaPaperStratifiedPCoA.plotList. Label and forced color lengths were not the same.") 
        return False
      lsLabelList = [ "".join([charForceColor[1][iLabelIndex], "_", lsLabelList[iLabelIndex]]) for iLabelIndex in xrange(0,len(charForceColor[1]))]
    #If the color is forced but the color does not vary, color all markers are the same.
    else:
      atupldLabelColors = charForceColor

    #Check to make sure unique colors are returned
    if(ValidateData.isValidList(atupldLabelColors)):
      if not len(acharUniqueValues)==len(list(set(atupldLabelColors))):
        logging.error("MicropitaPaperStratifiedPCoA.plotList. Non-uniques colors were generated for "+strName)
        logging.debug("Labels")
        logging.debug(lsLabelList)
        logging.debug("Colors")
        logging.debug(atupldLabelColors)
        return False

    logging.debug("lsLabelList")
    logging.debug(lsLabelList)
    objPCOA.plot(tempPlotName="".join([asFilePathPieces[0],"-metadata-",strName,asFilePathPieces[1]]), tempColorGrouping=atupldLabelColors, tempShape=alLabelShapes, tempLabels=lsLabelList, tempShapeSize = iSize, tempInvert = fInvert)


if __name__ == "__main__":
    _main( )

