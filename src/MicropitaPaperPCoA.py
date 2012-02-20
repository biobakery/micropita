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
import matplotlib.cm as cm
from MicroPITA import MicroPITA
import os
from PCoA import PCoA
from ValidateData import ValidateData

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperPCoA.py", 
    description = """Creates PCoA plots for MicroPITA results.""" )

#Arguments
#Abundance file
argp.add_argument( "strFileAbund", metavar = "Abundance_file", help = "An abundance table." )
#Select file
argp.add_argument( "strSelectionFile", metavar = "Select_file",
    help = "A file containing the samples selected which will be visualized." )
argp.add_argument( "sSVMPrediction", metavar = "SVM.predict", help = "The name of the input file specifying the SVm prediction." )

#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = "An optional output file" )

#Abundance table parameters
argp.add_argument("-n", dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering.")
argp.add_argument("-d", dest="iFirstDataRow", metavar= "FirstDataRow", default=1, 
                  help= "The row in the abundance file that is the first row to contain abundance data. This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata.")
argp.add_argument("-r", dest="fNormalize", metavar= "FirstDataRow", default=False, 
                  help= "Normalize the abundance data before working with it (default=False).")

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    c_fCheckFile = False
    c_NotSelected = "Not_Selected"
    c_shapeSize = 40

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable()
    abundance,metadata = rawData.textToStructuredArray(tempInputFile=args.strFileAbund, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=bool(args.fNormalize))
    sampleNames = abundance.dtype.names[1:]

    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(tempReadData=args.strFileAbund, tempIsRawData=True, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=bool(args.fNormalize), tempCheckFile=c_fCheckFile)
    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Draw known truths
    #Draw labeling from metadata
    iMetadataIndex = 0

    for asMetadata in metadata:
      plotList(objPCOA=analysis,lsLabelList=asMetadata,strName=str(iMetadataIndex),asFilePathPieces=asFilePathPieces,iSize=c_shapeSize,charForceColor='k')
      iMetadataIndex = iMetadataIndex + 1

    #Read in prediction file is supplied
    if(not args.sSVMPrediction == None):
        fHndlInput = open(args.sSVMPrediction,'r')
        strSVMSelection = fHndlInput.read()
        fHndlInput.close()

    #Make a list of predictions
    lsPredictions = list()
    for strSVMSelectionLine in filter(None,strSVMSelection.split(Constants.ENDLINE)):
        lsPredictElements = strSVMSelectionLine.split(Constants.WHITE_SPACE)
        lsPredictions.append(lsPredictElements[0])
    plotList(analysis,lsPredictions[1:],"SVMPredictions",asFilePathPieces,c_shapeSize)

    #Draw selections
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))
    for strSelectionMethod in lstrSelection:
        #Colors
        acharColors = []
        #Labels
        acharSelection = []

        #Color for selected samples in PCoA based on selection method
        charSelectedColor = ""
        astrSelectionMethod = strSelectionMethod.split(Constants.COLON)
        if astrSelectionMethod[0] == MicroPITA.c_DIVERSITY_1:
            charSelectedColor = Constants_Figures.c_charPCOADiversityInvS
        elif astrSelectionMethod[0] == MicroPITA.c_DIVERSITY_2:
            charSelectedColor = Constants_Figures.c_charPCOADiversityChao1
        elif astrSelectionMethod[0] == MicroPITA.c_EXTREME_DISSIMILARITY_1:
            charSelectedColor = Constants_Figures.c_charPCOAExtreme
        elif astrSelectionMethod[0] == MicroPITA.c_SVM_CLOSE:
            charSelectedColor = Constants_Figures.c_charPCOADiscriminant
        elif astrSelectionMethod[0] == MicroPITA.c_SVM_FAR:
            charSelectedColor = Constants_Figures.c_charPCOADistinctColor
        elif astrSelectionMethod[0] == MicroPITA.c_RANDOM:
            charSelectedColor = Constants_Figures.c_charPCOARandom
        elif astrSelectionMethod[0] == MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1:
            charSelectedColor = Constants_Figures.c_charPCOARepresentative
        elif astrSelectionMethod[0] == MicroPITA.c_USER_RANKED:
            charSelectedColor = Constants_Figures.c_charPCOATaxa

        #Parse samples
        astrSelectedSamples = astrSelectionMethod[1].split(Constants.COMMA)
        astrSelectedSamples = [strSelectedSample.strip() for strSelectedSample in astrSelectedSamples]

        #Generate colors
        for strSample in sampleNames:
            if(strSample in astrSelectedSamples):
                acharColors.append(charSelectedColor)
                acharSelection.append(astrSelectionMethod[0])
            else:
                acharColors.append(Constants_Figures.c_charPCOANoSelect)
                acharSelection.append(c_NotSelected)

        #Draw PCoA
        if astrSelectionMethod[0] in [MicroPITA.c_SVM_CLOSE, MicroPITA.c_SVM_FAR]:
          plotList(objPCOA=analysis,lsLabelList=lsPredictions[1:],strName=astrSelectionMethod[0],asFilePathPieces=asFilePathPieces,iSize=c_shapeSize, charForceColor=[acharColors,acharSelection])
        else:
          analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]), tempColorGrouping=acharColors, tempShape=acharShape, tempLabels=acharSelection, tempShapeSize=c_shapeSize, tempLegendLocation="lower left")

#charForceColor if set, automatic coloring will not occur but will occur based on the charForceColor
#CharForceColor should be a list of selection methods or not selected which will be automatically broken into colors
#charForceColor must contain the same elements as the lsLabeList
#Currently can be a list (1 color per marker in the order of the data), or 1 char to force all markers to
#charForceShapes if set, automatic shapes will not occur
#Currently can only be a char (forcing effects all markers equally)
def plotList(objPCOA,lsLabelList,strName,asFilePathPieces,iSize,charForceColor=None,charForceShape=None):
    #Get uniqueValues for labels
    acharUniqueValues = list(set(lsLabelList))
    iCountUniqueValues = len(acharUniqueValues)

    #Set colors
    atupldLabelColors = None

    #If the coloring is not forced, color so it is based on the labels
    if charForceColor == None:
      #Get colors based on labels
      atupldColors = [PCoA.RGBToHex(cm.jet(float(iUniqueValueIndex)/float(iCountUniqueValues))) for iUniqueValueIndex in xrange(0,iCountUniqueValues)]
      #Make label coloring
      atupldLabelColors = [ atupldColors[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
    #If the coloring is forced, color so it is based on the charForcedColor list
    elif(ValidateData.isValidList(charForceColor)):
      atupldLabelColors = charForceColor[0]
    #If the color is forced but the color does not vary, color all markers are the same.
    else:
      atupldLabelColors = charForceColor

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

    #Check to make sure unique colors are returned
    if(ValidateData.isValidList(atupldLabelColors)):
      if not len(acharUniqueValues)==len(list(set(atupldLabelColors))):
        print("Error, non-uniques colors were generated")
        print("Labels")
        print(lsLabelList)
        print("Colors")
        print(atupldLabelColors)
        return False
    objPCOA.plot(tempPlotName="".join([asFilePathPieces[0],"-metadata-",strName,asFilePathPieces[1]]), tempColorGrouping=atupldLabelColors, tempShape=alLabelShapes, tempLabels=lsLabelList, tempShapeSize = iSize)


if __name__ == "__main__":
    _main( )

