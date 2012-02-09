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

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperPCoA.py", 
    description = """Creates PCoA plots for MicroPITA results.""" )

#Arguments
#Abundance file
argp.add_argument( "strFileAbund", metavar = "Abundance_file", help = "An abundance table." )
#Select file
argp.add_argument( "strSelectionFile", metavar = "Select_file",
    help = "A file containing the samples selected which will be visualized." )
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = "An optional output file" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Analysis object
    analysis = PCoA()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    c_iNameRow = 0
    c_iFirstDataRow = 2
    c_fNormalize = True
    c_fCheckFile = False

    c_ColorScale = 100

    c_NotSelected = "Not_Selected"

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable()
    abundance,metadata = rawData.textToStructuredArray(tempInputFile=args.strFileAbund, tempDelimiter=Constants.TAB, tempNameRow=c_iNameRow, tempFirstDataRow=c_iFirstDataRow, tempNormalize=c_fNormalize)
    sampleNames = abundance.dtype.names[1:]

    #Shapes
    acharShape = Constants_Figures.c_charPCOAShape

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Generate PCoA
    #LoadData
    analysis.loadData(tempReadData=args.strFileAbund, tempIsRawData=True, tempDelimiter=Constants.TAB, tempNameRow=c_iNameRow, tempFirstDataRow=c_iFirstDataRow, tempNormalize=c_fNormalize, tempCheckFile=c_fCheckFile)
    #Make distance matrix
    pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

    #Draw known truths
    #Draw labeling from metadata
    iMetadataIndex = 0
    for asMetadata in metadata:
        #Get uniqueValues
        acharUniqueValues = list(set(asMetadata))
        iCountUniqueValues = len(acharUniqueValues)

        #Get colors
        atupldColors = [PCoA.RGBToHex(cm.jet(c_ColorScale*iUniqueValueIndex)) for iUniqueValueIndex in xrange(0,iCountUniqueValues)]

        #Make label coloring
        atupldLabelColors = [ atupldColors[acharUniqueValues.index(sMetadata)] for sMetadata in asMetadata ]

        #Plot
        iMetadataIndex = iMetadataIndex +1
        analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-metadata",str(iMetadataIndex),asFilePathPieces[1]]), tempColorGrouping=atupldLabelColors, tempShape=acharShape, tempColorLabels=asMetadata)

    #Draw selections
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))
    for strSelectionMethod in lstrSelection:
        #Colors
        acharColors = []
        #Labels
        acharLabels = []

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
                acharLabels.append(astrSelectionMethod[0])
            else:
                acharColors.append(Constants_Figures.c_charPCOANoSelect)
                acharLabels.append(c_NotSelected)

        #Draw PCoA
        analysis.plot(tempPlotName="".join([asFilePathPieces[0],"-",astrSelectionMethod[0],asFilePathPieces[1]]), tempColorGrouping=acharColors, tempShape=acharShape, tempColorLabels=acharLabels, tempLegendLocation="lower left")

if __name__ == "__main__":
    _main( )

