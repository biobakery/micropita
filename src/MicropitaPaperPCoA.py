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
from MicroPITA import MicroPITA
import os
from PCoA import PCoA

class MicropitaPaperPCoA:
    """
    Generate PCoA plots
    """

    def funcPCoAMethods(self, strInputFile = None, strOutputFile = None, acharColors = None, charShape = None, fRawData = True, charDelimiter = Constants.TAB, iNameRow = 0, firstDataRow = 1, fNormalize = True, fCheckFile = True):
        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=strInputFile, tempIsRawData=fRawData, tempDelimiter=charDelimiter, tempNameRow=iNameRow, tempFirstDataRow=firstDataRow, tempNormalize=fNormalize, tempCheckFile=fCheckFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)
        #Plot
        analysis.plot(tempPlotName=strOutputFile, tempColorGrouping=acharColors, tempShape=charShape)

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

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    #Read abundance file
    #Abundance table object to read in and manage data
    rawData = AbundanceTable()
    abundance,metadata = rawData.textToStructuredArray(tempInputFile=args.strFileAbund, tempDelimiter=Constants.TAB, tempNameRow=0, tempFirstDataRow=2, tempNormalize=True)
    sampleNames = abundance.dtype.names[1:]

    #File path components
    asFilePathPieces = os.path.splitext(args.strOutFile)

    #Draw known truths
    #Draw labeling
    for asMetadata in metadata:
        
#    print("metadata")
#    print(metadata)

    #Draw selections
    lstrSelection =  filter(None,strSelection.split(Constants.ENDLINE))
    for strSelectionMethod in lstrSelection:

        #Colors
        acharColors = []

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
            else:
                acharColors.append(Constants_Figures.c_charPCOANoSelect)

        #Shapes
        acharShape = Constants_Figures.c_charPCOAShape

        #Generate PCoA
        MicropitaPaperPCoA().funcPCoAMethods(strInputFile = args.strFileAbund, strOutputFile = "".join([asFilePathPieces[0],"-",astrSelectionMethod[0],"-",asFilePathPieces[1]]), acharColors = acharColors, charShape = acharShape, fRawData = True, charDelimiter = Constants.TAB, iNameRow = 0, firstDataRow = 2, fNormalize = False, fCheckFile = False)

if __name__ == "__main__":
    _main( )

