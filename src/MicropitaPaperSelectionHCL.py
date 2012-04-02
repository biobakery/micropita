#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
Generates an HCL of the selected samples (row) and selection method (columns).
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
from CommandLine import CommandLine
from Constants import Constants
from Constants_Arguments import Constants_Arguments
from Constants_Figures import Constants_Figures
from FileIO import FileIO
import logging
from MicroPITA import MicroPITA
import numpy as np
import os
from ValidateData import ValidateData

class MicropitaPaperSelectionHCL:
    """
    Generate HCL
    """

    def redefineMetricSelectionsAsAbsPres(self, tempSelectedSamplesDict=None, tempOutputFile=None):
        #Validate data
        if(not ValidateData.isValidDictionary(tempSelectedSamplesDict)):
          return False
        if(not ValidateData.isValidString(tempOutputFile)):
          return False

        #Get all unique samples/indexes
        totalItems = set()
        for itemsRow in tempSelectedSamplesDict:
            totalItems = totalItems.union(set(tempSelectedSamplesDict[itemsRow]))

        #Create a matrix of samples (r) by metrics (c) for presence or absence of sample in ranking
        #Row count
        sampleCount = len(totalItems)
        #Col Count
        metricsCount = len(tempSelectedSamplesDict.keys())
        #Create matrix
        selectedMatrix = np.zeros(shape=(sampleCount,metricsCount),dtype="a2")
        #Create id row (this is a column)
        rowIds = list(totalItems)
        #Create id col (this is a row)
        colIds = tempSelectedSamplesDict.keys()

        #Load the matrix based on the order of the rowIds and colIds
        for iMetric in xrange(0,len(colIds)):
            for iSample in xrange(0,len(rowIds)):
                if(rowIds[iSample] in tempSelectedSamplesDict[colIds[iMetric]]):
                    selectedMatrix[iSample][iMetric] = "1"
                else:
                    selectedMatrix[iSample][iMetric] = "0"

        #Write to file
        noWriteError = True
        output = FileIO(tempOutputFile,False,True,False)
        #Write id row
        nowriteError = noWriteError and output.writeToFile(Constants.TAB.join(["ID",Constants.TAB.join([cID for cID in colIds])])+"\n")
        #Write data
        for iRow in xrange(0,len(rowIds)):
            nowriteError = noWriteError and output.writeToFile(Constants.TAB.join([rowIds[iRow],Constants.TAB.join(selectedMatrix[iRow])])+"\n")
        #Close
        output.close()

        return noWriteError

    def generateColorFile(self, strColorFilePath, strLabelPath):
        #If the color file exists, delete
        if(os.path.exists(strColorFilePath)):
            os.remove(strColorFilePath)
        #If the label file exists, delete
        if(os.path.exists(strLabelPath)):
            os.remove(strLabelPath)
        #Create file handle to write to files
        colorFileWriter = FileIO(strColorFilePath,False,True,True)
        labelFileWriter = FileIO(strLabelPath,False,True,True)

        #Collect color and label data
        colorList = list()
        labelList = list()
        colorList.append("".join([MicroPITA.c_DIVERSITY_1+Constants.TAB+Constants_Figures.invSimpsonColor])) 
        colorList.append("".join([MicroPITA.c_DIVERSITY_2+Constants.TAB+Constants_Figures.chao1Color]))
        labelList.append("".join([MicroPITA.c_DIVERSITY_1+Constants.TAB+MicroPITA.c_DIVERSITY_1]))
        labelList.append("".join([MicroPITA.c_DIVERSITY_2+Constants.TAB+MicroPITA.c_DIVERSITY_2]))
        colorList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1+Constants.TAB+Constants_Figures.brayCurtisColor])) 
        colorList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_2+Constants.TAB+Constants_Figures.unifracColor])) 
        colorList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_3+Constants.TAB+Constants_Figures.weightedUnifracColor]))
        labelList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1+Constants.TAB+MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1]))
        labelList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_2+Constants.TAB+MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_2]))
        labelList.append("".join([MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_3+Constants.TAB+MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_3]))
        colorList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_1+Constants.TAB+Constants_Figures.inBrayCurtisColor]))
        colorList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_2+Constants.TAB+Constants_Figures.inUnifracColor]))
        colorList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_3+Constants.TAB+Constants_Figures.inWeightedUnifracColor]))
        labelList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_1+Constants.TAB+MicroPITA.c_EXTREME_DISSIMILARITY_1]))
        labelList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_2+Constants.TAB+MicroPITA.c_EXTREME_DISSIMILARITY_2]))
        labelList.append("".join([MicroPITA.c_EXTREME_DISSIMILARITY_3+Constants.TAB+MicroPITA.c_EXTREME_DISSIMILARITY_3]))
        colorList.append("".join([MicroPITA.c_USER_RANKED+Constants.TAB+Constants_Figures.userRanked]))
        labelList.append("".join([MicroPITA.c_USER_RANKED+Constants.TAB+MicroPITA.c_USER_RANKED]))
        colorList.append("".join([MicroPITA.c_RANDOM+Constants.TAB+Constants_Figures.randomColor]))
        labelList.append("".join([MicroPITA.c_RANDOM+Constants.TAB+MicroPITA.c_RANDOM]))
        colorList.append("".join([MicroPITA.c_SVM_CLOSE+Constants.TAB+Constants_Figures.svmClose]))
        colorList.append("".join([MicroPITA.c_SVM_FAR+Constants.TAB+Constants_Figures.svmFar]))
        labelList.append("".join([MicroPITA.c_SVM_CLOSE+Constants.TAB+MicroPITA.c_SVM_CLOSE]))
        labelList.append("".join([MicroPITA.c_SVM_FAR+Constants.TAB+MicroPITA.c_SVM_FAR]))

        #Close file handle
        colorFileWriter.writeToFile(Constants.ENDLINE.join(colorList))
        labelFileWriter.writeToFile(Constants.ENDLINE.join(labelList))
        colorFileWriter.close()
        labelFileWriter.close()


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperSelectionHCL.py", 
    description = """Generates an HCL of the selected samples (row) and selection method (columns).""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument( Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False",
	help = Constants_Arguments.c_strInvertHelp)
#Select file
argp.add_argument( "strSelectionFile", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp)
argp.add_argument( "strHCLLoc", metavar = "HClust_location", help = Constants_Arguments.c_strHCLLocation)
#Outputfile
argp.add_argument( "strOutHCLDataFile", metavar = "HCLData.txt", nargs = "?", help = Constants_Arguments.c_strHCLDataFile)
argp.add_argument( "strOutHCLColorFile", metavar = "HCLColor.txt", nargs = "?", help = Constants_Arguments.c_strHCLColorFile)
argp.add_argument( "strOutHCLLabelFile", metavar = "HCLLabel.txt", nargs = "?", help = Constants_Arguments.c_strHCLLabelFile)
argp.add_argument( "strOutFigure", metavar = "SelectionHCL.png", nargs = "?", help = Constants_Arguments.c_genericOutputFigureFileHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperSelectionHCL")

    mHCL = MicropitaPaperSelectionHCL()

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    #Dictionary to hold selection data
    dictSelection = dict()
    for strSelectionLine in filter(None,strSelection.split(Constants.ENDLINE)):
        astrSelectionMethod = strSelectionLine.split(Constants.COLON)
        dictSelection[astrSelectionMethod[0]] = filter(None,astrSelectionMethod[1].split(Constants.COMMA))

    if(len(dictSelection.keys())>1):
        #Figure 1b: Top ranked samples and how they are different
        mHCL.redefineMetricSelectionsAsAbsPres(tempSelectedSamplesDict=dictSelection, tempOutputFile=args.strOutHCLDataFile)
        #Generate the color and label files for HClust
        mHCL.generateColorFile(args.strOutHCLColorFile, args.strOutHCLLabelFile)
        #Hierarchical cluster absence presence metrix matrix
        CommandLine().runCommandLine([args.strHCLLoc, "--in", args.strOutHCLDataFile, "--out", args.strOutFigure, "--label2cols", args.strOutHCLColorFile, "-l", args.strOutHCLLabelFile, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1"])

    logging.info("Stop MicropitaPaperSelectionHCL")

if __name__ == "__main__":
    _main( )

