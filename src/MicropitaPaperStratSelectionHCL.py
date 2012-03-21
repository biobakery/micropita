#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
Generates an HCL of the study (samples = columns; taxa/otu = rows)
Highlights sample selection at the top / dendrograph
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
from Constants_Figures import Constants_Figures
import logging
from MicroPITA import MicroPITA
import os
import re

class MicropitaPaperStratSelectionHCL:
    """
    Generate HCL
    """

    def generateColorFile(self, dictSelection, strColorFilePath, strLabelPath):
        #If the color file exists, delete
        if(os.path.exists(strColorFilePath)):
            os.remove(strColorFilePath)
        #If the label file exists, delete
        if(os.path.exists(strLabelPath)):
            os.remove(strLabelPath)

        #Collect color and label data
        #Color samples (column) based on selection method
        colorList = list()
        labelList = list()

        dictSwitch = {MicroPITA.c_DIVERSITY_1:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.invSimpsonColor])),
                      MicroPITA.c_DIVERSITY_2:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.chao1Color])),
                      MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.brayCurtisColor])),
                      MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_2:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.unifracColor])),
                      MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_3:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.weightedUnifracColor])),
                      MicroPITA.c_EXTREME_DISSIMILARITY_1:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.inBrayCurtisColor])),
                      MicroPITA.c_EXTREME_DISSIMILARITY_2:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.inUnifracColor])),
                      MicroPITA.c_EXTREME_DISSIMILARITY_3:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.inWeightedUnifracColor])),
                      MicroPITA.c_USER_RANKED:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.userRanked])),
                      MicroPITA.c_RANDOM:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.randomColor])),
                      MicroPITA.c_SVM_CLOSE:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.svmClose])),
                      MicroPITA.c_SVM_FAR:lambda lsToAdd,strSample: lsToAdd.append(Constants.TAB.join([strSample,Constants_Figures.svmFar]))}

        #Indicate color for samples based on selection
        #What about multiple selection of the same sample????????? #TODO
        for sSelectionKey in dictSelection:
          for sSample in dictSelection[sSelectionKey]:
            dictSwitch[sSelectionKey](colorList,sSample)
            dictSwitch[sSelectionKey](labelList,sSample)

        #Close file handle
        #Create file handle to write colors and labels to file
        with open(strColorFilePath, 'w') as f:
          f.write(Constants.ENDLINE.join(colorList))
          f.close()
        with open(strLabelPath, 'w') as f:
          f.write(Constants.ENDLINE.join(labelList))
          f.close()

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperStratSelectionHCL.py", 
    description = """Generates an HCL of the study (samples = columns; taxa/otu = rows). Highlights sample selection at the top / dendrograph.""" )

#Arguments
#Logging
argp.add_argument("-l", dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], 
                  help= "Logging level which will be logged to a .log file with the same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
argp.add_argument( "-d", dest = "strFirstDataRow", action = "store", default="2",
	help = "The first Row that has data (starting from the top with 0)." )
argp.add_argument( "-id", dest = "strIDCol", action = "store", default="0",
	help = "The column with Taxa / OTU ids (starting with 0)." )
argp.add_argument( "-i", dest = "fInvert", action = "store", default="False",
	help = "Invert the image to a black background (default=False)." )
#Select file
argp.add_argument( "strSelectionFile", metavar = "Select_file",
    help = "A file containing the samples selected which will be visualized." )
argp.add_argument( "strAbundanceFile", metavar = "Abundance_file",
    help = "A file containing the sample abundances." )
#HCLUST
argp.add_argument( "strHCLLoc", metavar = "HClust_location",
    help = "The location to HClust (for example ./external/hclust/hclust.py)." )
#Outputfile
argp.add_argument( "strOutHCLColorFile", metavar = "HCLColor.txt", nargs = "?", help = "An output file that is used by HClust. This is the color file." )
argp.add_argument( "strOutHCLLabelFile", metavar = "HCLLabel.txt", nargs = "?", help = "An output file that is used by HClust. This is the label file." )
argp.add_argument( "strOutFigure", metavar = "StratSelectionHCL.png", nargs = "?", help = "The output HCL figure." )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperStratSelectionHCL")
    logging.info(str(args))
    mStratHCL = MicropitaPaperStratSelectionHCL()

    #Read in selection file
    strSelection = ""
    with open(args.strSelectionFile,'r') as fHndlInput:
      strSelection = fHndlInput.read()
      fHndlInput.close()

    #Dictionary to hold selection data
    dictSelection = dict()
    for strSelectionLine in filter(None,strSelection.split(Constants.ENDLINE)):
        astrSelectionMethod = strSelectionLine.split(Constants.COLON)
        dictSelection[astrSelectionMethod[0].split()[0]] = [ strSample.split()[0] for strSample in filter(None,astrSelectionMethod[1].split(Constants.COMMA))]

    if(len(dictSelection.keys())>1):
        #Generate the color and label files for HClust
        mStratHCL.generateColorFile(dictSelection, args.strOutHCLColorFile, args.strOutHCLLabelFile)
        #Hierarchical cluster the sample
        CommandLine().runCommandLine([args.strHCLLoc, "--in", args.strAbundanceFile, "--out", args.strOutFigure, "-d", "euclidean", "--label2cols", args.strOutHCLColorFile, "-l", args.strOutHCLLabelFile, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1", "-s", "log", "--ystart", args.strFirstDataRow, "--xstart", str(int(args.strIDCol)+1)])

    logging.info("Stop MicropitaPaperStratSelectionHCL")

if __name__ == "__main__":
    _main( )

