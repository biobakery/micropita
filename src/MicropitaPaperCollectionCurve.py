#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
Generates a collection curve of the Diversity / Top ranked taxa of a sample set (y)
given the number of samples in the set.
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
from pylab import *
from MicroPITA import MicroPITA
import os
import random
import re

class MicropitaPaperCollectionCurve:

    #TODO complete last
    def funcPlotCollectionCurve(self, strPlotName, dictMethods, dictdMaxDiversity, dictdDiversityBaseline, strMetric, fInvert):
        font = {'family':'arial',
           'color':'k',
           'weight':'normal',
           'size':12}

        #Get plot colors
        objColors = Constants_Figures()
        objColors.invertColors(fInvert=fInvert)

        #Get plot object
        imgFigure = plt.figure()
        imgFigure.set_facecolor(objColors.c_strBackgroundColorWord)
        imgSubplot = imgFigure.add_subplot(111,axisbg=objColors.c_strBackgroundColorLetter)
        imgSubplot.set_xlabel('Number of samples sampled from study (n)')
        imgSubplot.set_ylabel('Sampled Diversity (Inverse Simpson)')
        imgSubplot.spines['top'].set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.spines['bottom'].set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.spines['left'].set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.spines['right'].set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.xaxis.label.set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.yaxis.label.set_color(objColors.c_strDetailsColorLetter)
        imgSubplot.tick_params(axis='x', colors=objColors.c_strDetailsColorLetter)
        imgSubplot.tick_params(axis='y', colors=objColors.c_strDetailsColorLetter)
        charMarkerEdgeColor = objColors.c_strDetailsColorLetter

        iYMin = 0
        iYMax = 0
        iXMin = 0
        iXMax = 0

        liX = list()
        liY = list()
        #Plot max diversity
        for iCount in dictdMaxDiversity:
            liX.append(iCount)
            liY.append(dictdMaxDiversity[iCount])
        plot(liX, liY, color=objColors.c_strDetailsColorLetter, marker="*", linestyle='--', label = "Maximum "+strMetric)

        #Update range
        iYMin = min(liY)
        iYMax = max(liY)
        iXMin = min(liX)
        iXMax = max(liX)

        #Plot bootstrapped diversity
        liX = list()
        liY = list()
        for iCount in dictdDiversityBaseline:
            liX.append(iCount)
            liY.append(dictdDiversityBaseline[iCount])
        plot(liX, liY, color=objColors.c_strDetailsColorLetter, marker=".", linestyle='--', label = "Bootstrapped "+strMetric)

        #Update range
        iYMin = min(liY+[iYMin])
        iYMax = max(liY+[iYMax])
        iXMin = min(liX+[iXMin])
        iXMax = max(liX+[iXMax])

        #Plot methods
        for strMethod in dictMethods:
            dictCurMethod = dictMethods[strMethod]
            dictCurMetric = dictCurMethod[strMetric]
            liX = list()
            liY = list()
            for iCount in dictCurMetric:
                liX.append(iCount)
                liY.append(dictCurMetric[iCount])

            plot(liX, liY, color = objColors.dictConvertMethodToHEXColor[strMethod], marker = "o", linestyle='-', label = strMethod)

            #Update range
            iYMin = min(liY+[iYMin])
            iYMax = max(liY+[iYMax])
            iXMin = min(liX+[iXMin])
            iXMax = max(liX+[iXMax])

        xlim(iXMin*.9,iXMax*1.1)
        ylim(iYMin*.9,iYMax*1.1)
        objLegend = imgSubplot.legend(loc="upper right", scatterpoints=1, prop={'size':10})

        #Invert legend
        if(fInvert):
            objLegend.legendPatch.set_fc(objColors.c_strBackgroundColorWord)
            objLegend.legendPatch.set_ec(objColors.c_strDetailsColorLetter)
            plt.setp(objLegend.get_texts(),color=objColors.c_strDetailsColorLetter)

        imgFigure.savefig(strPlotName, facecolor=imgFigure.get_facecolor())

    def getAverageBootstrappedMetric(self, lsMetrics, iSelectSampleCount, iBootStrappingItr):
        ldAveragePerIteration = list()
        for iItr in xrange(iBootStrappingItr):
            ldSelection = [random.choice(lsMetrics) for iSampleCount in xrange(iSelectSampleCount)]
            ldAveragePerIteration.append(sum(ldSelection)/float(len(ldSelection)))
        return sum(ldAveragePerIteration)/float(len(ldAveragePerIteration))

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperCollectionCurve.py", 
    description = """Generates a collection curve of diversity and top ranked abudance Taxa/OTU given method and number of samples selected in studies.""" )

#Arguments
#Logging
argp.add_argument("-l", dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], 
                  help= "Logging level which will be logged to a .log file with the same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
argp.add_argument("-n", dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering.")
argp.add_argument("-d", dest="iFirstDataRow", metavar= "FirstDataRow", default=1, 
                  help= "The row in the abundance file that is the first row to contain abundance data. This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata.")
argp.add_argument( "-r", dest = "fNormalize", action = "store", default="False",
	help = "Normalize the abundance data before working with it (default=False)." )
argp.add_argument( "-i", dest = "fInvert", action = "store", default="False",
	help = "Invert the image to a black background (default=False)." )
#Select file
argp.add_argument( "strAbundanceFile", metavar = "Abundance_file",
    help = "A file containing the sample abundances." )

#Outputfile
argp.add_argument( "strOutFigure", metavar = "CollectionCurveOutputFile", help = "The output collection curve figure." )
argp.add_argument( "strOutText", metavar = "CollectionCurveOutputTextFile", help = "The data from collection curve figure." )
argp.add_argument( "strSelectionFiles", metavar = "InputSelectionFiles", nargs = "+", help = "A list of files that contain samples selection based on a specific, common abundance table." )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )
    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperCollectionCurve")
    logging.info("MicropitaPaperCollectionCurve. The following arguments were passed.")
    logging.info(str(args))

    #Normalize Abundance data
    fNormalize = (args.fNormalize.lower() == "true")

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Bootstrap interations
    c_BootstrapItr = 1000

    #Micropita instance
    microPITA = MicroPITA()

    #TODO maybe set up for multiple dviersity metrics at once
    #Diversity metric
    lsDiversityMetric = [microPITA.c_INV_SIMPSON_A_DIVERSITY]

    #Instance of plot collection curve script
    mCC = MicropitaPaperCollectionCurve()

    #Manage the cases where there are no or 1 selection file given
    if args.strSelectionFiles == None:
      logging.warning("MicropitaPaperCollectionCurve. No files were provided indicating sample selection so the collection curve was not made.")
      return
    if isinstance(args.strSelectionFiles, basestring):
      args.strSelectionFiles = [args.strSelectionFiles]

    #Get diversity of each sample in the abundance table
    totalData = AbundanceTable()
    rawAbundance,metadata = totalData.textToStructuredArray(tempInputFile=args.strAbundanceFile, tempDelimiter=Constants.TAB, 
                                                            tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow),
                                                            tempNormalize=fNormalize)
    #List of lists, one list per diversity metric for all samples
    lsSampleNames = rawAbundance.dtype.names[1:]
    lStudyDiversityMetrics = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = rawAbundance, tempSampleNames = lsSampleNames, tempDiversityMetricAlpha = lsDiversityMetric)
    lStudyDiversityMetrics = lStudyDiversityMetrics[0]

    #Get top ranked of each sample

    #Read through selection files and place contents in dictionary
    dictGrandSelection = dict()
    for strFile in args.strSelectionFiles:
        #Read in selection file
        strSelection = ""
        with open(strFile,'r') as fHndlInput:
            strSelection = fHndlInput.read()
            fHndlInput.close()

        #Dictionary to hold selection data
        dictSelection = dict()
        for strSelectionLine in filter(None,strSelection.split(Constants.ENDLINE)):
            astrSelectionMethod = strSelectionLine.split(Constants.COLON)
            dictSelection[astrSelectionMethod[0].split()[0]] = [strSample.split()[0] for strSample in filter(None,astrSelectionMethod[1].split(Constants.COMMA))]

        #Put sample selection into grand study selection
        if strFile in dictGrandSelection:
          logging.error("".join(["MicroPitaPaperCollectionCurve. File for study was already collected in the collection curve. The following file was ignored:",strFile,"."]))
        dictGrandSelection[strFile] = dictSelection

    #Break up the GrandSelection of subsampling by sample selection size
    #For each selection size
    #Get the average y for each sample subset (y={diversity, top ranked abundance})
    dictMetricsBySampleN = dict()
    iSelectionCount = -1
    fError = False

    print("dictGrandSelection")
    print(dictGrandSelection)

    #For each subsample (N) evaluate selection methods
    for dictStudy in dictGrandSelection:
        dictCurStudy = dictGrandSelection[dictStudy]
        iSelectionCount = -1
        #For each selection method make sure the other selection in the study selected the same N
        lsMethodStats = list()
        for strMethod in dictCurStudy:
            lsCurSampleSelections = dictCurStudy[strMethod]
            if iSelectionCount == -1:
                iSelectionCount = len(lsCurSampleSelections)
            elif not iSelectionCount == len(lsCurSampleSelections):
                logging.error("".join(["MicroPitaPaperCollectionCurve. Selection methods selected an uneven number of samples. Did not include this study in the collection curve.",dictStudy,"."]))
                fError = True
            #If there was an error do not add the study's methods to plot data
            if fError: break
            else:
                #Calculate diversity
                lsCurSampleDiversity = list()
                for strSample in lsCurSampleSelections:
                    lsCurSampleDiversity.append(lStudyDiversityMetrics[lsSampleNames.index(strSample)])
                dDiversity = sum(lsCurSampleDiversity)/float(len(lsCurSampleDiversity))
                #Caluclate top ranked
                dTopRanked = -1
                lsMethodStats.append([strMethod,["Diversity",dDiversity],["TopRank",dTopRanked]])

        #Build data structure to hold plotting data
        if not fError:
            for ldStats in lsMethodStats:
                if not ldStats[0] in dictMetricsBySampleN:
                    dictMetricsBySampleN[ldStats[0]]=dict()
                dictCurMethod = dictMetricsBySampleN[ldStats[0]]
                for lMetric in ldStats[1:]:
                    if not lMetric[0] in dictCurMethod:
                        dictCurMethod[lMetric[0]] = dict()
                    dictCurMetric = dictCurMethod[lMetric[0]]
                    if not iSelectionCount in dictCurMethod:
                        dictCurMetric[iSelectionCount] = list()
                    dictCurMetric[iSelectionCount].append(lMetric[1])

    print("lsMethodStats")
    print(lsMethodStats)

    #Most diverse set of samples  at each N
    dictdDiverseSamplesAtN = dict()
    #Booststrapped diversity level at N
    dictdBootstrappedDiversityAtN = dict()

    #Average multiple method metric instances at sample size N (if they exist)
    #N level
    for strMethod in dictMetricsBySampleN:
        dictCurMethod = dictMetricsBySampleN[strMethod]
        #Method level
        for strMethodMetric in dictCurMethod:
            dictMethodMetric = dictCurMethod[strMethodMetric]
            #Metric level
            for strN in dictMethodMetric:
                lsMethodMetricByN = dictMethodMetric[strN]
                #Average the list of metrics
                dictMethodMetric[strN] = sum(lsMethodMetricByN)/float(len(lsMethodMetricByN))

                #Get Most diverse sample per N
                if strN not in dictdDiverseSamplesAtN:
                    ldTopSamples = sorted(lStudyDiversityMetrics)[-1*int(strN):]
                    dictdDiverseSamplesAtN[strN] = sum(ldTopSamples)/float(len(ldTopSamples))

                #Get bootstrapped diversity at a level N
                if strN not in dictdBootstrappedDiversityAtN:
                    dictdBootstrappedDiversityAtN[strN] = mCC.getAverageBootstrappedMetric(lsMetrics = lStudyDiversityMetrics, iSelectSampleCount = int(strN), iBootStrappingItr = c_BootstrapItr)

    print(dictMetricsBySampleN)
    print(dictMetricsBySampleN)

    #Plot line graph
    mCC.funcPlotCollectionCurve(strPlotName=args.strOutFigure, dictMethods=dictMetricsBySampleN, dictdMaxDiversity=dictdDiverseSamplesAtN, dictdDiversityBaseline=dictdBootstrappedDiversityAtN, strMetric="Diversity", fInvert = fInvert)

    logging.info("Stop MicropitaPaperCollectionCurve")

if __name__ == "__main__":
    _main( )

