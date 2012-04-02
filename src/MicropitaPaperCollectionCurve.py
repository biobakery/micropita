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
from Diversity import Diversity
import logging
from pylab import *
from MicroPITA import MicroPITA
import numpy as np
import os
import random
import re

class MicropitaPaperCollectionCurve:

    c_UseDiversityMetrics = "Diversity"
    c_UseTopRankedMethods = "TopRanked"

    def funcPlotCollectionCurve(self, strPlotName, dictMethods, dictdMaxYMetric, dictdDiversityBaseline, strMetric, fInvert):
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
        imgSubplot.set_ylabel("".join(["Sampled Diversity (",strMetric,")"]))
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
        #Plot max y metric
        for iCount in dictdMaxYMetric:
            liX.append(iCount)
            liY.append(dictdMaxYMetric[iCount])
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

            plot(liX, liY, color = objColors.dictConvertMethodToHEXColor[strMethod], marker = "o", alpha=objColors.c_dAlpha, linestyle='-', label = strMethod)

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
        #Make legend background transparent
        objLegendFrame = objLegend.get_frame()
        objLegendFrame.set_alpha(objColors.c_dAlpha)

        imgFigure.savefig(strPlotName, facecolor=imgFigure.get_facecolor())

    #First a subselection occurs and pool. 
    def getMedianBootstrappedMetric(self, npaAbundance, lsSampleNames, iSelectSampleCount, sMetric, iBootStrappingItr):
        ldMeasurePerIteration = list()
        iSampleCount = len(lsSampleNames)
        for iItr in xrange(iBootStrappingItr):
            #Select population
            liSelection = random.sample(xrange(iSampleCount),iSelectSampleCount)
            lfSelection = [(iSelectIndex in liSelection) for iSelectIndex in xrange(iSampleCount)]

            #When combining combine counts (Chao1) by summing and (Relative abundance) by average
            if(sMetric == Diversity.c_INV_SIMPSON_A_DIVERSITY):
                ldPooledSample = list()
                for featureAbundance in npaAbundance:
                    ldSelectedAbundance = np.compress(lfSelection,list(featureAbundance)[1:])
                    ldPooledSample.append(sum(ldSelectedAbundance)/(len(ldSelectedAbundance)*1.0))
                if sum(ldPooledSample)==0.0:
                    ldMeasurePerIteration.append(0)
                else:
                    ldMeasurePerIteration.append(Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=np.array(ldPooledSample)))
            #Do not normalize
            elif(sMetric == Diversity.c_CHAO1_A_DIVERSITY):
                ldPooledSample = list()
                for featureAbundance in npaAbundance:
                    ldPooledSample.append(sum(np.compress(lfSelection, list(featureAbundance)[1:])))
                ldMeasurePerIteration.append(Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies=np.array(ldPooledSample)))
        return np.median(ldMeasurePerIteration)

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

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Bootstrap interations
    c_BootstrapItr = 1000

    #Micropita instance
    microPITA = MicroPITA()

    #Instance of plot collection curve script
    mCC = MicropitaPaperCollectionCurve()

    #Metric
    c_strMetricCategory = mCC.c_UseDiversityMetrics

    #Diversity metric
    lsDiversityRunMetrics = [microPITA.c_INV_SIMPSON_A_DIVERSITY,microPITA.c_CHAO1_A_DIVERSITY]

    #Manage the cases where there are no or 1 selection file given
    if args.strSelectionFiles == None:
      logging.warning("MicropitaPaperCollectionCurve. No files were provided indicating sample selection so the collection curve was not made.")
      return
    if isinstance(args.strSelectionFiles, basestring):
      args.strSelectionFiles = [args.strSelectionFiles]

    #The metrics within the metric category to run
    #For example, the category could be Diversity and the runMetrics could be Inverse Simpson and Choa1
    lsRunMetrics = []
    if c_strMetricCategory == mCC.c_UseDiversityMetrics:
        lsRunMetrics = lsDiversityRunMetrics
    elif c_strMetricCategory == mCC.c_UseTopRankedMethods:
        lsRunMetrics.append(mCC.c_UseTopRankedMethods)

    #For each run metric build a collection curve. The run metric will be the y axis.
    for strRunMetric in lsRunMetrics:

        #Do not normalize Choa1
        fNormalize = (not strRunMetric == Diversity.c_CHAO1_A_DIVERSITY)

        #Get diversity of each sample in the abundance table
        totalData = AbundanceTable()
        rawAbundance,metadata = totalData.textToStructuredArray(tempInputFile=args.strAbundanceFile, tempDelimiter=Constants.TAB, 
                                                            tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow),
                                                            tempNormalize=fNormalize)
        #Get sample names
        lsSampleNames = rawAbundance.dtype.names[1:]

        #Generate the metric measurements of samples for the y axis
        lStudyMetrics = None
        if c_strMetricCategory == mCC.c_UseDiversityMetrics:
            #List of lists, one list per diversity metric for all samples
            lStudyMetrics = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = rawAbundance, tempSampleNames = lsSampleNames, tempDiversityMetricAlpha = [strRunMetric])
            lStudyMetrics = lStudyMetrics[0]

        #Get top ranked of each sample
        elif c_strMetricCategory == mCC.c_UseTopRankedMethods:
            pass

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

        #Break up the grandSelection of subsampling by sample selection size
        #For each selection size
        #Get the average y for each sample subset (y={diversity, top ranked abundance})
        dictMetricsBySampleN = dict()
        iSelectionCount = -1
        fError = False

        logging.debug("dictGrandSelection")
        logging.debug(dictGrandSelection)

        #For each subsample (N) evaluate selection methods
        for dictStudy in dictGrandSelection:
            dictCurStudy = dictGrandSelection[dictStudy]

            logging.debug("dictCurStudy")
            logging.debug(dictCurStudy)

            iSelectionCount = -1
            #For each selection method make sure the other selection in the study selected the same N
            lsMethodStats = list()
            for strMethod in dictCurStudy:
                lsCurSampleSelections = dictCurStudy[strMethod]

                logging.debug("lsCurSampleSelections")
                logging.debug(lsCurSampleSelections)

                if iSelectionCount == -1:
                    iSelectionCount = len(lsCurSampleSelections)
                elif not iSelectionCount == len(lsCurSampleSelections):
                    print("".join(["MicroPitaPaperCollectionCurve. Selection methods selected an uneven number of samples. Did not include this study in the collection curve.",dictStudy,"."]))
                    logging.error("".join(["MicroPitaPaperCollectionCurve. Selection methods selected an uneven number of samples. Did not include this study in the collection curve.",dictStudy,"."]))
                    fError = True
                #If there was an error do not add the study's methods to plot data
                if fError:
                  lsMethodStats = list()
                  break
                else:
                    if c_strMetricCategory == mCC.c_UseDiversityMetrics:
                        #Calculate diversity
                        lsCurSampleDiversity = list()
                        for strSample in lsCurSampleSelections:
                            lsCurSampleDiversity.append(lStudyMetrics[lsSampleNames.index(strSample)])
                        dDiversity = sum(lsCurSampleDiversity)/float(len(lsCurSampleDiversity))
                        lsMethodStats.append([strMethod,[strRunMetric,dDiversity]])

                    elif c_strMetricCategory == mCC.c_UseTopRankedMethods:
                        #Caluclate top ranked
                        dTopRanked = -1
                        lsMethodStats.append([strMethod,[c_strMetricCategory,dTopRanked]])

            logging.debug("lsMethodStats")
            logging.debug(lsMethodStats)

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

        #Most diverse set of samples at each N
        dictdDiverseSamplesAtN = dict()
        #Booststrapped diversity level at N
        dictdBootstrappedDiversityAtN = dict()

        #Average multiple method metric instances at sample size N (if they exist)
        #N level
        logging.debug("dictMetricsBySampleN")
        logging.debug(dictMetricsBySampleN)

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
                        ldTopSamples = sorted(lStudyMetrics)[-1*int(strN):]
                        dictdDiverseSamplesAtN[strN] = sum(ldTopSamples)/float(len(ldTopSamples))

                    #Get bootstrapped diversity at a level N
                    if strN not in dictdBootstrappedDiversityAtN:
                        dictdBootstrappedDiversityAtN[strN] = mCC.getMedianBootstrappedMetric(npaAbundance=rawAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount = int(strN), sMetric=strRunMetric, iBootStrappingItr = c_BootstrapItr)

        #Update the plot name with the run metric
        strPlotNamePieces = filter(None,re.split(Constants.PATH_SEP,args.strOutFigure))
        strPlotName = Constants.PATH_SEP.join(strPlotNamePieces[0:-1])
        strPlotName = Constants.PATH_SEP.join([strPlotName,strRunMetric+"-"+strPlotNamePieces[-1:][0]])

        #Plot line graph
        mCC.funcPlotCollectionCurve(strPlotName=strPlotName, dictMethods=dictMetricsBySampleN, dictdMaxYMetric=dictdDiverseSamplesAtN, dictdDiversityBaseline=dictdBootstrappedDiversityAtN, strMetric=strRunMetric, fInvert = fInvert)

    logging.info("Stop MicropitaPaperCollectionCurve")

if __name__ == "__main__":
    _main( )

