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
from Constants_Arguments import Constants_Arguments
from Constants_Figures import Constants_Figures
from Diversity import Diversity
import logging
from pylab import *
from MicroPITA import MicroPITA
import numpy as np
import operator
import os
import random
import re

class MicropitaPaperCollectionCurve:

    #Bootstrap interations
    c_BootstrapItr = 100

    #Do not normalize if Choa1
    fNormalize = True

    #Metric
    c_strMetricCategory = MicroPITA.c_INV_SIMPSON_A_DIVERSITY

    #Allow switching to Choa1 for testing
    if True:
        fNormalize = False
        c_strMetricCategory = MicroPITA.c_CHAO1_A_DIVERSITY

    #Plot the collection curve
    def funcPlotCollectionCurve(self, strPlotName, dictMethods, dictdMaxYMetric, dictdMinYMetric, dictdDiversityBaseline, strMetric, fInvert):
        logging.info("Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve")
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strPlotName=",str(strPlotName)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictMethods=",str(dictMethods)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdMaxYMetric=",str(dictdMaxYMetric)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdMinYMetric=",str(dictdMinYMetric)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdDiversityBaseline=",str(dictdDiversityBaseline)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strMetric=",str(strMetric)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve fInvert=",str(fInvert)]))

        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strPlotName=",str(strPlotName)]))
        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictMethods=",str(dictMethods)]))
        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdMaxYMetric=",str(dictdMaxYMetric)]))
        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdMinYMetric=",str(dictdMinYMetric)]))
        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictdDiversityBaseline=",str(dictdDiversityBaseline)]))
        print("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strMetric=",str(strMetric)]))

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
        for iCount in sorted(dictdMaxYMetric.keys()):
            liX.append(iCount)
            liY.append(dictdMaxYMetric[iCount])
        plot(liX, liY, color=objColors.c_strDetailsColorLetter, marker="*", linestyle='--', label = "Average Max "+strMetric)

        #Update range
        iYMin = min(liY)
        iYMax = max(liY)
        iXMin = min(liX)
        iXMax = max(liX)

        #Plot min y metric
        liX = list()
        liY = list()
        for iCount in sorted(dictdMinYMetric.keys()):
            liX.append(iCount)
            liY.append(dictdMinYMetric[iCount])
        plot(liX, liY, color=objColors.c_strDetailsColorLetter, marker="-", linestyle='--', label = "Average Min "+strMetric)

        #Update range
        iYMin = min(liY)
        iYMax = max(liY)
        iXMin = min(liX)
        iXMax = max(liX)

        #Plot bootstrapped diversity
        liX = list()
        liY = list()
        for iCount in sorted(dictdDiversityBaseline.keys()):
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
            liX = list()
            liY = list()
            for iCount in sorted(dictCurMethod.keys()):
                liX.append(iCount)
                liY.append(dictCurMethod[iCount])
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
        logging.info("Stop MicropitaPaperCollectionCurve.funcPlotCollectionCurve")

    #First a subselection occurs and pool. 
    # Pool by sum then normalize then measure diversity
    def getMedianBootstrappedMetric(self, npaAbundance, lsSampleNames, iSelectSampleCount, iBootStrappingItr, microPITA):
        logging.info("Start MicropitaPaperCollectionCurve.getMedianBootstrappedMetric")
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedMetric lsSampleNames=",str(lsSampleNames)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedMetric iSelectSampleCount=",str(iSelectSampleCount)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedMetric iBootStrappingItr=",str(iBootStrappingItr)]))

        ldMeasurePerIteration = list()
        for iItr in xrange(iBootStrappingItr):
            #Select population
            lsSelectedSampleNames = random.sample(lsSampleNames,iSelectSampleCount)
            #When combining combine counts by summing
            ldPooledSample = np.array(self.funcSumRowsOfColumns(npaAbundance,lsSelectedSampleNames))
            if float(sum(ldPooledSample))==0.0:
                ldMeasurePerIteration.append(0.0)
            else:
                if self.fNormalize:
                    ldPooledSample = ldPooledSample/float(sum(ldPooledSample))
                ldMeasurePerIteration.append(microPITA.getAlphaMetric(tempAbundancies=ldPooledSample, tempMetric=self.c_strMetricCategory))
        logging.info("Stop MicropitaPaperCollectionCurve.getMedianBootstrappedMetric")
        return np.median(ldMeasurePerIteration)

    #Takes the column indices of a npArray and sums the rows into one column
    #Returns a list which is the row sums of the column
    def funcSumRowsOfColumns(self, npaAbundance, lsSampleNames):
        #Compress by data name
        npPooledSample = npaAbundance[lsSampleNames[0]]
        for strSampleName in lsSampleNames[1:]:
            #When combining, combine counts by summing
            npPooledSample = npPooledSample + npaAbundance[strSampleName]
        return list(npPooledSample)

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperCollectionCurve.py", 
    description = """Generates a collection curve of diversity and top ranked abudance Taxa/OTU given method and number of samples selected in studies.""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strSampleNameRowArgument, dest="iSampleNameRow", metavar= "SampleNameRow", default=0, help= Constants_Arguments.c_strSampleNameRowHelp)
argp.add_argument(Constants_Arguments.c_strFirstDataRow, dest="iFirstDataRow", metavar= "FirstDataRow", default=1, help= Constants_Arguments.c_strFirstDataRowHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strLoggingHelp)
#Select file
argp.add_argument( "strAbundanceFile", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument( "strOutFigure", metavar = "CollectionCurveOutputFile", help = Constants_Arguments.c_genericOutputFigureFileHelp)
argp.add_argument( "strOutText", metavar = "CollectionCurveOutputTextFile", help = Constants_Arguments.c_genericOutputDataFileHelp)
argp.add_argument( "strSelectionFiles", metavar = "InputSelectionFiles", nargs = "+", help = Constants_Arguments.c_strSelectionMethodsHelp)
#Selection parameter
argp.add_argument(Constants_Arguments.c_strPlotSelectedArgument, metavar = "Selection_Methods", nargs = "+", help = Constants_Arguments.c_strPlotSelectedHelp)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    print("args")
    print(args)

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start MicropitaPaperCollectionCurve")
    logging.info("MicropitaPaperCollectionCurve. The following arguments were passed.")
    logging.info(str(args))

    #Invert figure
    fInvert = (args.fInvert.lower() == "true")

    #Instance of plot collection curve script
    mCC = MicropitaPaperCollectionCurve()

    #Instance of microPITA; used to generate diversity matrices
    microPITA = MicroPITA()

    #Manage the cases where there are no or 1 selection file given
    if args.strSelectionFiles == None:
      logging.warning("MicropitaPaperCollectionCurve. No files were provided indicating sample selection so the collection curve was not made.")
      return
    if isinstance(args.strSelectionFiles, basestring):
      args.strSelectionFiles = [args.strSelectionFiles]

    #Get abundance table data
    totalData = AbundanceTable()
    rawAbundance,metadata = totalData.textToStructuredArray(tempInputFile=args.strAbundanceFile, tempDelimiter=Constants.TAB, 
                                                            tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow),
                                                            tempNormalize=False)
    #Get sample names
    lsSampleNames = rawAbundance.dtype.names[1:]

    #Calculate individual diversities
    #List of lists, one list per diversity metric for all samples
    #Sort the sample names by their diversity and store the names (Lowest diversity first)
    #This needs to be normalized for the diversity metric but only normalize a copy,
    #the abundances will be used for pooling later and so will need to stay unnormalized.
    lStudyMetrics = None
    if mCC.fNormalize:
        lStudyMetrics = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = totalData.normalizeColumns(tempStructuredArray = rawAbundance.copy(), tempColumns = list(lsSampleNames)),
                                                      tempSampleNames = lsSampleNames,
                                                      tempDiversityMetricAlpha = [mCC.c_strMetricCategory])
    else:
        lStudyMetrics = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = rawAbundance,
                                                      tempSampleNames = lsSampleNames,
                                                      tempDiversityMetricAlpha = [mCC.c_strMetricCategory])
    lStudyMetrics = zip(lStudyMetrics[0],lsSampleNames)
    lStudyMetrics.sort(key=operator.itemgetter(0))

    lStudyMetrics = [tplSample[1] for tplSample in lStudyMetrics]

    #Read through selection files and place contents in dictionary
    dictAllSelectionStudies = dict()
    for strFile in args.strSelectionFiles:

        #Put sample selection dictionaries into a larger dictionary by file name
        if strFile in dictAllSelectionStudies:
            logging.error("".join(["MicroPitaPaperCollectionCurve. File for study was already collected in the collection curve. The following file was ignored:",strFile,"."]))
        dictAllSelectionStudies[strFile] = MicroPITA.funcReadSelectionFileToDictionary(strFile)

    logging.debug("dictAllSelectionStudies")
    logging.debug(dictAllSelectionStudies)

    #Store the different sampling counts in all the studies
    setSampleCounts = set()

    #For each selection/study evaluate the methods by diversity or top ranked
    #Creating the following structure {"MethodName":{iSampleCount:metric}}
    dictMetricsBySampleN = dict()
    for dictStudy in dictAllSelectionStudies:
        dictCurStudy = dictAllSelectionStudies[dictStudy]
        logging.debug("dictCurStudy")
        logging.debug(dictCurStudy)

        for strCurrentMethod in dictCurStudy:
            #If the method parsed from the selection file is a method that is passed in as an argument and indicated as a method to plot
            if strCurrentMethod in args.pltSel:
                print("strCurrentMethod")
                print(strCurrentMethod)
                #If method is not in dictMetricsBySampleN add it
                if not strCurrentMethod in dictMetricsBySampleN:
                    dictMetricsBySampleN[strCurrentMethod] = dict()
                dictCurStudyMethod = dictMetricsBySampleN[strCurrentMethod]

                #Current selection of samples
                lsCurSampleSelections = dictCurStudy[strCurrentMethod]

                iSampleCount = len(lsCurSampleSelections)
                setSampleCounts.add(iSampleCount)

                #Calculate diversity
                #This assumes that a method is not ran multiple times in the same study at the same count
                #And if so that the same method at the same count will give the same results which is currently true
                ldSummedSubSet = np.array(mCC.funcSumRowsOfColumns(rawAbundance,lsCurSampleSelections))

                #Normalize
                if mCC.fNormalize:
                    ldSummedSubSet = ldSummedSubSet/float(sum(ldSummedSubSet))

                #Get diversity
                dictCurStudyMethod[iSampleCount]=microPITA.getAlphaMetric(tempAbundancies=ldSummedSubSet, tempMetric=mCC.c_strMetricCategory)

    logging.debug("dictMetricsBySampleN")
    logging.debug(dictMetricsBySampleN)

    #Most diverse set of samples at each N
    dictdMostMetricSamplesAtN = dict()
    #Least diverse set of samples at each N
    dictdLeastMetricSamplesAtN = dict()
    #Booststrapped diversity level at N
    dictdBootstrappedMetricAtN = dict()

    #Per sample selection count get the metric of the most, least, and bootstrapped samples given the metric
    for strSelectionCount in setSampleCounts:
        #Get Most and least diverse sample per N
        ldMostSamples = lStudyMetrics[-1*int(strSelectionCount):]
        ldLeastSamples = lStudyMetrics[0:int(strSelectionCount)]

        #Select most diverse samples, pool, and optionally normalize
        ldSummedSubSet = np.array(mCC.funcSumRowsOfColumns(rawAbundance, ldMostSamples))
        if mCC.fNormalize:
            ldSummedSubSet = ldSummedSubSet/float(sum(ldSummedSubSet))
        dictdMostMetricSamplesAtN[strSelectionCount] = microPITA.getAlphaMetric(tempAbundancies=ldSummedSubSet, tempMetric=mCC.c_strMetricCategory)

        #Select least diverse samples, pool, and optionally normalize
        ldSummedSubSet = np.array(mCC.funcSumRowsOfColumns(rawAbundance, ldLeastSamples))
        if mCC.fNormalize:
            ldSummedSubSet = ldSummedSubSet/float(sum(ldSummedSubSet))
        dictdLeastMetricSamplesAtN[strSelectionCount] = microPITA.getAlphaMetric(tempAbundancies=ldSummedSubSet, tempMetric=mCC.c_strMetricCategory)
        #Get bootstrapped metric
        dictdBootstrappedMetricAtN[strSelectionCount] = mCC.getMedianBootstrappedMetric(npaAbundance=rawAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount = int(strSelectionCount), iBootStrappingItr = mCC.c_BootstrapItr, microPITA = microPITA)

    #Update the plot name with the run metric
    strPlotNamePieces = filter(None,re.split(Constants.PATH_SEP,args.strOutFigure))
    strPlotName = Constants.PATH_SEP.join(strPlotNamePieces[0:-1])
    strPlotName = Constants.PATH_SEP.join([strPlotName,mCC.c_strMetricCategory+"-"+strPlotNamePieces[-1:][0]])

    #Plot line graph
    mCC.funcPlotCollectionCurve(strPlotName=strPlotName, dictMethods=dictMetricsBySampleN, dictdMaxYMetric=dictdMostMetricSamplesAtN, dictdMinYMetric=dictdLeastMetricSamplesAtN, dictdDiversityBaseline=dictdBootstrappedMetricAtN, strMetric=mCC.c_strMetricCategory, fInvert = fInvert)

    logging.info("Stop MicropitaPaperCollectionCurve")

if __name__ == "__main__":
    _main( )

