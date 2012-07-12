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
import itertools
import logging
from pylab import *
from MicroPITA import MicroPITA
import numpy as np
import operator
import os
import random
import re
from Utility_Math import Utility_Math

class MicropitaPaperCollectionCurve:

    #Measurement modes to use
    c_sRichness = "RICHNESS"
    c_sDiversity = "DIVERSITY"
    c_sEvenness = "EVENNESS"

    #Which measurement method to use
    c_sMeasurement = c_sRichness
    
    #Bootstrap interations
    c_BootstrapItr = 1000

    #Metric
    c_strMetricCategory = "Observed Counts"

    #Plot the collection curve
    def funcPlotCollectionCurve(self, abndData, strPlotName, dictMethods, setiSelectionCounts, strMetric, fInvert):
        logging.info("Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve")
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strPlotName=",str(strPlotName)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve dictMethods=",str(dictMethods)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve setiSelectionCounts=",str(setiSelectionCounts)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve strMetric=",str(strMetric)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.funcPlotCollectionCurve fInvert=",str(fInvert)]))

        font = {'family':'arial',
           'color':'k',
           'weight':'normal',
           'size':12}

        #Get plot colors
        objColors = Constants_Figures()
        objColors.invertColors(fInvert=fInvert)

        #Get plot object and set colors for inversion
        imgFigure = plt.figure()
        imgFigure.set_facecolor(objColors.c_strBackgroundColorWord)
        imgSubplot = imgFigure.add_subplot(111,axisbg=objColors.c_strBackgroundColorLetter)
        imgSubplot.set_xlabel('Number of samples sampled from study (n)')
        if self.c_sMeasurement == self.c_sRichness:
            imgSubplot.set_ylabel("".join(["Sampled Richness (",strMetric,")"]))
        elif self.c_sMeasurement == self.c_sEvenness:
            imgSubplot.set_ylabel("".join(["Sampled Evenness (",strMetric,")"]))
        elif self.c_sMeasurement == self.c_sDiversity:
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

        #Get sample names
        lsSampleNames = abndData.funcGetSampleNames()
        rawAbundance = abndData.funcGetAbundanceCopy()

        #Used to size the domain and range of the plots
        iYMin = 0
        iYMax = 0
        iXMin = 0
        iXMax = 0
        liX = list()
        liY = list()

        #Plot bootstrapped metric median values as a line
        liX = list()
        liY = list()

        #Booststrapped diversity level at N
        dictdBootstrappedMetricAtN = dict()

        #Size of boxplots
        dBoxPlotwidth = max([.5*((max(setiSelectionCounts)-min(setiSelectionCounts))/20.0),.5])

        #Per sample selection count get the bootstrapped metric
        for iSelectionCount in setiSelectionCounts:
            #Get bootstrapped metric
            ldMeasurements = self.getMedianBootstrappedObservedCount(npaAbundance=rawAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount = iSelectionCount, iBootStrappingItr = self.c_BootstrapItr)

            #Plot box plot
            bp = imgSubplot.boxplot(x=[ldMeasurements], positions=[iSelectionCount], notch=1, sym="", widths=dBoxPlotwidth, patch_artist=True)
            plt.setp(bp['boxes'], color=objColors.c_strDetailsColorLetter, facecolor=objColors.c_strGridLineColor, alpha=objColors.c_dAlpha)
            plt.setp(bp['whiskers'], color=objColors.c_strDetailsColorLetter)

            liX.append(iSelectionCount)
            liY.append(np.median(ldMeasurements))

            #Update range and domain
            iYMin = min([iYMin]+ldMeasurements)
            iYMax = max([iYMax]+ldMeasurements)
            iXMin = min([iXMin]+[iSelectionCount])
            iXMax = max([iXMax]+[iSelectionCount])

        #Plot median values as a line
        plot(liX, liY, color=objColors.c_strDetailsColorLetter, marker="*", linestyle='--', label = "Permuted "+strMetric)

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

        ###Aesthetics
        #Make the plot a little bigger than needed
        xlim(iXMin*.9,iXMax*1.1)
        ylim(iYMin*.9,iYMax*1.1)
        plt.xticks(setiSelectionCounts)
        objLegend = imgSubplot.legend(loc="lower right", scatterpoints=1, prop={'size':10})
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
    def getMedianBootstrappedObservedCount(self, npaAbundance, lsSampleNames, iSelectSampleCount, iBootStrappingItr):
        logging.info("Start MicropitaPaperCollectionCurve.getMedianBootstrappedObservedCount")
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedObservedCount lsSampleNames=",str(lsSampleNames)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedObservedCount iSelectSampleCount=",str(iSelectSampleCount)]))
        logging.debug("".join(["Start MicropitaPaperCollectionCurve.getMedianBootstrappedObservedCount iBootStrappingItr=",str(iBootStrappingItr)]))

        ldMeasurePerIteration = list()
        for iItr in xrange(iBootStrappingItr):
            #Select population
            lsSelectedSampleNames = random.sample(lsSampleNames,iSelectSampleCount)
#            lsSelectedSampleNames = Utility_Math.funcSampleWithReplacement(lsSampleNames,iSelectSampleCount)
            #When combining combine counts by summing
            ldPooledSample = np.array(Utility_Math.funcSumRowsOfColumns(npaAbundance,lsSelectedSampleNames))
            if float(sum(ldPooledSample))==0.0:
                ldMeasurePerIteration.append(0.0)
            else:
                if self.c_sMeasurement == self.c_sRichness:
                    ldMeasurePerIteration.append(Diversity.funcGetObservedCount(ldSampleAbundances=ldPooledSample))
                elif self.c_sMeasurement == self.c_sEvenness:
                    ldMeasurePerIteration.append(Diversity.funcGetPielouEvenness(ldSampleTaxaAbundancies=ldPooledSample))
                elif self.c_sMeasurement == self.c_sDiversity:
                    ldMeasurePerIteration.append(Diversity.funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldPooledSample))
        logging.info("Stop MicropitaPaperCollectionCurve.getMedianBootstrappedObservedCount")
        return ldMeasurePerIteration

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperCollectionCurve.py", 
    description = """Generates a collection curve of diversity and top ranked abudance Taxa/OTU given method and number of samples selected in studies.""" )

#Arguments
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strLoggingHelp)
argp.add_argument(Constants_Arguments.c_strEcologicalMeasurementArgument, dest = "sEcologicalMeasurement", action = "store", metavar="ecologicalMeasurementToMeasureSamples", help = Constants_Arguments.c_strEcologicalMeasurementHelp)

#Select file
argp.add_argument( "strAbundanceFile", metavar = "Abundance_file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument( "strOutFigure", metavar = "CollectionCurveOutputFile", help = Constants_Arguments.c_strGenericOutputFigureFileHelp)
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

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")

    #Instance of plot collection curve script
    mCC = MicropitaPaperCollectionCurve()
    mCC.c_sMeasurement = args.sEcologicalMeasurement
    if mCC.c_sMeasurement == mCC.c_sRichness:
        mCC.c_strMetricCategory = "Observed Occurence"
    elif mCC.c_sMeasurement == mCC.c_sEvenness:
        mCC.c_strMetricCategory = "Pielou"
    elif mCC.c_sMeasurement == mCC.c_sDiversity:
        mCC.c_strMetricCategory = "Inverse Simson"
    
    #Instance of microPITA; used to generate diversity matrices
    microPITA = MicroPITA()

    #Manage the cases where there are no or 1 selection file given
    if args.strSelectionFiles == None:
      logging.error("MicropitaPaperCollectionCurve. No files were provided indicating sample selection so the collection curve was not made.")
      return False
    if isinstance(args.strSelectionFiles, basestring):
      args.strSelectionFiles = [args.strSelectionFiles]

    #Read abundance file
    #Abundance table object to read in and manage data
    totalData = AbundanceTable.funcMakeFromFile(strInputFile=args.strAbundanceFile, fIsNormalized=fIsNormalized,
                                            fIsSummed=fIsSummed, sMetadataID=args.sIDName, sLastMetadata=args.sLastMetadataName)

    ##Certain metrics need different data states, check for them.
    #Do not produce a plot for summed data for any metric.
    if totalData.funcIsSummed():
        logging.error("MicropitaPaperCollectionCurve. Will not produce a refraction curve on summed data.")
        return False
    #Do not produce a plot for richness for normalized data
    if mCC.c_sMeasurement == mCC.c_sRichness:
        if totalData.funcIsNormalized():
            logging.error("MicropitaPaperCollectionCurve. Will not produce a refraction curve on normalized.")
            return False
    #If diversity or evennness is measured, normalize
    if (mCC.c_sMeasurement == mCC.c_sDiversity) or (mCC.c_sMeasurement == mCC.c_sEvenness):
        totalData.funcNormalize()

    rawAbundance = totalData.funcGetAbundanceCopy()

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
    setiSampleCounts = set()

    #For each selection/study evaluate the methods by diversity or top ranked
    #Creating the following structure {"MethodName":{iSampleCount:metric measurement}}
    # example {"Representative":{5:34}}
    dictMetricsBySampleN = dict()
    for dictStudy in dictAllSelectionStudies:
        dictCurStudy = dictAllSelectionStudies[dictStudy]
        logging.debug("dictCurStudy")
        logging.debug(dictCurStudy)

        for strCurrentMethod in dictCurStudy:
            #If the method parsed from the selection file is a method that is passed in as an argument and indicated as a method to plot
            if strCurrentMethod in args.pltSel:
                #If method is not in dictMetricsBySampleN add it
                if not strCurrentMethod in dictMetricsBySampleN:
                    dictMetricsBySampleN[strCurrentMethod] = dict()
                dictCurStudyMethod = dictMetricsBySampleN[strCurrentMethod]

                #Current selection of samples
                lsCurSampleSelections = dictCurStudy[strCurrentMethod]

                iSampleCount = len(lsCurSampleSelections)
                if iSampleCount > 0:
                    setiSampleCounts.add(iSampleCount)

                    #Calculate measurement
                    #This assumes that a method is not ran multiple times in the same study at the same count
                    #And if so that the same method at the same count will give the same results which is currently true
                    ldSummedSubSet = np.array(Utility_Math.funcSumRowsOfColumns(rawAbundance,lsCurSampleSelections))

                    #Get measurement
                    if mCC.c_sMeasurement == mCC.c_sRichness:
                        dictCurStudyMethod[iSampleCount]=Diversity.funcGetObservedCount(ldSampleAbundances=ldSummedSubSet)
                    elif mCC.c_sMeasurement == mCC.c_sEvenness:
                        dictCurStudyMethod[iSampleCount]=Diversity.funcGetPielouEvenness(ldSampleTaxaAbundancies=ldSummedSubSet)
                    elif mCC.c_sMeasurement == mCC.c_sDiversity:
                        dictCurStudyMethod[iSampleCount]=Diversity.funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldSummedSubSet)

    logging.debug("dictMetricsBySampleN")
    logging.debug(dictMetricsBySampleN)

    #Update the plot name with the run metric
    strPlotNamePieces = filter(None,re.split(Constants.PATH_SEP,args.strOutFigure))
    strPlotName = Constants.PATH_SEP.join(strPlotNamePieces[0:-1])
    strPlotName = Constants.PATH_SEP.join([strPlotName,mCC.c_sMeasurement+"-"+strPlotNamePieces[-1:][0]])

    #Plot line graph
    sortedCounts = list(setiSampleCounts)
    sortedCounts.sort()
    mCC.funcPlotCollectionCurve(abndData = totalData, strPlotName=strPlotName, dictMethods=dictMetricsBySampleN, setiSelectionCounts = sortedCounts, strMetric=mCC.c_strMetricCategory, fInvert = fInvert)

    logging.info("Stop MicropitaPaperCollectionCurve")

if __name__ == "__main__":
    _main( )

