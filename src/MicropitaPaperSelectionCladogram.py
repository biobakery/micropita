#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
Generates a cladogram of the selected samples.
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
from AbundanceTable import AbundanceTable
from Cladogram import Cladogram
from Constants import Constants
from Constants_Arguments import Constants_Arguments
from Constants_Figures import Constants_Figures
import logging
from MicroPITA import MicroPITA
import numpy as np
import os
import scipy.stats.stats as stats
from Utility_Math import Utility_Math

class MicropitaPaperSelectionCladogram:
    """
    Generate Cladogram
    """

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicropitaPaperSelectionCladogram.py", 
    description = """Generates a cladogram of the selected samples.""" )

#Arguments
#Select file
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, 
                  help= Constants_Arguments.c_strLoggingHelp)
argp.add_argument( "strSelectionFile", metavar = "Select_file", help = Constants_Arguments.c_strMicropitaSelectFileHelp )
argp.add_argument( "strInputFile", metavar = "Input_Abundance_File", help = Constants_Arguments.c_strAbundanceFileHelp )
argp.add_argument( "strStyleFile", metavar = "Style_file", help = Constants_Arguments.c_strCircladerStyleFileHelp )
argp.add_argument(Constants_Arguments.c_strTargetedSelectionFileArgument, dest="strTargetedTaxaFile", metavar= "TaxaFilePath", default=None, help= Constants_Arguments.c_strTargetedSelectionFileHelp)
argp.add_argument(Constants_Arguments.c_strHighlightCladeFileArgument, dest="iHighlightCladeFile", metavar= "HighlightFilePath", default=None, help= Constants_Arguments.c_strHighlightCladeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
argp.add_argument(Constants_Arguments.c_strRootArgument, dest = "sRoot", metavar= "root", default="None", help = Constants_Arguments.c_strRootHelp )
argp.add_argument(Constants_Arguments.c_strEnrichmentMethodArgument, dest="strEnrichmentIndicatorMethod", metavar= "enrichmentIndicatorMethod", default="FDR", 
                  choices=Constants_Arguments.c_strEnrichmentChoices, help= Constants_Arguments.c_strEnrichmentMethodHelp)
argp.add_argument(Constants_Arguments.c_strCladeFilterLevelArgument, dest = "iCladeFilterLevel", metavar= "CladeFilterLevel", default="None",help = Constants_Arguments.c_strCladeFilterLevelHelp)
argp.add_argument(Constants_Arguments.c_strCladeMeasureLevelArgument, dest = "iCladeFilterMeasure", metavar= "CladeFilterLevelMeasure", default="None", help = Constants_Arguments.c_strCladeMeasureLevelHelp)
argp.add_argument(Constants_Arguments.c_strCladeFilteringMinLevelArgument, dest = "iCladeFilterMinNumber", metavar= "CladeFilterMinNumber", default="None",
	help = Constants_Arguments.c_strCladeFilteringMinLevelHelp)
argp.add_argument(Constants_Arguments.c_strAbundanceFilterPercentileArgument, dest = "iAbundanceFilterPercentile", metavar= "AbundanceFilterPercentile", default=0.0,
	help = Constants_Arguments.c_strAbundanceFilterPercentileHelp)
argp.add_argument(Constants_Arguments.c_strAbundanceFilterCutoffArgument, dest = "iAbundanceFilterPercentCuttoff", metavar= "AbundanceFilterPercentCuttoff", default=0.0,
	help = Constants_Arguments.c_strAbundanceFilterCutoffHelp)
argp.add_argument(Constants_Arguments.c_strRingOrderArgument, dest = "iRingOrder", metavar= "Ring Order", default=None, help = Constants_Arguments.c_strRingOrderHelp)
argp.add_argument(Constants_Arguments.c_strCircladerTicksArgument, dest = "iTicks", metavar= "Internal Dendrogram Ticks", default=None, help = Constants_Arguments.c_strCircladerTicksHelp)
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "SampleRowName", default=None, help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "FirstDataRow", default=None, help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False", help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strEnrichmentThresholdArgument, dest = "dAlpha", action = "store", default = 0.05, help = Constants_Arguments.c_strEnrichmentThresholdHelp)
argp.add_argument(Constants_Arguments.c_strOccurenceFilterSequenceCountArgument, dest ="iMinSequenceCount", action = "store", default=0.0, help = Constants_Arguments.c_strOccurenceFilterSequenceHelp)
argp.add_argument(Constants_Arguments.c_strOccurenceFilterSampleCountArgument, dest ="iMinSampleCount", action = "store", default=0.0, help = Constants_Arguments.c_strOccurenceFilterSampleHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store", metavar= "flagIndicatingNormalization", 
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store", metavar= "flagIndicatingSummation", help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strSumDataArgument, dest="fSumData", action = "store", metavar= "WouldlikeDataSummed", help= Constants_Arguments.c_strSumDataHelp)
argp.add_argument(Constants_Arguments.c_strTerminalLevelArgument, dest = "iTerminalCladeLevel", metavar= "Terminal Clade Level", default=10, type=int, help = Constants_Arguments.c_strTerminalLevelHelp)

#Outputfile
argp.add_argument( "sTaxaFileName", metavar = "TaxaFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerTaxaFileHelp )
argp.add_argument( "sColorFileName", metavar = "ColorFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerColorFileHelp )
argp.add_argument( "sTickFileName", metavar = "TickFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerTickFileHelp )
argp.add_argument( "sHighlightFileName", metavar = "HighlightFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerHighlightFileHelp )
argp.add_argument( "sSizeFileName", metavar = "SizeFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerSizeFileHelp )
argp.add_argument( "sCircleFileName", metavar = "CircleFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerCircleFileHelp )
argp.add_argument( "strOutFigure", metavar = "SelectionCladogram.png", nargs = "?", help = Constants_Arguments.c_strCircladerOutputFigureHelp )
argp.add_argument( "strDetailOutputFile", metavar = "SelectionDetail.png", nargs = "?", help = Constants_Arguments.c_strCircladerOutputDetailsHelp )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start Circlader")
    logging.info(args)

    #Invert
    c_fInvert = (args.fInvert.lower() == "true")
    #Normalize
    c_Normalize = (args.fNormalize.lower() == "true")

    #Is summed and normalized
    fIsSummed = (args.fIsSummed.lower() == "true")
    fIsNormalized = (args.fIsNormalized.lower() == "true")
    fSumData = (args.fSumData.lower() == "true")

    #Position holders in t-stats data list
    c_IDINDEX = 0
    c_TSCOREINDEX = 1
    c_PVALUEINDEX = 2
    c_QVALUEINDEX = 3

    #Enrichment values
    c_IncreasedEnrichment = "^"
    c_DecreasedEnrichment = "v"
    c_NoChangeEnrichment = "R"

    #Get colors
    objColors = Constants_Figures()
    objColors.invertColors(c_fInvert)

    #Root tree
    c_fBacteriaRooted = (not args.sRoot.lower() == "none")

    #Delimiter for highlight and label information in the highlight clade files
    c_strLabelDelim = "="

    #Delimiter for taxa/otu lineages
    c_strLineageDelim = "|"

    #Contains details on the process that will be output
    sCladogramDetails = "".join(["Input File:",args.strInputFile,Constants.ENDLINE])

    #Get Abundance table data
    rawData = AbundanceTable.funcMakeFromFile(strInputFile=args.strInputFile, fIsNormalized=fIsNormalized,
                                          fIsSummed=fIsSummed, sMetadataID=args.sIDName, 
                                          sLastMetadata=args.sLastMetadataName,cFeatureNameDelimiter=c_strLineageDelim)

    #Update detail: original feature count
    sCladogramDetails = "".join([sCladogramDetails,"Original Feature Count=",str(rawData.funcGetFeatureCount()),Constants.ENDLINE,Constants.ENDLINE])

    #Sum clades before normalization and filtering
    if fSumData:
        rawData.funcSumClades()

    #Filtering
    #Indicate Filtering
    fFilterAbundance = (not args.iAbundanceFilterPercentile.lower() == "none")
    fFilterOccurence = (not args.iMinSequenceCount.lower() == "none")

    if fFilterOccurence:
      if rawData.funcIsNormalized():
          logging.error("MicropitaPaperSelectionCladogram::Can not filter on occurence on a normalized file.")
          return False 
      logging.debug("MicropitaPaperSelectionCladogram::Before occurence filtering the feature count is "+str(rawData.funcGetFeatureCount()))
      rawData.funcFilterAbundanceBySequenceOccurence(iMinSequence = int(args.iMinSequenceCount), iMinSamples = int(args.iMinSampleCount))
      logging.debug("MicropitaPaperSelectionCladogram::After occurence filtering the feature count is "+str(rawData.funcGetFeatureCount()))
      #Update detail: Occurence filter details
      sCladogramDetails = "".join([sCladogramDetails,"Occurence Filtering:",Constants.ENDLINE,"Minimum Sequence count=",args.iMinSequenceCount,Constants.ENDLINE,
          "Minimum sample count=",args.iMinSampleCount,Constants.ENDLINE,"Filtered Feature Count (After Occurence Filtering)=",str(rawData.funcGetFeatureCount()),Constants.ENDLINE,Constants.ENDLINE])

    if fFilterAbundance:
      logging.debug("MicropitaPaperSelectionCladogram::Before abundance filtering the feature count is "+str(rawData.funcGetFeatureCount()))
      rawData.funcFilterAbundanceByPercentile(dPercentileCutOff = float(args.iAbundanceFilterPercentile),
                                          dPercentageAbovePercentile = float(args.iAbundanceFilterPercentCuttoff))
      logging.debug("MicropitaPaperSelectionCladogram::After abundance filtering the feature count is "+str(rawData.funcGetFeatureCount()))
      #Update detail: abundance feature details
      sCladogramDetails = "".join([sCladogramDetails,"Abundance Filtering:",Constants.ENDLINE,"Percentile Abundance Threshold=",args.iAbundanceFilterPercentile,
          Constants.ENDLINE,"Percent Samples Above Percentile=",args.iAbundanceFilterPercentCuttoff,Constants.ENDLINE,
          "Filtered Feature Count (After Abundance Filtering)=",str(rawData.funcGetFeatureCount()),Constants.ENDLINE,Constants.ENDLINE])

    #Reduce the clades down to a certain clade level
    sCladogramDetails = "".join([sCladogramDetails," Feature count before reducing to a clade level:",str(rawData.funcGetFeatureCount()),Constants.ENDLINE])
    rawData.funcReduceFeaturesToCladeLevel(args.iTerminalCladeLevel)
    sCladogramDetails = "".join([sCladogramDetails," Reducing clades to the following level:",str(args.iTerminalCladeLevel),Constants.ENDLINE])
    sCladogramDetails = "".join([sCladogramDetails," Feature count AFTER reducing to a clade level:",str(rawData.funcGetFeatureCount()),Constants.ENDLINE])
    sCladogramDetails = "".join([sCladogramDetails,Constants.ENDLINE.join(rawData.funcGetFeatureNames())])+Constants.ENDLINE

    if c_Normalize:
      rawData.funcNormalize()

    #####
    ## Note all manipulations to the abundance table including feature filtering should occur before this point
    #####

    #Abundance data
    abundance = rawData.funcGetAbundanceCopy()
    #All sample names
    lsAllSampleNames = rawData.funcGetSampleNames()

    #All taxa in the study
    lsAllTaxa = rawData.funcGetFeatureNames()
    sCladogramDetails = "".join([sCladogramDetails," Total Taxa:",str(len(lsAllTaxa)),Constants.ENDLINE,Constants.ENDLINE.join(lsAllTaxa),Constants.ENDLINE])

    #Create a cladogram object
    cladogram = Cladogram()

    #Set the Abundance data
    cladogram.setAbundanceData(rawData)

    #Terminal taxa in the study
    lsTerminalTaxa = rawData.funcGetTerminalNodes()

    #Genus Level of the ancestry to filter on
    fFilterClades = (not args.iCladeFilterLevel.lower() == "none")
    iCladeLevelForFiltering = 0
    iCladeLevelForReducing = 0
    iCladeLevelMinimumCount = 0

    if(fFilterClades):
      iCladeLevelForFiltering = int(args.iCladeFilterMeasure)
      iCladeLevelForReducing = int(args.iCladeFilterLevel)
      iCladeLevelMinimumCount = int(args.iCladeFilterMinNumber)
    cladogram.setFilterByCladeSize(fCladeSizeFilter = fFilterClades, iCladeLevelToMeasure = iCladeLevelForFiltering, iCladeLevelToReduce = iCladeLevelForReducing, iMinimumCladeSize = iCladeLevelMinimumCount, cFeatureDelimiter=c_strLineageDelim)

    #Update detail: Clade filtering
    sCladogramDetails = "".join([sCladogramDetails,"Clade filtering:",Constants.ENDLINE,"Clade level for measuring=",str(iCladeLevelForFiltering),Constants.ENDLINE,
                                                "Clade level for filtering=",str(iCladeLevelForReducing),Constants.ENDLINE,"Minimum features needed in the measured clade=",
                                                str(iCladeLevelMinimumCount),Constants.ENDLINE,Constants.ENDLINE])

    #Dictionary to hold selection data
    dictSelection = MicroPITA.funcReadSelectionFileToDictionary(args.strSelectionFile)

    #Update detail: Sampling selection
    sCladogramDetails = "".join([sCladogramDetails,"Sample Selection:",Constants.ENDLINE,
                        Constants.ENDLINE.join(["".join([sKey,"=",str([dictSelection[sKey]])]) for sKey in dictSelection.keys()]),Constants.ENDLINE,Constants.ENDLINE])

    #Set Color Data
    dictColors = {MicroPITA.c_strDiversity1:objColors.invSimpsonColorN,
                  MicroPITA.c_strDiversity2:objColors.chao1ColorN,
                  MicroPITA.c_strRepresentativeDissimilarity1:objColors.brayCurtisColorN,
                  MicroPITA.c_strExtremeDissimiarity1:objColors.inBrayCurtisColorN,
                  MicroPITA.c_strUserRanked:objColors.userRankedN,
                  MicroPITA.c_strRandom:objColors.randomColorN,
                  MicroPITA.c_strSVMClose:objColors.svmCloseN,
                  MicroPITA.c_strSVMFar:objColors.svmFarN,
                  objColors.c_strBackgroundColorName:objColors.c_strBackgroundColor}
    cladogram.setColorData(dictColors)

    #Set tick data
    llsTicks = None
    if not args.iTicks == None:
      llsTicks = []
      #Parse string of ticks
      sFilteredTick = filter(None,[sTick.strip() for sTick in args.iTicks.split(Constants.COMMA)])

      #Make parsed data into the format [["#","tick"],...]
      iTickCount = 0
      for sTick in sFilteredTick:
        llsTicks.append([str(iTickCount),sTick])
        iTickCount = iTickCount + 1
    #Set ticks
    cladogram.setTicks(llsTicks)

    #Set Root
    if(c_fBacteriaRooted):
      cladogram.forceRoot(args.sRoot)

    #Add clade highlighting
    #Force the highlighting of taxa defined inputs (with highlight)
    lsUserDefinedTaxa = list()
    if(not args.strTargetedTaxaFile == "None"):
      fhndlTaxaInput = open(args.strTargetedTaxaFile,'r')
      lsUserDefinedTaxa = filter(None,fhndlTaxaInput.read().split(Constants.ENDLINE))

    #Update detail: Targeted taxa
    sCladogramDetails = "".join([sCladogramDetails,Constants.ENDLINE,"Targeted Taxa:",Constants.ENDLINE,"".join([sTaxa for sTaxa in lsUserDefinedTaxa]),Constants.ENDLINE,Constants.ENDLINE])

    lsUserDefinedHighlighted = list()
    if(not args.iHighlightCladeFile == "None"):
      fhndlHighlightInput = open(args.iHighlightCladeFile,'r')
      lsUserDefinedHighlighted = filter(None,fhndlHighlightInput.read().split(Constants.ENDLINE))

      #Highlight and relabel data
      dictTaxaHighlights = dict()
      dictRelabel = dict()
      iLabelCount = 1
      for sHighlight in lsUserDefinedHighlighted:
        if c_strLabelDelim in sHighlight:
          lsHighlightElements = filter(None,sHighlight.split(c_strLabelDelim))
          dictTaxaHighlights[lsHighlightElements[0]] = MicroPITA.c_strDiversity2
          dictRelabel[lsHighlightElements[0]] = "".join([str(iLabelCount),":",lsHighlightElements[1]])
          iLabelCount = iLabelCount + 1
        else:
          dictTaxaHighlights[sHighlight] = MicroPITA.c_strDiversity2
      cladogram.addHighLights(dictTaxaHighlights,True)
      if len(dictRelabel) > 0:
        cladogram.relabelIDs(dictRelabel)

    #Allows one to set the order of the rings if needed
    lsSelectedSampleMethod = dictSelection.keys()
    if not args.iRingOrder == None:
      lsSelectedSampleMethod = [filter(None,strMethod) for strMethod in (args.iRingOrder.split(Constants.COMMA))]

    #Update detail: Enrichment details
    sCladogramDetails = "".join([sCladogramDetails,"Enrichment threshold method:",str(args.strEnrichmentIndicatorMethod),Constants.ENDLINE,
                                              "Enrichment threshold:",str(args.dAlpha),Constants.ENDLINE,Constants.ENDLINE])

    #Measure enrichment
    if((args.strEnrichmentIndicatorMethod=="PVALUE") or (args.strEnrichmentIndicatorMethod=="FDR")):#P-value or qvalue
      #Are we using qvalues or pvalues
      fIsQValue = (args.strEnrichmentIndicatorMethod=="FDR")

      #This is performing t-test with p-values
      for selectedSampleMethod in lsSelectedSampleMethod:
        if selectedSampleMethod in dictSelection:
          sCladogramDetails = "".join([sCladogramDetails,"Selection Method:",selectedSampleMethod])
          #All samples that were selected by the method
          lsSelectedSamples = dictSelection[selectedSampleMethod]
          if len(lsSelectedSamples) > 0:
            #Samples selected
            lfSelectedSamplesUTest = []
            lfNotSelectedSamplesUTest = []
            #Set up boolean list to compress array
            for sample in lsAllSampleNames:
              wasSelected = sample in lsSelectedSamples
              lfSelectedSamplesUTest.append(wasSelected)
              lfNotSelectedSamplesUTest.append(not wasSelected)
            sCladogramDetails = "".join([sCladogramDetails,"\nSamples:",str(lsAllSampleNames)])
            sCladogramDetails = "".join([sCladogramDetails,"\nSelected:",str(lsSelectedSamples)])
            sCladogramDetails = "".join([sCladogramDetails,"\nSelected samples in total list: ",str(lfSelectedSamplesUTest)])
            sCladogramDetails = "".join([sCladogramDetails,"\nSamples not selected in total list: ",str(lfNotSelectedSamplesUTest)])

            #Holds t-tests,pvalues,qvalues as needed
            lsTaxaTScores = list()

            #Compress arrays to one or the other distribution
            #Conduct wilcoxon tests on all taxa
            for iTaxonIndex, strTaxaId in enumerate(lsAllTaxa):
              npaDistribution = np.array(list(abundance[iTaxonIndex,])[1:])
              #Hold info about the taxa to store in the cladogram info file
              sTaxaData = "".join([Constants.ENDLINE,strTaxaId])

              if(strTaxaId in lsTerminalTaxa):
                npaSelectedDistribution = np.compress(lfSelectedSamplesUTest,npaDistribution)
                npaNotSelectedDistribution = np.compress(lfNotSelectedSamplesUTest,npaDistribution)
                sCladogramDetails = "".join([sCladogramDetails,"\nTerminal Feature: ",strTaxaId])
                sCladogramDetails = "".join([sCladogramDetails,"\nnpaSelectedDistribution: ",str(npaSelectedDistribution)])
                sCladogramDetails = "".join([sCladogramDetails,"\nnpaNotSelectedDistribution: ",str(npaNotSelectedDistribution)])
                dScore, dPvalue = stats.ranksums(npaSelectedDistribution,npaNotSelectedDistribution)
                sTaxaData = " ".join(["Score",str(dScore),"P-value",str(dPvalue),Constants.ENDLINE])
                if(sum(npaSelectedDistribution)==0):
                  sTaxaData = "".join([sTaxaData,"Selected Average: 0"])
                else:
                  sTaxaData = "".join([sTaxaData,"Selected Average: ",str(sum(npaSelectedDistribution)/float(len(npaSelectedDistribution)))])
                if(sum(npaNotSelectedDistribution)==0):
                  sTaxaData = "".join([sTaxaData,", Not Selected Average: 0\n"])
                else:
                  sTaxaData = "".join([sTaxaData,", Not Selected Average: ",str(sum(npaNotSelectedDistribution)/float(len(npaNotSelectedDistribution))),"\n"])

                #[ID,TScore,PValue,QValue,SortOrder]
                lsTaxaTScores.append([strTaxaId,dScore,dPvalue,-1])
              else:
                sTaxaData = " ".join([sTaxaData,"Not terminal, not measured"])
              sCladogramDetails = "".join([sCladogramDetails,sTaxaData])

            #Get a list of pvalues preserving order
            #Holds terminal feature scores only at this point
            ldOrderedPValues = [iScoreDataIndex[c_PVALUEINDEX] for iScoreDataIndex in lsTaxaTScores]

            #If using qvalues generate them with FDR BH
            if fIsQValue:
              #Convert pvalues to qvalues
              ldOrderedQValues = Utility_Math.funcConvertToBHQValue(ldOrderedPValues)

              #Update the score data with qvalues
              for iQIndex in xrange(0,len(ldOrderedQValues)):
                lsTaxaTScores[iQIndex][c_QVALUEINDEX] = ldOrderedQValues[iQIndex]

            sCladogramDetails = "".join([sCladogramDetails,"\nAll Scores (Id, Score, pvalue) : ",str(lsTaxaTScores)])

            lsAlpha = list()
            lsShapes = list()
            lsTaxa = list()

            for lsCur in lsTaxaTScores:
              lsTaxa.append(lsCur[c_IDINDEX])
              dCurScore = lsCur[c_TSCOREINDEX]
              dValue = -1
              if fIsQValue:
                dValue = lsCur[c_QVALUEINDEX]
              else:
                dValue = lsCur[c_PVALUEINDEX]

              if(dValue <= float(args.dAlpha)):
                if dCurScore > 0.0:
                  lsAlpha.append(str(1-dValue))
                  lsShapes.append(c_IncreasedEnrichment)
                elif(dCurScore < 0.0):
                  lsAlpha.append(str(1-dValue))
                  lsShapes.append(c_DecreasedEnrichment)
                elif(dCurScore == 0.0):
                  logging.error("MicropitaPaperSelectionCladogram::Recieved a significant feature with no change in t-score...?")
                  return
              else:
                lsAlpha.append("0.0")
                lsShapes.append(c_NoChangeEnrichment)

            sCladogramDetails = "".join([sCladogramDetails,"\nlsAlpha : ",str(lsAlpha),"\n"])
            sCladogramDetails = "".join([sCladogramDetails,"\nlsShapes : ",str(lsShapes),"\n"])
            sCladogramDetails = "".join([sCladogramDetails,"\nlsTaxa : ",str(lsTaxa),"\n"])

            #Add circle for this data
            cladogram.addCircle(lsTaxa=lsTaxa, strShape=lsShapes, dAlpha=lsAlpha, strCircle=selectedSampleMethod, fForced=True)
            sCladogramDetails = "".join([sCladogramDetails,Constants.ENDLINE,Constants.ENDLINE,"Circle Data:",Constants.ENDLINE,Constants.ENDLINE.join([Constants.COMMA.join([str(lsTaxa[iIndex[0]]),str(lsShapes[iIndex[0]]), str(lsAlpha[iIndex[0]])]) for iIndex in enumerate(lsTaxa)])])

            #Update detail: Enrichment details
            #Build method enrichment detail
            iIncreasedEnrichmentCounts = sum([1 if cEnrichment == c_IncreasedEnrichment else 0 for cEnrichment in lsShapes])
            iDecreasedEnrichmentCounts = sum([1 if cEnrichment == c_DecreasedEnrichment else 0 for cEnrichment in lsShapes])
            iNoChangeEnrichmentCounts = sum([1 if cEnrichment == c_NoChangeEnrichment else 0 for cEnrichment in lsShapes])
            sMethodEnrichment = ""
            if iIncreasedEnrichmentCounts:
              sMethodEnrichment = "".join([sMethodEnrichment,"Terminal Taxa with Increased Enrichment= ",str(iIncreasedEnrichmentCounts),Constants.ENDLINE])
            if iDecreasedEnrichmentCounts:
              sMethodEnrichment = "".join([sMethodEnrichment,"Terminal Taxa with Decreased Enrichment= ",str(iDecreasedEnrichmentCounts),Constants.ENDLINE])
            if iNoChangeEnrichmentCounts:
              sMethodEnrichment = "".join([sMethodEnrichment,"Terminal Taxa with no change= ",str(iNoChangeEnrichmentCounts),Constants.ENDLINE])
            #Add method enrichment to total enrichment details
            if sMethodEnrichment:
              sCladogramDetails = "".join([sCladogramDetails, "\n", selectedSampleMethod,": ",Constants.ENDLINE,sMethodEnrichment,Constants.ENDLINE,Constants.ENDLINE])

    #Generate cladogram PDF
    cladogram.generate(strImageName=args.strOutFigure, strStyleFile=args.strStyleFile, sTaxaFileName=args.sTaxaFileName, iTerminalCladeLevel=args.iTerminalCladeLevel, sColorFileName=args.sColorFileName, sTickFileName=args.sTickFileName, sHighlightFileName=args.sHighlightFileName, sSizeFileName=args.sSizeFileName, sCircleFileName=args.sCircleFileName)

    #Write detail file out
    with open( args.strDetailOutputFile, 'w') as f:
      f.write(sCladogramDetails)

if __name__ == "__main__":
    _main( )
