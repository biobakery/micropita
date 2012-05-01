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
from ValidateData import ValidateData

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
argp.add_argument( "strStyleFile", metavar = "Style_file", help = Constants_Arguments.c_strCircladerStyleFile )
argp.add_argument(Constants_Arguments.c_strTaxaFilePath, dest="iTargetedTaxaFile", metavar= "TaxaFilePath", default=None, help= Constants_Arguments.c_strTaxaFileHelp)
argp.add_argument(Constants_Arguments.c_strHighlightCladeFile, dest="iHighlightCladeFile", metavar= "HighlightFilePath", default=None, help= Constants_Arguments.c_strHighlightCladeHelp)
argp.add_argument(Constants_Arguments.c_strInvertArgument, dest = "fInvert", action = "store", default="False", help = Constants_Arguments.c_strInvertHelp )
argp.add_argument(Constants_Arguments.c_strRoot, dest = "sRoot", metavar= "root", default="None", help = Constants_Arguments.c_strRootHelp )
argp.add_argument(Constants_Arguments.c_strEnrichmentMethod, dest="strEnrichmentIndicatorMethod", metavar= "enrichmentIndicatorMethod", default="FDR", 
                  choices=Constants_Arguments.c_strEnrichmentChoices, help= Constants_Arguments.c_strEnrichmentMethodHelp)
argp.add_argument(Constants_Arguments.c_strCladeFilterLevel, dest = "iCladeFilterLevel", metavar= "CladeFilterLevel", default="None",help = Constants_Arguments.c_strCladeFilterLevelHelp)
argp.add_argument(Constants_Arguments.c_strCladeMeasureLevel, dest = "iCladeFilterMeasure", metavar= "CladeFilterLevelMeasure", default="None", help = Constants_Arguments.c_strCladeMeasureLevelHelp)
argp.add_argument(Constants_Arguments.c_strCladeFilteringMinLevel, dest = "iCladeFilterMinNumber", metavar= "CladeFilterMinNumber", default="None",
	help = Constants_Arguments.c_strCladeFilteringMinLevelHelp)
argp.add_argument(Constants_Arguments.c_strAbundanceFilterPercentile, dest = "iAbundanceFilterPercentile", metavar= "AbundanceFilterPercentile", default="None",
	help = Constants_Arguments.c_strAbundanceFilterPercentileHelp)
argp.add_argument(Constants_Arguments.c_strAbundanceFilterCutoff, dest = "iAbundanceFilterPercentCuttoff", metavar= "AbundanceFilterPercentCuttoff", default="None",
	help = Constants_Arguments.c_strAbundanceFilterCutoffHelp)
argp.add_argument(Constants_Arguments.c_strRingOrder, dest = "iRingOrder", metavar= "Ring Order", default=None, help = Constants_Arguments.c_strRingOrder)
argp.add_argument(Constants_Arguments.c_strCircladerTicks, dest = "iTicks", metavar= "Internal Dendrogram Ticks", default=None, help = Constants_Arguments.c_strCircladerTicksHelp)
argp.add_argument(Constants_Arguments.c_strSampleNameRowArgument, dest="iSampleNameRow", metavar= "SampleNameRow", default=0, help= Constants_Arguments.c_strSampleNameRowHelp)
argp.add_argument(Constants_Arguments.c_strFirstDataRow, dest="iFirstDataRow", metavar= "FirstDataRow", default=1, help= Constants_Arguments.c_strFirstDataRowHelp)
argp.add_argument(Constants_Arguments.c_strNormalizeArgument, dest = "fNormalize", action = "store", default="False", help = Constants_Arguments.c_strNormalizeHelp)
argp.add_argument(Constants_Arguments.c_strEnrichmentThreshold, dest = "dAlpha", action = "store", default = 0.05, help = Constants_Arguments.c_strEnrichmentThresholdHelp)

#Outputfile
argp.add_argument( "sTaxaFileName", metavar = "TaxaFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerTaxaFile )
argp.add_argument( "sColorFileName", metavar = "ColorFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerColorFile )
argp.add_argument( "sTickFileName", metavar = "TickFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerTickFile )
argp.add_argument( "sHighlightFileName", metavar = "HighlightFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerHighlightFile )
argp.add_argument( "sSizeFileName", metavar = "SizeFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerSizeFile )
argp.add_argument( "sCircleFileName", metavar = "CircleFile.txt", nargs = "?", help = Constants_Arguments.c_strCircladerCircleFile )
argp.add_argument( "strOutFigure", metavar = "SelectionCladogram.png", nargs = "?", help = Constants_Arguments.c_strCircladerOutputFigure )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError("".join(["Invalid log level: ",strLogLevel," Try one of the following: "]+Constants_Arguments.c_lsLoggingChoices))
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFigure)[0],".log"]), filemode = 'w', level=iLogLevel)

    #Invert
    c_fInvert = (args.fInvert.lower() == "true")
    #Normalize
    c_Normalize = (args.fNormalize.lower() == "true")

    #Position holders in t-stats data list
    c_IDINDEX = 0
    c_TSCOREINDEX = 1
    c_PVALUEINDEX = 2
    c_QVALUEINDEX = 3

    #Get colors
    objColors = Constants_Figures()
    objColors.invertColors(c_fInvert)

    #Root tree
    c_fBacteriaRooted = (not args.sRoot.lower() == "none")

    #Delimiter for highlight and label information in the highlight clade files
    c_strLabelDelim = "="

    #Delimiter for taxa/otu lineages
    c_strLineageDelim = "|"

    #Read in abundance data
    rawData = AbundanceTable()
    abundance,metadata = rawData.textToStructuredArray(tempInputFile=args.strInputFile, tempDelimiter=Constants.TAB, tempNameRow=int(args.iSampleNameRow), tempFirstDataRow=int(args.iFirstDataRow), tempNormalize=c_Normalize)
    strSampleID = abundance.dtype.names[0]
    lsAllSampleNames = abundance.dtype.names[1:]

    #All taxa in the study
    lsAllTaxa = [strTaxaId for strTaxaId in list(abundance[strSampleID])]

    #Terminal taxa in the study
    lsTerminalTaxa = funcGetTerminalNodes(lsAllTaxa, c_strLineageDelim)

    #Create a cladogram object
    cladogram = Cladogram()

    #Genus Level of the ancestry to filter on
    fFilterClades = (not args.iCladeFilterLevel.lower() == "none")
    iCladeLevelForFiltering = 0
    iCladeLevelForReducing = 0
    iCladeLevelMinimumCount = 0
    if(fFilterClades):
      iCladeLevelForFiltering = int(args.iCladeFilterMeasure)
      iCladeLevelForReducing = int(args.iCladeFilterLevel)
      iCladeLevelMinimumCount = int(args.iCladeFilterMinNumber)
    cladogram.setFilterByCladeSize(fCladeSizeFilter = fFilterClades, iCladeLevelToMeasure = iCladeLevelForFiltering, iCladeLevelToReduce = iCladeLevelForReducing, iMinimumCladeSize = iCladeLevelMinimumCount)

    #Abundance filtering
    fFilterAbundance = (not args.iAbundanceFilterPercentile.lower() == "none")
    dPercentileCutOff = 0.0
    dPercentageAbovePercentile = 0.0
    if(fFilterAbundance):
      dPercentileCutOff = float(args.iAbundanceFilterPercentile)
      dPercentageAbovePercentile = float(args.iAbundanceFilterPercentCuttoff)
    cladogram.setFilterByAbundance(fAbundanceFilter = fFilterAbundance, dPercentileCutOff = dPercentileCutOff,  dPercentageAbovePercentile = dPercentageAbovePercentile)

    #Read in selection file
    fHndlInput = open(args.strSelectionFile,'r')
    strSelection = fHndlInput.read()
    fHndlInput.close()

    #Dictionary to hold selection data
    dictSelection = dict()
    for strSelectionLine in filter(None,strSelection.split(Constants.ENDLINE)):
        astrSelectionMethod = strSelectionLine.split(Constants.COLON)
        dictSelection[astrSelectionMethod[0]] = filter(None,[lsSample.strip() for lsSample in astrSelectionMethod[1].split(Constants.COMMA)])

    #Set Color Data
    dictColors = {MicroPITA.c_DIVERSITY_1:objColors.invSimpsonColorN,
                  MicroPITA.c_DIVERSITY_2:objColors.chao1ColorN,
                  MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1:objColors.brayCurtisColorN,
                  MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_2:objColors.unifracColorN,
                  MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_3:objColors.weightedUnifracColorN,
                  MicroPITA.c_EXTREME_DISSIMILARITY_1:objColors.inBrayCurtisColorN,
                  MicroPITA.c_EXTREME_DISSIMILARITY_2:objColors.inUnifracColorN,
                  MicroPITA.c_EXTREME_DISSIMILARITY_3:objColors.inWeightedUnifracColorN,
                  MicroPITA.c_USER_RANKED:objColors.userRankedN,
                  MicroPITA.c_RANDOM:objColors.randomColorN,
                  MicroPITA.c_SVM_CLOSE:objColors.svmCloseN,
                  MicroPITA.c_SVM_FAR:objColors.svmFarN,
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
    if(not args.iTargetedTaxaFile == "None"):
      fhndlTaxaInput = open(args.iTargetedTaxaFile,'r')
      lsUserDefinedTaxa = filter(None,fhndlTaxaInput.read().split(Constants.ENDLINE))
      fhndlTaxaInput.close()

    lsUserDefinedHighlighted = list()
    if(not args.iHighlightCladeFile == "None"):
      fhndlHighlightInput = open(args.iHighlightCladeFile,'r')
      lsUserDefinedHighlighted = filter(None,fhndlHighlightInput.read().split(Constants.ENDLINE))
      fhndlHighlightInput.close()

      #Highlight and relabel data
      dictTaxaHighlights = dict()
      dictRelabel = dict()
      iLabelCount = 1
      for sHighlight in lsUserDefinedHighlighted:
        if c_strLabelDelim in sHighlight:
          lsHighlightElements = filter(None,sHighlight.split(c_strLabelDelim))
          dictTaxaHighlights[lsHighlightElements[0]] = MicroPITA.c_DIVERSITY_2
          dictRelabel[lsHighlightElements[0]] = "".join([str(iLabelCount),":",lsHighlightElements[1]])
          iLabelCount = iLabelCount + 1
        else:
          dictTaxaHighlights[sHighlight] = MicroPITA.c_DIVERSITY_2
      cladogram.addHighLights(dictTaxaHighlights,True)
      if len(dictRelabel) > 0:
        cladogram.relabelIDs(dictRelabel)

#    if(args.strEnrichmentIndicatorMethod=="ABS"):#Abscense and Presence
#      #This is perfroming absense and presence
#      #Add circles for each method
#      dictTaxaInAllSelection = dict()
#
#      #Allows one to set the order of the rings if needed
#      lsSelectedSampleMethod = dictSelection.keys()
#      if not args.iRingOrder == None:
#        lsSelectedSampleMethod = [filter(None,strMethod) for strMethod in (args.iRingOrder.split(Constants.COMMA))]
#
#      for selectedSampleMethod in lsSelectedSampleMethod:
#        if selectedSampleMethod in dictSelection:
#          #All samples that were selected by the method
#          lsSelectedSamples = dictSelection[selectedSampleMethod]
#          #Get taxa interesting to this selection technique
#          setMethodTaxaTotal = set()
#          #For each selected sample get absence or presence of taxa
#          for strSelectedSample in lsSelectedSamples:
#            #Get taxa of each sample
#            lsSampleTaxaAbundance = abundance[strSelectedSample]
#            #Will hold taxa that have more than 0 abundance
#            lsSelectedTaxaPerSample = list()
#            #Go through each taxa abundance using the index to match it to the taxa name
#            for iAbundanceIndex in xrange(0,len(lsSampleTaxaAbundance)):
#              scurTaxon = abundance[strSampleID][iAbundanceIndex]
#              #If the taxa has more than 0.0 abundance add
#              if(float(lsSampleTaxaAbundance[iAbundanceIndex]) > 0.0):
#                lsSelectedTaxaPerSample.append(scurTaxon)
#            #Combine (union) taxa to a set
#            setMethodTaxaTotal = setMethodTaxaTotal | set(lsSelectedTaxaPerSample)
#
#          #Add circle for this data
#          #cladogram.addCircle(lsTaxa=lsTaxa, strShape=lsShapes, dAlpha=lsAlpha, strCircle=selectedSampleMethod, fForced=True)

    if((args.strEnrichmentIndicatorMethod=="PVALUE") or (args.strEnrichmentIndicatorMethod=="FDR")):#P-value or qvalue
      #Are we using qvalues or pvalues
      fIsQValue = (args.strEnrichmentIndicatorMethod=="FDR")

      #Allows one to set the order of the rings if needed
      lsSelectedSampleMethod = dictSelection.keys()
      if not args.iRingOrder == None:
        lsSelectedSampleMethod = [filter(None,strMethod) for strMethod in (args.iRingOrder.split(Constants.COMMA))]

      #This is performing t-test with p-values
      for selectedSampleMethod in lsSelectedSampleMethod:
        if selectedSampleMethod in dictSelection:
          #All samples that were selected by the method
          lsSelectedSamples = dictSelection[selectedSampleMethod]
          #Samples selected
          lfSelectedSamplesUTest = []
          lfNotSelectedSamplesUTest = []
          #Set up boolean list to compress array
          for sample in lsAllSampleNames:
            wasSelected = sample in lsSelectedSamples
            lfSelectedSamplesUTest.append(wasSelected)
            lfNotSelectedSamplesUTest.append(not wasSelected)

          #Holds t-tests,pvalues,qvalues as needed
          lsTaxaTScores = list()

          #Compress arrays to one or the other distribution
          #Conduct wilcoxon tests on all taxa
          for iTaxonIndex in xrange(0,len(lsAllTaxa)):
            npaTaxaData = list(abundance[iTaxonIndex,])
            strTaxaId = npaTaxaData[0]
            if(strTaxaId in lsTerminalTaxa):
              npaDistribution = np.array(npaTaxaData[1:])
              npaSelectedDistribution = np.compress(lfSelectedSamplesUTest,npaDistribution)
              npaNotSelectedDistribution = np.compress(lfNotSelectedSamplesUTest,npaDistribution)
              dScore, dPvalue = stats.ranksums(npaSelectedDistribution,npaNotSelectedDistribution)
              #[ID,TScore,PValue,QValue,SortOrder]
              lsTaxaTScores.append([strTaxaId,dScore,dPvalue,-1])

          #Get a list of pvalues preserving order
          ldOrderedPValues = list()
          for iScoreDataIndex in lsTaxaTScores:
            ldOrderedPValues.append(iScoreDataIndex[c_PVALUEINDEX])

          #If using qvalues generate them with FDR BH
          if fIsQValue:
            #Convert pvalues to qvalues
            ldOrderedQValues = Utility_Math.convertToBHQValue(ldOrderedPValues)

            #Update the score data with qvalues
            for iQIndex in xrange(0,len(ldOrderedQValues)):
              lsTaxaTScores[iQIndex][c_QVALUEINDEX] = ldOrderedQValues[iQIndex]

          lsAlpha = list()
          lsShapes = list()
          lsTaxa = list()
          for iTaxa in xrange(0,len(lsTaxaTScores)):
            lsCur = lsTaxaTScores[iTaxa]
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
                lsShapes.append("^")
              elif(dCurScore < 0.0):
                lsAlpha.append(str(1-dValue))
                lsShapes.append("v")
              elif(dCurScore == 0.0):
                lsAlpha.append("0.0")
                lsShapes.append("R")
            else:
              lsAlpha.append("0.0")
              lsShapes.append("R")

          #Add circle for this data
          cladogram.addCircle(lsTaxa=lsTaxa, strShape=lsShapes, dAlpha=lsAlpha, strCircle=selectedSampleMethod, fForced=True)

    #Generate cladogram
    cladogram.generate(strInputFile=args.strInputFile, strImageName=args.strOutFigure, strStyleFile=args.strStyleFile, sTaxaFileName=args.sTaxaFileName, charDelimiter=Constants.TAB, iNameRow=int(args.iSampleNameRow), iFirstDataRow=int(args.iFirstDataRow), fNormalize=False, sColorFileName=args.sColorFileName, sTickFileName=args.sTickFileName, sHighlightFileName=args.sHighlightFileName, sSizeFileName=args.sSizeFileName, sCircleFileName=args.sCircleFileName)

##Returns only terminal nodes given the list's structure
def funcGetTerminalNodes(lsTaxa,cDelim):
  #Return list
  lsRetList = list()

  #Build hash
  dictCounts = dict()
  for strTaxaName in lsTaxa:
    #Split into the elements of the clades
    lsClades = filter(None,strTaxaName.split(cDelim))
    #Count clade levels
    iCladeLength = len(lsClades)
    
    #Make sure there is data to work with
    if iCladeLength < 0:
      pass

    #Evaluate first element
    sClade = lsClades[0]
    if sClade in dictCounts:
      dictCounts[sClade] = False
    else:
      dictCounts[sClade] = True

    #Evaluate the rest of the elements
    if iCladeLength > 1:
      for iIndex in xrange(1,iCladeLength):
        prevClade = sClade
        sClade = cDelim.join([sClade,lsClades[iIndex]])
        if sClade in dictCounts:
          dictCounts[sClade] = False
          dictCounts[prevClade] = False
        else:
          dictCounts[sClade] = True
          dictCounts[prevClade] = False

  #Return only the elements that were of count 1
  for sName in dictCounts:
    if dictCounts[sName]==True:
      lsRetList.append(sName)
  return lsRetList

if __name__ == "__main__":
    _main( )
