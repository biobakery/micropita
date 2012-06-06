import re
import sfle
import sys
import os

Import( "*" )

#TODO make this dynamic
sys.path.append("/home/ttickle/Desktop/ttickle/sfle/input/micropita/src")
from Constants_Arguments import Constants_Arguments
#from MicroPITA import MicroPITA

pE = DefaultEnvironment( )

#Process flags
c_fGenerateInSilicoDataSets = False
c_fRunCollectionCurve = True
fMakeMovies = False

#Normalize
c_NormalizePCOA = "False"

#Extentions
c_strSufTable = ".pcl"
c_strPCOAMovieEnding = "/PCOA.avi"
c_strStratPCOAMovieEnding = "/Strat-PCOA.avi"
c_strCladogramMovieEnding = "/Cladogram.avi"
c_strSufActualData = "-actual.txt"
c_strSufCircCircle = ".CCircle"
c_strSufCircColor = ".CColor"
c_strSufCircDetail = ".CDetail"
c_strSufCircHighlight = ".CHLight"
c_strSufCircSize = ".CSize"
c_strSufCircStyle = ".CStyle"
c_strSufCircStyleInv = "-invert.CStyle"
c_strSufCircTaxa = ".CTaxa"
c_strSufCircTick = ".CTick"
c_strSufCheckedTable = "".join(["-checked",c_strSufTable])
c_strSufCollectionCurveFigure = "-ColCurv.pdf"
c_strSufCollectionCurveText = "-ColCurv.txt"
c_strSufCombinedPCOA = "-Combined.png"
c_strSufCombinedStratPCOA = "-Combined-StratPCoA.pdf"
c_strSufConfig = ".config"
c_strSufConfusionMatrix = "-Confusion.png"
c_strSufFig2 = "-Fig2.pdf"
c_strSufHCLUSTColor = ".HCLColor"
c_strSufHCLUSTData = ".HCLData"
c_strSufHCLUSTFig = "-SHCL.pdf"
c_strSufHCLUSTLabel = ".HCLLabel"
c_strSufInsilicoData = c_strSufTable
c_strSufMetaMatrix = ".meta"
c_sufMetaConfusionMatrix = "-metaConfusion.png"
c_sufMetaOverlapMatrix = "-metaOverlap.pdf"
c_strSufMicropita = ".txt"
c_strSufOverlapMatrix = "-overlap.pdf"
c_strSufPCOA = "-PCoA.pdf"
c_strSufPng = ".pdf"
c_strSufPredict = "-SVM.predict"
c_strSufSelectedTaxa = ".taxa"
c_strSufStratHCLUSTColor = "-Strat.HCLColor"
c_strSufStratHCLUSTFig = "-StratHCL.pdf"
c_strSufStratHCLUSTLabel = "-Strat.HCLLabel"
c_strSufStratPCOA = "-StratPCoA.pdf"
c_strSufUncheckedTable = c_strSufTable
c_strSufValidatedDiversity = "-ValidatedMaxDiv.pdf"
c_strSufValidatedFeature = "-ValidatedFeature.pdf"
c_strSufValidatedFeatureHistogram = "-ValidatedFeatureHist.pdf"
c_strSufValidatedExtreme = "-ValidatedMaxDis.pdf"
c_strSufValidatedRepresentative = "-ValidatedMaxRep.pdf"
c_strSufValidatedDistinct = "-ValidatedDistinct.pdf"
c_strSufValidatedDiscriminant = "-ValidatedDisrciminant.pdf"

#Dict keys
c_strSelectionFiles = "SelectionFiles"

#Special characters
c_strConfigFileHeaderChar = "["
c_strConfigFileCommentChar = "#"
c_strExtDelim = "."
c_strPathDelim = "/"
c_strComma = ","

#Metrics
c_DIVERSITY = "Diversity"
c_TARGETED_TAXA = "Taxa_Defined"
c_EXTREME = "Extreme"
c_EXTREME_DISSIMILARITY_1 = "Extreme_B"
c_REPRESENTATIVE = "Representative"
c_REPRESENTATIVE_DISSIMILARITY_1 = "Representative_B"
c_DISCRIMINANT = "Discriminant"
c_DISTINCT = "Distinct"

#Directories
strOutputSummaryFolder = "".join([c_strPathDelim,"Metasummary"])
strDataFolder = "Data"
fileDirDataName = c_strPathDelim.join([fileDirInput.get_abspath(),strDataFolder])

#Data files to generate
c_strInsilicoDataDiversity = c_strPathDelim.join([strDataFolder,"DiversityTest"+c_strSufInsilicoData])
c_strInsilicoDataDiversityKey = "Diversity"
c_strInsilicoDataRepresentative = c_strPathDelim.join([strDataFolder,"RepresentativeTest"+c_strSufInsilicoData])
c_strInsilicoDataUnbalanced = c_strPathDelim.join([strDataFolder,"UnbalancedTest"+c_strSufInsilicoData])
c_strInsilicoDataUnbalancedKey = "Unbalanced"
lsInsilicoDataNames = [c_strInsilicoDataDiversity,c_strInsilicoDataRepresentative,c_strInsilicoDataUnbalanced]

#The Column that holds the TAXA / OTU IDS
c_strAbundanceIDCol = "0"

#Scons Configuration file headers
c_strConfigAbundanceFilter = "[Abundance Filter]"
c_strConfigAbundanceFilterPercentile = "[Abundance Filter Percentile]"
c_strConfigAbundanceFilterPercent = "[Abundance Filter Percent Above Percentile]"
c_strConfigActualFile = "[File Containing Actual Classes of Samples if Known]"
c_strConfigCladeFilter = "[Clade Filter]"
c_strConfigCladeFilterMeasure = "[Clade Level to Measure]"
c_strConfigCladeFilterLevel = "[Clade Level to Filter]"
c_strConfigCladeFilterMinSize = "[Minimum Clade Size]"
c_strConfigCladogramAlpha = "[Enrichment Threshold]"
c_strConfigCladogramRingOrder = "[Cladogram Ring Order]"
c_strConfigCladogramTicks = "[Cladogram Dendrogram Ticks]"
c_strConfigLastMetadataRow = "[Last Metadata Name]"
c_strConfigEnrichmentMeasurement = "[Taxa or OTU Enrichment Measurement]"
c_strConfigFileIsNormalized = "[File is Already Normalized]"
c_strConfigFileIsSummed = "[File is Already Summed]"
c_strConfigHighlightClades = "[Higlight Clades]"
c_strConfigInputFile = "[Input File]"
c_strConfigInvertImage = "[Invert Image]"
c_strConfigLabels = "[Labels]"
c_strConfigLogging = "[Logging]"
c_strConfigNormalizeAbundance = "[Show Plots as Normalized Abundance]"
c_strOccurenceFilter = "[Occurence Filter]"
c_strOccurenceFilterMinSequence = "[Occurence Filter Minimum Sequence]"
c_strOccurenceFilterMinSample = "[Occurence Filter Minimum Sample]"
c_strConfigPlotCombinedSelectionTechniques = "[Unsupervised Selection Methods to plot in a combined graph]"
c_strConfigPlotCombinedSupervisedSelectionTechniques = "[Supervised Selection Methods to plot in a combined graph]"
c_strConfigPairingMetadata = "[Pairing Metadata]"
c_strConfigProjects = "[Projects]"
c_strConfigRoot = "[Root]"
c_strConfigSampleRow = "[Sample ID Name]"
c_strConfigSelection = "[Selection]"
c_strConfigSelectedTaxa = "[Selected Taxa]"
c_strConfigSelectionTechniques = "[Selection Techniques]"
c_strConfigSelectionTechniquesCollectorCurve = "[Methods to plot in collection curve]"
c_strConfigSumData = "[Sum Clades for Analysis/Plotting]"
c_strConfigSupervisedLabel = "[Supervised Label]"
c_strConfigSupervisedCount = "[Supervised Selection Count]"
c_strConfigTargetedSelection = "[Targeted Feature Selection Method]"
c_strConfigTerminalLevel = "[Show Clade Level]"
c_strConfigUnsupervisedCount = "[Unsupervised Selection Count]"
c_strConfigUnsupervisedStratify = "[Stratify by Metadata]"
c_strConfigValidationFile = "[Validation Input File]"
c_strConfigValidationIsNormalized = "[Validation File is Already Normalized]"
c_strConfigValidationIsSummed = "[Validation File is Already Summed]"
c_strConfigValidationIDName = "[Validation Sample ID Name]"
c_strConfigValidationLastMetadataName = "[Validation Last Metadata Name]"
c_strConfigValidationMethodologies = "[Validation Selection Techniques]"

#Cladogram style file for micropita
c_fileCladogramStyleFile = File(sfle.d( fileDirInput, c_strPathDelim.join([strDataFolder,"microPITA"])+c_strSufCircStyle ))

#External programs
#c_progHCLPath = "./external/hclust/hclust.py"
c_progHCL = "../../external/hclust/hclust.py"

#Src Code
c_fileProgAbundanceTable = File( sfle.d( fileDirSrc, "AbundanceTable.py" ) )
c_fileProgCladogram = File( sfle.d( fileDirSrc, "Cladogram.py" ) )
c_fileProgCheckFile = File( sfle.d( fileDirSrc, "CheckAbundanceTable.py"))
c_fileProgCombinedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCombinedPCoA.py") )
c_fileProgCombinedStratifiedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCombinedStratifiedPCoA.py" ) )
c_fileProgCommandLine = File( sfle.d( fileDirSrc, "CommandLine.py" ) )
c_fileProgCollectionCurveFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCollectionCurve.py") )
c_fileProgConfusionMatrixFigure = File( sfle.d( fileDirSrc, "MicropitaPaperConfusionMatrix.py") )
c_fileProgConstants = File( sfle.d( fileDirSrc, "Constants.py" ) )
c_fileProgConstantsFigures = File( sfle.d( fileDirSrc, "Constants_Figures.py" ) )
c_fileProgDiversity = File( sfle.d( fileDirSrc, "Diversity.py" ) )
c_fileProgMencoder = File("/usr/bin/mencoder")
c_fileProgConfusionMatrixMetaFigure = File( sfle.d( fileDirSrc, "MicropitaPaperConfusionMetaMatrix.py" ) )
c_fileProgOverlapMatrixMetaFigure = File( sfle.d( fileDirSrc, "MicropitaPaperOverlapMetaMatrix.py" ) )
c_fileProgMicroPITA = File( sfle.d( fileDirSrc, "MicroPITA.py" ) )
c_fileProgMLPYDistanceAdaptor = File( sfle.d( fileDirSrc, "MLPYDistanceAdaptor.py" ) )
c_fileProgOverlapMatrixFigure = File( sfle.d( fileDirSrc, "MicropitaPaperOverlapMatrix.py" ) )
c_fileProgPCOA = File( sfle.d( fileDirSrc, "PCoA.py" ) )
c_fileProgPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperPCoA.py" ) )
c_fileProgSelectionCladogramFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionCladogram.py") )
c_fileProgSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionHCL.py" ) )
c_fileProgStratifiedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperStratifiedPCoA.py" ) )
c_fileProgStratSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperStratSelectionHCL.py" ) )
c_fileProgSVM = File( sfle.d( fileDirSrc, "SVM.py" ) )
c_fileProgUtilityData = File( sfle.d( fileDirSrc, "MicropitaPaperConstructDataSets.py" ) )
c_fileProgValidateData = File( sfle.d( fileDirSrc, "ValidateData.py" ) )
c_fileProgValidateDiversity = File( sfle.d( fileDirSrc, "MicropitaPaperValidateDiversity.py" ))
c_fileProgValidateFeature = File( sfle.d( fileDirSrc, "MicropitaPaperValidateFeatureAbundance.py" ))
c_fileProgValidateFeatureHistogram = File( sfle.d( fileDirSrc, "MicropitaPaperValidateFeatureAbundanceHistogram.py" ))
c_fileProgValidatePCoA = File( sfle.d( fileDirSrc, "MicropitaPaperValidatePCoA.py" ))

#Lists of sources needed for different python scripts that are ran so that they can be used as libraries
c_filesSecondarySrc = [c_fileProgAbundanceTable, c_fileProgCommandLine, c_fileProgConstants, c_fileProgDiversity,
	               c_fileProgMLPYDistanceAdaptor, c_fileProgPCOA, c_fileProgSVM, c_fileProgUtilityData]
c_filePrimarySrc = c_filesSecondarySrc + [c_fileProgMicroPITA]
ls_srcFig1 = [c_fileProgCommandLine, c_fileProgConstants, c_fileProgConstantsFigures, c_fileProgMicroPITA, c_fileProgValidateData]
ls_srcFig2 = [c_fileProgConstants, c_fileProgConstantsFigures, c_fileProgMicroPITA, c_fileProgValidateData, c_fileProgCladogram, c_fileProgAbundanceTable]

#Selection parameters
c_strDiversity = "Diversity"
c_strExtremeDissimilarity = "Extreme"
c_strDiscriminant = "Discriminant"
c_strDistinct = "Distinct"
c_strRandom = "Random"
c_strRepresentativeDissimilarity = "Representative"
c_lsSupervisedMethods = set([c_strDiscriminant,c_strDistinct])
c_strTaxa = "Taxa_Defined"
c_lsUnsupervisedMethods = set([c_strDiversity,c_strExtremeDissimilarity,c_strRandom,c_strRepresentativeDissimilarity,c_strTaxa])
#c_lstrAllSelectionMethods = [c_strRandom,c_strTaxa,c_strDiversity,c_strRepresentativeDissimilarity,c_strExtremeDissimilarity,c_strDiscriminant,c_strDistinct]
#c_lstrAllSelectionMethods = [c_strTaxa,c_strDiversity,c_strRepresentativeDissimilarity,c_strExtremeDissimilarity,c_strDiscriminant,c_strDistinct]

##Start flow
#Scripts to call
##Generate insilico data sets
def funcGenerateUnbalancedData(strDataSetKey):
  def funcGenerateUnbalancedDataRet( target, source, env, strDataSetKey=strDataSetKey ):
    strT, astrSs = sfle.ts( target, source )
    strProg = astrSs[0]
    return sfle.ex( [strProg,strT,strDataSetKey, target[1].get_abspath()] )
  return funcGenerateUnbalancedDataRet

def funcGenerateDiversityData(strDataSetKey):
  def funcGenerateDiversityDataRet( target, source, env, strDataSetKey=strDataSetKey ):
    strT, astrSs = sfle.ts( target, source )
    strProg = astrSs[0]
    return sfle.ex( [strProg,strT,strDataSetKey] )
  return funcGenerateDiversityDataRet

#Visualize insilico data
def funcVisualizeInsilicoData( strXStart, strYStart ):
  def funcVisualizeInsilicoDataRet( target, source, env, strXStart=strXStart, strYStart=strYStart ):
    strT, astrSs, = sfle.ts( target, source )
    strProg, strInputDataFile = astrSs[0], astrSs[1]
    return sfle.ex([strProg, "--in", strInputDataFile, "--out", strT, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1","--flabel","1","--font_size","5","--xstart",strXStart,"--ystart",strYStart])
  return funcVisualizeInsilicoDataRet

#Check abundance table data
def funcCheckAbundanceData( sLogging, sMetadataName ):
  def funcCheckAbundanceDataRet( target, source, env, sLogging=sLogging, sMetadataName=sMetadataName ):
    strT, astrSs, = sfle.ts( target, source )
    strProg, strInputDataFile = astrSs[0], astrSs[1]
    return sfle.ex([strProg, sLogging, sMetadataName, strInputDataFile, strT])
  return funcCheckAbundanceDataRet

##Call micropita to run
def funcMicroPita( strLoggingLevel, sIDName, sLastMetadata, fNormalized, fSummed, fSumData, sFeatureSelection, iUnsupervisedCount, strTaxaFile, strUnsupervisedStratify, strSupervisedLabel, iSupervisedCount, strDirTmp, lsSelectionMethods):
  def funcMicroPitaRet( target, source, env, strLoggingLevel=strLoggingLevel, sIDName=sIDName, sLastMetadata=sLastMetadata, fNormalized=fNormalized, fSummed=fSummed, fSumData=fSumData, sFeatureSelection=sFeatureSelection, iUnsupervisedCount=iUnsupervisedCount, strTaxaFile=strTaxaFile, strUnsupervisedStratify=strUnsupervisedStratify, strSupervisedLabel=strSupervisedLabel, iSupervisedCount=iSupervisedCount, strDirTmp=strDirTmp, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd = astrSs[0], astrSs[1]
    lsCommandline = [strProg,strLoggingLevel, sIDName, sLastMetadata]
    if fNormalized:
        lsCommandline = lsCommandline + [Constants_Arguments.c_strIsNormalizedArgument]
    if fSummed:
        lsCommandline = lsCommandline + [Constants_Arguments.c_strIsSummedArgument]
    if not fSumData:
        lsCommandline = lsCommandline + [Constants_Arguments.c_strSumDataArgument]
    return sfle.ex( lsCommandline + [sFeatureSelection, iUnsupervisedCount, strTaxaFile, strUnsupervisedStratify, strSupervisedLabel,
                    iSupervisedCount, strDirTmp]+[strAbnd, strT] + lsSelectionMethods)
  return funcMicroPitaRet

##Create figures
#Visualize output with PCoA (Figure 1A)
def funcPCoASelectionMethods( strLoggingLevel, sIDName, sLastMetadata, iNormalize, strInvert,
                              fNormalized, fSummed, fSumData, strPredictFileArgument ):
  def funcPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, sIDName=sIDName, sLastMetadata=sLastMetadata,
                   iNormalize=iNormalize, strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed, fSumData=fSumData,
                   strPredictFileArgument=strPredictFileArgument):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[strLoggingLevel, sIDName,sLastMetadata,iNormalize,strInvert,fNormalized,
                   fSummed, fSumData, strPredictFileArgument]+[ strAbnd, strSelection, strT])
  return funcPCoARet

#Visualize output with Combined PCoA (Figure 1A)
def funcCombinedPCoASelectionMethods( strLoggingLevel, sIDName, sLastMetadata, iNormalize, strInvert,
                                      fNormalized, fSummed, fSumData, lsSelectionMethods):
  def funcCombinedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, sIDName=sIDName,
                           sLastMetadata=sLastMetadata, iNormalize=iNormalize, fNormalized=fNormalized,
                           fSummed=fSummed, fSumData=fSumData, strInvert=strInvert, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[strLoggingLevel,sIDName,sLastMetadata,iNormalize,fNormalized,fSummed,fSumData,strInvert]+[ strAbnd, strSelection, strT]+lsSelectionMethods)
  return funcCombinedPCoARet

#Visualize selected output with HCL (Figure 1B)
#def funcHCLSelectionMethods( strLoggingLevel, strInvert ):
#  def funcHCLSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert ):
#    strT, astrSs = sfle.ts( target, source )
#    strProg, strSelection = astrSs[0], astrSs[1]
#    strData, strColor, strLabel = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath()
#    return sfle.ex([strProg, strLoggingLevel, strInvert, strSelection, c_progHCLPath, strData, strColor, strLabel, strT])
#  return funcHCLSelectionRet

#Visualize selected output with Confusion matrix (Figure 1B Alt 1)
def funcConfusionMatrix( strLoggingLevel, strInvert, lsSelectionMethods ):
  def funcConfusionMatrixRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, lsSelectionMethods=lsSelectionMethods ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection, strActualFile = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg, strLoggingLevel, strInvert, strSelection, strActualFile, strT]+lsSelectionMethods)
  return funcConfusionMatrixRet

#Visualize selected output with Overap Matrix (Figure 1B Alt 2)
def funcOverlapMatrix( strLoggingLevel, strInvert, lsSelectionMethods ):
  def funcOverlapMatrixRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, lsSelectionMethods=lsSelectionMethods ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection = astrSs[0], astrSs[1]
    return sfle.ex([strProg, strLoggingLevel, strInvert, strSelection, strT]+lsSelectionMethods)
  return funcOverlapMatrixRet

#Visualize selected output with Cladogram (Figure 2)
def funcCladogramSelectionMethods( strLoggingLevel, strTargetedTaxaFile, strHighlightCladeFile, sIDName, sLastMetadata, iNormalize, strInvert,
                                   fNormalized, fSummed, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                                   strAbundanceFilterPercentile, strAbundanceFilterPercent, strRingOrder, strTicks, strAlpha, strMinSequence, strMinSample, fSumData, iTerminalClade):

  def funcCladogramSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strTargetedTaxaFile=strTargetedTaxaFile, strHighlightCladeFile=strHighlightCladeFile, 
                                 sIDName=sIDName, sLastMetadata=sLastMetadata, iNormalize=iNormalize, strInvert=strInvert,
                                 fNormalized=fNormalized, fSummed=fSummed, strRoot=strRoot, strEnrichment=strEnrichment,
                                 strCladeFilterLevel=strCladeFilterLevel, strCladeFilterMeasure=strCladeFilterMeasure, strCladeFilterMin=strCladeFilterMin, strAbundanceFilterPercentile=strAbundanceFilterPercentile,
                                 strAbundanceFitlerPercent=strAbundanceFilterPercent, strRingOrder=strRingOrder, strTicks=strTicks, strAlpha=strAlpha, strMinSequence=strMinSequence, strMinSample=strMinSample, fSumData=fSumData, iTerminalClade=iTerminalClade):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection, strAbundance, strStyleFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
    strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile, strDetailFile = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath(), target[4].get_abspath(), target[5].get_abspath(), target[6].get_abspath(), target[7].get_abspath()
    return sfle.ex([strProg] + [strTargetedTaxaFile, strHighlightCladeFile, sIDName, sLastMetadata, iNormalize, strInvert, fNormalized, fSummed, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                               strAbundanceFilterPercentile, strAbundanceFilterPercent, strAlpha, strMinSequence, strMinSample, fSumData, iTerminalClade, strRingOrder, strTicks] + [strSelection, strAbundance, strStyleFile, strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile, strT, strDetailFile])
  return funcCladogramSelectionRet

#Visualize with HCL selected samples with stratification (Figure 3)
#def funcHCLStratSelectionMethods( strLoggingLevel, strInvert, strDataRow, strIdCol ):
#  def funcHCLStratSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, strDataRow=strDataRow, strIdCol=strIdCol ):
#    strT, astrSs = sfle.ts( target, source )
#    strProg, strSelection, strAbundance = astrSs[0], astrSs[1], astrSs[2]
#    strColor, strLabel = target[1].get_abspath(), target[2].get_abspath()
#    return sfle.ex([strProg, strLoggingLevel, strDataRow, strIdCol, strInvert, strSelection, strAbundance, c_progHCLPath, strColor, strLabel, strT])
#  return funcHCLStratSelectionRet

#Visualize output with Stratified PCoAs (Figure 4)
def funcStratifiedPCoASelectionMethods( strLoggingLevel, sIDName, sLastMetadata, strStratifyMetadata, iNormalize, strInvert,
                                        fNormalized, fSummed, fSumData, strPredictFileArgument ):
  def funcStratifiedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, sIDName=sIDName,
                             sLastMetadata=sLastMetadata, strStratifyMetadata=strStratifyMetadata, iNormalize=iNormalize,
                             strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed, fSumData=fSumData,
                             strPredictFileArgument=strPredictFileArgument):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[sIDName,sLastMetadata,strStratifyMetadata,iNormalize,strInvert,
                              fNormalized,fSummed,fSumData,strPredictFileArgument]+[ strAbnd, strSelection, strT])
  return funcStratifiedPCoARet

#Visualize output with Stratified PCoAs (Figure 3/4 Combined)
def funcCombinedStratifiedPCoASelectionMethods( strLoggingLevel, sIDName, sLastMetadata, strStratifyMetadata, iNormalize, strInvert,
                                                fNormalized, fSummed, fSumData, lsSelectionMethods ):
  def funcCombinedStratifiedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, sIDName=sIDName,
                                     sLastMetadata=sLastMetadata, strStratifyMetadata=strStratifyMetadata, iNormalize=iNormalize,
                                     strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed, fSumData=fSumData,
                                     lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[sIDName,sLastMetadata,strStratifyMetadata,iNormalize,strInvert, fNormalized,
                   fSummed, fSumData]+[ strAbnd, strSelection, strT] + lsSelectionMethods)
  return funcCombinedStratifiedPCoARet

#Create a summary chart with collection curves
def funcCollectionCurveSummary( strLoggingLevel, strSampleNameRow, strLastMetadataName, strInvert, fNormalized, fSummed, lsSelectionMethods):
  def funcCollectionCurveRet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow,
                              strLastMetadataName=strLastMetadataName, strInvert=strInvert, fNormalized=fNormalized,
                              fSummed=fSummed, lsSelectionMethods=lsSelectionMethods):
    strFigureT, astrSs = sfle.ts( target, source )
    strProg, strAbundance = astrSs[0], astrSs[1]
    #Get input selection files (discriminate from the src files also passed as a list)
    #Check to see where the input selection files end and the python src begin
    iFileCount = 0
    for strFile in astrSs[1:]:
      if strFile[-3:] == ".py":
        break
      iFileCount = iFileCount + 1
    lsSelectionFiles = astrSs[2:iFileCount+1]
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strLastMetadataName, fNormalized,
                    fSummed, strInvert, strAbundance, strFigureT]+lsSelectionFiles+[lsSelectionMethods])
  return funcCollectionCurveRet

#Validation steps
#ValidateDiversity
def funcMeasureDiversityByGroup( strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName,
                                 strInvert, fNormalized, fSummed, fValidateNormalized, fValidateSummed, sPair):
  def funcMeasureDiversityByGroupRet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow,
                                      strValidateSampleNameRow=strValidateSampleNameRow, strValidateLastMetadataName=strValidateLastMetadataName,
                                      strLastMetadataName=strLastMetadataName, strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed,
                                      fValidateNormalized=fValidateNormalized, fValidateSummed=fValidateSummed, sPair=sPair):
    strFigureT, astrSs = sfle.ts( target, source )
    strProg, strValidationAbundance, strAbundance, strSelectionFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName, strInvert, fNormalized,
                    fSummed, fValidateNormalized, fValidateSummed, sPair, strValidationAbundance, strAbundance, strSelectionFile, strFigureT])
  return funcMeasureDiversityByGroupRet

#ValidateFeature
def funcMeasureFeatureByGroup( strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName,
                               strInvert, fNormalized, fSummed, fValidateNormalized, fValidateSummed, sPair, sMeasureMethod):
  def funcMeasureFeatureByGroupRet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow,
                                      strValidateSampleNameRow=strValidateSampleNameRow, strValidateLastMetadataName=strValidateLastMetadataName,
                                      strLastMetadataName=strLastMetadataName, strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed,
                                      fValidateNormalized=fValidateNormalized, fValidateSummed=fValidateSummed, sPair=sPair, sMeasureMethod=sMeasureMethod):
    strFigureT, astrSs = sfle.ts( target, source )
    strProg, strValidationAbundance, strAbundance, strSelectionFile, strFeatureFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3], astrSs[4]
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName, strInvert, fNormalized,
                    fSummed, fValidateNormalized, fValidateSummed, sPair, sMeasureMethod, strValidationAbundance, strAbundance, strSelectionFile, strFeatureFile, strFigureT])
  return funcMeasureFeatureByGroupRet

def funcMeasureFeatureByGroupHistogram( strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName,
                               strInvert, fNormalized, fSummed, fValidateNormalized, fValidateSummed, sPair):
  def funcMeasureFeatureByGroupHistogramRet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow,
                                      strValidateSampleNameRow=strValidateSampleNameRow, strValidateLastMetadataName=strValidateLastMetadataName,
                                      strLastMetadataName=strLastMetadataName, strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed,
                                      fValidateNormalized=fValidateNormalized, fValidateSummed=fValidateSummed, sPair=sPair):
    strFigureT, astrSs = sfle.ts( target, source )
    strProg, strValidationAbundance, strAbundance, strSelectionFile, strFeatureFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3], astrSs[4]
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName, strInvert, fNormalized,
                    fSummed, fValidateNormalized, fValidateSummed, sPair, strValidationAbundance, strAbundance, strSelectionFile, strFeatureFile, strFigureT])
  return funcMeasureFeatureByGroupHistogramRet

#Validate Extreme, Representative
def funcValidateInPCoA( strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName,
                        strInvert, fNormalized, fSummed, fValidateNormalized, fValidateSummed, sPair, sMetric, sStratLabel):
  def funcValidateInPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow,
                                      strValidateSampleNameRow=strValidateSampleNameRow, strValidateLastMetadataName=strValidateLastMetadataName,
                                      strLastMetadataName=strLastMetadataName, strInvert=strInvert, fNormalized=fNormalized, fSummed=fSummed,
                                      fValidateNormalized=fValidateNormalized, fValidateSummed=fValidateSummed, sPair=sPair, sMetric=sMetric, sStratLabel=sStratLabel):
    strFigureT, astrSs = sfle.ts( target, source )
    strProg, strValidationAbundance, strAbundance, strSelectionFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strLastMetadataName, strValidateSampleNameRow, strValidateLastMetadataName, strInvert, fNormalized,
                    fSummed, fValidateNormalized, fValidateSummed, sPair, sMetric, sStratLabel, strValidationAbundance, strAbundance, strSelectionFile, strFigureT])
  return funcValidateInPCoARet

#Meta matrix creation
def funcOverlapMetaMatrix( strLoggingLevel, strInvert, sSelectionMethods ):
  def funcOverlapMetaMatrixRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, sSelectionMethods=sSelectionMethods ):
    strT, astrSs = sfle.ts( target, source )
    strProg, lsSelectionFiles = astrSs[0], astrSs[1:]
    return sfle.ex([strProg, strLoggingLevel, strInvert, strT, sSelectionMethods]+lsSelectionFiles)
  return funcOverlapMetaMatrixRet

def funcConfusionMetaMatrix( strLoggingLevel, strInvert, sSelectionMethods ):
  def funcConfusionMetaMatrixRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, sSelectionMethods=sSelectionMethods ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strActualFile, lsSelectionFiles = astrSs[0], astrSs[1], astrSs[2:]
    return sfle.ex([strProg, strLoggingLevel, strInvert, strActualFile, strT, sSelectionMethods]+lsSelectionFiles)
  return funcConfusionMetaMatrixRet

###################################### General methods
#Make movies from a list of images
def funcMakeMovie( lsImageFiles ):
  def funcMakeMovieRet( target, source, env, lsImageFiles=lsImageFiles ):
    strMovie, astrSs = sfle.ts( target, source )
    strProg = astrSs[0]
    return sfle.ex([strProg,"mf://"+",".join([s.get_abspath() for s in lsImageFiles]),'-mf',"type="+lsImageFiles[0].get_abspath()[-3:]+":w=800:h=600:fps=2",'-ovc','lavc','-lavcopts','vcodec=mpeg4','-oac','copy','-o',strMovie])
  return funcMakeMovieRet

def globFilesByExtension(strDirectory,strExtension,strDelim = c_strExtDelim):
  #Read in all input files in the input directory
  lFileInputFiles = Glob( sfle.d( strDirectory, "*" ) )

  #Make sure they are the correct input files for the extension
  lCleanedFiles = []
  for fileName in lFileInputFiles:
    strPathPieces = [filter(None,strPathPiece) for strPathPiece in (fileName.get_abspath().split(c_strExtDelim))]
    if c_strExtDelim+strPathPieces[-1] == strExtension:
      lCleanedFiles.append(fileName)
  return lCleanedFiles

def funcReadConfigFile( strConfigFile ):
  #Dictionary to hold config data
  dictFileConfiguration = dict()

  #Read in contents into a list of lines to read
  sFileContents = list()
  with open(strConfigFile) as f:
    sFileContents = f.read()
    f.close()
  sFileContents = filter(None,re.split("\n",sFileContents))

  #Parse config file
  iIndex = 0
  iLengthContents = len(sFileContents)
  #All lines
  while iIndex < iLengthContents:
    sLine = sFileContents[iIndex]
    #Header lines
    if(sLine[0] == c_strConfigFileHeaderChar):
      strConfigKey = sLine
      iIndex = iIndex + 1
      if(iIndex < iLengthContents):
        sLine = sFileContents[iIndex]
        lsConfigData = list()
        #Contents of headers which are not headers
        while((iIndex < iLengthContents) and (not sLine[0] == c_strConfigFileHeaderChar)):
          #Ignore comments
          if(not sLine[0] == c_strConfigFileCommentChar):
            lsConfigData.append(sLine)
          iIndex = iIndex + 1
          if iIndex < iLengthContents:
            sLine = sFileContents[iIndex]
        if(len(lsConfigData) > 0):
          dictFileConfiguration[strConfigKey]=lsConfigData[0]
    else:
      iIndex = iIndex + 1
  return dictFileConfiguration

###Start process
#Generate insilico data sets
#Remember when creating insilico data sets, the Unbalanced test set's labels are off given built in randomness.
#When you remake this data set you will have to adjust the labels
#lDataFiles = globFilesByExtension(strDirectory=fileDirInput,strExtension=c_strSufInsilicoData)
lDataFiles = Glob( sfle.d( fileDirInput, "".join(["*",c_strSufInsilicoData]) ) )
lsInsilicoDataFiles = []

if c_fGenerateInSilicoDataSets:
  for strInsilicoData in lsInsilicoDataNames:
    if not strInsilicoData in lDataFiles:
      if strInsilicoData == c_strInsilicoDataDiversity:
        lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,strInsilicoData )).get_abspath())
        lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )).get_abspath())
        Command(File( sfle.d( fileDirInput,strInsilicoData )), [c_fileProgUtilityData], funcGenerateDiversityData(strDataSetKey=c_strInsilicoDataDiversityKey))
        Command(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )),[c_progHCL,File( sfle.d( fileDirInput,strInsilicoData ))],funcVisualizeInsilicoData(strXStart="1", strYStart="2"))
      elif strInsilicoData == c_strInsilicoDataUnbalanced:
        lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,strInsilicoData )).get_abspath())
        lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )).get_abspath())
        Command([File( sfle.d( fileDirInput,strInsilicoData )),File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufActualData) ))], [c_fileProgUtilityData], funcGenerateUnbalancedData(strDataSetKey=c_strInsilicoDataUnbalancedKey))
        Command(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )),[c_progHCL,File( sfle.d( fileDirInput,strInsilicoData ))],funcVisualizeInsilicoData(strXStart="1", strYStart="3"))

#Get micropita config files
#lMicropitaFiles = globFilesByExtension(strDirectory=fileDirInput,strExtension=c_strSufConfig)
lMicropitaFiles = Glob( sfle.d( fileDirInput, "".join(["*",c_strSufConfig]) ) )

#Collect input files later used for multistudy summary graphics
#{"InputFile:[SelectionFile1],[SelectionFile2],[SelectionFile3]...]}
dictSelectionFiles = dict()

#Read in each config file
#An use their contents to perform analysis
for fileConfigMicropita in lMicropitaFiles:
#  print("fileConfigMicropita")
#  print(fileConfigMicropita)
  #Get Contents of config file
  sFileConfiguration = funcReadConfigFile(fileConfigMicropita.get_abspath())
  #Set defaults
  if not c_strConfigSelectedTaxa in sFileConfiguration:
    sFileConfiguration[c_strConfigSelectedTaxa]="None"

  #Parse the supervised runs selection into multiple runs if needed
  if(c_strConfigSupervisedCount in sFileConfiguration):
    sSupervisedCounts = sFileConfiguration[c_strConfigSupervisedCount]
    sFileConfiguration[c_strConfigSupervisedCount] = filter(None,re.split(",",sSupervisedCounts))

  #Parse the unsupervised runs selection into multipley runs if needed
  if(c_strConfigUnsupervisedCount in sFileConfiguration):
    sSupervisedCounts = sFileConfiguration[c_strConfigUnsupervisedCount]
    sFileConfiguration[c_strConfigUnsupervisedCount] = filter(None,re.split(",",sSupervisedCounts))

  #Common configurations
  lsSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigSelectionTechniques]))
  lsPlotCombinedSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigPlotCombinedSelectionTechniques]))
  lsPlotCombinedSupervisedSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigPlotCombinedSupervisedSelectionTechniques]))
  strInvertImage = sFileConfiguration[c_strConfigInvertImage]
  if strInvertImage.lower()=="true":
    c_fileCladogramStyleFile = File(sfle.d( fileDirInput, c_strPathDelim.join([strDataFolder,"microPITA"])+c_strSufCircStyleInv ))

  #Indicator of if this is a supervised or unsupervised (or both run)
  fUnsupervisedRun = (len(c_lsUnsupervisedMethods.intersection(set(lsSelectionMethods))) > 0)
  fSupervisedRun = (len(c_lsSupervisedMethods.intersection(set(lsSelectionMethods))) > 0)
  fIsStratified = False
  if c_strConfigUnsupervisedStratify in sFileConfiguration:
    fIsStratified = not sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none"

  #Check to make sure the stratified/unsupervised count matches the supervised count
  if(fUnsupervisedRun and fSupervisedRun):
    #If stratified selection is not enabled for unsupervised
    if(not len(sFileConfiguration[c_strConfigUnsupervisedCount]) == len(sFileConfiguration[c_strConfigSupervisedCount])):
      print("Unsupervised selection count length was not equal to the supervised selection count. Did not run study.")
      print("=".join(["Config file",fileConfigMicropita.get_abspath()]))
      print("".join(["Unsupervised selection =",str(len(sFileConfiguration[c_strConfigUnsupervisedCount])),". Supervised selection =",str(len(sFileConfiguration[c_strConfigSupervisedCount]))]))
      pass

  #Update flag variables to pass
  fIsNormalized = sFileConfiguration[c_strConfigFileIsNormalized].lower() == "true"
  fIsSummed = sFileConfiguration[c_strConfigFileIsSummed].lower() == "true"
  fSumData = sFileConfiguration[c_strConfigSumData].lower() == "true"

  #Update the configurations that can be toggled on or off
  if sFileConfiguration[c_strConfigCladeFilter].lower() == "false":
    sFileConfiguration[c_strConfigCladeFilterLevel] = "none"
    sFileConfiguration[c_strConfigCladeFilterMeasure] = "none"
    sFileConfiguration[c_strConfigCladeFilterMinSize] = "none"

  if sFileConfiguration[c_strConfigAbundanceFilter].lower() == "false":
    sFileConfiguration[c_strConfigAbundanceFilterPercentile] = "none"
    sFileConfiguration[c_strConfigAbundanceFilterPercent] = "none"

  if sFileConfiguration[c_strOccurenceFilter].lower() == "false":
    sFileConfiguration[c_strOccurenceFilterMinSequence] = "none"
    sFileConfiguration[c_strOccurenceFilterMinSample] = "none"

  #Get supervised and unsupervised counts
  lsUnsupervisedCounts = sFileConfiguration[c_strConfigUnsupervisedCount]
  lsSupervisedCounts = sFileConfiguration[c_strConfigSupervisedCount]
  #Indicate if the unsupervised or supervised index are going to be the drivers of the loop
  #This is done because there may not be supervised or unsupervised selection performed
  lsIndexCounts = lsUnsupervisedCounts
  if(not fUnsupervisedRun):
    lsIndexCounts = lsSupervisedCounts

  #Current working directory
  sFileBase = os.path.basename(fileConfigMicropita.get_abspath()).split(c_strExtDelim)[0]
  sOutputDir = sfle.d(fileDirOutput,sFileBase)

  #Holds figures to make into movies
  lsCladogramFigures = []
  lsPCOAFigures = []
  lsStratPCOAFigures = []

  #If there is a validation file
  fileValidationFile = None
  fileCheckedValidationFile = None

  #Check the validation file
  if(not sFileConfiguration[c_strConfigValidationFile].lower() == "none"):
    fileValidationFile = File(sfle.d(fileDirInput.get_abspath(),sFileConfiguration[c_strConfigValidationFile]))
    fileCheckedValidationFile = File(sfle.d(fileDirDataName,sfle.rebase(sFileConfiguration[c_strConfigValidationFile], c_strSufUncheckedTable, c_strSufCheckedTable)))
    Command(fileCheckedValidationFile,[c_fileProgCheckFile, fileValidationFile], 
            funcCheckAbundanceData( " ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                    " ".join([Constants_Arguments.c_strLastMetadataNameArgument, sFileConfiguration[c_strConfigValidationLastMetadataName]])))

  #Loop through the count selection using the unsupervised counts length for indexing
  for iCountIndex in xrange(0,len(lsIndexCounts)):

    sAbundanceFileName = sFileConfiguration[c_strConfigInputFile]
    sCheckedAbundanceFileName = sfle.d(fileDirDataName, sfle.rebase(sAbundanceFileName, c_strSufUncheckedTable, c_strSufCheckedTable))
    sCheckedAbundanceFile = File(sCheckedAbundanceFileName)

    #Build the micropita output name
    sPrefix = "".join(["S",lsSupervisedCounts[iCountIndex],"-"])
    if fUnsupervisedRun:
      if fIsStratified:
        sPrefix = "".join(["Strat",lsUnsupervisedCounts[iCountIndex],"-"])
      else:
        sPrefix = "".join(["U",lsUnsupervisedCounts[iCountIndex],"-"])

    sMicropitaOutput = sfle.rebase(sAbundanceFileName, c_strSufUncheckedTable, c_strSufMicropita)
    sMicropitaOutput = "".join([sPrefix,sMicropitaOutput])
    sMicropitaOutput = File(sfle.d(fileDirOutput,c_strPathDelim.join([sFileBase,sMicropitaOutput])))

    #Make file names from the input micropita file in the temp directory
    sMicropitaPredictFile = File(sfle.d( fileDirTmp, sfle.rebase( sCheckedAbundanceFile, c_strSufTable, c_strSufPredict )))

    #Make files names from the micropita output file
    sOutputFigure1APCoA,sOutputCombinedFigure1APCoA,sOutputFigure4PCoA,sOutputFigure4CombinedPCoA,sOutputFigure1BHCL,sOutputFigure1BConfusion,sOutputFigure1BOverlap,sHCLSelectData,sHCLSelectColor,sHCLSelectLabel,sCladogramSelectFig,sCladogramTaxaFile,\
    sCladogramColorFile,sCladogramSizeFile,sCladogramTickFile,sCladogramHighlightFile,sCladogramCircleFile, sCladogramDetailFile =  [File(sfle.d( sOutputDir,
    sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufPCOA,c_strSufCombinedPCOA,c_strSufStratPCOA,c_strSufCombinedStratPCOA,c_strSufHCLUSTFig,c_strSufConfusionMatrix,c_strSufOverlapMatrix,c_strSufHCLUSTData,c_strSufHCLUSTColor,
    c_strSufHCLUSTLabel,c_strSufFig2,c_strSufCircTaxa,c_strSufCircColor,c_strSufCircSize,c_strSufCircTick,c_strSufCircHighlight, c_strSufCircCircle, c_strSufCircDetail)]

    #Make more file names, these for cladogram options
    cCladogramSelectedTaxa = None
    if c_strConfigSelectedTaxa in sFileConfiguration:
      if(not sFileConfiguration[c_strConfigSelectedTaxa].lower() == "none"):
        cCladogramSelectedTaxa = File(c_strPathDelim.join(["input",sFileConfiguration[c_strConfigSelectedTaxa]]))

    cCladogramHighlightedTaxa = None
    if c_strConfigHighlightClades in sFileConfiguration:
      if(not sFileConfiguration[c_strConfigHighlightClades].lower() == "none"):
        cCladogramHighlightedTaxa = File(c_strPathDelim.join(["input",sFileConfiguration[c_strConfigHighlightClades]]))

    #Make more files, now for the figure 3 HCL
    sOutputFigure3HCL, sStratHCLSelectColor, sStratHCLSelectLabel =  [File(sfle.d( sOutputDir,
    sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufStratHCLUSTFig,c_strSufStratHCLUSTColor,c_strSufStratHCLUSTLabel)]

    #Run micropita analysis
    #Note, micropita will check the input file so no need to send in the checked abundance file,
    #Send in the raw abundance file and micropita will make a checked version to use.
    #The lsInsilicoDataFiles is added to this command to order the actions correctly.
    #The insilico file creation must go before micropita
    sAbundanceFileName = File(c_strPathDelim.join(["input",sAbundanceFileName]))
    lFileMicropitaConditionalDependencies = []+lsInsilicoDataFiles
    if os.path.exists(sCheckedAbundanceFile.get_abspath()):
      lFileMicropitaConditionalDependencies.append(sCheckedAbundanceFile.get_abspath())
    Command(sMicropitaOutput, [c_fileProgMicroPITA, sAbundanceFileName] + c_filesSecondarySrc + [fileConfigMicropita]+lFileMicropitaConditionalDependencies,
        funcMicroPita(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
        " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
        " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
        fIsNormalized, fIsSummed, fSumData,
        " ".join([Constants_Arguments.c_strTargetedFeatureMethodArgument,sFileConfiguration[c_strConfigTargetedSelection]]),
        " ".join([Constants_Arguments.c_strUnsupervisedCountArgument,sFileConfiguration[c_strConfigUnsupervisedCount][iCountIndex]]),
        " ".join([Constants_Arguments.c_strTargetedSelectionFileArgument,cCladogramSelectedTaxa.get_abspath()]),
        " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadataArgument,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
        " ".join([Constants_Arguments.c_strSupervisedLabelArgument,sFileConfiguration[c_strConfigSupervisedLabel]]),
        " ".join([Constants_Arguments.c_strSupervisedLabelCountArgument,sFileConfiguration[c_strConfigSupervisedCount][iCountIndex]]),
        " ".join([Constants_Arguments.c_strTemporaryDirectoryArgument,fileDirTmp.get_abspath()]),
        lsSelectionMethods))

    #Create figure 1A PCoA
    #Unstratified PCoA of selection
    #Manage optional files
    if fSupervisedRun:
      strPredictFileArgument = " ".join(["-p",sMicropitaPredictFile.get_abspath()])
    else:
      strPredictFileArgument = " ".join(["-p","None"])
    if sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none":
      lsPCOAFigures.append(sOutputFigure1APCoA)
      Command(sOutputFigure1APCoA, [c_fileProgPCoAFigure, sCheckedAbundanceFile, sMicropitaOutput] + c_filePrimarySrc, funcPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strSumDataArgument,sFileConfiguration[c_strConfigSumData]]),
                                                                                                                                                                strPredictFileArgument))
    #Create figure 1A PCoA combined
    #Unstratified PCoA of selection
    #Manage optional files
    if sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none":
      Command(sOutputCombinedFigure1APCoA, [c_fileProgCombinedPCoAFigure, sCheckedAbundanceFile, sMicropitaOutput] + c_filePrimarySrc, funcCombinedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strSumDataArgument,sFileConfiguration[c_strConfigSumData]]),
                                                                                                                                                                 lsPlotCombinedSelectionMethods))

    #Create figure 1B HCL
    #Selection of samples by selection method
#    Command([sOutputFigure1BHCL,sHCLSelectData, sHCLSelectColor, sHCLSelectLabel], [c_fileProgSelectionHCLFigure, sMicropitaOutput] + ls_srcFig1, 
#        funcHCLSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]), " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage])))

    #Create alt figure 1B overlapp matrix
    Command(sOutputFigure1BOverlap, [c_fileProgOverlapMatrixFigure, sMicropitaOutput] + ls_srcFig1, 
        funcOverlapMatrix(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]), " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),lsPlotCombinedSelectionMethods+lsPlotCombinedSupervisedSelectionMethods))

    #Create alt figure 1B Confusion matrix
    #Do not generate the confusion on data without a file defining the actual classes.
    if not sFileConfiguration[c_strConfigActualFile].lower() == "none":

        #Create alt figure 1B confusion matrix
        #Selection of samples by selection method
        Command(sOutputFigure1BConfusion, [c_fileProgConfusionMatrixFigure, sMicropitaOutput, File( sfle.d( fileDirInput,sFileConfiguration[c_strConfigActualFile] ))] + ls_srcFig1, 
            funcConfusionMatrix(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                lsPlotCombinedSelectionMethods+lsPlotCombinedSupervisedSelectionMethods))

    #Create a cladogram figure 2
    #Manage files which are optional
    strTaxaArgForCladogram = "-t"
    strHighlightArgForCladogram = "-c"
    lsFig2Sources = [c_fileProgSelectionCladogramFigure, sMicropitaOutput, sCheckedAbundanceFile, c_fileCladogramStyleFile]
    if not cCladogramSelectedTaxa == None:
      lsFig2Sources.append(cCladogramSelectedTaxa)
      strTaxaArgForCladogram = " ".join([strTaxaArgForCladogram,cCladogramSelectedTaxa.get_abspath()])
    else:
      strTaxaArgForCladogram = " ".join([strTaxaArgForCladogram,"None"])
    if not cCladogramHighlightedTaxa == None:
      lsFig2Sources.append(cCladogramHighlightedTaxa)
      strHighlightArgForCladogram = " ".join([strHighlightArgForCladogram,cCladogramHighlightedTaxa.get_abspath()])
    else:
      strHighlightArgForCladogram = " ".join([strHighlightArgForCladogram,"None"])
    #Taxonomic Enrichment
    lsCladogramFigures.append(sCladogramSelectFig)

#    print("sFileConfiguration")
#    print(sFileConfiguration)

    Command([sCladogramSelectFig, sCladogramTaxaFile, sCladogramColorFile, sCladogramTickFile, sCladogramHighlightFile,
         sCladogramSizeFile, sCladogramCircleFile, sCladogramDetailFile], lsFig2Sources + ls_srcFig2, 
         funcCladogramSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                       strTaxaArgForCladogram,
                                       strHighlightArgForCladogram,
                                       " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                       " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                       " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                       " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                       " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                       " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                       " ".join([Constants_Arguments.c_strRootArgument,sFileConfiguration[c_strConfigRoot]]),
                                       " ".join([Constants_Arguments.c_strEnrichmentMethodArgument,sFileConfiguration[c_strConfigEnrichmentMeasurement]]),
                                       " ".join([Constants_Arguments.c_strCladeFilterLevelArgument,sFileConfiguration[c_strConfigCladeFilterLevel]]),
                                       " ".join([Constants_Arguments.c_strCladeMeasureLevelArgument,sFileConfiguration[c_strConfigCladeFilterMeasure]]),
                                       " ".join([Constants_Arguments.c_strCladeFilteringMinLevelArgument,sFileConfiguration[c_strConfigCladeFilterMinSize]]),
                                       " ".join([Constants_Arguments.c_strAbundanceFilterPercentileArgument,sFileConfiguration[c_strConfigAbundanceFilterPercentile]]),
                                       " ".join([Constants_Arguments.c_strAbundanceFilterCutoffArgument,sFileConfiguration[c_strConfigAbundanceFilterPercent]]),
                                       " ".join([Constants_Arguments.c_strRingOrderArgument,sFileConfiguration[c_strConfigCladogramRingOrder]]),
                                       " ".join([Constants_Arguments.c_strCircladerTicksArgument,sFileConfiguration[c_strConfigCladogramTicks]]),
                                       " ".join([Constants_Arguments.c_strEnrichmentThresholdArgument,sFileConfiguration[c_strConfigCladogramAlpha]]),
                                       " ".join([Constants_Arguments.c_strOccurenceFilterSequenceCountArgument,sFileConfiguration[c_strOccurenceFilterMinSequence]]),
                                       " ".join([Constants_Arguments.c_strOccurenceFilterSampleCountArgument,sFileConfiguration[c_strOccurenceFilterMinSample]]),
                                       " ".join([Constants_Arguments.c_strSumDataArgument,sFileConfiguration[c_strConfigSumData]]),
                                       " ".join([Constants_Arguments.c_strTerminalLevelArgument,sFileConfiguration[c_strConfigTerminalLevel]])))

    #Create figure 3 HCL version
#    #Selection in Stratification
#    Command([sOutputFigure3HCL, sStratHCLSelectColor, sStratHCLSelectLabel], [c_fileProgStratSelectionHCLFigure, sMicropitaOutput, sCheckedAbundanceFile] + ls_srcFig1, 
#        funcHCLStratSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
#                                " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
#                                " ".join(["-id", c_strAbundanceIDCol]),
#                                " ".join([Constants_Arguments.c_strInvertArgument, strInvertImage])))

    #Create figure 4 PCoA
    #PCoA of stratified selection
    if not sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none":
      lsStratPCOAFigures.append(sOutputFigure4PCoA)
      Command(sOutputFigure4PCoA, [c_fileProgStratifiedPCoAFigure, sCheckedAbundanceFile, sMicropitaOutput] + c_filePrimarySrc, funcStratifiedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadataArgument,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strSumDataArgument,sFileConfiguration[c_strConfigSumData]]),
                                                                                                                                                   strPredictFileArgument))

      Command(sOutputFigure4CombinedPCoA, [c_fileProgCombinedStratifiedPCoAFigure, sCheckedAbundanceFile, sMicropitaOutput] + c_filePrimarySrc, funcCombinedStratifiedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadataArgument,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strSumDataArgument,sFileConfiguration[c_strConfigSumData]]),
                                                                                                                                                   lsPlotCombinedSelectionMethods+lsPlotCombinedSupervisedSelectionMethods))

    #Validate data
    #Validate diversityfuncMeasureDiversityByGroup
    if fileCheckedValidationFile:
      if c_DIVERSITY in lsSelectionMethods:
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedDiversity)]) )),
              [c_fileProgValidateDiversity, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput],
              funcMeasureDiversityByGroup(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]])))

      if c_TARGETED_TAXA in lsSelectionMethods:
        #Validate Feature
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedFeature)]) )),
              [c_fileProgValidateFeature, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput, cCladogramSelectedTaxa],
              funcMeasureFeatureByGroup(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]]),
                                      " ".join([Constants_Arguments.c_strTargetedFeatureMethodArgument,sFileConfiguration[c_strConfigTargetedSelection]])))

        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedFeatureHistogram)]) )),
              [c_fileProgValidateFeatureHistogram, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput, cCladogramSelectedTaxa],
              funcMeasureFeatureByGroupHistogram(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]])))

      if c_EXTREME in lsSelectionMethods:
        #Validate Extreme selection
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedExtreme)]) )),
              [c_fileProgValidatePCoA, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput],
              funcValidateInPCoA(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]]),
                                      " ".join([Constants_Arguments.c_strMetricArgument,c_EXTREME_DISSIMILARITY_1]),
                                      " ".join([Constants_Arguments.c_strSupervisedLabelArgument,"None"])))

      if c_REPRESENTATIVE in lsSelectionMethods:
        #Validate Representative selection
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedRepresentative)]) )),
              [c_fileProgValidatePCoA, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput],
              funcValidateInPCoA(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]]),
                                      " ".join([Constants_Arguments.c_strMetricArgument,c_REPRESENTATIVE_DISSIMILARITY_1]),
                                      " ".join([Constants_Arguments.c_strSupervisedLabelArgument,"None"])))

      if c_DISTINCT in lsSelectionMethods:
        #Validate Distinct selection
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedDistinct)]) )),
              [c_fileProgValidatePCoA, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput],
              funcValidateInPCoA(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]]),
                                      " ".join([Constants_Arguments.c_strMetricArgument,c_DISTINCT]),
                                      " ".join([Constants_Arguments.c_strSupervisedLabelArgument,sFileConfiguration[c_strConfigSupervisedLabel]])))

      if c_DISCRIMINANT in lsSelectionMethods:
        #Validate Discriminant selection
        Command(File(sfle.d( sOutputDir,"".join([sPrefix,sfle.rebase(fileCheckedValidationFile.get_abspath(), c_strSufCheckedTable, c_strSufValidatedDiscriminant)]) )),
              [c_fileProgValidatePCoA, fileCheckedValidationFile, sCheckedAbundanceFile, sMicropitaOutput],
              funcValidateInPCoA(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                      " ".join([Constants_Arguments.c_strIDNameArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                      " ".join([Constants_Arguments.c_strLastMetadataNameArgument,sFileConfiguration[c_strConfigLastMetadataRow]]),
                                      " ".join([Constants_Arguments.c_strValidationIDNameArgument,sFileConfiguration[c_strConfigValidationIDName]]),
                                      " ".join([Constants_Arguments.c_strValidationLastMetadataNameArgument,sFileConfiguration[c_strConfigValidationLastMetadataName]]),
                                      " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                      " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                                      " ".join([Constants_Arguments.c_strValidationIsNormalizedArgument,sFileConfiguration[c_strConfigValidationIsNormalized]]),
                                      " ".join([Constants_Arguments.c_strValidationIsSummedArgument,sFileConfiguration[c_strConfigValidationIsSummed]]),
                                      " ".join([Constants_Arguments.c_strPairingMetadataArgument,sFileConfiguration[c_strConfigPairingMetadata]]),
                                      " ".join([Constants_Arguments.c_strMetricArgument,c_DISCRIMINANT]),
                                      " ".join([Constants_Arguments.c_strSupervisedLabelArgument,sFileConfiguration[c_strConfigSupervisedLabel]])))

    #Add input files to be later summarized in the metaplots
    #Create key combining abundance file and stratification status
    sKey = "".join(["STRAT-",sFileConfiguration[c_strConfigUnsupervisedStratify],"-INPUT-",os.path.basename(sCheckedAbundanceFile.get_abspath())])
    if sKey not in dictSelectionFiles:
      dictSelectionFiles[sKey] = dict()
      curSummaryFileDict = dictSelectionFiles[sKey]
      curSummaryFileDict[c_strSelectionFiles] = list()

    dictSelectionFiles[sKey][c_strConfigInputFile] = sCheckedAbundanceFile
    dictSelectionFiles[sKey][c_strSelectionFiles].append(sMicropitaOutput)

    dictSelectionFiles[sKey][c_strConfigSelectionTechniquesCollectorCurve] = sFileConfiguration[c_strConfigSelectionTechniquesCollectorCurve]
    dictSelectionFiles[sKey][c_strConfigLogging] = sFileConfiguration[c_strConfigLogging]
    dictSelectionFiles[sKey][c_strConfigSampleRow] = sFileConfiguration[c_strConfigSampleRow]
    dictSelectionFiles[sKey][c_strConfigLastMetadataRow] = sFileConfiguration[c_strConfigLastMetadataRow]
    dictSelectionFiles[sKey][c_strConfigInvertImage] = sFileConfiguration[c_strConfigInvertImage]
    dictSelectionFiles[sKey][c_strConfigFileIsNormalized] = sFileConfiguration[c_strConfigFileIsNormalized]
    dictSelectionFiles[sKey][c_strConfigFileIsSummed] = sFileConfiguration[c_strConfigFileIsSummed]

  #If this configuration file had multiple sampling
  if fMakeMovies:
    if lsIndexCounts > 1:
#      if len(lsPCOAFigures) > 1:
#        Command(sOutputDir+c_strPCOAMovieEnding,[c_fileProgMencoder]+lsPCOAFigures, funcMakeMovie(lsImageFiles=lsPCOAFigures))
#      if len(lsStratPCOAFigures) > 1:
#        Command(sOutputDir+c_strStratPCOAMovieEnding,[c_fileProgMencoder]+lsStratPCOAFigures, funcMakeMovie(lsImageFiles=lsStratPCOAFigures))
      if len(lsCladogramFigures) > 1:
        Command(sOutputDir+c_strCladogramMovieEnding,[c_fileProgMencoder]+lsCladogramFigures, funcMakeMovie(lsImageFiles=lsCladogramFigures))



if c_fRunCollectionCurve:
  #For each input file
  #Create summary files
  #Create collection curves
  for strInputSummaryKey in dictSelectionFiles:
    curSummaryDict = dictSelectionFiles[strInputSummaryKey]

    if (curSummaryDict[c_strConfigFileIsNormalized].lower()=="false") and (curSummaryDict[c_strConfigFileIsSummed].lower()=="false"):

      curInputFile = curSummaryDict[c_strConfigInputFile]
      curListofSelection = curSummaryDict[c_strSelectionFiles]
      curLogging = curSummaryDict[c_strConfigLogging]
      curSampleNameRow = curSummaryDict[c_strConfigSampleRow]
      curLastMetadataRow = curSummaryDict[c_strConfigLastMetadataRow]
      curInvert = curSummaryDict[c_strConfigInvertImage]

      lsPlotCollectorSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigSelectionTechniquesCollectorCurve]))

      #Make more files, now for figure 5 Collection Curve
      sOutputFigure5CC =  File(sfle.d( fileDirOutput.get_abspath()+strOutputSummaryFolder, sfle.rebase( strInputSummaryKey, c_strSufTable, c_strSufCollectionCurveFigure)))

      #Create Figure 5
      #Collection Curve
      Command([sOutputFigure5CC], [c_fileProgCollectionCurveFigure, curInputFile] + curListofSelection + ls_srcFig1, 
          funcCollectionCurveSummary(" ".join([Constants_Arguments.c_strLoggingArgument, curLogging]),
                              " ".join([Constants_Arguments.c_strIDNameArgument, curSampleNameRow]),
                              " ".join([Constants_Arguments.c_strLastMetadataNameArgument, curLastMetadataRow]),
                              " ".join([Constants_Arguments.c_strInvertArgument, curInvert]),
                              " ".join([Constants_Arguments.c_strIsNormalizedArgument,sFileConfiguration[c_strConfigFileIsNormalized]]),
                              " ".join([Constants_Arguments.c_strIsSummedArgument,sFileConfiguration[c_strConfigFileIsSummed]]),
                              " ".join([Constants_Arguments.c_strPlotSelectedArgument]+lsPlotCollectorSelectionMethods)))

#Generate the Confusion and Overlap matrices based on multiple projects
lMetamatrixConfigFiles = Glob( sfle.d( fileDirInput, "".join(["*",c_strSufMetaMatrix]) ) )
for fileMetamatrixConfig in lMetamatrixConfigFiles:
  sMetaFileName = fileMetamatrixConfig.get_abspath()
  dictMetadetails = funcReadConfigFile(sMetaFileName)
  
  #Actual file
  strActualFile = File(sfle.d(fileDirInput,dictMetadetails[c_strConfigActualFile]))
  #Invert
  strInvert = dictMetadetails[c_strConfigInvertImage]
  #Logging
  strLogging = dictMetadetails[c_strConfigLogging]
  #Selection techniques
  strSelectionTechniques = dictMetadetails[c_strConfigSelectionTechniques]
  #Projects of whose output will be combined in the plot
  strProjects = dictMetadetails[c_strConfigProjects]
  #Projects as a list of files
  sProjectFiles = [File(sfle.d(fileDirOutput,sFile)) for sFile in strProjects.split(c_strComma)]

  #Create output file name for overlap matrix
  sOutputFigure1BOverlapMeta = File(sfle.d(fileDirOutput.get_abspath()+strOutputSummaryFolder,"".join([os.path.splitext(os.path.split(sMetaFileName)[1])[0],c_sufMetaOverlapMatrix])))
  sOutputFigure1BConfusionMeta = File(sfle.d(fileDirOutput.get_abspath()+strOutputSummaryFolder,"".join([os.path.splitext(os.path.split(sMetaFileName)[1])[0],c_sufMetaConfusionMatrix])))

  #TODO add in import dependencies from inside script
  #Create Figure 1B- Overlap matrix
  Command(sOutputFigure1BOverlapMeta, [c_fileProgOverlapMatrixMetaFigure] + sProjectFiles, 
    funcOverlapMetaMatrix(" ".join([Constants_Arguments.c_strLoggingArgument, strLogging]),
                          " ".join([Constants_Arguments.c_strInvertArgument, strInvert]),
                          strSelectionTechniques))

  #TODO add in import dependencies from inside script
  #Create Figure 1B- Confusion matrix
  Command(sOutputFigure1BConfusionMeta, [c_fileProgConfusionMatrixMetaFigure, strActualFile] + sProjectFiles, 
    funcConfusionMetaMatrix(" ".join([Constants_Arguments.c_strLoggingArgument, strLogging]),
                            " ".join([Constants_Arguments.c_strInvertArgument, strInvert]),
                            strSelectionTechniques))

