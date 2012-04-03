import re
import sfle
import sys
import os

Import( "*" )

#TODO make this dynamic
sys.path.append("/home/ttickle/Desktop/ttickle/sfle/input/micropita/src")
from Constants_Arguments import Constants_Arguments

c_strDataFile = sfle.d(fileDirInput, "data")

pE = DefaultEnvironment( )

#Process flags
c_fRunCollectionCurve = False

#Normalize
c_NormalizePCOA = "False"

#Extentions
c_strPCOAMovieEnding = "/PCOA.avi"
c_strStratPCOAMovieEnding = "/Strat-PCOA.avi"
c_strCladogramMovieEnding = "/Cladogram.avi"
c_strSufCollectionCurveFigure = "-ColCurv.png"
c_strSufCollectionCurveText = "-ColCurv.txt"
c_strSufCombinedPCOA = "-Combined.png"
c_strSufSelectedTaxa = ".taxa"
c_strSufPng = ".png"
c_strSufPCOA = "-PCoA.png"
c_strSufCombinedStratPCOA = "-Combined-StratPCoA.png"
c_strSufStratPCOA = "-StratPCoA.png"
c_strSufHCLUSTData = ".HCLData"
c_strSufStratHCLUSTColor = "-Strat.HCLColor"
c_strSufHCLUSTColor = ".HCLColor"
c_strSufStratHCLUSTLabel = "-Strat.HCLLabel"
c_strSufHCLUSTLabel = ".HCLLabel"
c_strSufStratHCLUSTFig = "-StratHCL.png"
c_strSufHCLUSTFig = "-SHCL.png"
c_strSufInsilicoData = ".txt"
c_strSufFig2 = "-Fig2.png"
c_strSufCircTaxa = ".CTaxa"
c_strSufCircColor = ".CColor"
c_strSufCircSize = ".CSize"
c_strSufCircStyle = ".CStyle"
c_strSufCircTick = ".CTick"
c_strSufCircHighlight = ".CHLight"
c_strSufCircCircle = ".CCircle"
c_strSufConfig = ".config"
c_strSufCircStyleInv = "-invert.CStyle"
c_strSufMicropita = ".txt"
c_strSufPredict = "-SVM.predict"

#Data files to generate
c_strInsilicoDataDiversity = "DiversityTest"+c_strSufInsilicoData
c_strInsilicoDataDiversityKey = "Diversity"
c_strInsilicoDataRepresentative = "RepresentativeTest"+c_strSufInsilicoData
c_strInsilicoDataUnbalanced = "UnbalancedTest"+c_strSufInsilicoData
c_strInsilicoDataUnbalancedKey = "Unbalanced"
lsInsilicoDataNames = [c_strInsilicoDataDiversity,c_strInsilicoDataRepresentative,c_strInsilicoDataUnbalanced]

#Dict keys
c_strSelectionFiles = "SelectionFiles"

#Special characters
c_strConfigFileHeaderChar = "["
c_strConfigFileCommentChar = "#"
c_strExtDelim = "."
c_strPathDelim = "/"

#Directories
strOutputSummaryFolder = "".join([c_strPathDelim,"Metasummary"])

#The Column that holds the TAXA / OTU IDS
c_strAbundanceIDCol = "0"

#Scons Configuration file headers
c_strConfigAbundanceFilter = "[Abundance Filter]"
c_strConfigAbundanceFilterPercentile = "[Abundance Filter Percentile]"
c_strConfigAbundanceFilterPercent = "[Abundance Filter Percent Above Percentile]"
c_strConfigCladeFilter = "[Clade Filter]"
c_strConfigCladeFilterMeasure = "[Clade Level to Measure]"
c_strConfigCladeFilterLevel = "[Clade Level to Filter]"
c_strConfigCladeFilterMinSize = "[Minimum Clade Size]"
c_strConfigCladogramRingOrder = "[Cladogram Ring Order]"
c_strConfigCladogramTicks = "[Cladogram Dendrogram Ticks]"
c_strConfigDataRow = "[Data Name Row Index]"
c_strConfigEnrichmentMeasurement = "[Taxa or OTU Enrichment Measurement]"
c_strConfigHighlightClades = "[Higlight Clades]"
c_strConfigInputFile = "[Input File]"
c_strConfigInvertImage = "[Invert Image]"
c_strConfigLabels = "[Labels]"
c_strConfigLogging = "[Logging]"
c_strConfigNormalizeAbundance = "[Normalize Abundance]"
c_strConfigPlotCombinedSelectionTechniques = "[Unsupervised Selection Methods to plot in a combined graph]"
c_strConfigPlotCombinedSupervisedSelectionTechniques = "[Supervised Selection Methods to plot in a combined graph]"
c_strConfigRoot = "[Root]"
c_strConfigSampleRow = "[Sample Name Row Index]"
c_strConfigSelection = "[Selection]"
c_strConfigSelectedTaxa = "[Selected Taxa]"
c_strConfigSelectionTechniques = "[Selection Techniques]"
c_strConfigSupervisedLabel = "[Supervised Label]"
c_strConfigSupervisedCount = "[Supervised Selection Count]"
c_strConfigUnsupervisedCount = "[Unsupervised Selection Count]"
c_strConfigUnsupervisedStratify = "[Stratify by Metadata]"

#None generated output micropita files
#c_fileIBDWGSSelectionManual = File(sfle.d(fileDirInput,"IBDWGSSelection.txt"))
#c_fileIBDWGSSelectionOutput = File(sfle.d(fileDirOutput,"IBDWGSSelection/IBDWGSSelection.txt"))

#Cladogram style file for micropita
c_fileCladogramStyleFile = File(sfle.d( fileDirInput, "microPITA"+c_strSufCircStyle ))

#External programs
c_progHCLPath = "./external/hclust/hclust.py"
c_progHCL = "../../external/hclust/hclust.py"

#Src Code
c_fileProgAbundanceTable = File( sfle.d( fileDirSrc, "AbundanceTable.py" ) )
c_fileProgCladogram = File( sfle.d( fileDirSrc, "Cladogram.py" ) )
c_fileProgCommandLine = File( sfle.d( fileDirSrc, "CommandLine.py" ) )
c_fileProgConstants = File( sfle.d( fileDirSrc, "Constants.py" ) )
c_fileProgConstantsFigures = File( sfle.d( fileDirSrc, "Constants_Figures.py" ) )
c_fileProgDiversity = File( sfle.d( fileDirSrc, "Diversity.py" ) )
c_fileProgFileIO = File( sfle.d( fileDirSrc, "FileIO.py" ) )
c_fileProgMencoder = File("/usr/bin/mencoder")
c_fileProgMicroPITA = File( sfle.d( fileDirSrc, "MicroPITA.py" ) )
c_fileProgCollectionCurveFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCollectionCurve.py") )
c_fileProgPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperPCoA.py" ) )
c_fileProgCombinedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCombinedPCoA.py") )
c_fileProgCombinedStratifiedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCombinedStratifiedPCoA.py" ) )
c_fileProgStratifiedPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperStratifiedPCoA.py" ) )
c_fileProgSelectionCladogramFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionCladogram.py") )
c_fileProgSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionHCL.py" ) )
c_fileProgStratSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperStratSelectionHCL.py" ) )
c_fileProgMLPYDistanceAdaptor = File( sfle.d( fileDirSrc, "MLPYDistanceAdaptor.py" ) )
c_fileProgPCOA = File( sfle.d( fileDirSrc, "PCoA.py" ) )
c_fileProgSVM = File( sfle.d( fileDirSrc, "SVM.py" ) )
c_fileProgUtilityData = File( sfle.d( fileDirSrc, "MicropitaPaperConstructDataSets.py" ) )
c_fileProgValidateData = File( sfle.d( fileDirSrc, "ValidateData.py" ) )

#Lists of sources needed for different python scripts that are ran so that they can be used as libraries
c_filesSecondarySrc = [c_fileProgAbundanceTable, c_fileProgCommandLine, c_fileProgConstants, c_fileProgDiversity,
	c_fileProgFileIO, c_fileProgMLPYDistanceAdaptor, c_fileProgPCOA, c_fileProgSVM, c_fileProgUtilityData]
c_filePrimarySrc = c_filesSecondarySrc + [c_fileProgMicroPITA]
ls_srcFig1 = [c_fileProgCommandLine, c_fileProgConstants, c_fileProgConstantsFigures, c_fileProgMicroPITA, c_fileProgValidateData, c_fileProgFileIO]
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
def funcGenerateData(strDataSetKey):
  def funcGenerateDataRet( target, source, env, strDataSetKey=strDataSetKey ):
    strT, astrSs = sfle.ts( target, source )
    strProg = astrSs[0]
    return sfle.ex( [strProg,strT,strDataSetKey] )
  return funcGenerateDataRet

#Visualize insilico data
def funcVisualizeInsilicoData( target, source, env ):
  strT, astrSs, = sfle.ts( target, source )
  strProg, strInputDataFile = astrSs[0], astrSs[1]
  return sfle.ex([strProg, "--in", strInputDataFile, "--out", strT, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1",'--flabel',"1",'--font_size',"5"])

##Call micropita to run
def funcMicroPita( strLoggingLevel, iUnsupervisedCount, iSampleNameRow, iFirstDataRow, iSupervisedCount, strSupervisedLabel, strUnsupervisedStratify, strTaxaFile, lsSelectionMethods):
  def funcMicroPitaRet( target, source, env, iUnsupervisedCount=iUnsupervisedCount, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, 
                        iSupervisedCount=iSupervisedCount, strSupervisedLabel=strSupervisedLabel, strUnsupervisedStratify=strUnsupervisedStratify, strTaxaFile=strTaxaFile, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd = astrSs[0], astrSs[1]
    aArgs = [iUnsupervisedCount]+lsSelectionMethods
    return sfle.ex( [strProg]+[strLoggingLevel, iSampleNameRow, iFirstDataRow, iSupervisedCount, strSupervisedLabel, strUnsupervisedStratify]+[strT, strAbnd, strTaxaFile, fileDirTmp] + aArgs)
  return funcMicroPitaRet

##Create figures
#Visualize output with PCoA (Figure 1A)
def funcPCoASelectionMethods( strLoggingLevel, iSampleNameRow, iFirstDataRow, iNormalize, strInvert, strPredictFileArgument ):
  def funcPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, iNormalize=iNormalize, strInvert=strInvert, strPredictFileArgument=strPredictFileArgument):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[strLoggingLevel, iSampleNameRow,iFirstDataRow,iNormalize,strInvert,strPredictFileArgument]+[ strAbnd, strSelection, strT])
  return funcPCoARet

#Visualize output with Combined PCoA (Figure 1A)
def funcCombinedPCoASelectionMethods( strLoggingLevel, iSampleNameRow, iFirstDataRow, iNormalize, strInvert, lsSelectionMethods):
  def funcCombinedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, iNormalize=iNormalize, strInvert=strInvert, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[strLoggingLevel,iSampleNameRow,iFirstDataRow,iNormalize,strInvert]+[ strAbnd, strSelection, strT]+lsSelectionMethods)
  return funcCombinedPCoARet

#Visualize selected output with HCL (Figure 1B)
def funcHCLSelectionMethods( strLoggingLevel, strInvert ):
  def funcHCLSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection = astrSs[0], astrSs[1]
    strData, strColor, strLabel = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath()
    return sfle.ex([strProg, strLoggingLevel, strInvert, strSelection, c_progHCLPath, strData, strColor, strLabel, strT])
  return funcHCLSelectionRet

#Visualize selected output with Cladogram (Figure 2)
def funcCladogramSelectionMethods( strLoggingLevel, strTargetedTaxaFile, strHighlightCladeFile, strInvert, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                                   strAbundanceFilterPercentile, strAbundanceFilterPercent, strRingOrder, strTicks):

  def funcCladogramSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strTargetedTaxaFile=strTargetedTaxaFile, strHighlightCladeFile=strHighlightCladeFile, strInvert=strInvert, strRoot=strRoot, strEnrichment=strEnrichment,
                                 strCladeFilterLevel=strCladeFilterLevel, strCladeFilterMeasure=strCladeFilterMeasure, strCladeFilterMin=strCladeFilterMin, strAbundanceFilterPercentile=strAbundanceFilterPercentile,
                                 strAbundanceFitlerPercent=strAbundanceFilterPercent, strRingOrder=strRingOrder, strTicks=strTicks):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection, strAbundance, strStyleFile = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
    strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath(), target[4].get_abspath(), target[5].get_abspath(), target[6].get_abspath()
    return sfle.ex([strProg] + [strTargetedTaxaFile, strHighlightCladeFile, strInvert, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                               strAbundanceFilterPercentile, strAbundanceFilterPercent, strRingOrder, strTicks] + [strSelection, strAbundance, strStyleFile, strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile, strT])
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
def funcStratifiedPCoASelectionMethods( strLoggingLevel, iSampleNameRow, iFirstDataRow, strStratifyMetadata, iNormalize, strInvert, strPredictFileArgument ):
  def funcStratifiedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, strStratifyMetadata=strStratifyMetadata, iNormalize=iNormalize, strInvert=strInvert, strPredictFileArgument=strPredictFileArgument):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[iSampleNameRow,iFirstDataRow,strStratifyMetadata,iNormalize,strInvert,strPredictFileArgument]+[ strAbnd, strSelection, strT])
  return funcStratifiedPCoARet

#Visualize output with Stratified PCoAs (Figure 3/4 Combined)
def funcCombinedStratifiedPCoASelectionMethods( strLoggingLevel, iSampleNameRow, iFirstDataRow, strStratifyMetadata, iNormalize, strInvert, lsSelectionMethods ):
  def funcCombinedStratifiedPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, strStratifyMetadata=strStratifyMetadata, iNormalize=iNormalize, strInvert=strInvert, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection = astrSs[0], astrSs[1], astrSs[2]
    return sfle.ex([strProg]+[iSampleNameRow,iFirstDataRow,strStratifyMetadata,iNormalize,strInvert]+[ strAbnd, strSelection, strT] + lsSelectionMethods)
  return funcCombinedStratifiedPCoARet

#Create a summary chart with collection curves
def funcCollectionCurveSummary( strLoggingLevel, strSampleNameRow, strFirstDataRow, strInvert):
  def funcCollectionCurveRet( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow, strFirstDataRow=strFirstDataRow, strInvert=strInvert):
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
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strFirstDataRow, strInvert, strAbundance, strFigureT, target[1].get_abspath()]+lsSelectionFiles)
  return funcCollectionCurveRet

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

###Start process
#Generate insilico data sets
lDataFiles = globFilesByExtension(strDirectory=fileDirInput,strExtension=c_strSufInsilicoData)
lsInsilicoDataFiles = []
for strInsilicoData in lsInsilicoDataNames:
  if not strInsilicoData in lDataFiles:
    if strInsilicoData == c_strInsilicoDataDiversity:
      lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,strInsilicoData )).get_abspath())
      lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )).get_abspath())
      Command(File( sfle.d( fileDirInput,strInsilicoData )), [c_fileProgUtilityData], funcGenerateData(strDataSetKey=c_strInsilicoDataDiversityKey))
      Command(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )),[c_progHCL,File( sfle.d( fileDirInput,strInsilicoData ))],funcVisualizeInsilicoData)
    elif strInsilicoData == c_strInsilicoDataUnbalanced:
      lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,strInsilicoData )).get_abspath())
      lsInsilicoDataFiles.append(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )).get_abspath())
      Command(File( sfle.d( fileDirInput,strInsilicoData )), [c_fileProgUtilityData], funcGenerateData(strDataSetKey=c_strInsilicoDataUnbalancedKey))
      Command(File( sfle.d( fileDirInput,sfle.rebase(strInsilicoData,c_strSufInsilicoData,c_strSufHCLUSTFig) )),[c_progHCL,File( sfle.d( fileDirInput,strInsilicoData ))],funcVisualizeInsilicoData)

#Get micropita config files
lMicropitaFiles = globFilesByExtension(strDirectory=fileDirInput,strExtension=c_strSufConfig)

#Collect input files later used for multistudy summary graphics
#{"InputFile:[SelectionFile1],[SelectionFile2],[SelectionFile3]...]}
dictSelectionFiles = dict()

#Read in each config file
#An use their contents to perform analysis
for fileConfigMicropita in lMicropitaFiles:
  print("fileConfigMicropita")
  print(fileConfigMicropita)
  #Get Contents of config file
  sFileConfiguration = dict()
  #Set defaults
  sFileConfiguration[c_strConfigSelectedTaxa]="None"
  sFileContents = list()
  with open(fileConfigMicropita.get_abspath()) as f:
    sFileContents = f.read()
    sFileContents = filter(None,re.split("\n",sFileContents))
    f.close()

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
          sFileConfiguration[strConfigKey]=lsConfigData[0]
    else:
      iIndex = iIndex + 1

  #Parse the supervised runs selection into multiple runs if needed
  if(c_strConfigSupervisedCount in sFileConfiguration):
    sSupervisedCounts = sFileConfiguration[c_strConfigSupervisedCount]
    sFileConfiguration[c_strConfigSupervisedCount] = filter(None,re.split(",",sSupervisedCounts))

  #Parse the unsupervised runs selection into multiple runs if needed
  if(c_strConfigUnsupervisedCount in sFileConfiguration):
    sSupervisedCounts = sFileConfiguration[c_strConfigUnsupervisedCount]
    sFileConfiguration[c_strConfigUnsupervisedCount] = filter(None,re.split(",",sSupervisedCounts))

  #Common configurations
  lsSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigSelectionTechniques]))
  lsPlotCombinedSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigPlotCombinedSelectionTechniques]))
  lsPlotCombinedSupervisedSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigPlotCombinedSupervisedSelectionTechniques]))
  strInvertImage = sFileConfiguration[c_strConfigInvertImage]
  if strInvertImage.lower()=="true":
    c_fileCladogramStyleFile = File(sfle.d( fileDirInput, "microPITA"+c_strSufCircStyleInv ))

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

  #Update the configurations that can be toggled on or off
  if sFileConfiguration[c_strConfigCladeFilter].lower() == "false":
    sFileConfiguration[c_strConfigCladeFilterLevel] = "none"
    sFileConfiguration[c_strConfigCladeFilterMeasure] = "none"
    sFileConfiguration[c_strConfigCladeFilterMinSize] = "none"

  if sFileConfiguration[c_strConfigAbundanceFilter].lower() == "false":
    sFileConfiguration[c_strConfigAbundanceFilterPercentile] = "none"
    sFileConfiguration[c_strConfigAbundanceFilterPercent] = "none"

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

  #Loop through the count selection using the unsupervised counts length for indexing
  for iCountIndex in xrange(0,len(lsIndexCounts)):

    sAbundanceFileName = sFileConfiguration[c_strConfigInputFile]

    #Build the micropita output name
    sMicropitaOutput = sAbundanceFileName
    if fSupervisedRun:
      sMicropitaOutput = "".join(["S",lsSupervisedCounts[iCountIndex],"-",sMicropitaOutput])
    if fUnsupervisedRun:
      if fIsStratified:
        sMicropitaOutput = "".join(["Strat",lsUnsupervisedCounts[iCountIndex],"-",sMicropitaOutput])
      else:
        sMicropitaOutput = "".join(["U",lsUnsupervisedCounts[iCountIndex],"-",sMicropitaOutput])
    sMicropitaOutput = File(sfle.d(fileDirOutput,c_strPathDelim.join([sFileBase,sMicropitaOutput])))

    #Make file names from the input micropita file in the temp directory
    sMicropitaPredictFile = File(sfle.d( fileDirTmp, sfle.rebase( File(sAbundanceFileName), c_strSufMicropita, c_strSufPredict )))

    #Make files from the micropita output file
    sOutputFigure1APCoA,sOutputCombinedFigure1APCoA,sOutputFigure4PCoA,sOutputFigure4CombinedPCoA,sOutputFigure1BHCL,sHCLSelectData,sHCLSelectColor,sHCLSelectLabel,sCladogramSelectFig,sCladogramTaxaFile,\
    sCladogramColorFile,sCladogramSizeFile,sCladogramTickFile,sCladogramHighlightFile,sCladogramCircleFile =  [File(sfle.d( sOutputDir,
    sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufPCOA,c_strSufCombinedPCOA,c_strSufStratPCOA,c_strSufCombinedStratPCOA,c_strSufHCLUSTFig,c_strSufHCLUSTData,c_strSufHCLUSTColor,
    c_strSufHCLUSTLabel,c_strSufFig2,c_strSufCircTaxa,c_strSufCircColor,c_strSufCircSize,c_strSufCircTick,c_strSufCircHighlight, c_strSufCircCircle)]

    #Make more files, these for cladogram options
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
    sAbundanceFileName = File(c_strPathDelim.join(["input",sAbundanceFileName]))
    Command(sMicropitaOutput, [c_fileProgMicroPITA, sAbundanceFileName] + c_filesSecondarySrc + [fileConfigMicropita] + lsInsilicoDataFiles, funcMicroPita(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                    sFileConfiguration[c_strConfigUnsupervisedCount][iCountIndex],
                                                                                                                                    " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                    " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                    " ".join([Constants_Arguments.c_strSupervisedLabelCount,sFileConfiguration[c_strConfigSupervisedCount][iCountIndex]]),
                                                                                                                                    " ".join([Constants_Arguments.c_strSupervisedLabel,sFileConfiguration[c_strConfigSupervisedLabel]]),
                                                                                                                                    " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadata,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                    cCladogramSelectedTaxa,
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
      Command(sOutputFigure1APCoA, [c_fileProgPCoAFigure, sAbundanceFileName, sMicropitaOutput] + c_filePrimarySrc, funcPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                                strPredictFileArgument))
    #Create figure 1A PCoA combined
    #Unstratified PCoA of selection
    #Manage optional files
    if sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none":
      Command(sOutputCombinedFigure1APCoA, [c_fileProgCombinedPCoAFigure, sAbundanceFileName, sMicropitaOutput] + c_filePrimarySrc, funcCombinedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                                " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                                 lsPlotCombinedSelectionMethods))

    #Create figure 1B HCL
    #Selection of samples by selection method
    Command([sOutputFigure1BHCL,sHCLSelectData, sHCLSelectColor, sHCLSelectLabel], [c_fileProgSelectionHCLFigure, sMicropitaOutput] + ls_srcFig1, 
        funcHCLSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]), " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage])))

    #Create a cladogram figure 2
    #Manage files which are optional
    strTaxaArgForCladogram = "-t"
    strHighlightArgForCladogram = "-c"
    lsFig2Sources = [c_fileProgSelectionCladogramFigure, sMicropitaOutput, sAbundanceFileName, c_fileCladogramStyleFile]
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
    #Taxanomic Enrichment
    lsCladogramFigures.append(sCladogramSelectFig)
    Command([sCladogramSelectFig, sCladogramTaxaFile, sCladogramColorFile, sCladogramTickFile, sCladogramHighlightFile,
         sCladogramSizeFile, sCladogramCircleFile], lsFig2Sources + ls_srcFig2, 
         funcCladogramSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                       strTaxaArgForCladogram,
                                       strHighlightArgForCladogram,
                                       " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                       " ".join([Constants_Arguments.c_strRoot,sFileConfiguration[c_strConfigRoot]]),
                                       " ".join([Constants_Arguments.c_strEnrichmentMethod,sFileConfiguration[c_strConfigEnrichmentMeasurement]]),
                                       " ".join([Constants_Arguments.c_strCladeFilterLevel,sFileConfiguration[c_strConfigCladeFilterLevel]]),
                                       " ".join([Constants_Arguments.c_strCladeMeasureLevel,sFileConfiguration[c_strConfigCladeFilterMeasure]]),
                                       " ".join([Constants_Arguments.c_strCladeFilteringMinLevel,sFileConfiguration[c_strConfigCladeFilterMinSize]]),
                                       " ".join([Constants_Arguments.c_strAbundanceFilterPercentile,sFileConfiguration[c_strConfigAbundanceFilterPercentile]]),
                                       " ".join([Constants_Arguments.c_strAbundanceFilterCutoff,sFileConfiguration[c_strConfigAbundanceFilterPercent]]),
                                       " ".join([Constants_Arguments.c_strRingOrder,sFileConfiguration[c_strConfigCladogramRingOrder]]),
                                       " ".join([Constants_Arguments.c_strCircladerTicks,sFileConfiguration[c_strConfigCladogramTicks]])))


    #Create figure 3
#    #Selection in Stratification
#    Command([sOutputFigure3HCL, sStratHCLSelectColor, sStratHCLSelectLabel], [c_fileProgStratSelectionHCLFigure, sMicropitaOutput, sAbundanceFileName] + ls_srcFig1, 
#        funcHCLStratSelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
#                                " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
#                                " ".join(["-id", c_strAbundanceIDCol]),
#                                " ".join([Constants_Arguments.c_strInvertArgument, strInvertImage])))

    #Create figure 4 PCoA
    #PCoA of stratified selection
    if not sFileConfiguration[c_strConfigUnsupervisedStratify].lower() == "none":
      lsStratPCOAFigures.append(sOutputFigure4PCoA)
      Command(sOutputFigure4PCoA, [c_fileProgStratifiedPCoAFigure, sAbundanceFileName, sMicropitaOutput] + c_filePrimarySrc, funcStratifiedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadata,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                   strPredictFileArgument))

      Command(sOutputFigure4CombinedPCoA, [c_fileProgCombinedStratifiedPCoAFigure, sAbundanceFileName, sMicropitaOutput] + c_filePrimarySrc, funcCombinedStratifiedPCoASelectionMethods(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strUnsupervisedStratifyMetadata,sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strNormalizeArgument,sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                   " ".join([Constants_Arguments.c_strInvertArgument,strInvertImage]),
                                                                                                                                                   lsPlotCombinedSelectionMethods+lsPlotCombinedSupervisedSelectionMethods))

    #Add input files to be later summarized in the metaplots
    #Create key combining abundance file and stratification status
    sKey = "".join(["STRAT-",sFileConfiguration[c_strConfigUnsupervisedStratify],"-INPUT-",os.path.basename(sAbundanceFileName.get_abspath())])
    if sKey not in dictSelectionFiles:
      dictSelectionFiles[sKey] = dict()
      curSummaryFileDict = dictSelectionFiles[sKey]
      curSummaryFileDict[c_strSelectionFiles] = list()

    dictSelectionFiles[sKey][c_strConfigInputFile] = sAbundanceFileName
    dictSelectionFiles[sKey][c_strSelectionFiles].append(sMicropitaOutput)

  #If this configuration file had multiple sampling
  if lsIndexCounts > 1:
#    if len(lsPCOAFigures) > 1:
#      Command(sOutputDir+c_strPCOAMovieEnding,[c_fileProgMencoder]+lsPCOAFigures, funcMakeMovie(lsImageFiles=lsPCOAFigures))
#    if len(lsStratPCOAFigures) > 1:
#      Command(sOutputDir+c_strStratPCOAMovieEnding,[c_fileProgMencoder]+lsStratPCOAFigures, funcMakeMovie(lsImageFiles=lsStratPCOAFigures))
    if len(lsCladogramFigures) > 1:
      Command(sOutputDir+c_strCladogramMovieEnding,[c_fileProgMencoder]+lsCladogramFigures, funcMakeMovie(lsImageFiles=lsCladogramFigures))

if c_fRunCollectionCurve:
  #For each input file
  #Create summary files
  #Create collection curves
  for strInputSummaryKey in dictSelectionFiles:
    curSummaryDict = dictSelectionFiles[strInputSummaryKey]
    curInputFile = curSummaryDict[c_strConfigInputFile]
    curListofSelection = curSummaryDict[c_strSelectionFiles]

    #Make more files, now for figure 5 Collection Curve
    sOutputFigure5CC, sOutputFigureText5CC =  [File(sfle.d( fileDirOutput.get_abspath()+strOutputSummaryFolder,
    sfle.rebase( strInputSummaryKey, c_strSufMicropita, s ))) for s in (c_strSufCollectionCurveFigure, c_strSufCollectionCurveText)]

    #Create Figure 5
    #Collection Curve
    Command([sOutputFigure5CC, sOutputFigureText5CC], [c_fileProgCollectionCurveFigure, curInputFile] + curListofSelection + ls_srcFig1, 
        funcCollectionCurveSummary(" ".join([Constants_Arguments.c_strLoggingArgument, sFileConfiguration[c_strConfigLogging]]),
                              " ".join([Constants_Arguments.c_strSampleNameRowArgument,sFileConfiguration[c_strConfigSampleRow]]),
                              " ".join([Constants_Arguments.c_strFirstDataRow,sFileConfiguration[c_strConfigDataRow]]),
                              " ".join([Constants_Arguments.c_strInvertArgument, strInvertImage])))

