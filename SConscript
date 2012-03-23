import re
import sfle
import sys

Import( "*" )

c_strDataFile = sfle.d(fileDirInput, "data")

pE = DefaultEnvironment( )

#Normalize
c_NormalizePCOA = "False"

#Extentions
c_strSufCollectionCurveFigure = "-ColCurv.png"
c_strSufCollectionCurveText = "-ColCurv.txt"
c_strSufSelectedTaxa = ".taxa"
c_strSufPng = ".png"
c_strSufPCOA = "-PCoA.png"
c_strSufHCLUSTData = ".HCLData"
c_strSufStratHCLUSTColor = "-Strat.HCLColor"
c_strSufHCLUSTColor = ".HCLColor"
c_strSufStratHCLUSTLabel = "-Strat.HCLLabel"
c_strSufHCLUSTLabel = ".HCLLabel"
c_strSufStratHCLUSTFig = "-StratHCL.png"
c_strSufHCLUSTFig = "-SHCL.png"
c_strSufFig2 = "-Fig2.svg"
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
c_strConfigDataRow = "[Data Name Row Index]"
c_strConfigEnrichmentMeasurement = "[Taxa or OTU Enrichment Measurement]"
c_strConfigHighlightClades = "[Higlight Clades]"
c_strConfigInputFile = "[Input File]"
c_strConfigInvertImage = "[Invert Image]"
c_strConfigLabels = "[Labels]"
c_strConfigLogging = "[Logging]"
c_strConfigNormalizeAbundance = "[Normalize Abundance]"
c_strConfigRoot = "[Root]"
c_strConfigSampleRow = "[Sample Name Row Index]"
c_strConfigSelection = "[Selection]"
c_strConfigSelectedTaxa = "[Selected Taxa]"
c_strConfigSelectionTechniques = "[Selection Techniques]"
c_strConfigStratifiedCount = "[Stratified Count]"
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
c_progHC = "./external/hclust/hclust.py"

#Src Code
c_fileProgAbundanceTable = File( sfle.d( fileDirSrc, "AbundanceTable.py" ) )
c_fileProgCladogram = File( sfle.d( fileDirSrc, "Cladogram.py" ) )
c_fileProgCommandLine = File( sfle.d( fileDirSrc, "CommandLine.py" ) )
c_fileProgConstants = File( sfle.d( fileDirSrc, "Constants.py" ) )
c_fileProgConstantsFigures = File( sfle.d( fileDirSrc, "Constants_Figures.py" ) )
c_fileProgDiversity = File( sfle.d( fileDirSrc, "Diversity.py" ) )
c_fileProgFileIO = File( sfle.d( fileDirSrc, "FileIO.py" ) )
c_fileProgMicroPITA = File( sfle.d( fileDirSrc, "MicroPITA.py" ) )
c_fileProgCollectionCurveFigure = File( sfle.d( fileDirSrc, "MicropitaPaperCollectionCurve.py") )
c_fileProgPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperPCoA.py" ) )
c_fileProgSelectionCladogramFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionCladogram.py") )
c_fileProgSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionHCL.py" ) )
c_fileProgStratSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperStratSelectionHCL.py" ) )
c_fileProgMLPYDistanceAdaptor = File( sfle.d( fileDirSrc, "MLPYDistanceAdaptor.py" ) )
c_fileProgPCOA = File( sfle.d( fileDirSrc, "PCoA.py" ) )
c_fileProgSVM = File( sfle.d( fileDirSrc, "SVM.py" ) )
c_fileProgUtilityData = File( sfle.d( fileDirSrc, "Utility_Data.py" ) )
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
c_strSVM = [c_strDiscriminant,c_strDistinct]
c_strTaxa = "Taxa_Defined"
#c_lstrAllSelectionMethods = [c_strRandom,c_strTaxa,c_strDiversity,c_strRepresentativeDissimilarity,c_strExtremeDissimilarity,c_strDiscriminant,c_strDistinct]
#c_lstrAllSelectionMethods = [c_strTaxa,c_strDiversity,c_strRepresentativeDissimilarity,c_strExtremeDissimilarity,c_strDiscriminant,c_strDistinct]

##Start flow
#Scripts to call
##Call micropita to run
def funcMicroPita( strLoggingLevel, iUnsupervisedCount, iSampleNameRow, iFirstDataRow, iSupervisedCount, strSupervisedLabel, strUnsupervisedStratify, strTaxaFile, lsSelectionMethods):
  def funcMicroPitaRet( target, source, env, iUnsupervisedCount=iUnsupervisedCount, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, 
                        iSupervisedCount=iSupervisedCount, strSupervisedLabel=strSupervisedLabel, strUnsupervisedStratify=strUnsupervisedStratify, strTaxaFile=strTaxaFile, lsSelectionMethods=lsSelectionMethods):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strTaxa = astrSs[0], astrSs[1], astrSs[2]
    aArgs = [iUnsupervisedCount]+lsSelectionMethods
    return sfle.ex( [strProg]+[strLoggingLevel, iSampleNameRow, iFirstDataRow, iSupervisedCount, strSupervisedLabel, strUnsupervisedStratify]+[strT, strAbnd, strTaxaFile, fileDirTmp] + aArgs)
  return funcMicroPitaRet

##Create figures
#Visualize output with PCoA (Figure 1A)
def funcPCoASelectionMethods( strLoggingLevel, iSampleNameRow, iFirstDataRow, iNormalize, strInvert ):
  def funcPCoARet( target, source, env, strLoggingLevel=strLoggingLevel, iSampleNameRow=iSampleNameRow, iFirstDataRow=iFirstDataRow, iNormalize=iNormalize, strInvert=strInvert):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strSelection, strPrediction = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
    return sfle.ex([strProg]+[iSampleNameRow,iFirstDataRow,iNormalize,strInvert]+[ strAbnd, strSelection, strPrediction, strT])
  return funcPCoARet

#Visualize selected output with HCL (Figure 1B)
def funcHCLSelectionMethods( strLoggingLevel, strInvert ):
  def funcHCLSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection = astrSs[0], astrSs[1]
    strData, strColor, strLabel = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath()
    return sfle.ex([strProg, strLoggingLevel, strInvert, strSelection, c_progHC, strData, strColor, strLabel, strT])
  return funcHCLSelectionRet

#Visualize selected output with Cladogram (Figure 2)
def funcCladogramSelectionMethods( strTargetedTaxaFile, strHighlightCladeFile, strInvert, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                                   strAbundanceFilterPercentile, strAbundanceFilterPercent):

  def funcCladogramSelectionRet( target, source, env, strTargetedTaxaFile=strTargetedTaxaFile, strHighlightCladeFile=strHighlightCladeFile, strInvert=strInvert, strRoot=strRoot, strEnrichment=strEnrichment,
                                 strCladeFilterLevel=strCladeFilterLevel, strCladeFilterMeasure=strCladeFilterMeasure, strCladeFilterMin=strCladeFilterMin, strAbundanceFilterPercentile=strAbundanceFilterPercentile,
                                 strAbundanceFitlerPercent=strAbundanceFilterPercent):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection, strAbundance, strStyleFile, strTargetedTaxa = astrSs[0], astrSs[1], astrSs[2], astrSs[3], astrSs[4]
    strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath(), target[4].get_abspath(), target[5].get_abspath(), target[6].get_abspath()
    return sfle.ex([strProg] + [strTargetedTaxaFile, strHighlightCladeFile, strInvert, strRoot, strEnrichment, strCladeFilterLevel, strCladeFilterMeasure, strCladeFilterMin,
                               strAbundanceFilterPercentile, strAbundanceFilterPercent] + [strSelection, strAbundance, strStyleFile, strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile, strT])
  return funcCladogramSelectionRet

#Visualize with HCL selected samples with stratification (Figure 3)
def funcHCLStratSelectionMethods( strLoggingLevel, strInvert, strDataRow, strIdCol ):
  def funcHCLStratSelectionRet( target, source, env, strLoggingLevel=strLoggingLevel, strInvert=strInvert, strDataRow=strDataRow, strIdCol=strIdCol ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strSelection, strAbundance = astrSs[0], astrSs[1], astrSs[2]
    strColor, strLabel = target[1].get_abspath(), target[2].get_abspath()
    return sfle.ex([strProg, strLoggingLevel, strDataRow, strIdCol, strInvert, strSelection, strAbundance, c_progHC, strColor, strLabel, strT])
  return funcHCLStratSelectionRet

#Create a summary chart with collection curves
def funcCollectionCurveSummary( strLoggingLevel, strSampleNameRow, strFirstDataRow, strNormalize, strInvert):
  def funcCollectionCurveSummary( target, source, env, strLoggingLevel=strLoggingLevel, strSampleNameRow=strSampleNameRow, strFirstDataRow=strFirstDataRow, strNormalize=strNormalize, strInvert=strInvert):
    print("START")
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
    return sfle.ex([strProg, strLoggingLevel, strSampleNameRow, strFirstDataRow, strNormalize, strInvert, strAbundance, strFigureT, target[1].get_abspath()]+lsSelectionFiles)
  return funcCollectionCurveSummary

###Start process
#Read in all input files in the input directory
lFileInputFiles = Glob( sfle.d( fileDirInput, "*" ) )

#Make sure they are the correct input files for micropita
lMicropitaFiles = []
for fileName in lFileInputFiles:
  strPathPieces = [filter(None,strPathPiece) for strPathPiece in (fileName.get_abspath().split(c_strExtDelim))]
  if c_strExtDelim+strPathPieces[-1] == c_strSufConfig:
    lMicropitaFiles.append(fileName)

print("lMicropitaFiles")
print([strFile.get_abspath() for strFile in lMicropitaFiles])

#Collect input files later used for multistudy summary graphics
#{"InputFile:[SelectionFile1],[SelectionFile2],[SelectionFile3]...]}
dictSelectionFiles = dict()

#Read in each config file
#An use thier contents to perform analysis
for fileConfigMicropita in lMicropitaFiles:
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

  #Common configurations
  lsSelectionMethods = filter(None,re.split(",",sFileConfiguration[c_strConfigSelectionTechniques]))
  strInvertImage = sFileConfiguration[c_strConfigInvertImage]
  if strInvertImage.lower()=="true":
    c_fileCladogramStyleFile = File(sfle.d( fileDirInput, "microPITA"+c_strSufCircStyleInv ))

  #Update the configurations that can be toggled on or off
  if sFileConfiguration[c_strConfigCladeFilter].lower() == "false":
    sFileConfiguration[c_strConfigCladeFilterLevel] = "None"
    sFileConfiguration[c_strConfigCladeFilterMeasure] = "None"
    sFileConfiguration[c_strConfigCladeFilterMinSize] = "None"

  if sFileConfiguration[c_strConfigAbundanceFilter].lower() == "false":
    sFileConfiguration[c_strConfigAbundanceFilterPercentile] = "None"
    sFileConfiguration[c_strConfigAbundanceFilterPercent] = "None"

  sAbundanceFileName = sFileConfiguration[c_strConfigInputFile]
  sFileBase = os.path.basename(fileConfigMicropita.get_abspath()).split(c_strExtDelim)[0]
  sOutputDir = sfle.d(fileDirOutput,sFileBase)
  sMicropitaOutput = File(sfle.d(fileDirOutput,c_strPathDelim.join([sFileBase,sAbundanceFileName])))

  #Make file names from the input micropita file in the temp directory
  sMicropitaPredictFile = File(sfle.d( fileDirTmp, sfle.rebase( File(sAbundanceFileName), c_strSufMicropita, c_strSufPredict )))

  #Make files from the micropita output file
  sOutputFigure1APCoA,sOutputFigure1BHCL,sHCLSelectData,sHCLSelectColor,sHCLSelectLabel,sCladogramSelectFig,sCladogramTaxaFile,\
  sCladogramColorFile,sCladogramSizeFile,sCladogramTickFile,sCladogramHighlightFile,sCladogramCircleFile =  [File(sfle.d( sOutputDir,
  sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufPCOA,c_strSufHCLUSTFig,c_strSufHCLUSTData,c_strSufHCLUSTColor,
  c_strSufHCLUSTLabel,c_strSufFig2,c_strSufCircTaxa,c_strSufCircColor,c_strSufCircSize,c_strSufCircTick,c_strSufCircHighlight, c_strSufCircCircle)]

  #Make more files, these for cladogram options
  cCladogramSelectedTaxa,cCladogramHighlightedTaxa = [File(sfle.d( fileDirInput, sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) 
  for s in (c_strSufSelectedTaxa,c_strSufSelectedTaxa)]

  #Make more files, now for the figure 3 HCL
  sOutputFigure3HCL, sStratHCLSelectColor, sStratHCLSelectLabel =  [File(sfle.d( sOutputDir,
  sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufStratHCLUSTFig,c_strSufStratHCLUSTColor,c_strSufStratHCLUSTLabel)]

  #Run micropita analysis
  sAbundanceFileName = File(c_strPathDelim.join(["input",sAbundanceFileName]))
  Command(sMicropitaOutput, [c_fileProgMicroPITA, sAbundanceFileName, cCladogramSelectedTaxa] + c_filesSecondarySrc + [fileConfigMicropita], funcMicroPita(" ".join(["-l", sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                    sFileConfiguration[c_strConfigUnsupervisedCount],
                                                                                                                                    " ".join(["-n",sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                    " ".join(["-d",sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                    " ".join(["-s",sFileConfiguration[c_strConfigSupervisedCount]]),
                                                                                                                                    " ".join(["-p",sFileConfiguration[c_strConfigSupervisedLabel]]),
                                                                                                                                    " ".join(["-u",sFileConfiguration[c_strConfigUnsupervisedStratify]]),
                                                                                                                                    sFileConfiguration[c_strConfigSelectedTaxa],
                                                                                                                                    lsSelectionMethods))
  #Create figure 1A PCoA
  #PCoA of selection
  Command(sOutputFigure1APCoA, [c_fileProgPCoAFigure, sAbundanceFileName, sMicropitaOutput, sMicropitaPredictFile] + c_filePrimarySrc, funcPCoASelectionMethods(" ".join(["-l", sFileConfiguration[c_strConfigLogging]]),
                                                                                                                                                                " ".join(["-n",sFileConfiguration[c_strConfigSampleRow]]),
                                                                                                                                                                " ".join(["-d",sFileConfiguration[c_strConfigDataRow]]),
                                                                                                                                                                " ".join(["-r",sFileConfiguration[c_strConfigNormalizeAbundance]]),
                                                                                                                                                                " ".join(["-i",strInvertImage])))
  #Create figure 1B HCL
  #Selection of samples by selection method
  Command([sOutputFigure1BHCL,sHCLSelectData, sHCLSelectColor, sHCLSelectLabel], [c_fileProgSelectionHCLFigure, sMicropitaOutput] + ls_srcFig1, 
        funcHCLSelectionMethods(" ".join(["-l", sFileConfiguration[c_strConfigLogging]]), " ".join(["-i",strInvertImage])))

  #Create a cladogram figure 2
  #Taxanomic Enrichment
  Command([sCladogramSelectFig, sCladogramTaxaFile, sCladogramColorFile, sCladogramTickFile, sCladogramHighlightFile,
         sCladogramSizeFile, sCladogramCircleFile], [c_fileProgSelectionCladogramFigure, sMicropitaOutput, sAbundanceFileName, c_fileCladogramStyleFile, cCladogramSelectedTaxa] + ls_srcFig2, 
         funcCladogramSelectionMethods(" ".join(["-t",sFileConfiguration[c_strConfigSelectedTaxa]]),
                                       " ".join(["-c",sFileConfiguration[c_strConfigHighlightClades]]),
                                       " ".join(["-i",strInvertImage]),
                                       " ".join(["-rt",sFileConfiguration[c_strConfigRoot]]),
                                       " ".join(["-e",sFileConfiguration[c_strConfigEnrichmentMeasurement]]),
                                       " ".join(["-cfl",sFileConfiguration[c_strConfigCladeFilterLevel]]),
                                       " ".join(["-cfm",sFileConfiguration[c_strConfigCladeFilterMeasure]]),
                                       " ".join(["-cfn",sFileConfiguration[c_strConfigCladeFilterMinSize]]),
                                       " ".join(["-afp",sFileConfiguration[c_strConfigAbundanceFilterPercentile]]),
                                       " ".join(["-afc",sFileConfiguration[c_strConfigAbundanceFilterPercent]])))

  #Create figure 3
  #Selection in Stratification
  Command([sOutputFigure3HCL, sStratHCLSelectColor, sStratHCLSelectLabel], [c_fileProgStratSelectionHCLFigure, sMicropitaOutput, sAbundanceFileName] + ls_srcFig1, 
        funcHCLStratSelectionMethods(" ".join(["-l", sFileConfiguration[c_strConfigLogging]]),
                                " ".join(["-d",sFileConfiguration[c_strConfigDataRow]]),
                                " ".join(["-id", c_strAbundanceIDCol]),
                                " ".join(["-i", strInvertImage])))

  #Add input files to be later summarized in the metaplots
  #Create key combining abundance file and stratification status
  sKey = "".join(["STRAT-",sFileConfiguration[c_strConfigUnsupervisedStratify],"-INPUT-",os.path.basename(sAbundanceFileName.get_abspath())])
  if sKey not in dictSelectionFiles:
    dictSelectionFiles[sKey] = dict()
    curSummaryFileDict = dictSelectionFiles[sKey]
    curSummaryFileDict[c_strSelectionFiles] = list()
  dictSelectionFiles[sKey][c_strConfigInputFile] = sAbundanceFileName
  dictSelectionFiles[sKey][c_strSelectionFiles].append(sMicropitaOutput)

#For each input file
print("dictSelectionFiles")
print(dictSelectionFiles)
for strInputSummaryKey in dictSelectionFiles:
  curSummaryDict = dictSelectionFiles[strInputSummaryKey]
  curInputFile = curSummaryDict[c_strConfigInputFile]
  curListofSelection = curSummaryDict[c_strSelectionFiles]

  #Make more files, now for figure 5 Collection Curve
  sOutputFigure5CC, sOutputFigureText5CC =  [File(sfle.d( fileDirOutput.get_abspath()+strOutputSummaryFolder,
  sfle.rebase( strInputSummaryKey, c_strSufMicropita, s ))) for s in (c_strSufCollectionCurveFigure, c_strSufCollectionCurveText)]

  print([sOutputFigure5CC.get_abspath(), sOutputFigureText5CC.get_abspath()])

  #Create Figure 5
  #Collection Curve
  print("curInputFile")
  print(curInputFile.get_abspath())
  Command([sOutputFigure5CC, sOutputFigureText5CC], [c_fileProgCollectionCurveFigure, curInputFile] + curListofSelection + ls_srcFig1, 
      funcCollectionCurveSummary(" ".join(["-l", sFileConfiguration[c_strConfigLogging]]),
                              " ".join(["-n",sFileConfiguration[c_strConfigSampleRow]]),
                              " ".join(["-d",sFileConfiguration[c_strConfigDataRow]]),
                              " ".join(["-r",sFileConfiguration[c_strConfigNormalizeAbundance]]),
                              " ".join(["-i", strInvertImage])))

