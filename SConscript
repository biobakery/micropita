import sfle
import sys

Import( "*" )

c_strDataFile = sfle.d(fileDirInput, "data")

pE = DefaultEnvironment( )

#Invert
c_fInvert = "False"
#Normalize
c_NormalizePCOA = "False"

#Extentions
c_strSufSelectedTaxa = ".taxa"
c_strSufPng = ".png"
c_strSufPCOA = "-PCoA.png"
c_strSufHCLUSTData = ".HCLData"
c_strSufHCLUSTColor = ".HCLColor"
c_strSufHCLUSTLabel = ".HCLLabel"
c_strSufHCLUSTFig = "-SHCL.png"
c_strSufFig2 = "-Fig2.svg"
c_strSufCircTaxa = ".CTaxa"
c_strSufCircColor = ".CColor"
c_strSufCircSize = ".CSize"
c_strSufCircStyle = ".CStyle"
c_strSufCircTick = ".CTick"
c_strSufCircHighlight = ".CHLight"
c_strSufCircCircle = ".CCircle"
c_strSufMicropita = ".txt"
c_strSufPredict = "-SVM.predict"

c_strExtDelim = "."
c_strPathDelim = "/"

#The label used to determine the phenotype
#TODO make so I dont have to do this
strLabel = "dx"
#strLabel = "BMI"

#None generated output micropita files
#c_fileIBDWGSSelectionManual = File(sfle.d(fileDirInput,"IBDWGSSelection.txt"))
#c_fileIBDWGSSelectionOutput = File(sfle.d(fileDirOutput,"IBDWGSSelection/IBDWGSSelection.txt"))

#Cladogram style file for micropita
c_fileCladogramStyleFile = File(sfle.d( fileDirInput, "microPITA"+c_strSufCircStyle ))

#External programs
c_progHC = "./external/hclust/hclust.py"

#SRC Code
c_fileProgPCoAFigure = File( sfle.d( fileDirSrc, "MicropitaPaperPCoA.py" ) )
c_fileProgSelectionHCLFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionHCL.py" ) )
c_fileProgSelectionCladogramFigure = File( sfle.d( fileDirSrc, "MicropitaPaperSelectionCladogram.py") )
c_fileProgAbundanceTable = File( sfle.d( fileDirSrc, "AbundanceTable.py" ) )
c_fileProgCladogram = File( sfle.d( fileDirSrc, "Cladogram.py" ) )
c_fileProgCommandLine = File( sfle.d( fileDirSrc, "CommandLine.py" ) )
c_fileProgConstants = File( sfle.d( fileDirSrc, "Constants.py" ) )
c_fileProgConstantsFigures = File( sfle.d( fileDirSrc, "Constants_Figures.py" ) )
c_fileProgDiversity = File( sfle.d( fileDirSrc, "Diversity.py" ) )
c_fileProgFileIO = File( sfle.d( fileDirSrc, "FileIO.py" ) )
c_fileProgMicroPITA = File( sfle.d( fileDirSrc, "MicroPITA.py" ) )
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
c_lstrAllSelectionMethods = [c_strTaxa,c_strDiversity,c_strRepresentativeDissimilarity,c_strExtremeDissimilarity,c_strDiscriminant,c_strDistinct]

#Logging level
c_strLoggingLevel = "DEBUG"

##Start flow
#Scripts to call
##Call micropita to run
def funcMicroPita( iN, strParameters, strLabel ):
  def funcRet( target, source, env, iN = iN, strParameters = strParameters, strLabel = strLabel ):
    strT, astrSs = sfle.ts( target, source )
    strProg, strAbnd, strTaxa = astrSs[0], astrSs[1], astrSs[2]
    aArgs = [iN]+c_lstrAllSelectionMethods
    return sfle.ex( [strProg]+["-l", c_strLoggingLevel, strParameters,strLabel]+[strT, strAbnd, strTaxa, fileDirTmp] + aArgs)
  return funcRet

##Create figures
#Visualize output with PCoA (Figure 1A)
def funcPCoASelectionMethods( target, source, env ):
	strT, astrSs = sfle.ts( target, source )
	strProg, strAbnd, strSelection, strPrediction = astrSs[0], astrSs[1], astrSs[2], astrSs[3]
	return sfle.ex([strProg]+["-n","0","-d","2","-r",c_NormalizePCOA,"-i",c_fInvert]+[ strAbnd, strSelection, strPrediction, strT])

#Visualize selected output with HCL (Figure 1B)
def funcHCLSelectionMethods( target, source, env ):
	strT, astrSs = sfle.ts( target, source )
	strProg, strSelection = astrSs[0], astrSs[1]
        strData, strColor, strLabel = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath()
	return sfle.ex([strProg, strSelection, c_progHC, strData, strColor, strLabel, strT])

#Visualize selected output with Cladogram (Figure 2)
def funcCladogramSelectionMethods( target, source, env ):
	strT, astrSs = sfle.ts( target, source )
	strProg, strSelection, strAbundance, strStyleFile, strTargetedTaxa = astrSs[0], astrSs[1], astrSs[2], astrSs[3], astrSs[4]
        strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile = target[1].get_abspath(), target[2].get_abspath(), target[3].get_abspath(), target[4].get_abspath(), target[5].get_abspath(), target[6].get_abspath()
	return sfle.ex([strProg] + ["-t",strTargetedTaxa, "-i",c_fInvert] + [strSelection, strAbundance, strStyleFile, strTaxaFile, strColorFile, strTickFile, strHighlightFile, strSizeFile, strCircleFile, strT])

###Start process
#Read in all input files in the input directory
lFileInputFiles = Glob( sfle.d( fileDirInput, "*" ) )

#Make sure they are the correct input files for micropita
lMicropitaFiles = []
for fileName in lFileInputFiles:
  strPathPieces = [filter(None,strPathPiece) for strPathPiece in (fileName.get_abspath().split(c_strExtDelim))]
  if c_strExtDelim+strPathPieces[-1] == c_strSufMicropita:
    lMicropitaFiles.append(fileName)

for fileMicropita in lMicropitaFiles:
  sAbundanceFileName = fileMicropita.get_abspath().split(c_strPathDelim)[-1]
  sFileBase = sAbundanceFileName.split(c_strExtDelim)[0]
  sOutputDir = sfle.d(fileDirOutput,sFileBase)
  sMicropitaOutput = File(sfle.d(fileDirOutput,c_strPathDelim.join([sFileBase,sAbundanceFileName])))

  #Make file names from the input micropita file in the temp directory
  sMicropitaPredictFile = File(sfle.d( fileDirTmp, sfle.rebase( fileMicropita, c_strSufMicropita, c_strSufPredict )))

  #Make files from the micropita output file
  sOutputFigure1APCoA,sOutputFigure1BHCL,sHCLSelectData,sHCLSelectColor,sHCLSelectLabel,sCladogramSelectFig,sCladogramTaxaFile,\
  sCladogramColorFile,sCladogramSizeFile,sCladogramTickFile,sCladogramHighlightFile,sCladogramCircleFile =  [File(sfle.d( sOutputDir,
  sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) for s in (c_strSufPCOA,c_strSufHCLUSTFig,c_strSufHCLUSTData,c_strSufHCLUSTColor,
  c_strSufHCLUSTLabel,c_strSufFig2,c_strSufCircTaxa,c_strSufCircColor,c_strSufCircSize,c_strSufCircTick,c_strSufCircHighlight, c_strSufCircCircle)]

  cCladogramSelectedTaxa,cCladogramHighlightedTaxa = [File(sfle.d( fileDirInput, sfle.rebase( sMicropitaOutput, c_strSufMicropita, s ))) 
  for s in (c_strSufSelectedTaxa,c_strSufSelectedTaxa)]

  #Run micropita analysis
  sAbundanceFileName = c_strPathDelim.join(["input",sAbundanceFileName])

  Command(sMicropitaOutput, [c_fileProgMicroPITA, sAbundanceFileName, cCladogramSelectedTaxa] + c_filesSecondarySrc, funcMicroPita(8, "-n 0 -d 2 -s 4 -p", strLabel))
  #Create figure 1A PCoA
  Command(sOutputFigure1APCoA, [c_fileProgPCoAFigure, fileMicropita, sMicropitaOutput, sMicropitaPredictFile] + c_filePrimarySrc, funcPCoASelectionMethods)
  #Create figure 1B HCL
  Command([sOutputFigure1BHCL,sHCLSelectData, sHCLSelectColor, sHCLSelectLabel], [c_fileProgSelectionHCLFigure, sMicropitaOutput] + ls_srcFig1, 
        funcHCLSelectionMethods)
  #Create a cladogram figure 2
  Command([sCladogramSelectFig, sCladogramTaxaFile, sCladogramColorFile, sCladogramTickFile, sCladogramHighlightFile,
         sCladogramSizeFile, sCladogramCircleFile], [c_fileProgSelectionCladogramFigure, sMicropitaOutput, fileMicropita, c_fileCladogramStyleFile, cCladogramSelectedTaxa] + ls_srcFig2, 
        funcCladogramSelectionMethods)
