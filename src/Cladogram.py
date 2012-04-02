#######################################################
# Author: Timothy Tickle
# Description: This manages the basic infrastructure
# Needed for Circlader. Inherit this class and write over
# over methods you want to customize.
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from AbundanceTable import AbundanceTable
from CommandLine import CommandLine
from Constants import Constants
from Constants_Figures import Constants_Figures
from FileIO import FileIO
import math
import numpy as np
import os
import re
import scipy.stats
from ValidateData import ValidateData
#import scipy.stats.stats as stats

class Cladogram:
  """
  This class manages creating files for Circlader and calling circulator.
  """
  #Script name
  circladerScript = "/home/ttickle/Desktop/ttickle/PythonProjects/external/circlader/circlader.py"

  #Constants
  c_sTaxa="Taxa"
  c_sCircle="Circle"
  c_sBorder="Border"
  c_sShape="Shape"
  c_sAlpha="Alpha"
  c_sForced="Forced"

  #Numpy array (structured array) holding data
  #Should be SampleID, Sample Abundances/Data (samples = columns).....
  npaAbundance = None
  #List of sample names
  lsSampleNames = None
  #Name of output image
  strImageName = "Cladogram.png"
  #String used to call the sample id column
  strSampleID = "ID"

  #Minimum size of clade (terminal node count for clade)
  iMinCladeSize = -1
  #Level of ancestry to filter at (starts with 0 and based on the input file)
  iCladeLevelToMeasure = -1
  iCladeLevelToReduce = -1

  #Flags
  #Turns on (True) or off (False) abundance-based filtering
  fAbundanceFilter = True
  #Turns on (True) or off (False) clade size-based filtering
  fCladeSizeFilter = True
  #Indicate if the following files were made
  fSizeFileMade=False
  fCircleFileMade=False
  fColorFileMade=False
  fTickFileMade=False
  fHighlightFileMade=False

  #Circlader files
  strTreeFilePath="_Taxa.txt"
  strCircleFilePath = "_Circle.txt"
  strColorFilePath="_Color.txt"
  strTickFilePath="_Tick.txt"
  strHighLightFilePath="_HighLight.txt"
  strSizeFilePath="_Size.txt"
  strStyleFilePath=""

  #Thresholds
  #Controls the showing of taxa
  c_dPercentileCutOff = 90.0
  c_dPercentageAbovePercentile = 1.0

  #Minimum average abundance score when using log scale
  c_dMinLogSize = 0.0000000001
  #Constant used to maginfy the size difference in the taxa (log only)
  c_dLogScale = 1000000
  #When after log10, an addition scaling adjustment (use this)
  c_dCircleScale = 3

  #Data for circular files
  #Used to change IDs to proper labels
  dictConvertIDs = None
  #Labels to be relabeled
  dictRelabels = None
  #Colors
  dictColors = None
  #Elements that are forced to be highlighted
  dictForcedHighLights = None
  #Ticks
  llsTicks = None
  #Forced root of the tree, discarding data as needed.
  strRoot = None
  #Holds circle data as a list of dictionaries
  #One dictionary per circle
  ldictCircleData = None

  def __init__(self):
    self.dictForcedHighLights = dict()

  #Happy Path Tested
  def addHighLights(self, dictClades,fOverwrite):
    """
    This methods allows highlighting to be added.
    When an element is added in this manner it will not be filtered out.
    These elements, if existing in the tree will be highlighted the named color given.
    This color name should be supplied in the set Color Data method
    {strName1:strColorName1,strName2:strColorName2,...}

    :param dictClades Names of elements, if found in the tree which should be highlighted
    :type dictClades Dictionary of element name (string) and element color (string) 
    :param fOverwrite If element is already indicated to be highlighted, overwrite the color to the one provided here.
    :type fOverwrite boolean (True == overwrite color)
    """
    if ValidateData.isValidDictionary(dictClades):
        if ValidateData.isValidBoolean(fOverwrite):
            for strElement in dictClades:
                if(strElement in self.dictForcedHighLights):
                    if(fOverwrite):
                        self.dictForcedHighLights[strElement] = dictClades[strElement]
                else:
                    self.dictForcedHighLights[strElement] = dictClades[strElement]

  #Not tested
  def getHighLights(self):
    return self.dictForcedHighLights

  #Not tested
  def forceRoot(self, strRoot):
    """
    This method allows one to root the tree at a certain level and value
    Only taxa that contain this value in their ancestry will be plotted
    The root will be the value given, any previous heirachy will be ignored
    This will remove highlighted data if inidicated to do so

    :params strRoot Where to root the tree
    :type strRoot String
    """
    self.strRoot = strRoot

  def generate(self, strInputFile, strImageName, strStyleFile, sTaxaFileName, charDelimiter=Constants.TAB, iNameRow=0, iFirstDataRow=1, fNormalize=False, sColorFileName=None, sTickFileName=None, sHighlightFileName=None, sSizeFileName=None, sCircleFileName=None):
    """
    This is the method to call to generate a cladogram using circlader.
    The default data file is an abundance table unless the readData function is overwritten.

    :param strInputFile File to read in data
    :type strInputFile File path (string)
    :param strImageName File name to save the output cladogram image
    :type strImageName File name (string)
    :param strStyleFile File path indicating the style file to use
    :type strStyleFile File path (string)
    :param sTaxaFileName File path indicating the taxa file to use
    :type sTaxaFileName File path (string)
    :param strColorFile File path indicating the color file to use
    :type strColorFile File path (string)
    :param strTickFile File path indicating the tick file to use
    :type strTickFile File path (string)
    :param strHighlightFile File path indicating the highlight file to use
    :type strHighlightFile File path (string)
    :param strSizeFile File path indicating the size file to use
    :type strSizeFile File path (string)
    """
    #Set output file name
    self.strImageName = strImageName

    #Check files exist and remove files which will be written
    self.manageFilePaths(sTaxaFileName, strStyleFile, sColorFileName, sTickFileName, sHighlightFileName, sSizeFileName, sCircleFileName)

    #Read in the data
    self.npaAbundance, self.strSampleID = self.readData(strInputFile, charDelimiter=charDelimiter, iNameRow=iNameRow, iFirstDataRow=iFirstDataRow, fNormalize=fNormalize)
    self.lsSampleNames = self.npaAbundance.dtype.names[1:]

    #Get IDs
    lsIDs = [strId for strId in list(self.npaAbundance[self.strSampleID])]

    #Generate a dictionary to convert the ids to correct format
    #Fix unclassified names
    #Make numeric labels as indicated
    self.dictConvertIDs = self.generateLabels(lsIDs)

    print("lsIDs: baseline")
    print(len(lsIDs))

    #Filter by abundance
    if(self.fAbundanceFilter):
      lsIDs = self.filterByAbundance(lsIDs)

    print("lsIDs: filter by abundance")
    print(len(lsIDs))

    #Update to the correct root
    lsIDs = self.updateToRoot(lsIDs)

    #Set highlights to root for consistency
    if(not self.strRoot == None):
      dictRootedHighLights = dict()
      if not self.dictForcedHighLights == None:
        for sKey in self.dictForcedHighLights.keys():
          strUpdatedKey = self.updateToRoot([sKey])
          dictRootedHighLights[strUpdatedKey[0]]=self.dictForcedHighLights[sKey]
        self.dictForcedHighLights = dictRootedHighLights

    #Set relabels to root for consistency
    if(not self.strRoot == None):
      dictRootedLabels = dict()
      if not self.dictRelabels == None:
        for sKey in self.dictRelabels.keys():
          strUpdatedKey = self.updateToRoot([sKey])
          dictRootedLabels[strUpdatedKey[0]]=self.dictRelabels[sKey]
        self.dictRelabels = dictRootedLabels

    #Filter by clade size
    if(self.fCladeSizeFilter):
      lsIDs = self.filterByCladeSize(lsIDs)

    print("lsIDs: filter by clade size")
    print(len(lsIDs))

    #Add in forced highlighting
    lsIDs.extend(self.dictForcedHighLights.keys())
    lsIDs = list(set(lsIDs))

    print("lsIDs: forced highlights")
    print(len(lsIDs))

    #Create circle files (needs to be after any filtering because it has a forcing option).
    self.createCircleFile(lsIDs)

    print("lsIDs: create circle file")
    print(len(lsIDs))

    #Generate / Write Tree file
    self.createTreeFile(lsIDs)

    print("lsIDs: write tree file")
    print(len(lsIDs))

    #Generate / Write Highlight file
    self.createHighlightFile(lsIDs)

    print("lsIDs: write highlight file")
    print(len(lsIDs))

    #Generate / write color file
    if(self.dictColors is not None):
        lsColorData = [Constants.TAB.join([sColorKey,self.dictColors[sColorKey]]) for sColorKey in self.dictColors]
        self.writeToFile(self.strColorFilePath, Constants.ENDLINE.join(lsColorData), False)
        self.fColorFileMade=True

    #Generate / write tick file
    if(self.llsTicks is not None):
        lsTickData = [Constants.TAB.join(lsTicks) for lsTicks in self.llsTicks]
        self.writeToFile(self.strTickFilePath, Constants.ENDLINE.join(lsTickData), False)
        self.fTickFileMade=True

    #Generate / Write size data
    self.createSizeFile(lsIDs)

    print("lsIDs: write size file")
    print(len(lsIDs))

    #Call commandline
    lsCommand = [self.circladerScript, self.strTreeFilePath, self.strImageName, "--style_file", self.strStyleFilePath, "--tree_format", "tabular"]
    if(self.fSizeFileMade):
      lsCommand.extend(["--size_file", self.strSizeFilePath])
    if(self.fColorFileMade):
      lsCommand.extend(["--color_file", self.strColorFilePath])
    if(self.fTickFileMade):
      lsCommand.extend(["--tick_file", self.strTickFilePath])
    if(self.fHighlightFileMade):
      lsCommand.extend(["--highlight_file", self.strHighLightFilePath])
    if(self.fCircleFileMade):
      lsCommand.extend(["--circle_file", self.strCircleFilePath])
    CommandLine().runCommandLine(lsCommand)

  #Happy path tested
  def setColorData(self, dictColors):
    """
    This methods allows color information to be specified.
    Need to give a dictionary having a name (key)(string) and color (value)(string RGB)data
    {strName1:Color,strName2:Color...}
    Name will be a string name that references what needs to be this color
    Color data should be a string in the RGB format 0-255,0-255,0-255

    :param dictColors Color Name and RGB specification
    :type dictColorsDictionary strings 
    """
    if ValidateData.isValidDictionary(dictColors):
      self.dictColors = dictColors
      if not Constants_Figures.c_strBackgroundColorName in self.dictColors:
        self.dictColors[Constants_Figures.c_strBackgroundColorName]=Constants_Figures.c_strBackgroundColor

  #Not tested
  def setFilterByAbundance(self, fAbundanceFilter, dPercentileCutOff = 90.0,  dPercentageAbovePercentile = 1.0):
    """
    Switch filtering by abundance on and off.
    fAbundanceFilter == True indicates filtering is on

    :param fAbundanceFilter Switch to turn on (true) and off (false) abundance-based filtering
    :type fAbundanceFilter boolean
    """
    self.fAbundanceFilter = fAbundanceFilter
    self.c_dPercentileCutOff = dPercentileCutOff
    self.c_dPercentageAbovePercentile = dPercentageAbovePercentile

  def setCircleScale(self, iScale):
    """
    Is a scale used to increase or decrease node sizes in the the cladogram to make more visible
    iScale default is 3

    :param iScale Integer to increase the relative sizes of nodes
    :type iScale integer
    """
    self.c_dCircleScale = iScale

  #Not tested
  def setFilterByCladeSize(self, fCladeSizeFilter, iCladeLevelToMeasure = 3, iCladeLevelToReduce = 1, iMinimumCladeSize = 5):
    """
    Switch filtering by clade size on and off.
    fAbundanceFilter == True indicates filtering is on

    :param fAbundanceFilter Switch to turn on (true) and off (false) clade size-based filtering
    :type fAbundanceFilter boolean
    """
    self.fCladeSizeFilter = fCladeSizeFilter
    self.iCladeLevelToMeasure = iCladeLevelToMeasure
    self.iCladeLevelToReduce = iCladeLevelToReduce
    self.iMinCladeSize = iMinimumCladeSize

  #Not tested
  def setTicks(self, llsTicks):
    """
    This methods allows tick information to be specified.
    Need to generate a list of lists each having a tick level (number starting at 0 as a string), and tick name
    #Lowest numbers are closest to the center of the tree
    [[#,Name1],[#,Name2]...]

    :param llsTicks Level # and Name of level
    :type llsTicks List of lists of strings 
    """
    self.llsTicks = llsTicks

  def addCircle(self, lsTaxa, strCircle, dBorder=0.0, strShape="R", dAlpha=1.0, fForced=False):
    """
    This methods allows one to add a circle to the outside of the cladogram.
    Note, this will automatically add any valid circle to the tick marks if tickmarks exist
    One can add pieces of the circle at a time, the itck marks will reflect the order of the first
    Instance of a specific circle being added to the figure.

    :param lsTaxa Taxa to highlight with this circle
    :type lsTaxa List of strings (taxa names)
    :param strCircle Circle the elements will be in, indicates color and circle level.
    :type strCircle String circle 
    :param dBorder Border size for the circle element border (between 0.0 and 1.0)
      can also be a list of dBorders.  If list, position must match lsTaxa.
    :type dBorder Float of border size (or list of floats).
    :param strShape String Indicator of shape or method to determine shape.
      Can also be a list of shapes.  If list, position must match lsTaxa.
    :type strShape String to indicate the shape (may also be a list of strings).
        Default value is square.
        Valid shapes are R(Square), v(inward pointing triangle), ^(outward pointing triangle)
    :param dAlpha The transparency of the circle element (between 0.0[clear] and 1.0[solid]).
      Can also be a list of floats.  If list, position must match lsTaxa.
    :type dAlpha Float to indicate the transparency of the shape (may also be a list of strings).
    """
    if(self.ldictCircleData == None):
      self.ldictCircleData = list()
    dictCircleData = dict()
    dictCircleData[self.c_sTaxa]=lsTaxa
    dictCircleData[self.c_sCircle]=strCircle
    dictCircleData[self.c_sBorder]=dBorder
    dictCircleData[self.c_sShape]=strShape
    dictCircleData[self.c_sAlpha]=dAlpha
    dictCircleData[self.c_sForced]=fForced
    #TODO validate data, if not valid, return false and no add
    self.ldictCircleData.append(dictCircleData)
    return True

  def createCircleFile(self, lsIDs):
    """
    Write circle data to file.
    """
    #If there is circle data
    if(not self.ldictCircleData == None):
      if self.strCircleFilePath == None:
        print("Error, there is no circle file specified to write to.")
        return False
      #Holds circle data {Taxaname:string updates correctly for output to file}
      dictCircleDataMethods = dict()
      lsCircleData = list()

      for dictCircleData in self.ldictCircleData:
        lsTaxa = dictCircleData[self.c_sTaxa]
        #Shape/s for taxa
        datShape = dictCircleData[self.c_sShape]
        fShapeIsList = (str(type(datShape)) == "<type 'list'>")
        #Border/s for taxa
        datBorder = dictCircleData[self.c_sBorder]
        fBorderIsList = (str(type(datBorder)) == "<type 'list'>")
        #Alpha/s for taxa
        datAlpha = dictCircleData[self.c_sAlpha]
        fAlphaIsList = (str(type(datAlpha)) == "<type 'list'>")
        #Circle name
        sCircleMethod = dictCircleData[self.c_sCircle]

        #Check to make sure the lengths of the array match up
        if(fShapeIsList):
          if not len(datShape) == len(lsTaxa):
            print("".join(["Error, Shapes were given as an list not of the size of the taxa list. Shape list length: ",str(len(datShape)),". Taxa list length: ",str(len(lsTaxa)),"."]))
            return False
        if(fBorderIsList):
          if not len(datBorder) == len(lsTaxa):
            print("".join(["Error, Border sizes were given as an list not of the size of the taxa list. Border list length: ",str(len(datBorder)),". Taxa list length: ",str(len(lsTaxa)),"."]))
            return False
        if(fAlphaIsList):
          if not len(datAlpha) == len(lsTaxa):
            print("".join(["Error, Alpha sizes were given as an list not of the size of the taxa list. Alpha list length: ",str(len(datAlpha)),". Taxa list length: ",str(len(lsTaxa)),"."]))
            return False

        #Update taxa to root if needed
        #When doing this if any of the other data is an array we have to edit them
        #as the taxa are edited for updating root
        if((not fShapeIsList) and (not fBorderIsList) and (not fAlphaIsList)):
          lsTaxa = self.updateToRoot(dictCircleData[self.c_sTaxa])
        else:
          #Initilize as lists or as the string value they already are
          lsUpdatedTaxa = list()
          datUpdatedShapes=list()
          if(not fShapeIsList):
            datUpdatedShapes = datShape
          datUpdatedBorders=list()
          if(not fBorderIsList):
            datUpdatedBorders = datBorder
          datUpdatedAlphas=list()
          if(not fAlphaIsList):
            datUpdatedAlphas = datAlpha

          #If a taxa is kept, keep associated list information
          #If not a list data, leave alone, it will be used globally for all taxa.
          iTaxaIndex = -1
          for sTaxa in lsTaxa:
            iTaxaIndex = iTaxaIndex + 1
            sUpdatedTaxa=self.updateToRoot([sTaxa])

            if len(sUpdatedTaxa)==1:
              lsUpdatedTaxa.append(sUpdatedTaxa[0])
              if(fShapeIsList):
                datUpdatedShapes.append(datShape[iTaxaIndex])
              if(fBorderIsList):
                datUpdatedBorders.append(datBorder[iTaxaIndex])
              if(fAlphaIsList):
                datUpdatedAlphas.append(datAlpha[iTaxaIndex])

          #Reset data to rooted data
          lsTaxa=lsUpdatedTaxa
          datShape=datUpdatedShapes
          datBorder=datUpdatedBorders
          datAlpha=datUpdatedAlphas

        #QC passes so we will add the circle to the figure and the ticks.
        #If there are ticks and if the circle is not already in the ticks.
        if(not self.llsTicks == None):
          strCircleName = dictCircleData[self.c_sCircle]
          fFound = False
          iHighestNumber = -1
          for tick in self.llsTicks:
            #Look for name
            if tick[1] == strCircleName:
              fFound = True
            #Find highest count
            if int(tick[0]) > iHighestNumber:
              iHighestNumber = int(tick[0])
          if not fFound:
            self.llsTicks.append([str(iHighestNumber+1),strCircleName])

        #If the circle is forced, add the taxa to the lsIDs
        #Otherwise we will only plot those that are matching 
        #the lsIDs and the circle taxa list.
        if dictCircleData[self.c_sForced]:
          for iAlpha in xrange(0,len(datAlpha)):
            if(not datAlpha[iAlpha] == "0.0"):
              lsIDs.append(lsTaxa[iAlpha])
          lsIDs = list(set(lsIDs))

        #For all taxa in the cladogram
        for sTaxa in lsTaxa:
          #Store circle content name in dictionary
          if not sTaxa in dictCircleDataMethods:
            #Reset name to . delimited
            asNameElements = filter(None,re.split("\|",sTaxa))

            sCurTaxaName = asNameElements[len(asNameElements)-1]
            if(len(asNameElements)>1):
              if(sCurTaxaName=="unclassified"):
                sCurTaxaName = ".".join([asNameElements[len(asNameElements)-2],sCurTaxaName])
            sCurTaxa = ".".join(asNameElements)
            #Add to dictionary
            dictCircleDataMethods[sTaxa] = sCurTaxa

          #If the taxa is in the selected method
          if sTaxa in lsTaxa:
            #Index of the id in the circle data
            iTaxaIndex = lsTaxa.index(sTaxa)
            #Get border
            sBorder = ""
            if(fBorderIsList):
              sBorder = str(datBorder[iTaxaIndex])
            else:
              sBorder = str(datBorder)
            #Get shape
            sShape = ""
            if(fShapeIsList):
              sShape = datShape[iTaxaIndex]
            else:
              sShape = datShape
            #Get alpha
            sAlpha = ""
            if(fAlphaIsList):
              sAlpha = str(datAlpha[iTaxaIndex])
            else:
              sAlpha = str(datAlpha)
            dictCircleDataMethods[sTaxa]=dictCircleDataMethods[sTaxa]+"".join([Constants.TAB,sCircleMethod,":",sAlpha,"!",sShape,"#",sBorder])
          else:
            dictCircleDataMethods[sTaxa]=dictCircleDataMethods[sTaxa]+"".join([Constants.TAB,sCircleMethod,":0.0!R#0.0"])

      if len(dictCircleDataMethods)>0:
        lsTaxaKeys = dictCircleDataMethods.keys()
        sCircleContent = dictCircleDataMethods[lsTaxaKeys[0]]
        for sTaxaKey in lsTaxaKeys[1:len(lsTaxaKeys)]:
          sCircleContent = Constants.ENDLINE.join([sCircleContent,dictCircleDataMethods[sTaxaKey]])
        self.writeToFile(self.strCircleFilePath, sCircleContent, False)
        self.fCircleFileMade=True

        return True
    self.fCircleFileMade=False
    return False

############################################################

  def createHighlightFile(self, lsIDs):
    """
    Write highlight data to file

    :param lsIDs Ids to include in the highlight file
    :type lsIDs List of strings
    """
    lsHighLightData = list()
    #Each taxa name
    for sID in lsIDs:
      sCurColor = ""
      #Rename taxa to be consisten with the . delimit format
      asNameElements = filter(None,re.split("\|",sID))
      sCurTaxaName = asNameElements[len(asNameElements)-1]
      if(len(asNameElements)>1):
        if(sCurTaxaName=="unclassified"):
          sCurTaxaName = ".".join([asNameElements[len(asNameElements)-2],sCurTaxaName])
      sCurTaxa = ".".join(asNameElements)

      sCurLabel = ""
      #Get color
      sColorKey = ""
      if(sID in self.dictForcedHighLights):
        sColorKey = self.dictForcedHighLights[sID]
        if(sColorKey in self.dictColors):
          sCurColor = self.formatRGB(self.dictColors[sColorKey])
        #Get label
        if(self.dictRelabels is not None):
          if(sID in self.dictRelabels):
            sCurLabel = self.dictRelabels[sID]
        if(sCurLabel == ""):
          lsHighLightData.append(Constants.TAB.join([sCurTaxa,sCurTaxaName,sCurLabel,sCurColor]))
        else:
          lsHighLightData.append(Constants.TAB.join([sCurTaxa,sCurLabel,sCurLabel,sCurColor]))

    if len(lsHighLightData)>0:
      self.writeToFile(self.strHighLightFilePath, Constants.ENDLINE.join(lsHighLightData), False)
      self.fHighlightFileMade=True
    return True

  #Happy path tested
  def createSizeFile(self, lsIDs):
    """
    Write size data to file

    :param lsIDs Ids to include in the size file
    :type lsIDs List of strings
    """
    if self.npaAbundance is not None:
      dMinimumValue = (self.c_dMinLogSize*self.c_dLogScale)+1
      lsWriteData = list()
      for rowData in self.npaAbundance:
        strCurrentId = rowData[0]
        #Reset to root if needed to match current data
        if(not self.strRoot == None):
          strCurrentId = self.updateToRoot([strCurrentId])
          if(len(strCurrentId) > 0):
            strCurrentId = strCurrentId[0]
        if(strCurrentId in lsIDs):
          dAverage = np.average(list(rowData)[1:])
          dSize = max([dMinimumValue,(dAverage*self.c_dLogScale)+1])
          lsWriteData.append(".".join(re.split("\|",strCurrentId))+Constants.TAB+str(math.log10(dSize)*self.c_dCircleScale))
      if len(lsWriteData)>0:
        self.writeToFile(self.strSizeFilePath, Constants.ENDLINE.join(lsWriteData), False)
        self.fSizeFileMade=True
    return True

  def createTreeFile(self, lsIDs):
    """
    Write tree data to file. The tree file defines the internal cladogram and all it's points.

    :param lsIDs Ids to include in the tree file as well as their ancestors
    :type lsIDs List of strings
    """
    lsFullTree = list()
    for sID in lsIDs:
      lsIDElements = filter(None,re.split("\|",sID))
      sElementCur = lsIDElements[0]
      if(not sElementCur in lsFullTree):
        lsFullTree.append(sElementCur)
      if(len(lsIDElements) > 1):
        sNodePath = ""
        for iEndLevel in xrange(1,len(lsIDElements)+1):
          sCurAncestry = lsIDElements[0:iEndLevel]
          sNodePath = ".".join(sCurAncestry)
          if(not sNodePath in lsFullTree):
            lsFullTree.append(sNodePath)

    if len(lsFullTree)>0:
      self.writeToFile(self.strTreeFilePath, Constants.ENDLINE.join(lsFullTree), False)
    return True

  #Happy Path tested
  def filterByAbundance(self, lsIDs):
    """
    Filter by abundance. Specifically this version requires elements of
    the tree to have a certain percentage of a certain percentile in samples.

    :param lsIDs Ids to filter
    :type lsIDs List of strings
    """
    #list of ids to return that survived the filtering
    retls = list()
    if(self.npaAbundance is not None):
      #Hold the cuttoff score (threshold) for the percentile of interest {SampleName(string):score(double)}
      dictPercentiles = dict()
      for index in xrange(1,len(self.npaAbundance.dtype.names)):
        dScore = scipy.stats.scoreatpercentile(self.npaAbundance[self.npaAbundance.dtype.names[index]],self.c_dPercentileCutOff)
        dictPercentiles[self.npaAbundance.dtype.names[index]] = dScore

      #Sample count (Ignore sample id [position 0] which is not a name)
      dSampleCount = float(len(self.npaAbundance.dtype.names[1:]))

      #Check each taxa
      for rowTaxaData in self.npaAbundance:
        sCurTaxaName = rowTaxaData[0]
        #Only look at the IDs given
        if(sCurTaxaName in lsIDs):
          dCountAbovePercentile = 0.0
          ldAbundanceMeasures = list(rowTaxaData)[1:]
          #Check to see if the abundance score meets the threshold and count if it does
          for iScoreIndex in xrange(0,len(ldAbundanceMeasures)):
            if(ldAbundanceMeasures[iScoreIndex] >= dictPercentiles[self.lsSampleNames[iScoreIndex]]):
              dCountAbovePercentile = dCountAbovePercentile + 1.0
          dPercentOverPercentile = dCountAbovePercentile / dSampleCount
          if(dPercentOverPercentile >= (self.c_dPercentageAbovePercentile/100.0)):
            retls.append(sCurTaxaName)
    return retls

  #TODO need to test, definition of filtering redefined
  def filterByCladeSize(self, lsIDs):
    """
    Filter by the count of individuals in the clade.

    :param lsIDs Ids to filter
    :type lsIDs List of strings
    """
    #Count how many of a given level are in the clade
    dictCladeCount = dict()
    for sID in lsIDs:
        sIDElements = filter(None,re.split("\|",sID))
        if((len(sIDElements) >= self.iCladeLevelToMeasure) or ((sIDElements[-1] == "unclassified") and (len(sIDElements) >= self.iCladeLevelToMeasure))):
          if not sIDElements[self.iCladeLevelToReduce] in dictCladeCount:
              dictCladeCount[sIDElements[self.iCladeLevelToReduce]] = 0
          dictCladeCount[sIDElements[self.iCladeLevelToReduce]] = dictCladeCount[sIDElements[self.iCladeLevelToReduce]] + 1

    retls = list()
    for strID in lsIDs:
        strIDElements = filter(None,re.split("\|",strID))
        if(len(strIDElements) >= self.iCladeLevelToMeasure):
            if( strIDElements[self.iCladeLevelToReduce] in dictCladeCount):
                if dictCladeCount[strIDElements[self.iCladeLevelToReduce]] >= self.iMinCladeSize:
                    retls.append(strID)
    return retls

  #Happy path tested
  def formatRGB(self, sColor):
    """
    Takes a string that is of the format 0-255,0-255,0-255 and converts it to the 
    color format of circlader _c_[0-1,0-1,0-1]

    :param sColor String RGB format
    :type sColor String
    """
    sCircladerColor = "_c_[1,1,1]"
    if(sColor is not None):
      sColorElements = filter(None,re.split(",",sColor))
      if(len(sColorElements)==3):
        iR = int(sColorElements[0])/255.0
        iG = int(sColorElements[1])/255.0
        iB = int(sColorElements[2])/255.0
        sCircladerColor = "".join(["_c_[",str(iR),",",str(iG),",",str(iB),"]"])
    return sCircladerColor

  def generateLabels(self, lsIDs):
    """
    Labels for visualization. 
    Changes unclassified to one_level_higher.unclassified and enables numeric labeling / relabeling.
    Will only rename, will not add the label. The key must exist for the value to be used in replacing.
    """
    dictRet = dict()
    for sID in lsIDs:
        lsIDElements = filter(None,re.split("\|",sID))
        iIDElementsCount = len(lsIDElements)
        sLabel = lsIDElements[iIDElementsCount-1]
        #Fix unclassified
        if((sLabel == "unclassified") and (iIDElementsCount > 1)):
            sLabel = ".".join([lsIDElements[iIDElementsCount-2],sLabel])
        #Change to relabels if given
        if(self.dictRelabels is not None):
            if(sLabel in self.dictRelabels):
                sLabel = self.dictRelabels[sLabel]
        #Store lable
        dictRet[sID] = sLabel
    return dictRet

  def manageFilePaths(self, sTaxaFileName, strStyleFile, sColorFileName=None, sTickFileName=None, sHighlightFileName=None, sSizeFileName=None, sCircleFileName=None):
    """
    This method sets the naming to the files generated that Circlader acts on.
    These files include the tree, color, highlight, tick, and size files.
    Checks to make sure the file path to the syle file provided is an existing file.
    Deletes any existing files with these generated names (except for the style files).
    Private method

    :param strStyleFile File path indicating the style file to use
    :type strStyleFile File path (string)
    :param strTaxaFile File path indicating the taxa file to use
    :type strTaxaFile File path (string)
    :param strColorFile File path indicating the color file to use
    :type strColorFile File path (string)
    :param strTickFile File path indicating the tick file to use
    :type strTickFile File path (string)
    :param strHighlightFile File path indicating the highlight file to use
    :type strHighlightFile File path (string)
    :param strSizeFile File path indicating the size file to use
    :type strSizeFile File path (string)
    """
    #Do not remove the style file, it is static
    if(not os.path.exists(strStyleFile)):
      print("Error, no log file")
      return(-1)
    else:
      self.strStyleFilePath = strStyleFile

    #Set file names and remove if they exist
    #Tree file
    self.strTreeFilePath = sTaxaFileName
    if not self.strTreeFilePath == None:
      if(os.path.exists(self.strTreeFilePath)):
        os.remove(self.strTreeFilePath)
    #Color file
    self.strColorFilePath = sColorFileName
    if not self.strColorFilePath == None:
      if(os.path.exists(self.strColorFilePath)):
        os.remove(self.strColorFilePath)
    #Tick file
    self.strTickFilePath = sTickFileName
    if not self.strTickFilePath == None:
      if(os.path.exists(self.strTickFilePath)):
        os.remove(self.strTickFilePath)
    #Highlight file
    self.strHighLightFilePath = sHighlightFileName
    if not self.strHighLightFilePath == None:
      if(os.path.exists(self.strHighLightFilePath)):
        os.remove(self.strHighLightFilePath)
    #Size file
    self.strSizeFilePath = sSizeFileName
    if not self.strSizeFilePath == None:
      if(os.path.exists(self.strSizeFilePath)):
        os.remove(self.strSizeFilePath)
    #Circle file
    self.strCircleFilePath = sCircleFileName
    if not self.strCircleFilePath == None:
      if(os.path.exists(self.strCircleFilePath)):
        os.remove(self.strCircleFilePath)

  def readData(self, strInputFile, charDelimiter=Constants.TAB, iNameRow=0, iFirstDataRow=1, fNormalize=False):
    """
    This reads data from a file to build a cladogram.

    :param strInputFile File to read in data
    :type strInputFile File path 
    """
    #Abundance table object to read in and manage data
    rawData = AbundanceTable()

    #Read in abundance data
    #Abundance is a structured array. Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
    npaData,metadata = rawData.textToStructuredArray(tempInputFile=strInputFile, tempDelimiter=charDelimiter, tempNameRow=iNameRow, tempFirstDataRow=iFirstDataRow, tempNormalize=fNormalize)
    self.strSampleID = npaData.dtype.names[0]

    #Return the data array and the sample ID
    return([npaData, self.strSampleID])

  def relabelIDs(self, dictLabels):
    """
    Allows the relabeling of ids. Can be used to make numeric labeling of ids or renaming

    :param dictLabels Should label (key) (after unclassified is modified) and new label (value)
    :type dictLabels Dictionary of string (key) string (value)
    """
    self.dictRelabels = dictLabels

  #Happy path tested
  def updateToRoot(self, lsIDs):
    """ Updates the clade to the root given. The clade must contain the root and the level of the 
    root in the clade will be rest to it's first level, ignoring the previous levels of the clade.

    :param lsIDs List of Clades that will be reset to the root specified by setRoot
    :type lsIDs List of strings. Each string representing a clade.
    """

    if(self.strRoot is None):
      return lsIDs
    #Force root tree if indicated to do so
    lsRootedIDs = list()
    for sID in lsIDs:
      sIDElements = filter(None,re.split("\|",sID))
      if(self.strRoot in sIDElements):
        iRootIndex = sIDElements.index(self.strRoot)
        #If multiple levels of the clade exist after the new root merge them.
        if(len(sIDElements)>iRootIndex+2):
          lsRootedIDs.append("|".join(sIDElements[iRootIndex+1:]))
        #If only one level of the clade exists after the new root, return it.
        elif(len(sIDElements)>iRootIndex+1):
          lsRootedIDs.append(sIDElements[iRootIndex+1])
    return(lsRootedIDs)

  def writeToFile(self, strFileName, strDataToWrite, fAppend):
    """
    Helper function that writes a string to a file

    :param strFileName File to write to
    :type strFileName File path (string)
    :param strDataToWrite Data to write to file
    :type strDataToWrite String
    :param fAppend Indicates if an append should occur (True == Append)
    :type fAppend boolean
    """
    fileHandle = FileIO(strFileName,False,True,fAppend)
    fileHandle.writeToFile(strDataToWrite)
    fileHandle.close()
