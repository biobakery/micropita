#######################################################
# Author: Timothy Tickle
# Description: Class to Run Principle Coordinates Analysis
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from AbundanceTable import AbundanceTable
from Constants import Constants
from Constants_Figures import Constants_Figures
from Diversity import Diversity
from cogent.cluster.nmds import NMDS
from scipy.spatial.distance import squareform
from ValidateData import ValidateData
from matplotlib import pyplot as plt
import re

#To run PCoA first load the raw data or distance matrix using the "load" method, 
#  then use the "run" method to derive points, and then use "plot" to plot the graph.
#The process is structured in this way so that data is read once but can be transformed to different
#distance matricies and after analysis can be plotted with multiple sample highlighting.
#One can always reload or rerun data by calling
class PCoA:

    #Supported distance metrics
    c_BRAY_CURTIS="BRAY_CURTIS"

    #Holds the data Matrix
    dataMatrix=None
    #Indicates if the data matrix is raw data (True) or a distance matrix (False)
    isRawData=None

    #Current pcoa object
    pcoa = None

    #Loads data into PCoA (given the matrix or a valid file path)
    #Data can be passed or read in
    #Data can be the original data to be converted to a distance matrix or a distance matrix
    #If it is the orignal data, indicate that it is rawData (tempIsRawData=True)
     #If it is the distance matrix already generated indicate (tempIsRawData=False)
    #  and no conversion will occur in subsequent methods
    #Raw data file is assumed to have a column at index 0 that is the id for the rows, this is ignored.
    #Otherwise assumed to have samples as columns and features as rows.
    #If there are metadata rows at the top of the file, use tempNameRow to identify 
    #  the id row and tempFirstDataRow to skip the rest of the metadata rows.
    #Indices are 0 based.
    #Distance matrices are expected to be ######TODO 
    #@params tempReadData Either a Structured matrix of data or a valid file path to read from
    #@params tempDelimiter Delimiter for reading in the file
    #@params tempNameRow The index of the row that identifies the columns
    #@params tempFirstDataRow The index of teh first row to contain actual data
    #@params tempNormalize True normalizes each column by the sum of the column (columnElement=columnElement/sum(column))
    #@params tempCheckFile True indicates the files should be check. 
    #@return Return boolean indicator of success (True=Was able to load data)
    def loadData(self,tempReadData, tempIsRawData, tempDelimiter=Constants.TAB, tempNameRow=0, tempFirstDataRow=1, tempNormalize=True, tempCheckFile=True):

        #Indicates if data needs to be read or if a structured matrix is given
        readData=None

        #Validate parameters
        if(ValidateData.isValidFileName(tempReadData)):
            readData=True
        elif(ValidateData.isValidStructuredArray(tempReadData)):
            readData=False
        else:
            print("".join(["PCoA:loadData::Error tempReadData was not an existing file or an np array."]))
            return False

        #If a file path is given, read data.
        if(ValidateData.isTrue(readData)):
            #Object to read in and manipulate raw data
            rawData = AbundanceTable()
            #If indicated, check the input file
            if(ValidateData.isTrue(tempCheckFile)):
                tempReadData = rawData.checkRawDataFile(tempReadData)
                if(ValidateData.isFalse(tempReadData)):
                    print("".join(["PCoA:loadData::Error when checking raw data file, did not perform PCoA. File name:",str(tempReadData)]))
                    return False
            #Read in the file data to a numpy array.
            #Samples (column) by Taxa (rows)(lists) without the column
            data = rawData.textToArray(tempInputFile=tempReadData, tempDelimiter=tempDelimiter, tempNameRow=tempNameRow, tempFirstDataRow=tempFirstDataRow, tempNormalize=tempNormalize)
            if(ValidateData.isFalse(data)):
                print("PCoA:loadData::Error when reading checked raw data file, did not perform PCoA.")
                return False

            #Transpose data to be Taxa (columns) by samples (rows)(lists)
            data = rawData.transposeDataMatrix(tempMatrix=data[0], tempRemoveAdornments=False)
            if(ValidateData.isFalse(data)):
                print("PCoA:loadData::Error when transposing read raw data file, did not perform PCoA.")
                return False
            else:
                self.dataMatrix=data
                self.isRawData=tempIsRawData
        #Otherwise load the data directly as passed.
        else:
            self.dataMatrix=tempReadData
            self.isRawData=tempIsRawData
        return True

    #Runs analysis on loaded data
    #@params tempReadData Either a Structured matrix of data from an abundance table or a valid file path
    #@return Return string file path
    def run(self,tempDistanceMetric=None):

        #If distance metric is none, check to see if the matrix is a distance matrix
        #If so, run NMDS on the distance matrix
        #Otherwise return a false and do not run
        if(tempDistanceMetric==None):
            if(ValidateData.isTrue(self.isRawData)):
                print("PCoA:run::Error, no distance metric was specified but the previous load was not of a distance matrix.")
                return False
            elif(ValidateData.isFalse(self.isRawData)):
                self.pcoa = NMDS(dataMatrix, verbosity=0)
                return self.pcoa
        
        #Make sure the distance metric was a valid string type
        if(not ValidateData.isValidString(tempDistanceMetric)):
            print("PCoA:run::Error, distance metric was not a valid string type.")
            return False

        #Supported distances
        if(tempDistanceMetric==self.c_BRAY_CURTIS):
            distanceMatrix=Diversity().getBrayCurtisDissimilarity(tempSampleTaxaAbundancies=self.dataMatrix)
            if(ValidateData.isFalse(distanceMatrix)):
                print "ERROR"
                return False
            self.pcoa = NMDS(squareform(distanceMatrix), verbosity=0)
            return self.pcoa
        else:
            print("PCoA:run::Error, not a supported distance metric. Please generate the distance matrix and load.")
            return False

        return False

    #@params tempPlotName A valid file path to save the image of the plot
    def plot(self,tempPlotName="PCOA.png", tempColorGrouping='g', tempShape='o', tempLabels=["Green"], tempShapeLabels=["Circle"], tempShapeSize = 20, tempLegendLocation="upper right", tempInvert=False):
        if(not self.pcoa == None):
            #Get point count
            adPoints = self.pcoa.getPoints()
            iPointCount = len(adPoints[:,0])

            #Check shapes
            if(ValidateData.isValidList(tempShape)):
              if not len(tempShape) == iPointCount:
                print("Error, the list of shapes was given but was not the same size as the points so nothing was plotted.")
                return

            #Check colors
            if(ValidateData.isValidList(tempColorGrouping)):
              if not len(tempColorGrouping) == iPointCount:
                print("Error, the list of colors was given but was not the same size as the points so nothing was plotted.")
                return

            #Check sizes
            if(ValidateData.isValidList(tempShapeSize)):
              if not len(tempShapeSize) == iPointCount:
                print("Error, the list of sizes was given but was not the same size as the points so nothing was plotted.")
                return

            #Get plot object
            imgFigure = plt.figure()

            #Get plot colors
            objColors = Constants_Figures()
            objColors.invertColors(fInvert=tempInvert)

            #Color/Invert figure
            imgFigure.set_facecolor(objColors.c_strBackgroundColorWord)
            imgSubplot = imgFigure.add_subplot(111,axisbg=objColors.c_strBackgroundColorLetter)
            imgSubplot.set_xlabel("Dimension 1")
            imgSubplot.set_ylabel("Dimension 2")
            imgSubplot.spines['top'].set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.spines['bottom'].set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.spines['left'].set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.spines['right'].set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.xaxis.label.set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.yaxis.label.set_color(objColors.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='x', colors=objColors.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='y', colors=objColors.c_strDetailsColorLetter)
            charMarkerEdgeColor = objColors.c_strDetailsColorLetter

            #Plot colors seperately so the legend will pick up on the labels and make a legend
            if(ValidateData.isValidList(tempColorGrouping)):
                if len(tempColorGrouping) == iPointCount:
                    #Get unique colors and plot each individually
                    acharUniqueColors = list(set(tempColorGrouping))
                    for iColorIndex in xrange(0,len(acharUniqueColors)):
                        #Get the color
                        charColor = acharUniqueColors[iColorIndex]
                        #Get indices of colors
                        aiColorPointPositions = self.getIndices(tempColorGrouping,charColor)

                        #Reduce the labels by color
                        acharLabelsByColor = self.reduceList(tempLabels,aiColorPointPositions)

                        #Reduces sizes to indices if a list
                        reducedSizes = tempShapeSize
                        #Reduce sizes if a list
                        if(ValidateData.isValidList(reducedSizes)):
                          reducedSizes = self.reduceList(reducedSizes,aiColorPointPositions)

                        #Reduce to the current color grouping
                        aiXPoints = self.reduceList(adPoints[:,0],aiColorPointPositions)
                        aiYPoints = self.reduceList(adPoints[:,1],aiColorPointPositions)

                        #If the shapes are not a list plot
                        #Otherwise plot per shape per color (can not plot list of shapes in matplotlib)
                        reducedShapes = tempShape
                        if(not ValidateData.isValidList(reducedShapes)):
                          reducedShapes = reducedShapes[0]
                          imgSubplot.scatter(aiXPoints,aiYPoints, s=reducedSizes, c=[charColor], marker=reducedShapes, label=tempLabels[tempColorGrouping.index(charColor)], edgecolor=charMarkerEdgeColor)
                        #Shapes are supplied as a list so plot each shape
                        else:
                          #Reduce to shape sof the current colors
                          reducedShapes = self.reduceList(reducedShapes,aiColorPointPositions)
                          acharReducedShapesElements = list(set(reducedShapes))
                          #If there are multiple shapes, plot seperately because one is not allowed to plot them as a list
                          for aCharShapeElement in acharReducedShapesElements:
                            #Get indices
                            aiShapeIndices = self.getIndices(reducedShapes,aCharShapeElement)
                            #Reduce label to shapes
                            strShapeLabel = self.reduceList(acharLabelsByColor,aiShapeIndices)
                            #Get points per shape
                            aiXPointsPerShape = self.reduceList(aiXPoints,aiShapeIndices)
                            aiYPointsPerShape = self.reduceList(aiYPoints,aiShapeIndices)
                            #Get sizes per shape
                            #Reduce sizes if a list
                            reducedSizesPerShape = reducedSizes
                            if(ValidateData.isValidList(reducedSizes)):
                              reducedSizesPerShape = self.reduceList(reducedSizes,aiShapeIndices)
                            #Plot
                            imgSubplot.scatter(aiXPointsPerShape,aiYPointsPerShape, s=reducedSizesPerShape, c=[charColor], marker=aCharShapeElement, label=strShapeLabel[0], edgecolor=charMarkerEdgeColor)

            elif((not ValidateData.isValidList(tempColorGrouping)) and (ValidateData.isValidList(tempShape))):
                if len(tempShape) == iPointCount:
                    #Get unique shapes and plot each individually
                    acharUniqueShapes = list(set(tempShape))
                    for iShapeIndex in xrange(0,len(acharUniqueShapes)):
                        aiShapePointPositions = self.getIndices(tempShape,acharUniqueShapes[iShapeIndex])
                        #Reduce sizes if needed
                        reducedSizes = tempShapeSize
                        if(ValidateData.isValidList(reducedSizes)):
                          reducedSizes = self.reduceList(reducedSizes,aiShapePointPositions)
                        imgSubplot.scatter(self.reduceList(adPoints[:,0],aiShapePointPositions),self.reduceList(adPoints[:,1],aiShapePointPositions), s=reducedSizes, c=[tempColorGrouping], marker=acharUniqueShapes[iShapeIndex], label=tempLabels[tempShape.index(acharUniqueShapes[iShapeIndex])], edgecolor=charMarkerEdgeColor)
            else:
                imgSubplot.scatter(adPoints[:,0],adPoints[:,1], s=tempShapeSize, c=tempColorGrouping, marker=tempShape, label=tempLabels, edgecolor=charMarkerEdgeColor)

            objLegend = imgSubplot.legend(loc=tempLegendLocation, scatterpoints=1, prop={'size':10})

            #Invert legend
            if(tempInvert):
              objLegend.legendPatch.set_fc(objColors.c_strBackgroundColorWord)
              objLegend.legendPatch.set_ec(objColors.c_strDetailsColorLetter)
              plt.setp(objLegend.get_texts(),color=objColors.c_strDetailsColorLetter)

            imgFigure.savefig(tempPlotName, facecolor=imgFigure.get_facecolor())

    #TODO put in utilities
    #Returns the indicies of the element in the array as a list
    def getIndices(self, aList, dataElement):
        aretIndices = []
        for dataIndex in xrange(0,len(aList)):
            if(aList[dataIndex] == dataElement):
                aretIndices.append(dataIndex)
        return aretIndices

    #TODO put in utilities
    def reduceList(self, aList, dataIndicies):
        aretList = []
        for dataIndex in xrange(0,len(dataIndicies)):
            aretList.append(aList[dataIndicies[dataIndex]])
        return aretList

    #TODO put in utilities
    #The adColor should be in the range between 0 and 1
    @staticmethod
    def RGBToHex(adColor, charAlpha = "99"):
        charR = (hex(int(adColor[0]*255)))[2:]
        if(str(charR) == "0"):
            charR = "00"
        charG = (hex(int(adColor[1]*255)))[2:]
        if(str(charG) == "0"):
            charG = "00"
        charB = (hex(int(adColor[2]*255)))[2:]
        if(str(charB) == "0"):
            charB = "00"
        return "".join(["#",charR, charG, charB])

    @staticmethod
    def getShapes(intShapeCount):
      lsShapes = ['o','^','s','D','v','<','>','8','p','h']
      if intShapeCount > len(lsShapes):
        print("".join(["Error, PCoA.getShapes. Do not have enough shapes to give. Received request for ",str(intShapeCount)," shapes. Max available shape count is ",str(len(lsShapes)),"."]))
        return None
      return lsShapes[0:intShapeCount]
