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
import logging
import math
import matplotlib.cm as cm
import numpy as np
from scipy.spatial.distance import squareform
from Utility_Math import Utility_Math
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

    #Get plot colors
    objFigureControl = Constants_Figures()

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
#    def loadData(self,tempReadData, tempIsRawData, tempDelimiter=Constants.TAB, tempNameRow=0, tempFirstDataRow=1, tempNormalize=True, tempCheckFile=True):
#
#        #Indicates if data needs to be read or if a structured matrix is given
#        readData=None
#
#        #Validate parameters
#        if(ValidateData.isValidFileName(tempReadData)):
#            readData=True
#        elif(ValidateData.isValidStructuredArray(tempReadData)):
#            readData=False
#        else:
#            print("".join(["PCoA:loadData::Error tempReadData was not an existing file or an np array."]))
#            return False
#
#        #If a file path is given, read data.
#        if(ValidateData.isTrue(readData)):
#            #Object to read in and manipulate raw data
#            rawData = AbundanceTable()
#            #If indicated, check the input file
#            if(ValidateData.isTrue(tempCheckFile)):
#                tempReadData = rawData.checkRawDataFile(tempReadData)
#                if(ValidateData.isFalse(tempReadData)):
#                    print("".join(["PCoA:loadData::Error when checking raw data file, did not perform PCoA. File name:",str(tempReadData)]))
#                    return False
#            #Read in the file data to a numpy array.
#            #Samples (column) by Taxa (rows)(lists) without the column
#            data = rawData.textToArray(tempInputFile=tempReadData, tempDelimiter=tempDelimiter, tempNameRow=tempNameRow, tempFirstDataRow=tempFirstDataRow, tempNormalize=tempNormalize)
#            if(ValidateData.isFalse(data)):
#                print("PCoA:loadData::Error when reading checked raw data file, did not perform PCoA.")
#                return False
#
#            #Transpose data to be Taxa (columns) by samples (rows)(lists)
#            data = rawData.transposeDataMatrix(tempMatrix=data[0], tempRemoveAdornments=False)
#            if(ValidateData.isFalse(data)):
#                print("PCoA:loadData::Error when transposing read raw data file, did not perform PCoA.")
#                return False
#            else:
#                self.dataMatrix=data
#                self.isRawData=tempIsRawData
#        #Otherwise load the data directly as passed.
#        else:
#            self.dataMatrix=tempReadData
#            self.isRawData=tempIsRawData
#        return True

    #Loads data into PCoA (given the matrix or a valid file path)
    #Data can be the Abundance Table to be converted to a distance matrix or a distance matrix
    #If it is the AbundanceTable, indicate that it is rawData (tempIsRawData=True)
    #If it is the distance matrix already generated indicate (tempIsRawData=False)
    #and no conversion will occur in subsequent methods
    #@params tempReadData Either a Structured matrix of data or a valid file path to read from
    #@params xData Abundancetable or Distance matrix . Taxa (columns) by samples (rows)(lists)
    #@return Return boolean indicator of success (True=Was able to load data)
    def loadData(self, xData, fIsRawData):

        if fIsRawData:

            #Read in the file data to a numpy array.
            #Samples (column) by Taxa (rows)(lists) without the column
            data = xData.funcToArray()
            if data==None:
                print("PCoA:loadData::Error when converting AbundanceTable to Array, did not perform PCoA.")
                return False

            #Transpose data to be Taxa (columns) by samples (rows)(lists)
            data = Utility_Math.transposeDataMatrix(data,tempRemoveAdornments=False)
            if(ValidateData.isFalse(data)):
                print("PCoA:loadData::Error when transposing data file, did not perform PCoA.")
                return False
            else:
                self.dataMatrix=data
                self.isRawData=fIsRawData

        #Otherwise load the data directly as passed.
        else:
            self.dataMatrix=tempReadData
            self.isRawData=fIsRawData
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
    def plot(self,tempPlotName="PCOA.png", tempColorGrouping='g', tempShape='o', tempLabels=["Green"], tempShapeLabels=["Circle"], tempShapeSize = 20, tempAlpha = 1.0, tempLegendLocation="upper right", tempInvert=False):

        if(not self.pcoa == None):
            #Get point count
            adPoints = self.pcoa.getPoints()
            ldXPoints = list(adPoints[:,0])
            ldYPoints = list(adPoints[:,1])
            iPointCount = len(ldXPoints)

            #Check shapes
            if(ValidateData.isValidList(tempShape)):
              if not len(tempShape) == iPointCount:
                print("Error, the list of shapes was given but was not the same size as the points so nothing was plotted.")
                print("tempShape")
                print(tempShape)
                print("iPointCount")
                print(iPointCount)
                print("len adPoints")
                print(iPointCount)
                return

            #Check colors
            if(ValidateData.isValidList(tempColorGrouping)):
              if not len(tempColorGrouping) == iPointCount:
                print("Error, the list of colors was given but was not the same size as the points so nothing was plotted.")
                print("tempColorGrouping")
                print(tempColorGrouping)
                print("iPointCount")
                print(iPointCount)
                print("len adPoints")
                print(iPointCount)
                return

            #Check sizes
            if(ValidateData.isValidList(tempShapeSize)):
              if not len(tempShapeSize) == iPointCount:
                print("Error, the list of sizes was given but was not the same size as the points so nothing was plotted.")
                print("tempShapeSize")
                print(tempShapeSize)
                print("iPointCount")
                print(iPointCount)
                print("len adPoints")
                print(iPointCount)
                return

            #Get plot object
            imgFigure = plt.figure()

            self.objFigureControl.invertColors(fInvert=tempInvert)

            #Color/Invert figure
            imgFigure.set_facecolor(self.objFigureControl.c_strBackgroundColorWord)
            imgSubplot = imgFigure.add_subplot(111,axisbg=self.objFigureControl.c_strBackgroundColorLetter)
            imgSubplot.set_xlabel("Dimension 1")
            imgSubplot.set_ylabel("Dimension 2")
            imgSubplot.spines['top'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['bottom'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['left'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['right'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.xaxis.label.set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.yaxis.label.set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='x', colors=self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='y', colors=self.objFigureControl.c_strDetailsColorLetter)
            charMarkerEdgeColor = self.objFigureControl.c_strDetailsColorLetter

            #If given a list of colors, each color will be plotted individually stratified by shape
            #Plot colors seperately so the legend will pick up on the labels and make a legend
            if(ValidateData.isValidList(tempColorGrouping)):
                if len(tempColorGrouping) == iPointCount:
                    #Check for lists in the list which indicate the need to plot pie charts
                    lfAreLists = [ValidateData.isValidList(objColor) for objIndex, objColor in enumerate(tempColorGrouping)]

                    #Pie chart data seperated out
                    lsColorsPieCharts = None
                    lcShapesPieCharts = None
                    lsLabelsPieCharts = None
                    lsSizesPieCharts = None
                    ldXPointsPieCharts = None
                    ldYPointsPieCharts = None
                    #Split out piechart data
                    if sum(lfAreLists) > 0:
                        #Get lists of index that are and are not lists
                        liAreLists = []
                        liAreNotLists = []
                        curIndex = 0
                        for fIsList in lfAreLists:
                            if fIsList: liAreLists.append(curIndex)
                            else: liAreNotLists.append(curIndex)
                            curIndex = curIndex + 1

                        lsColorsPieCharts = self.reduceList(tempColorGrouping, liAreLists)
                        tempColorGrouping = self.reduceList(tempColorGrouping, liAreNotLists)
                        #Split out shapes
                        if ValidateData.isValidList(tempShape):
                            lcShapesPieCharts = self.reduceList(tempShape, liAreLists)
                            tempShape = self.reduceList(tempShape, liAreNotLists)
                        else:
                            lcShapesPieCharts = tempShape
                        #Split out labels
                        if ValidateData.isValidList(tempLabels):
                            lsLabelsPieCharts = self.reduceList(tempLabels, liAreLists)
                            tempLabels = self.reduceList(tempLabels, liAreNotLists)
                        else:
                            lsLabelsPieCharts = tempLabels
                        #Split out sizes
                        if ValidateData.isValidList(tempShapeSize):
                            lsSizesPieCharts = self.reduceList(tempShapeSize, liAreLists)
                            tempShapeSize = self.reduceList(tempShapeSize, liAreNotLists)
                        else:
                            lsSizesPieCharts = tempShapeSize
                        #Split out xpoints
                        if ValidateData.isValidList(ldXPoints):
                            ldXPointsPieCharts = self.reduceList(ldXPoints, liAreLists)
                            ldXPoints = self.reduceList(ldXPoints, liAreNotLists)
                        else:
                            ldXPointsPieCharts = ldXPoints
                        #Split out ypoints
                        if ValidateData.isValidList(ldYPoints):
                            ldYPointsPieCharts = self.reduceList(ldYPoints, liAreLists)
                            ldYPoints = self.reduceList(ldYPoints, liAreNotLists)
                        else:
                            ldYPointsPieCharts = ldYPoints

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
                        aiXPoints = self.reduceList(ldXPoints,aiColorPointPositions)
                        aiYPoints = self.reduceList(ldYPoints,aiColorPointPositions)

                        #There are 3 options for shapes which are checked in this order.
                        #1. 1 shape character is given which is used for all markers
                        #2. A list is given of marker characters or lists of decimals which will be used to make pie chart markers
                        #This is handled after the rest this block of code
                        #3. A list of char are given each indicating the marker for a sample
                        #If the shapes are not a list plot
                        #Otherwise plot per shape per color (can not plot list of shapes in matplotlib)
                        reducedShapes = tempShape
                        if(not ValidateData.isValidList(reducedShapes)):
                          reducedShapes = reducedShapes[0]
                          imgSubplot.scatter(aiXPoints,aiYPoints, c=[charColor], marker=reducedShapes, alpha=tempAlpha, label=tempLabels[tempColorGrouping.index(charColor)], s=reducedSizes, edgecolor=charMarkerEdgeColor)
                        #Shapes are supplied as a list so plot each shape
                        else:
                          #Reduce to shapes of the current colors
                          reducedShapes = self.reduceList(reducedShapes,aiColorPointPositions)
                          acharReducedShapesElements = list(set(reducedShapes))
                          #If there are multiple shapes, plot seperately because one is not allowed to plot them as a list
                          for aCharShapeElement in acharReducedShapesElements:
                            #Get indices
                            aiShapeIndices = self.getIndices(reducedShapes,aCharShapeElement)
                            #Reduce label by shapes
                            strShapeLabel = self.reduceList(acharLabelsByColor,aiShapeIndices)
                            #Reduce sizes by shapes
                            strShapeSizes = reducedSizes
                            if ValidateData.isValidList(reducedSizes):
                              strShapeSizes = self.reduceList(reducedSizes,aiShapeIndices)
                            #Get points per shape
                            aiXPointsPerShape = self.reduceList(aiXPoints,aiShapeIndices)
                            aiYPointsPerShape = self.reduceList(aiYPoints,aiShapeIndices)
                            #Get sizes per shape
                            #Reduce sizes if a list
                            reducedSizesPerShape = reducedSizes
                            if(ValidateData.isValidList(reducedSizes)):
                              reducedSizesPerShape = self.reduceList(reducedSizes,aiShapeIndices)
                            #Plot
                            imgSubplot.scatter(aiXPointsPerShape,aiYPointsPerShape, c=[charColor], marker=aCharShapeElement, alpha=tempAlpha, label=strShapeLabel[0], s=strShapeSizes, edgecolor=charMarkerEdgeColor)

                    #Plot pie charts
                    if not lsColorsPieCharts == None:
                        self.plotWithPieMarkers(imgSubplot=imgSubplot, aiXPoints=ldXPointsPieCharts, aiYPoints=ldYPointsPieCharts, dSize=lsSizesPieCharts, llColors=lsColorsPieCharts, lsLabels=lsLabelsPieCharts, lcShapes=lcShapesPieCharts, edgeColor=charMarkerEdgeColor, dAlpha=tempAlpha)

            #If the Color is not a list but shapes are a list then plot by each unique shape
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
                        imgSubplot.scatter(self.reduceList(ldXPoints,aiShapePointPositions),self.reduceList(ldYPoints,aiShapePointPositions), s=reducedSizes, c=[tempColorGrouping], marker=acharUniqueShapes[iShapeIndex], alpha=dAlpha, label=tempLabels[tempShape.index(acharUniqueShapes[iShapeIndex])], edgecolor=charMarkerEdgeColor)

            #This is the simple case where color and shape are constant and do not change.
            else:
                imgSubplot.scatter(ldXPoints,ldYPoints, c=tempColorGrouping, marker=tempShape, alpha=tempAlpha, label=tempLabels, s=self.objFigureControl.iMarkerSize, edgecolor=charMarkerEdgeColor)

            objLegend = imgSubplot.legend(loc=tempLegendLocation, scatterpoints=1, prop={'size':10})

            #Invert legend
            if(tempInvert):
              if objLegend:
                objLegend.legendPatch.set_fc(self.objFigureControl.c_strBackgroundColorWord)
                objLegend.legendPatch.set_ec(self.objFigureControl.c_strDetailsColorLetter)
                plt.setp(objLegend.get_texts(),color=self.objFigureControl.c_strDetailsColorLetter)

            #Make legend background transparent
            if objLegend:
              objLegendFrame = objLegend.get_frame()
              objLegendFrame.set_alpha(self.objFigureControl.c_dAlpha)

            imgFigure.savefig(tempPlotName, facecolor=imgFigure.get_facecolor())

    #The all lists should be in the same order
    #aiXPoints List of X axis points (one element per color list)
    #aiYPoints List of X axis points (one element per color list)
    #dSize double or List of doubles (one element per color list)
    #llColors List of Lists of colors, one list of colors is for 1 piechart/multiply highlighted feature
    #Example ["red","blue","green"] for a marker with 3 sections
    #lsLabels List of labels  (one element per color list)
    #lcShapes indicates whihc shape of a pie chart to use, currently supported 'o' and 's'  (one element per color list)
    #edgeColor One color entry for the edge of the piechart
    def plotWithPieMarkers(self, imgSubplot, aiXPoints, aiYPoints, dSize, llColors, lsLabels, lcShapes, edgeColor, dAlpha):
        #Zip up points to pairs
        xyPoints = zip(aiXPoints,aiYPoints)
        #For each pair of points
        for iIndex,dXY in enumerate(xyPoints):
            ldWedges = []
            #Get colors
            lcurColors = llColors[iIndex]
            #Get pie cut shape
            cPieChartType = lcShapes[iIndex]
            if cPieChartType == Constants_Figures().c_charPCOAPieChart:
                ldWedges = self.makePieWedges(len(lcurColors),20)
            elif cPieChartType == Constants_Figures().c_charPCOASquarePieChart:
                ldWedges = self.makeSquarePieWedges(len(lcurColors))
            for iWedgeIndex,dWedge in enumerate(ldWedges):
                imgSubplot.scatter(x=dXY[0], y=dXY[1], marker=(dWedge,0), s=dSize[iIndex], label=lsLabels[iIndex], facecolor=lcurColors[iWedgeIndex], edgecolor=edgeColor, alpha=dAlpha)

    #iWedgecount Number of equal wedge sizes to make
    #iSplineResolution The higher the number the better smoothing for the circle parameter
    def makePieWedges(self, iWedgeCount,iSplineResolution = 10):
        ldWedge = []
        dLastValue = 0.0

        #Create a list of equal percentages for all wedges
        #Do not include a last wedge it gets all the space from the 2nd to last wedge to the end
        #Which should still be equal to the others
        ldPercentages = [1.0/iWedgeCount]*(iWedgeCount-1)

        for dPercentage in ldPercentages:
            ldX = [0] + np.cos(np.linspace(2*math.pi*dLastValue,2*math.pi*(dLastValue+dPercentage),iSplineResolution)).tolist()
            ldY = [0] + np.sin(np.linspace(2*math.pi*dLastValue,2*math.pi*(dLastValue+dPercentage),iSplineResolution)).tolist()
            ldWedge.append(zip(ldX,ldY))
            dLastValue = dLastValue+dPercentage
        ldX = [0] + np.cos(np.linspace(2*math.pi*dLastValue,2*math.pi,iSplineResolution)).tolist()
        ldY = [0] + np.sin(np.linspace(2*math.pi*dLastValue,2*math.pi,iSplineResolution)).tolist()
        ldWedge.append(zip(ldX,ldY))
        return ldWedge

    #iWedgecount Number of equal wedge sizes to make
    #iSplineResolution The higher the number the better smoothing for the circle parameter
    def makeSquarePieWedges(self, iWedgeCount):
        ldWedge = []
        dLastPercentageValue = 0.0
        dLastSquareValue = 0.0
        dCumulativePercentageValue = 0.0
        dRadius = None
        fXYSwitched = False
        fAfterCorner = False
        iSwitchCounts = 0
        iMagicNumber =(1.0/4)

        #Create a list of equal percentages for all wedges
        #Do not include a last wedge it gets all the space from the 2nd to last wedge to the end
        #Which should still be equal to the others
        ldPercentages = [1.0/iWedgeCount]*(iWedgeCount)

        for dPercentage in ldPercentages:
          ldCircleXs = np.cos([2*math.pi*dLastPercentageValue,2*math.pi*(dLastPercentageValue+dPercentage)])
          ldCircleYs = np.sin([2*math.pi*dLastPercentageValue,2*math.pi*(dLastPercentageValue+dPercentage)])

          if dRadius == None:
            dRadius = ldCircleXs[0]

          #Check to see if at corner
          fAtCorner = False
          iDistance = int((dLastPercentageValue+dPercentage+(iMagicNumber/2))/iMagicNumber
                  ) - int((dLastPercentageValue+(iMagicNumber/2))/iMagicNumber)
          if(iDistance > 0):
            fAtCorner = True
            if iDistance > 1:
              fXYSwitched = not fXYSwitched
              iSwitchCounts = iSwitchCounts + 1

          #Check to see if at a side center
          fAtSide = False
          if (int((dLastPercentageValue+dPercentage)/iMagicNumber) > int(dLastPercentageValue/iMagicNumber)):
            fAtSide = True

          #Handle corner xy switching
          if fAtCorner:
            fXYSwitched = not fXYSwitched
            iSwitchCounts = iSwitchCounts + 1
          #Make sure the xy switching occurs to vary the slope at the corner.
          if fXYSwitched:
              ldCircleXs,ldCircleYs = ldCircleYs,ldCircleXs

          dSquarePoint = dRadius * (ldCircleYs[1]/float(ldCircleXs[1]))
          dRadiusSq1 = dRadius
          dRadiusSq2 = dRadius
          dLastSquareValueSq = dLastSquareValue
          dSquarePointSq = dSquarePoint

          #If in quadrants 2,3 make sign changes
          if iSwitchCounts in [2,3]:
            if iSwitchCounts == 2:
              dRadiusSq1 = dRadiusSq1 *-1
            elif iSwitchCounts == 3:
              dRadiusSq1 = dRadiusSq1 * -1
              dRadiusSq2 = dRadiusSq2 * -1
            dLastSquareValueSq = dLastSquareValueSq * -1.0
            dSquarePointSq = dSquarePointSq * -1.0

          if fAtCorner:
            #Corner 1
            if iSwitchCounts==1:
              ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,dLastSquareValueSq,dRadiusSq2,dRadiusSq2,0]))
            #Corner 2
            elif iSwitchCounts==2:
              if iDistance > 1:
                ldWedge.append(zip([0,-dRadiusSq1,-dRadiusSq1,dRadiusSq1,dRadiusSq1,0],[0,-dLastSquareValueSq,dRadiusSq2,dRadiusSq2,dSquarePointSq,0]))
              else:
                ldWedge.append(zip([0,-dLastSquareValueSq,dRadiusSq1,dRadiusSq1,0],[0,dRadiusSq2,dRadiusSq2,dSquarePointSq,0]))
            #Corner 3
            elif iSwitchCounts==3:
              if iDistance > 1:
                ldWedge.append(zip([0,-dLastSquareValueSq,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,-dRadiusSq2,-dRadiusSq2,dRadiusSq2,dRadiusSq2,0]))
              else:
                ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,dLastSquareValueSq,dRadiusSq2,dRadiusSq2,0]))
            #Corner 4
            elif iSwitchCounts==4:
              if iDistance > 1:
                ldWedge.append(zip([0,-dRadiusSq1,-dRadiusSq1,dRadiusSq1,dRadiusSq1,0],[0,-dLastSquareValueSq,-dRadiusSq2,-dRadiusSq2,dSquarePointSq,0]))
              else:
                ldWedge.append(zip([0,(-1.0*dLastSquareValueSq),dRadiusSq1,dRadiusSq1,0],[0,(-1.0*dRadiusSq2),(-1.0*dRadiusSq2),dSquarePointSq,0]))

            fAfterCorner = True
          else:
            if iSwitchCounts%2:
              ldWedge.append(zip([0,dLastSquareValueSq,dSquarePointSq,0],[0,dRadiusSq2,dRadiusSq2,0]))
            else:
              ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,0],[0,dLastSquareValueSq,dSquarePointSq,0]))

          dLastSquareValue = dSquarePoint
          dCumulativePercentageValue = dCumulativePercentageValue + dLastSquareValue
          dLastPercentageValue = dLastPercentageValue+dPercentage

        return ldWedge

    #charForceColor if set, automatic coloring will not occur but will occur based on the charForceColor
    #CharForceColor should be a list of selection methods or not selected which will be automatically broken into colors
    #charForceColor must contain the same elements as the lsLabeList
    #Currently can be a list (1 color per marker in the order of the data), or 1 char to force all markers to
    #charForceShapes if set, automatic shapes will not occur
    #Currently can only be a char (forcing effects all markers equally)
    def plotList(self,lsLabelList, strOutputFileName, iSize, dAlpha = 1.0, charForceColor=None, charForceShape=None, fInvert=False):

        #Get uniqueValues for labels
        acharUniqueValues = list(set(lsLabelList))
        iCountUniqueValues = len(acharUniqueValues)

        #Set colors
        atupldLabelColors = None

        #Set shapes
        alLabelShapes = None
        if charForceShape == None:
            #Get shapes
            acharShapes = PCoA.getShapes(iCountUniqueValues)
            if acharShapes == None:
                return False
            #Make label shapes
            alLabelShapes = [ acharShapes[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
        else:
            alLabelShapes = charForceShape

        #If the coloring is not forced, color so it is based on the labels
        if charForceColor == None:
            #Get colors based on labels
            atupldColors = [PCoA.RGBToHex(cm.jet(float(iUniqueValueIndex)/float(iCountUniqueValues))) for iUniqueValueIndex in xrange(0,iCountUniqueValues)]
            #Make sure generated colors are unique
            if not iCountUniqueValues == len(set(atupldColors)):
                logging.error("PCoA.plotList. Generated colors were not unique for each unique label value.")
                logging.error("Labels")
                logging.error(lsLabelList)
                logging.error(len(lsLabelList))
                logging.error("Unique Labels")
                logging.error(set(lsLabelList))
                logging.error(len(set(lsLabelList)))
                logging.error("Colors")
                logging.error(atupldColors)
                logging.error(len(atupldColors))
                logging.error("Unique Colors")
                logging.error(set(atupldColors))
                logging.error(len(set(atupldColors)))
                return False
            #Make label coloring
            atupldLabelColors = [ atupldColors[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
        #If the coloring is forced, color so it is based on the charForcedColor list
        elif(ValidateData.isValidList(charForceColor)):
            atupldLabelColors = charForceColor[0]
            if not len(lsLabelList) == len(atupldLabelColors):
                logging.error("PCoA.plotList. Label and forced color lengths were not the same.")
                logging.error("Labels")
                logging.error(lsLabelList)
                logging.error(len(lsLabelList))
                logging.error("Forced Colors")
                logging.error(charForceColor[0])
                logging.error(len(charForceColor[0]))
                return False
            lsLabelList = [ "".join([charForceColor[1][iLabelIndex], "_", lsLabelList[iLabelIndex]]) for iLabelIndex in xrange(0,len(charForceColor[1]))]
        #If the color is forced but the color does not vary, color all markers are the same.
        else:
            atupldLabelColors = charForceColor

        logging.debug("lsLabelList")
        logging.debug(lsLabelList)
        self.plot(tempPlotName=strOutputFileName, tempColorGrouping=atupldLabelColors, tempShape=alLabelShapes, tempLabels=lsLabelList, tempShapeSize = iSize, tempAlpha=dAlpha, tempInvert = fInvert)


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
      #Shapes for the data (typically used when stratifying)
      #Currently pie cut objects are just o and s
      #TODO extend pie cut objects to all convex polygons derived from straight lines
      lsPointShapes = ['o','s','^','v','<','>','8','p','h']
      if intShapeCount > len(lsPointShapes):
        print("".join(["Error, PCoA.getShapes. Do not have enough shapes to give. Received request for ",str(intShapeCount)," shapes. Max available shape count is ",str(len(lsPointShapes)),"."]))
        return None
      return lsPointShapes[0:intShapeCount]
