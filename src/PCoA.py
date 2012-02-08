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
    def plot(self,tempPlotName="PCOA.png", tempColorGrouping='g', tempShape='o', tempColorLabels=["Green"], tempShapeLabels=["Circle"], tempLegendLocation="upper right"):
        if(not self.pcoa == None):
            print("tempColorLabels")
            print(tempColorLabels)
            adPoints = self.pcoa.getPoints()
            iPointCount = len(adPoints[:,0])
            imgFigure = plt.figure()
            imgSubplot = imgFigure.add_subplot(111)
            #Plot colors seperately so the legend will pick up on the labels and make a legend
            if(ValidateData.isValidList(tempColorGrouping)):
                if len(tempColorGrouping) == iPointCount:
                    #Get unique colors and plot each individually
                    acharUniqueColors = list(set(tempColorGrouping))
                    for iColorIndex in xrange(0,len(acharUniqueColors)):
                        aiColorPointPositions = self.getIndices(tempColorGrouping,acharUniqueColors[iColorIndex])
                        imgSubplot.scatter(self.reduceList(adPoints[:,0],aiColorPointPositions),self.reduceList(adPoints[:,1],aiColorPointPositions), c=[acharUniqueColors[iColorIndex]], marker=tempShape, label=tempColorLabels[tempColorGrouping.index(acharUniqueColors[iColorIndex])])
            else:
                imgSubplot.scatter(adPoints[:,0],adPoints[:,1], c=tempColorGrouping, marker=tempShape, label=tempColorLabels)
#            if(ValidateData.isValidList(tempShape)):
#                if(len(tempShape) == len(adPoints[:,0])):
#                    for iShapeIndex in xrange(0,len(tempShape)):
#                        imgSubplot.scatter(adPoints[:,0][iShapeIndex],adPoints[:,1][iShapeIndex], c=tempColorGrouping[iShapeIndex], marker=tempShape[iShapeIndex])
#            else:
#                imgSubplot.scatter(adPoints[:,0],adPoints[:,1], c=tempColorGrouping, marker=tempShape)
            imgSubplot.legend(loc=tempLegendLocation, scatterpoints=1, prop={'size':10})
            imgFigure.savefig(tempPlotName)

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

