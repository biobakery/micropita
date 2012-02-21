#######################################################
# Author: Timothy Tickle
# Description: Class to test the PCoA class
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libraries
from Constants import Constants
from numpy import array
from PCoA import PCoA
import os
import unittest

##
#Tests the Blog object
class PCoATest(unittest.TestCase):

    def testPCoAFromArray(self):

        normalize = False
        checkFile = False

        #Inputs
        aData = array([[0,3,4,8],
                [3,0,1,27],
                [4,1,0,3.5],
                [8,27,3.5,0]],'d')

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=aData, tempIsRawData=True, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)
        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAFromArray.png")

        #Generate result
        result = ""

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))


    def testPCoAFromPorportional(self):

        #Inputs
        inputFile = "./testData/PCoA/TestDataGeneration_Porportional.txt"
        checkedFile = "./testData/PCoA/TestDataGeneration_Porportional-checked.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        checkFile = True

        acharColors = ["g","g","g","g","g","g","g","g","r","r","r","r","r","r","r","b","b","b","b","b","b","b","y","y"]
        cShape = "o"

        #Remove checked file if it exists before the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=inputFile, tempIsRawData=True, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)
        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAFromPorportional.png", tempColorGrouping=acharColors, tempShape=cShape)

        #Generate result
        result = ""

        #Remove checked file if it exists after the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testPCoAFromAbsolute(self):

        #Inputs
        inputFile = "./testData/PCoA/TestDataGeneration.txt"
        checkedFile = "./testData/PCoA/TestDataGeneration-checked.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        checkFile = True

        acharColors = ["g","g","g","g","g","g","g","g","r","r","r","r","r","r","r","b","b","b","b","b","b","b","y","y"]
        cShape = "o"

        #Remove checked file if it exists before the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=inputFile, tempIsRawData=True, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAFromAbsolute.png", tempColorGrouping=acharColors, tempShape=cShape)

        #Generate result
        result = ""

        #Remove checked file if it exists after the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testPCoAForListOfShapes(self):

        #Inputs
        inputFile = "./testData/PCoA/TestDataGeneration.txt"
        checkedFile = "./testData/PCoA/TestDataGeneration-checked.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        checkFile = True

        acharColors = ["g","g","g","g","g","g","g","g","r","r","r","r","r","r","r","b","b","b","b","b","b","b","y","y"]
        cShape = ['o','o','o','o','o','o','o','o','x','x','x','x','x','x','x','+','+','+','+','+','+','+','<','<']

        #Remove checked file if it exists before the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=inputFile, tempIsRawData=True, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAForListOfShapes.png", tempColorGrouping=acharColors, tempShape=cShape)

        #Generate result
        result = ""

        #Remove checked file if it exists after the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testPCoAForSelectTaxa(self):

        #Inputs
        inputFile = "./testData/PCoA/TestDataGeneration.txt"
        checkedFile = "./testData/PCoA/TestDataGeneration-checked.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        checkFile = True

        acharColors = ["g","g","g","g","g","g","g","g","y","g","y","g","g","y","g","y","g","y","g","g","y","g","y","y"]
        cShape = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']

        #Remove checked file if it exists before the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=inputFile, tempIsRawData=True, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAForSelectTaxa.png", tempColorGrouping=acharColors, tempShape=cShape)

        #Generate result
        result = ""

        #Remove checked file if it exists after the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testPCoAForSVMSamples(self):

        #Inputs
        inputFile = "./testData/PCoA/TestDataGeneration.txt"
        checkedFile = "./testData/PCoA/TestDataGeneration-checked.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        checkFile = True

        acharColors = ["g","y","y","g","g","y","g","y","y","g","y","g","g","y","g","y","g","y","g","g","y","g","y","y"]
        cShape = ['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']

        #Remove checked file if it exists before the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Correct Answer
        answer = ""

        #Analysis object
        analysis = PCoA()
        #LoadData
        analysis.loadData(tempReadData=inputFile, tempIsRawData=True, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize, tempCheckFile=checkFile)
        #Run Analysis
        pcoaResults = analysis.run(tempDistanceMetric=analysis.c_BRAY_CURTIS)

        #Plot
        analysis.plot(tempPlotName="./testData/PCoA/testPCoAForSVMSamples.png", tempColorGrouping=acharColors, tempShape=cShape)

        #Generate result
        result = ""

        #Remove checked file if it exists after the test
        if(os.path.exists(checkedFile)):
            os.remove(checkedFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(PCoATest)
    return suite
