#######################################################
#
#	Title:		CladogramTest
#	Author:		Timothy Tickle 
#	Date:		February 13, 2012
#	Purpose:	Test Cladogram class
#
#######################################################

#Import local code
from AbundanceTable import AbundanceTable
import os
import unittest
from Cladogram import Cladogram
from Constants_Figures import Constants_Figures
from Constants_Testing import Constants_Testing

class CladogramTest(unittest.TestCase):

    ##Set up for tests
    def setUp(self): pass

    def testFormatRGBForBlack(self):
        methodName = "testFormatRGBForBlack"

        tempRGB = "0,0,0"
        correctAnswer = "_c_[0.0,0.0,0.0]"
        result = Cladogram().formatRGB(tempRGB)
        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFormatRGBForWhite(self):
        methodName = "testFormatRGBForWhite"

        tempRGB = "255,255,255"
        correctAnswer = "_c_[1.0,1.0,1.0]"
        result = Cladogram().formatRGB(tempRGB)
        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFormatRGBForMid(self):
        methodName = "testFormatRGBForMid"

        tempRGB = "128,128,128"
        correctAnswer = "_c_[0.501960784314,0.501960784314,0.501960784314]"
        result = Cladogram().formatRGB(tempRGB)
        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFilterByCladeSizeFor0Removed(self):
        methodName = "testFilterByCladeSizeFor0Removed"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A"]
        correctAnswer = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A","1A|2A|3A|4A|5A","1A|2A|3A|4A|6A"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 2
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        result = cladogram.filterByCladeSize(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    # Test filter by size
    def testFilterByCladeSizeForAllFiltering(self):
        methodName = "testFilterByCladeSizeForAllFiltering"

        lsIds = ["A|B|C","B|C","B|D","A|B","B|E"]
        answerFileContents = []

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        lsRet = cladogram.filterByCladeSize(lsIds)

        self.assertEquals(answerFileContents,lsRet, methodName+"Did not get the correct file, expected:"+str(answerFileContents)+" but received:"+str(lsRet)+".")

    def testFilterByCladeSizeForNoFiltering(self):
        methodName = "testFilterByCladeSizeForNoFiltering"

        lsIds = ["A|B|C|1","A|B|2","A|B|C|3","A|B|4","A|B|C|5"]
        answerFileContents = lsIds

        #Filter settings
        fFilter = True
        iMeasureLevel = 2
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        lsRet = cladogram.filterByCladeSize(lsIds)

        self.assertEquals(answerFileContents,lsRet, methodName+"Did not get the correct file, expected:"+str(answerFileContents)+" but received:"+str(lsRet)+".")

    def testFilterByCladeSizeForSomeFiltering(self):
        methodName = "testFilterByCladeSizeSomeFiltering"

        lsIds = ["A|B|C|1","A2|B|C|6","A|B|2","A3|B|C|7","A3|B|C|6","A|B|C|3","A|B|4","A4|B|C|6","A|B|C|5"]
        answerFileContents = ["A|B|C|1","A|B|2","A|B|C|3","A|B|4","A|B|C|5"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 2
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        lsRet = cladogram.filterByCladeSize(lsIds)

        self.assertEquals(answerFileContents,lsRet, methodName+"Did not get the correct file, expected:"+str(answerFileContents)+" but received:"+str(lsRet)+".")

    def testFilterByCladeSizeForSomeFilteringUnclassified(self):
        methodName = "testFilterByCladeSizeSomeFilteringUnclassified"

        lsIds = ["A|B|C|1","A2|B|C|6","A|B|2","A3|B|C|7","A3|B|C|6","A|B|unclassified","A|B|4","A4|B|unclassified","A|B|C|5"]
        answerFileContents = ["A|B|C|1","A|B|2","A|B|unclassified","A|B|4","A|B|C|5"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 2
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize, strUnclassified="unclassified")

        #Run to get result
        lsRet = cladogram.filterByCladeSize(lsIds)

        self.assertEquals(answerFileContents,lsRet, methodName+"Did not get the correct file, expected:"+str(answerFileContents)+" but received:"+str(lsRet)+".")

    def testFilterByCladeSizeFor1Removed(self):
        methodName = "testFilterByCladeSizeFor1Removed"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B"]
        correctAnswer = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A","1A|2A|3A|4A|5A","1A|2A|3A|4A|6A"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        result = cladogram.filterByCladeSize(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFilterByCladeSizeFor2Removed(self):
        methodName = "testFilterByCladeSizeFor2Removed"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B"]
        correctAnswer = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A","1A|2A|3A|4A|5A","1A|2A|3A|4A|6A"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        result = cladogram.filterByCladeSize(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFilterByCladeSizeFor3Removed(self):
        methodName = "testFilterByCladeSizeFor3Removed"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B","1C|2C|3C|4C|3C"]
        correctAnswer = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A","1A|2A|3A|4A|5A","1A|2A|3A|4A|6A"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        result = cladogram.filterByCladeSize(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFilterByCladeSizeFor4Removed(self):
        methodName = "testFilterByCladeSizeFor3Removed"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                        "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]
        correctAnswer = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B"]

        #Filter settings
        fFilter = True
        iMeasureLevel = 5
        iReduceLevel = 1
        iCladeSize = 5
        cladogram = Cladogram()
        cladogram.setFilterByCladeSize(fCladeSizeFilter=fFilter, iCladeLevelToMeasure = iMeasureLevel,
                                       iCladeLevelToReduce = iReduceLevel, iMinimumCladeSize = iCladeSize)

        #Run to get result
        result = cladogram.filterByCladeSize(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testUpdateToRootForLevel0(self):
        methodName = "testUpdateToRootForLevel0"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                        "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]
        correctAnswer = ["2A|3A|4A|1A","2A|3A|4A|2A","2A|3A|4A|3A","2A|3A|4A|4A",
                      "2A|3A|4A|5A","2A|3A|4A|6A"]
        cladogram = Cladogram()
        cladogram.forceRoot("1A")
        result = cladogram.updateToRoot(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testUpdateToRootForLevel1(self):
        methodName = "testUpdateToRootForLevel1"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                        "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]
        correctAnswer = ["3B|4B|1B","3B|4B|2B","3B|4B|3B","3B|4B|4B","3B|4B|5B","3B|4B|6B","3B|4B|7B"]
        cladogram = Cladogram()
        cladogram.forceRoot("2B")
        result = cladogram.updateToRoot(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testUpdateToRootForLevel4(self):
        methodName = "testUpdateToRootForLevel4"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]
        correctAnswer = ["1C"]
        cladogram = Cladogram()
        cladogram.forceRoot("4C")
        result = cladogram.updateToRoot(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testUpdateToRootForNone(self):
        methodName = "testUpdateToRootForNone"

        tempClades = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                        "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]
        correctAnswer = []
        cladogram = Cladogram()
        cladogram.forceRoot("1F")
        result = cladogram.updateToRoot(tempClades)
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnTrue(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnTrue"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder = 0.45
        strTestShape = "v"
        dTestAlpha = 0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = True
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnFalseLengthBorder(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnFalseLengthBorder"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=[0.45,0.0]
        strTestShape="v"
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = False
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnFalseLengthShape(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnFalseLengthBorder"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape=["v","v","v"]
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = False
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnFalseLengthAlpha(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnFalseLengthAlpha"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape="v"
        dTestAlpha=[0.5,0.0,0.0,0.0,0.0]

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = False
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnTrueLengthBorder(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnTrueLengthBorder"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=[0.16,0.15,0.14,0.13,0.12,0.11,0.10,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
        strTestShape="v"
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = True
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnTrueLengthShape(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnTrueLengthBorder"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape=["v","v","v","v","v","v","v","v","v","v","v","v","v","v","v","v"]
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = True
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForReturnTrueLengthAlpha(self):
        methodName = "testAddingAndGeneratingCircleDataForReturnTrueLengthAlpha"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape="v"
        dTestAlpha=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.10,0.11,0.12,0.13,0.14,0.15]

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)

        result = cladogram.createCircleFile(lsTestTaxa)
        correctAnswer = True
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthAlpha(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthAlpha"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape="v"
        dTestAlpha=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.10,0.11,0.12,0.13,0.14,0.15]

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(lsTestTaxa)

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthBorder(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthBorder"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=[0.16,0.15,0.14,0.13,0.12,0.11,0.10,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
        strTestShape="v"
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(lsTestTaxa)

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthShape(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthShape"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape=["v","^","^","v","v","^","^","v","v","^","v","^","v","^","v","^"]
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(lsTestTaxa)

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthShapeForcedRoot(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthShapeForcedRoot"

        lsTestTaxa = ["1A|2A|3A|4A|1","1A|2A|3A|4A|2","1A|2A|3A|4A|3","1A|2A|3A|4A|4",
                      "1A|2A|3A|4A|5","1A|2A|3A|4A|6","1B|2B|3B|4B|1","1B|2B|3B|4B|2",
                      "1B|2B|3B|4B|3","1B|2B|3B|4B|4","1B|2B|3B|4B|5","1B|2B|3B|4B|6","1B|2B|3B|4B|7",
                      "1C|2C|3C|4C|1","1D|2D|3D|4D|1","1E|2E|3E|4E|1"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape=["v","^","^","v","v","^","^","v","v","^","v","^","v","^","v","^"]
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.forceRoot("2A")
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(cladogram.updateToRoot(lsTestTaxa))

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthShapeForcedRoot2(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthShapeForcedRoot2"

        lsTestTaxa = ["1A|2A|3A|4A|1","1A|2A|3A|4A|2","1A|2A|3A|4A|3","1A|2A|3A|4A|4",
                      "1A|2A|3A|4A|5","1A|2A|3A|4A|6","1B|2B|3B|4B|1","1B|2B|3B|4B|2",
                      "1B|2B|3B|4B|3","1B|2B|3B|4B|4","1B|2B|3B|4B|5","1B|2B|3B|4B|6","1B|2B|3B|4B|7",
                      "1C|2C|3C|4C|1","1D|2D|3D|4D|1","1E|2E|3E|4E|1"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape=["v","^","^","v","v","^","^","v","v","^","v","^","v","^","v","^"]
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.forceRoot("1B")
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(cladogram.updateToRoot(lsTestTaxa))

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthBorderForcedRoot(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthBorderForcedRoot"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=[0.16,0.15,0.14,0.13,0.12,0.11,0.10,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
        strTestShape="v"
        dTestAlpha=0.5

        cladogram = Cladogram()
        cladogram.forceRoot("1B")
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(cladogram.updateToRoot(lsTestTaxa))

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def donttestAddingAndGeneratingCircleDataForContentTrueLengthAlphaForcedRoot(self):
        methodName = "testAddingAndGeneratingCircleDataForContentTrueLengthAlphaForcedRoot"

        lsTestTaxa = ["1A|2A|3A|4A|1A","1A|2A|3A|4A|2A","1A|2A|3A|4A|3A","1A|2A|3A|4A|4A",
                      "1A|2A|3A|4A|5A","1A|2A|3A|4A|6A","1B|2B|3B|4B|1B","1B|2B|3B|4B|2B",
                      "1B|2B|3B|4B|3B","1B|2B|3B|4B|4B","1B|2B|3B|4B|5B","1B|2B|3B|4B|6B","1B|2B|3B|4B|7B",
                      "1C|2C|3C|4C|1C","1D|2D|3D|4D|1D","1E|2E|3E|4E|1E"]

        strCircleAnswerFile = Constants_Testing.c_strTestingTMP+methodName+"Correct.txt"
        strCircleOutputFile = Constants_Testing.c_strTestingTMP+methodName+".txt"
        strTestCircle = "Test1"
        dTestBorder=0.45
        strTestShape="v"
        dTestAlpha=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.10,0.11,0.12,0.13,0.14,0.15]

        cladogram = Cladogram()
        cladogram.forceRoot("1B")
        cladogram.strCircleFilePath = strCircleOutputFile
        cladogram.addCircle(lsTaxa=lsTestTaxa, strCircle=strTestCircle, dBorder=dTestBorder, strShape=strTestShape, dAlpha=dTestAlpha, fForced=False)
        cladogram.createCircleFile(cladogram.updateToRoot(lsTestTaxa))

        fhdlCorrectFile = open(strCircleAnswerFile,"r")
        correctAnswer = fhdlCorrectFile.read()
        fhdlCorrectFile.close()
        fhdlOutputFile = open(strCircleOutputFile,"r")
        result = fhdlOutputFile.read()
        fhdlOutputFile.close()
        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testAddHighLightsForAddingNothingToNothing(self):
        methodName = "testAddHighLightsForAddingNothingToNothing"

        cladogram = Cladogram()

        cladogram.addHighLights(dictClades=None,fOverwrite=False)

        result = len(cladogram.dictForcedHighLights.keys())
        correctAnswer = 0

        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testAddHighLightsForAdding(self):
        methodName = "testAddHighLightsForAdding"

        cladogram = Cladogram()
        cladogram.addHighLights(dictClades=None,fOverwrite=False)
        cladogram.addHighLights(dictClades={"1":"One","2":"Two"},fOverwrite=False)
        cladogram.addHighLights(dictClades=None,fOverwrite=False)
        cladogram.addHighLights(dictClades={"5":"Five","6":"Six"},fOverwrite=False)

        result = len(cladogram.dictForcedHighLights.keys())
        correctAnswer = 4

        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testAddHighLightsForAdding(self):
        methodName = "testAddHighLightsForAdding"

        cladogram = Cladogram()
        cladogram.addHighLights(dictClades=None,fOverwrite=False)
        cladogram.addHighLights(dictClades={"1":"One","2":"Two"},fOverwrite=False)
        cladogram.addHighLights(dictClades=None,fOverwrite=True)
        cladogram.addHighLights(dictClades={"5":"Five","6":"Six"},fOverwrite=False)


        result = len(cladogram.dictForcedHighLights.keys())
        correctAnswer = 4

        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testAddHighLightsForAddingOverwritting(self):
        methodName = "testAddHighLightsForAddingOverwritting"

        cladogram = Cladogram()
        cladogram.addHighLights(dictClades=None,fOverwrite=True)
        cladogram.addHighLights(dictClades={"1":"One","2":"Two"},fOverwrite=False)
        cladogram.addHighLights(dictClades=None,fOverwrite=False)
        cladogram.addHighLights(dictClades={"2":"OOPs","5":"Five","6":"EEEk"},fOverwrite=False)
        cladogram.addHighLights(dictClades={"6":"Six", "3":"Three"},fOverwrite=True)

        errorOccured=False
        errorMessage=" "

        correctLength = 5
        if not len(cladogram.dictForcedHighLights.keys()) == correctLength:
          errorMessage = errorMessage + " Did not get the correct length of keys, expected:"+str(correctLength)+". but received:"+str(len(cladogram.dictForcedHighLights.keys()))+"."
          errorOccured=True
        if not cladogram.dictForcedHighLights["1"] == "One":
          errorMessage = errorMessage + " Did not get the correct data for key 1, expected:One. but received:"+str(cladogram.dictForcedHighLights["1"])+"."
          errorOccured=True
        if not cladogram.dictForcedHighLights["2"] == "Two":
          errorMessage = errorMessage + " Did not get the correct data for key 2, expected:Two. but received:"+str(cladogram.dictForcedHighLights["2"])+"."
          errorOccured=True
        if not cladogram.dictForcedHighLights["3"] == "Three":
          errorMessage = errorMessage + " Did not get the correct data for key 3, expected:Three. but received:"+str(cladogram.dictForcedHighLights["3"])+"."
          errorOccured=True
        if not cladogram.dictForcedHighLights["5"] == "Five":
          errorMessage = errorMessage + " Did not get the correct data for key 5, expected:Five. but received:"+str(cladogram.dictForcedHighLights["5"])+"."
          errorOccured=True
        if not cladogram.dictForcedHighLights["6"] == "Six":
          errorMessage = errorMessage + " Did not get the correct data for key 6, expected:Six. but received:"+str(cladogram.dictForcedHighLights["6"])+"."
          errorOccured=True

        self.assertEquals(errorOccured, False, methodName+errorMessage)

    def testSetColorDataForSetting(self):
        methodName = "testSetColorDataForSetting"

        cladogram = Cladogram()

        dictColors = {"Color1":"Color1.1", "Color2":"Color2.2"}
        cladogram.setColorData(dictColors)
        cladogram.setColorData(None)

        result = len(cladogram.dictColors.keys())
        correctAnswer = 3 #One extra for the background color automatically added if you dont.

        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSetColorDataForBackgroundColor(self):
        methodName = "testSetColorDataForBackgroundColor"

        cladogram = Cladogram()

        dictColors = {Constants_Figures.c_strBackgroundColorName:"Color2.2"}
        cladogram.setColorData(dictColors)

        result = len(cladogram.dictColors.keys())
        correctAnswer = 1

        self.assertEquals(str(correctAnswer), str(result), methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSetColorDataForSetting2(self):
        methodName = "testSetColorDataForSetting"

        cladogram = Cladogram()

        dictColors = {"Color1":"Color1.1", "Color2":"Color2.2"}
        cladogram.setColorData(dictColors)
        cladogram.setColorData(None)

        result = len(cladogram.dictColors.keys())
        errorOccured=False
        errorMessage=" "

        correctLength = 3
        if not len(cladogram.dictColors.keys()) == correctLength:
          errorMessage = errorMessage + " Did not get the correct length of keys, expected:"+str(correctLength)+". but received:"+str(len(cladogram.dictForcedHighLights.keys()))+"."
          errorOccured=True
        if not cladogram.dictColors["Color1"] == "Color1.1":
          errorMessage = errorMessage + " Did not get the correct data for key Color1, expected:One. but received:"+str(cladogram.dictColors["Color1"])+"."
          errorOccured=True
        if not cladogram.dictColors["Color2"] == "Color2.2":
          errorMessage = errorMessage + " Did not get the correct data for key Color2, expected:Two. but received:"+str(cladogram.dictColors["Color2"])+"."
          errorOccured=True
        if not cladogram.dictColors[Constants_Figures.c_strBackgroundColorName] == Constants_Figures.c_strBackgroundColor:
          errorMessage = errorMessage + " Did not get the correct data for background, expected:"+str(Constants_Figures.c_strBackgroundColor)+". but received:"+str(cladogram.dictColors[Constants_Figures.c_strBackgroundColorName])+"."
          errorOccured=True

        self.assertEquals(errorOccured, False, methodName+errorMessage)

    #Test filter by abundance
    def testFilterByAbundanceForGoodCase(self):
        methodName = "testFilterByAbundanceForGoodCase"

        strTestFile = "".join([Constants_Testing.c_strTestingData,"FilterByAbundance.txt"])
        iCurNameRow = 0
        iCurFirstDataRow = 1
        lsIds = ["Taxa1","Taxa2","Taxa3","Taxa4","Taxa5","Taxa6","Taxa7","Taxa8","Taxa9","Taxa10"]
        lsCorrectIds = ["Taxa1","Taxa2","Taxa8","Taxa9","Taxa10"]

        #Get Abundance table data
        rawData = AbundanceTable.makeFromFile(strInputFile=strTestFile, fIsNormalized=False,
                                          fIsSummed=True, iNameRow = iCurNameRow,
                                          iFirstDataRow = iCurFirstDataRow)
        cladogram = Cladogram()
        cladogram.setAbundanceData(rawData)

        lsRet = cladogram.filterByAbundance(lsIds)

        self.assertEquals(str(lsCorrectIds),str(lsRet), methodName+"Did not get the correct file, expected:"+str(lsCorrectIds)+" but received:"+str(lsRet)+".")

    # Test create size file
    def testCreateSizeFileForGoodCase(self):
        methodName = "testCreateSizeFileForGoodCase"

        strTestFile = "".join([Constants_Testing.c_strTestingData,"SmallAbundance.txt"])
        iCurNameRow = 0
        iCurFirstDataRow = 1
        lsIds = ["Taxa1","Taxa2","Taxa3","Taxa4","Taxa5"]
        answerFileContents = "Taxa1\t19.1938205472\nTaxa2\t24.2073385653\nTaxa3\t0.000130281830588\nTaxa4\t20.7493620036\nTaxa5\t23.3284684971"

        #Get Abundance table data
        rawData = AbundanceTable.makeFromFile(strInputFile=strTestFile, fIsNormalized=False,
                                              fIsSummed=True, iNameRow = iCurNameRow,
                                              iFirstDataRow = iCurFirstDataRow)
        cladogram = Cladogram()
        cladogram.setAbundanceData(rawData)

        lsRet = cladogram.filterByAbundance(lsIds)

        if(os.path.exists(cladogram.strSizeFilePath)):
          os.remove(cladogram.strSizeFilePath)

        cladogram.createSizeFile(lsIds)
        fhdlResultFile = open(cladogram.strSizeFilePath,"r")
        resultFileContents =  fhdlResultFile.read()
        fhdlResultFile.close()

        self.assertEquals(answerFileContents,resultFileContents, methodName+"Did not get the correct file, expected:"+answerFileContents+" but received:"+resultFileContents+".")

##
#Create a suite to be called to test
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(CladogramTest)
    return suite

#if __name__=='__main__':
#    unittest.main()
