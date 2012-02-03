#######################################################
#
#	Title:		Utility_FileTest
#	Author:		Timothy Tickle 
#	Date:		June 28, 2010
#	Purpose:	Test Utility_File class
#
#######################################################

#Import local code
import Constants
import FileIO
import os
import unittest
import Utility_File

class Utility_FileTest(unittest.TestCase):

    ##Set up for tests
    def setUp(self): pass

    def testSplitFileForBadCaseWrongNumberOfLines(self):
        methodName = "testSplitFileForBadCaseWrongNumberOfLines"

        tempNumberOfLines = 99
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = True
        error = ""
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)
        if(not result):
            error = error + "Received a false in method call."

        readFileLineCount = Utility_File.Utility_File.fileLineCount(tempReadFileName)
        writeFileLineCount = Utility_File.Utility_File.fileLineCount("./Testing/TestSplitFileOutput_1.txt")
        if(not (readFileLineCount == writeFileLineCount)):
            result = False
            error = error + "Split created a file of line count :"+str(writeFileLineCount)+". From a file of line count :"+str(readFileLineCount)+"."
        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+"."+error)

    def testIncrementFileForNull(self):
        methodName = "testIncrementFileForNull"

        testString = None
        correctAnswer = False
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testIncrementFileForWrongType(self):
        methodName = "testIncrementFileForWrongType"

        testString = 456754
        correctAnswer = False
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testIncrementFileForBlank(self):
        methodName = "testIncrementFileForBlank"

        testString = "                                "
        correctAnswer = False
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testIncrementFileForEmpty(self):
        methodName = "testIncrementFileForEmpty"

        testString = ""
        correctAnswer = False
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")


    def testIncrementFileForNoPeriod(self):
        methodName = "testIncrementFileForNoPeriod"

        testString = "simple"
        correctAnswer = "simple_1"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForNoPeriodOneLineLetter(self):
        methodName = "testIncrementFileForNoPeriodOneLineLetter"

        testString = "simple_A"
        correctAnswer = "simple_A_1"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForNoPeriodOneLineNumber(self):
        methodName = "testIncrementFileForNoPeriodOneLineLetter"

        testString = "simple_9"
        correctAnswer = "simple_10"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForNoPeriodTwoLineLetter(self):
        methodName = "testIncrementFileForNoPeriodOneLineLetter"

        testString = "simple_B_A"
        correctAnswer = "simple_B_A_1"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForNoPeriodTwoLineNumber(self):
        methodName = "testIncrementFileForNoPeriodOneLineLetter"

        testString = "simple_D_9"
        correctAnswer = "simple_D_10"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriod(self):
        methodName = "testIncrementFileForOnePeriod"

        testString = "simple.txt"
        correctAnswer = "simple_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodManyLines(self):
        methodName = "testIncrementFileForOnePeriod"

        testString = "simple_____.txt"
        correctAnswer = "simple______1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodManyDots(self):
        methodName = "testIncrementFileForOnePeriod"

        testString = "simple......txt"
        correctAnswer = "simple....._1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodOneLineLetter(self):
        methodName = "testIncrementFileForOnePeriodOneLineLetter"

        testString = "simple_S.txt"
        correctAnswer = "simple_S_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodOneLineNumber(self):
        methodName = "testIncrementFileForOnePeriodOneLineNumber"

        testString = "simple_5.txt"
        correctAnswer = "simple_6.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodTwoLineLetter(self):
        methodName = "testIncrementFileForOnePeriodTwoLineLetter"

        testString = "simple_2_S.txt"
        correctAnswer = "simple_2_S_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForOnePeriodTwoLineNumber(self):
        methodName = "testIncrementFileForOnePeriodTwoLineNumber"

        testString = "simple_g_5.txt"
        correctAnswer = "simple_g_6.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseTwoPeriods(self):
        methodName = "testIncrementFileForTwoPeriod"

        testString = "simple.mostly.txt"
        correctAnswer = "simple.mostly_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseTwoPeriodsOneLineLetter(self):
        methodName = "testIncrementFileForTwoPeriodsOneLineLetter"

        testString = "simple.mostly_k.txt"
        correctAnswer = "simple.mostly_k_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseTwoPeriodsOneLineNumber(self):
        methodName = "testIncrementFileForTwoPeriodsOneLineNumber"

        testString = "simple.mostly_7.txt"
        correctAnswer = "simple.mostly_8.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseTwoPeriodsTwoLineLetter(self):
        methodName = "testIncrementFileForTwoPeriodsTwoLineLetter"

        testString = "simple.mostly_0_k.txt"
        correctAnswer = "simple.mostly_0_k_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseTwoPeriodsTwoLineNumber(self):
        methodName = "testIncrementFileForTwoPeriodsTwoLineNumber"

        testString = "simple.mostly_0_7.txt"
        correctAnswer = "simple.mostly_0_8.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseThreePeriods(self):
        methodName = "testIncrementFileForThreePeriod"

        testString = "simple.mostly.notreally.txt"
        correctAnswer = "simple.mostly.notreally_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseThreePeriodsOneLineLetter(self):
        methodName = "testIncrementFileForThreePeriodOneLineLetter"

        testString = "simple.mostly.notreally_w.txt"
        correctAnswer = "simple.mostly.notreally_w_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseThreePeriodsOneLineNumber(self):
        methodName = "testIncrementFileForThreePeriodOneLineNumber"

        testString = "simple.mostly.notreally_2.txt"
        correctAnswer = "simple.mostly.notreally_3.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseThreePeriodsTwoLineLetter(self):
        methodName = "testIncrementFileForThreePeriodTwoLineLetter"

        testString = "simple.mostly.notreally_0_w.txt"
        correctAnswer = "simple.mostly.notreally_0_w_1.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")

    def testIncrementFileForGoodCaseThreePeriodsTwoLineNumber(self):
        methodName = "testIncrementFileForThreePeriodTwoLineNumber"

        testString = "simple.mostly.notreally_tyfg_2.txt"
        correctAnswer = "simple.mostly.notreally_tyfg_3.txt"
        result = Utility_File.Utility_File.incrementFileName(testString)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+correctAnswer+". but received ."+result+".")


    def testSplitFileForBadCaseNegativeNumberOfLines(self):
        methodName = "testSplitFileForBadCaseNegativeNumberOfLines"

        tempNumberOfLines = -1
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileForBadCaseNegativeNumberOfLines.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseWrongTypeNumberOfLines(self):
        methodName = "testSplitFileForBadCaseWrongTypeNumberOfLines"

        tempNumberOfLines = "l"
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseNullNumberOfLines(self):
        methodName = "testSplitFileForBadCaseNullNumberOfLines"

        tempNumberOfLines = None
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseBlankReadFile(self):
        methodName = "testSplitFileForBadCaseBlankReadFile"

        tempNumberOfLines = 10
        tempReadFileName = "     "
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseWrongTypeReadFile(self):
        methodName = "testSplitFileForBadCaseWrongTypeReadFile"

        tempNumberOfLines = 10
        tempReadFileName = True
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseEmptyReadFile(self):
        methodName = "testSplitFileForBadCaseEmptyReadFile"

        tempNumberOfLines = 10
        tempReadFileName = ""
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseNoneReadFile(self):
        methodName = "testSplitFileForBadCaseNoneReadFile"

        tempNumberOfLines = 10
        tempReadFileName = None
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseWrongReadFile(self):
        methodName = "testSplitFileForBadCaseWrongReadFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFileOutputNottaFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseWrongOuputFile(self):
        methodName = "testSplitFileForBadCaseWrongReadFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/Notta a directory//TestSplitFileOutput.txt"

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseNoneOutputFile(self):
        methodName = "testSplitFileForBadCaseNoneOutputFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = None

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseBlankOutputFile(self):
        methodName = "testSplitFileForBadCaseBlankOutputFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "                                  "

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseEmptyOutputFile(self):
        methodName = "testSplitFileForBadCaseEmptyOutputFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = ""

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForBadCaseWrongTypeOutputFile(self):
        methodName = "testSplitFileForBadCaseEmptyOutputFile"

        tempNumberOfLines = 10
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = 876543

        correctAnswer = False
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testSplitFileForGoodCase1OutputFile(self):
        methodName = "testSplitFileForGoodCase1OutputFile"

        hasError = False
        errorMessage = ""

        tempNumberOfLines = 100
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"
        fileData = "1\n2\n3\n4\n5\n6\n7\n8\n9\n0\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n31\n32\n33"

        if os.path.isfile("./Testing/TestSplitFile.txt"):
            os.remove("./Testing/TestSplitFile.txt")
        fileWriter = FileIO.FileIO("./Testing/TestSplitFile.txt",False,True,True)
        fileWriter.writeToFile(fileData)
        fileWriter.close()

        correctAnswer = True
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)
        if(not correctAnswer == result):
            hasError = True
            errorMessage = errorMessage + " Did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+"."

        #Check for file 1
        if not os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_1.txt"

        #Check contents of file 1
        if os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_1.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == fileData:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+fileData+". but received: "+contents+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_1.txt")
        self.assertEquals(False, hasError, methodName+errorMessage)

    def testSplitFileForGoodCase2OutputFile(self):
        methodName = "testSplitFileForGoodCase2OutputFile"

        hasError = False
        errorMessage = ""

        tempNumberOfLines = 20
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"
        fileData = "1\n2\n3\n4\n5\n6\n7\n8\n9\n0\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n31\n32\n33"
        results1 = "1\n2\n3\n4\n5\n6\n7\n8\n9\n0\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n"
        results2 = "21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n31\n32\n33"

        if os.path.isfile("./Testing/TestSplitFile.txt"):
            os.remove("./Testing/TestSplitFile.txt")
        fileWriter = FileIO.FileIO("./Testing/TestSplitFile.txt",False,True,True)
        fileWriter.writeToFile(fileData)
        fileWriter.close()

        correctAnswer = True
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)
        if(not correctAnswer == result):
            hasError = True
            errorMessage = errorMessage + " Did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+"."

        #Check for file 1
        if not os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_1.txt"

        #Check contents of file 1
        if os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_1.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == results1:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+results1+". but received: "+str(contents)+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_1.txt")

        #Check for file 2
        if not os.path.isfile("./Testing/TestSplitFileOutput_2.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_2.txt"

        #Check contents of file 2
        if os.path.isfile("./Testing/TestSplitFileOutput_2.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_2.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == results2:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+results2+". but received: "+str(contents)+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_2.txt")

        self.assertEquals(False, hasError, methodName+errorMessage)

    def testSplitFileForGoodCase3OutputFile(self):
        methodName = "testSplitFileForGoodCase3OutputFile"

        hasError = False
        errorMessage = ""

        tempNumberOfLines = 15
        tempReadFileName = "./Testing/TestSplitFile.txt"
        tempOutputFileName = "./Testing/TestSplitFileOutput.txt"
        fileData = "1\n2\n3\n4\n5\n6\n7\n8\n9\n0\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n31\n32\n33"
        results1 = "1\n2\n3\n4\n5\n6\n7\n8\n9\n0\n11\n12\n13\n14\n15\n"
        results2 = "16\n17\n18\n19\n20\n21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n"
        results3 = "31\n32\n33"

        if os.path.isfile("./Testing/TestSplitFile.txt"):
            os.remove("./Testing/TestSplitFile.txt")
        fileWriter = FileIO.FileIO("./Testing/TestSplitFile.txt",False,True,True)
        fileWriter.writeToFile(fileData)
        fileWriter.close()

        correctAnswer = True
        result = Utility_File.Utility_File.splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName)
        if(not correctAnswer == result):
            hasError = True
            errorMessage = errorMessage + " Did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+"."

        #Check for file 1
        if not os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_1.txt"

        #Check contents of file 1
        if os.path.isfile("./Testing/TestSplitFileOutput_1.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_1.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == results1:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+results1+". but received: "+str(contents)+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_1.txt")

        #Check for file 2
        if not os.path.isfile("./Testing/TestSplitFileOutput_2.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_2.txt"

        #Check contents of file 2
        if os.path.isfile("./Testing/TestSplitFileOutput_2.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_2.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == results2:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+results2+". but received: "+str(contents)+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_2.txt")

        #Check for file 3
        if not os.path.isfile("./Testing/TestSplitFileOutput_3.txt"):
            hasError = True
            errorMessage = errorMessage + " Did not generate file:"+"./Testing/TestSplitFileOutput_3.txt"

        #Check contents of file 3
        if os.path.isfile("./Testing/TestSplitFileOutput_3.txt"):
            fileReader = FileIO.FileIO("./Testing/TestSplitFileOutput_3.txt",True, False, False)
            contents = fileReader.readFullFile()
            if not contents == results3:
                hasError = True
                errorMessage = errorMessage + " Contents were not correct. Expected :"+results3+". but received: "+str(contents)+"."
            fileReader.close()
            os.remove("./Testing/TestSplitFileOutput_3.txt")

        self.assertEquals(False, hasError, methodName+errorMessage)

    def testCombineFilesForNullOutputName(self):
        methodName = "testCombineFilesForNullOutputName"
        #String
        outputFileName = None
        #List
        listOfFiles = ["TestCombineFilesInoutFile1.txt","TestCombineFilesInoutFile2.txt","TestCombineFilesInoutFile3.txt"]
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForWrongTypeOutputName(self):
        methodName = "testCombineFilesForWrongTypeOutputName"
        #String
        outputFileName = 1
        #List
        listOfFiles = ["TestCombineFilesInoutFile1.txt","TestCombineFilesInoutFile2.txt","TestCombineFilesInoutFile3.txt"]
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForNullListOfFiles(self):
        methodName = "testCombineFilesForNullListOfFiles"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = None
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForWrongTypeListOfFiles(self):
        methodName = "testCombineFilesForWrongTypeListOfFiles"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = 1
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForNullInListOfFiles(self):
        methodName = "testCombineFilesForNullInListOfFiles"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["Test1.txt",None,"Test2.txt"]
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForNoneRemoveBoolean(self):
        methodName = "testCombineFilesForNoneRemoveBoolean"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInoutFile1.txt","TestCombineFilesInoutFile2.txt","TestCombineFilesInoutFile3.txt"]
        #Bool
        tempRemoveHeader = None
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForWrongTypeRemoveBoolean(self):
        methodName = "testCombineFilesForWrongTypeRemoveBoolean"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInoutFile1.txt","TestCombineFilesInoutFile2.txt","TestCombineFilesInoutFile3.txt"]
        #Bool
        tempRemoveHeader = 2
        #Bool
        correctReturn = False
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testCombineFilesForGood3Files(self):
        methodName = "testCombineFilesForGood3Files"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInputFile1.txt","TestCombineFilesInputFile2.txt","TestCombineFilesInputFile3.txt"]
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = True

        #Build file 1
        for fileName in listOfFiles:
            fileWriter = open(fileName, 'w')
            fileWriter.write("I am a header\n")
            fileWriter.write("Line 1\n")
            fileWriter.write("Line 2\n")
            fileWriter.write("Line 3\n")
            fileWriter.write("Line 4\n")
            fileWriter.write("Line 5\n")
            fileWriter.write("Line 6\n")
            fileWriter.write("Line 7\n")
            fileWriter.write("Line 8\n")
            fileWriter.write("Line 9\n")
            fileWriter.write("Line 10\n")
            fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = "I am a header\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\nI am a header\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\nI am a header\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\n"
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(outputFileName)
        for fileName in listOfFiles:
            os.remove(fileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(cumulativeResult)+". but received ."+str(result)+". Error found at"+errorSource)

    def testCombineFilesForGood1File(self):
        methodName = "testCombineFilesForGood1File"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInputFile1.txt"]
        #Bool
        tempRemoveHeader = False
        #Bool
        correctReturn = True

        #Build file 1
        for fileName in listOfFiles:
            fileWriter = open(fileName, 'w')
            fileWriter.write("I am a header\n")
            fileWriter.write("Line 1\n")
            fileWriter.write("Line 2\n")
            fileWriter.write("Line 3\n")
            fileWriter.write("Line 4\n")
            fileWriter.write("Line 5\n")
            fileWriter.write("Line 6\n")
            fileWriter.write("Line 7\n")
            fileWriter.write("Line 8\n")
            fileWriter.write("Line 9\n")
            fileWriter.write("Line 10\n")
            fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = "I am a header\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\n"
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(outputFileName)
        for fileName in listOfFiles:
            os.remove(fileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(cumulativeResult)+". but received ."+str(result)+". Error found at"+errorSource)

    def testCombineFilesForGood3FilesNoHeader(self):
        methodName = "testCombineFilesForGood3FilesNoHeader"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInputFile1.txt","TestCombineFilesInputFile2.txt","TestCombineFilesInputFile3.txt"]
        #Bool
        tempRemoveHeader = True
        #Bool
        correctReturn = True

        #Build file 1
        for fileName in listOfFiles:
            fileWriter = open(fileName, 'w')
            fileWriter.write("I am a header\n")
            fileWriter.write("Line 1\n")
            fileWriter.write("Line 2\n")
            fileWriter.write("Line 3\n")
            fileWriter.write("Line 4\n")
            fileWriter.write("Line 5\n")
            fileWriter.write("Line 6\n")
            fileWriter.write("Line 7\n")
            fileWriter.write("Line 8\n")
            fileWriter.write("Line 9\n")
            fileWriter.write("Line 10\n")
            fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = "Line 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\n"
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(outputFileName)
        for fileName in listOfFiles:
            os.remove(fileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(cumulativeResult)+". but received ."+str(result)+". Error found at"+errorSource)

    def testCombineFilesForGood1FileNoHeader(self):
        methodName = "testCombineFilesForGood1FileNoHeader"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #List
        listOfFiles = ["TestCombineFilesInputFile1.txt"]
        #Bool
        tempRemoveHeader = True
        #Bool
        correctReturn = True

        #Build file 1
        for fileName in listOfFiles:
            fileWriter = open(fileName, 'w')
            fileWriter.write("I am a header\n")
            fileWriter.write("Line 1\n")
            fileWriter.write("Line 2\n")
            fileWriter.write("Line 3\n")
            fileWriter.write("Line 4\n")
            fileWriter.write("Line 5\n")
            fileWriter.write("Line 6\n")
            fileWriter.write("Line 7\n")
            fileWriter.write("Line 8\n")
            fileWriter.write("Line 9\n")
            fileWriter.write("Line 10\n")
            fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.combineFiles(outputFileName,listOfFiles,tempRemoveHeader)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = "Line 1\nLine 2\nLine 3\nLine 4\nLine 5\nLine 6\nLine 7\nLine 8\nLine 9\nLine 10\n"
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(outputFileName)
        for fileName in listOfFiles:
            os.remove(fileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(cumulativeResult)+". but received ."+str(result)+". Error found at"+errorSource)

    def testFileNamePrefixForNone(self):
        methodName = "testFileNamePrefixForNone"

        fileName = None
        correctAnswer = False

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForBlank(self):
        methodName = "testFileNamePrefixForBlank"

        fileName = "   "
        correctAnswer = False

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForEmpty(self):
        methodName = "testFileNamePrefixForEmpty"

        fileName = ""
        correctAnswer = False

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForWrongType(self):
        methodName = "testFileNamePrefixForWrongType"

        fileName = list()
        correctAnswer = False

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForGoodCase(self):
        methodName = "testFileNamePrefixForGoodCase"

        fileName = "GoodCase.txt"
        correctAnswer = "GoodCase"

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForNoPeriod(self):
        methodName = "testFileNamePrefixForNoPeriod"

        fileName = "GoodCasetxt"
        correctAnswer = "GoodCasetxt"

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForSpecialCharacters(self):
        methodName = "testFileNamePrefixForSpecialCharacters"

        fileName = "!@#$%^&&*()_+.txt"
        correctAnswer = "!@#$%^&&*()_+"

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testFileNamePrefixForMultPeriod(self):
        methodName = "testFileNamePrefixForSpecialCharacters"

        fileName = "GoodCase..txt"
        correctAnswer = "GoodCase."

        #Check result return
        result = Utility_File.Utility_File.getFileNamePrefix(fileName)

        self.assertEquals(correctAnswer, result, methodName+" did not give correct result. Expected ."+str(correctAnswer)+". but received ."+str(result)+".")

    def testAbridgeFileForGoodCase(self):
        methodName = "testAbridgeFileForGoodCase"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 6
        #Indicator if the test is supposed to pass or fail
        correctReturn = True

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = "I am a header\nLine 1\nLine 2\nLine 3\nLine 4\nLine 5\n"
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(inputFileName)
        os.remove(outputFileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+". Error found at"+errorSource)

    def testAbridgeFileFor0CountCase(self):
        methodName = "testAbridgeFileFor0CountCase"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 0
        #Indicator if the test is supposed to pass or fail
        correctReturn = True

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        fileReader = open(outputFileName, 'r')
        outputFileString = fileReader.read()
        fileReader.close()
        correctContent = ""
        contentResult = (outputFileString == correctContent)

        #Delete files
        os.remove(inputFileName)
        os.remove(outputFileName)
        
        #Set up cumulative result
        cumulativeResult = (result and contentResult)
        errorSource = ""
        if(not result):
            errorSource = errorSource +" Method Return"
        if(not contentResult):
            errorSource = errorSource +" Content."+outputFileString+".outputvscorrect."+correctContent
    
        self.assertEquals(correctReturn, cumulativeResult, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+". Error found at"+errorSource)


    def testAbridgeFileForNoneCount(self):
        methodName = "testAbridgeFileForNoneCount"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = None
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForNegCount(self):
        methodName = "testAbridgeFileForNegCount"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = -7
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForStringCount(self):
        methodName = "testAbridgeFileForStringCount"
        #String
        outputFileName = "TestCombineFilesOutput.txt"
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = "10"
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForNoneOutputFile(self):
        methodName = "testAbridgeFileForNoneOutputFile"
        #String
        outputFileName = None
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForNoneInputFile(self):
        methodName = "testAbridgeFileForNoneInputFile"
        #String
        outputFileName =  "TestCombineFilesOutput.txt"
        #String
        inputFileName = None
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(inputFileName)
            os.remove(outputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForBlankOutputFile(self):
        methodName = "testAbridgeFileForBlankOutputFile"
        #String
        outputFileName = "    "
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForBlankInputFile(self):
        methodName = "testAbridgeFileForBlankInputFile"
        #String
        outputFileName =  "TestCombineFilesOutput.txt"
        #String
        inputFileName = "  "
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(inputFileName)
            os.remove(outputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForEmptyOutputFile(self):
        methodName = "testAbridgeFileForEmptyOutputFile"
        #String
        outputFileName = ""
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForEmptyInputFile(self):
        methodName = "testAbridgeFileForEmptyInputFile"
        #String
        outputFileName =  "TestCombineFilesOutput.txt"
        #String
        inputFileName = ""
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(inputFileName)
            os.remove(outputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForWrongTypeOutputFile(self):
        methodName = "testAbridgeFileForWrongTypeOutputFile"
        #String
        outputFileName = False
        #String
        inputFileName = Constants.Constants.TEST_ANSWER_DOCUMENTS+methodName+".txt"
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Build file 1
        fileWriter = open(inputFileName, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(outputFileName)

        #Delete files
        os.remove(inputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAbridgeFileForWrongTypeInputFile(self):
        methodName = "testAbridgeFileForWrongTypeInputFile"
        #String
        outputFileName =  "TestCombineFilesOutput.txt"
        #String
        inputFileName = dict()
	#Int
        count = 5
        #Indicator if the test is supposed to pass or fail
        correctReturn = False

        #Check result return
        result = Utility_File.Utility_File.abridgeFile(count, inputFileName, outputFileName)

        #Check result
        if(result):
            os.remove(inputFileName)
            os.remove(outputFileName)

        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForGoodCase(self):
        methodName = "testAreEqualForGoodCase"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = True

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseBlank1(self):
        methodName = "testAreEqualForBlank1"
        #String
        fileName1 = " "
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseBlank2(self):
        methodName = "testAreEqualForBlank2"
        #String
        fileName2 = " "
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseEmpty1(self):
        methodName = "testAreEqualForEmpty1"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual("", fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseEmpty2(self):
        methodName = "testAreEqualForEmpty2"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, "")

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseNone1(self):
        methodName = "testAreEqualForNone1"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(None, fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseNone2(self):
        methodName = "testAreEqualForNone2"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, None)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseType1(self):
        methodName = "testAreEqualForType1"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(1, fileName2)

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseType2(self):
        methodName = "testAreEqualForType2"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = False

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, list())

        #Delete files
        os.remove(fileName1)
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseNotExist1(self):
        methodName = "testAreEqualForNotExist1"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = None

        #Build file 1
        fileWriter = open(fileName2, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, fileName2)

        #Delete files
        os.remove(fileName2)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testAreEqualForBadCaseNotExist2(self):
        methodName = "testAreEqualForNotExist2"
        #String
        fileName2 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase1.txt"
        #String
        fileName1 = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"testAreEqualForGoodCase2.txt"

        correctReturn = None

        #Build file 1
        fileWriter = open(fileName1, 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("Line 1\n")
        fileWriter.write("Line 2\n")
        fileWriter.write("Line 3\n")
        fileWriter.write("Line 4\n")
        fileWriter.write("Line 5\n")
        fileWriter.write("Line 6\n")
        fileWriter.write("Line 7\n")
        fileWriter.write("Line 8\n")
        fileWriter.write("Line 9\n")
        fileWriter.write("Line 10\n")
        fileWriter.close()

        #Check result return
        result = Utility_File.Utility_File.areEqual(fileName1, fileName2)

        #Delete files
        os.remove(fileName1)
        
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForBadCaseEmptyFileName(self):
        methodName = "testClearDirectoryForBadCaseEmptyFileName"
        tempFileName = ""
        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForBadCaseBlankFileName(self):
        methodName = "testClearDirectoryForBadCaseBlankFileName"
        tempFileName = "      "
        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForBadCaseNoneFileName(self):
        methodName = "testClearDirectoryForBadCaseNoneFileName"
        tempFileName = None
        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForBadCaseWrongTypeFileName(self):
        methodName = "testClearDirectoryForBadCaseWrongTypeFileName"
        tempFileName = 1
        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForBadCaseWrongFileName(self):
        methodName = "testClearDirectoryForBadCaseWrongFileName"
        tempFileName = "..///ImNotaFile!.what"
        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testClearDirectoryForGoodCase(self):
        methodName = "testClearDirectoryForGoodCase"
        correctReturn = True
        tempFileName = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"

        if(not os.path.exists(tempFileName)):
            os.mkdir(tempFileName)
        fileWriter = open(Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt", 'w')
        fileWriter.write("I am a header\n")
        fileWriter.close()

        result = Utility_File.Utility_File.clearDirectory(tempDirectory = tempFileName)
        errorNotFound = os.path.exists(tempFileName)
        error = ""
        if(not errorNotFound):
            error = error + " Directory not remade."
        else:
            fileCount = len(os.listdir(tempFileName))
            if(not fileCount == 0):
                errorNotFound = False
                error = error + "Did not receive an empty file count, recieved file count of "+str(fileCount)+"."

        self.assertEquals(correctReturn, errorNotFound, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(errorNotFound)+". "+error)



    def testFileLineCountForBadCaseEmptyFileName(self):
        methodName = "testFileLineCountForBadCaseEmptyFileName"
        tempFileName = ""
        result = Utility_File.Utility_File.fileLineCount(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testFileLineCountForBadCaseBlankFileName(self):
        methodName = "testFileLineCountForBadCaseBlankFileName"
        tempFileName = "      "
        result = Utility_File.Utility_File.fileLineCount(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testFileLineCountForBadCaseNoneFileName(self):
        methodName = "testFileLineCountForBadCaseNoneFileName"
        tempFileName = None
        result = Utility_File.Utility_File.fileLineCount(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testFileLineCountForBadCaseWrongTypeFileName(self):
        methodName = "testFileLineCountForBadCaseWrongTypeFileName"
        tempFileName = 1
        result = Utility_File.Utility_File.fileLineCount(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testFileLineCountForBadCaseWrongFileName(self):
        methodName = "testFileLineCountForBadCaseWrongFileName"
        tempFileName = "..///ImNotaFile!.what"
        result = Utility_File.Utility_File.fileLineCount(tempDirectory = tempFileName)
        correctReturn = False
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received ."+str(result)+".")

    def testFileLineCountForGoodCase10Lines(self):
        methodName = "testFileLineCountForGoodCase10Lines"
        correctReturn = 10
        tempFileName = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"
        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)

        if(not os.path.exists(tempFileName)):
            os.mkdir(tempFileName)

        fileWriter = open(Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt", 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.close()

        result = Utility_File.Utility_File.fileLineCount(tempDirectory = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt")

        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received "+str(result)+".")

    def testFileLineCountForGoodCase30Lines(self):
        methodName = "testFileLineCountForGoodCase30Lines"
        correctReturn = 30
        tempFileName = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"
        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)

        if(not os.path.exists(tempFileName)):
            os.mkdir(tempFileName)

        fileWriter = open(Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt", 'w')
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.write("I am a header\n")
        fileWriter.close()

        result = Utility_File.Utility_File.fileLineCount(tempDirectory = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt")

        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received "+str(result)+".")

    def testFileLineCountForGoodCase0Lines(self):
        methodName = "testFileLineCountForGoodCase10Lines"
        correctReturn = 0
        tempFileName = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"
        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)

        if(not os.path.exists(tempFileName)):
            os.mkdir(tempFileName)

        fileWriter = open(Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt", 'w')
        fileWriter.close()

        result = Utility_File.Utility_File.fileLineCount(tempDirectory = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt")

        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received "+str(result)+".")

    def testFileLineCountForGoodCase1Lines(self):
        methodName = "testFileLineCountForGoodCase10Lines"
        correctReturn = 1
        tempFileName = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"
        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)

        if(not os.path.exists(tempFileName)):
            os.mkdir(tempFileName)

        fileWriter = open(Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt", 'w')
        fileWriter.write("I am a header\n")
        fileWriter.close()

        result = Utility_File.Utility_File.fileLineCount(tempDirectory = Constants.Constants.TEST_DATA_TEMP_DIRECTORY+"ClearMe/"+"DeleteMe.txt")

        Utility_File.Utility_File.clearDirectory(Constants.Constants.TEST_DATA_TEMP_DIRECTORY)
        self.assertEquals(correctReturn, result, methodName+" did not give correct result. Expected ."+str(correctReturn)+". but received "+str(result)+".")

##
#Create a suite to be called to test
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(Utility_FileTest)
    return suite

#if __name__=='__main__':
#    unittest.main()
