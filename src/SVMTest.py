#######################################################
# Author: Timothy Tickle
# Description: Class to test the SVM class
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
import Constants
import FileIO
import numpy as np
import os
import SVM
import unittest

##
#Tests the Blog object
class SVMTest(unittest.TestCase):

    def testConvertAbundanceFileToSVMFileForGoodCase(self):

        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        outputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
        delimiter = Constants.Constants.TAB
        labels = [-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]
        firstDataRow = 2
        normalize = False
        skipFirst = True

        #Correct Answer
        answer = "-1 1:1 2:0 3:3 4:0 5:5\n1 1:0 2:10 3:0 4:45 5:0\n-1 1:0 2:0 3:0 4:0 5:0\n1 1:12 2:43 3:29 4:34 5:2\n-1 1:0 2:6 3:0 4:3 5:0\n1 1:6 2:0 3:45 4:0 5:6\n-1 1:0 2:23 3:0 4:0 5:0\n1 1:2 2:0 3:1 4:0 5:1\n-1 1:1 2:1 3:1 4:1 5:1\n1 1:0 2:0 3:0 4:0 5:0\n"

        #Call method
        SVM.SVM.convertAbundanceFileToSVMFile(tempInputFile=inputFile, tempOutputSVMFile=outputFile, tempDelimiter=delimiter, tempLabels=labels, tempFirstDataRow=2, tempSkipFirstColumn=skipFirst, tempNormalize=normalize)
        read = FileIO.FileIO(outputFile, True,False,False)
        result = read.readFullFile()
        read.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testConvertAbundanceFileToSVMFileForGoodCaseNormalize(self):

        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        outputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
        delimiter = Constants.Constants.TAB
        labels = [-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]
        firstDataRow = 2
        normalize = True
        skipFirst = True

        #Correct Answer
        answer = "-1 1:0.111111111111 2:0.0 3:0.333333333333 4:0.0 5:0.555555555556\n1 1:0.0 2:0.181818181818 3:0.0 4:0.818181818182 5:0.0\n-1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n1 1:0.1 2:0.358333333333 3:0.241666666667 4:0.283333333333 5:0.0166666666667\n-1 1:0.0 2:0.666666666667 3:0.0 4:0.333333333333 5:0.0\n1 1:0.105263157895 2:0.0 3:0.789473684211 4:0.0 5:0.105263157895\n-1 1:0.0 2:1.0 3:0.0 4:0.0 5:0.0\n1 1:0.5 2:0.0 3:0.25 4:0.0 5:0.25\n-1 1:0.2 2:0.2 3:0.2 4:0.2 5:0.2\n1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n"

        #Call method
        SVM.SVM.convertAbundanceFileToSVMFile(tempInputFile=inputFile, tempOutputSVMFile=outputFile, tempDelimiter=delimiter, tempLabels=labels, tempFirstDataRow=2, tempSkipFirstColumn=skipFirst, tempNormalize=normalize)
        read = FileIO.FileIO(outputFile, True,False,False)
        result = read.readFullFile()
        read.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testCreateLinearModelForTestDataTrain1(self):

        #Inputs
        inputFile = "./testData/SVM/train.1"
        dataPredictionFile = inputFile
        scaling = -1
        costRange = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10"
        probabilisticOutput = False

        #Correct Answer
        answerCost = "0.0"
        answerAccuracy = "95.4354"

        #Call method
        svm = SVM.SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)

        if(not outputFiles == False):
            results = svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
            resultCost = outputFiles[svm.c_COST_VALUE]
            resultAccuracy = outputFiles[svm.c_ACCURACY]
        else:
            resultCost = "False"
            resultAccuracy = "False"

        #Check result against answer
        self.assertEqual(str(resultCost)+" "+str(resultAccuracy),str(answerCost)+" "+str(answerAccuracy),"".join([str(self),"::Expected=",str(answerCost)+" "+str(answerAccuracy),". Received=",str(resultCost)+" "+str(resultAccuracy),"."]))

    def testCreateLinearModelForTestDataSVMGuide1(self):

        #Inputs
        inputFile = "./testData/SVM/svmguide1.t"
        dataPredictionFile = inputFile
        scaling = -1
        costRange = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10"
        probabilisticOutput = False

        #Correct Answer
        answerCost = "10.0"
        answerAccuracy = "95.75"

        #Call method
        svm = SVM.SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)

        if(not outputFiles == False):
            results = svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
            resultCost = outputFiles[svm.c_COST_VALUE]
            resultAccuracy = outputFiles[svm.c_ACCURACY]
        else:
            resultCost = "False"
            resultAccuracy = "False"

        #Check result against answer
        self.assertEqual(str(resultCost)+" "+str(resultAccuracy),str(answerCost)+" "+str(answerAccuracy),"".join([str(self),"::Expected=",str(answerCost)+" "+str(answerAccuracy),". Received=",str(resultCost)+" "+str(resultAccuracy),"."]))

    def nottestCreateLinearModelNotProbalistic(self):

        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
        dataPredictionFile = inputFile
        scaling = 0
        costRange = "1,2,3,4,5"
        probabilisticOutput = False

        #Correct Answer
        answer = ""

        #Call method
        svm = SVM.SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)

        if(not outputFiles == False):
            svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
        result = "."

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def nottestCreateLinearModelProbalistic(self):

        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
        dataPredictionFile = inputFile
        scaling = 0
        costRange = "1,2,3,4,5"
        probabilisticOutput = True

        #Correct Answer
        answer = ""

        #Call method
        svm = SVM.SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)
        if(not outputFiles == False):
            svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
        result = "."

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(SVMTest)
    return suite
