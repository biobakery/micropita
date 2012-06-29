#######################################################
# Author: Timothy Tickle
# Description: Class to test the LIBSVM class
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
from AbundanceTable import AbundanceTable
from Constants import Constants
from Constants_Testing import Constants_Testing
import numpy as np
import os
from LIBSVM import LIBSVM
import unittest

##
#Tests the Blog object
class LIBSVMTest(unittest.TestCase):

    def holdtestCreateLinearModelForTestDataTrain1(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"SVM/train.1"])
        dataPredictionFile = inputFile
        scaling = -1
        costRange = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10"
        probabilisticOutput = False

        #Correct Answer
        answerCost = "0.0"
        answerAccuracy = "95.4354"

        #Call method
        svm = SVM()
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

    def holdtestCreateLinearModelForTestDataSVMGuide1(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"SVM/svmguide1.t"])
        dataPredictionFile = inputFile
        scaling = -1
        costRange = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10"
        probabilisticOutput = False

        #Correct Answer
        answerCost = "10.0"
        answerAccuracy = "95.75"

        #Call method
        svm = SVM()
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
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"])
        dataPredictionFile = inputFile
        scaling = 0
        costRange = "1,2,3,4,5"
        probabilisticOutput = False

        #Correct Answer
        answer = ""

        #Call method
        svm = SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)

        if(not outputFiles == False):
            svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
        result = "."

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def nottestCreateLinearModelProbalistic(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"])
        dataPredictionFile = inputFile
        scaling = 0
        costRange = "1,2,3,4,5"
        probabilisticOutput = True

        #Correct Answer
        answer = ""

        #Call method
        svm = SVM()
        outputFiles = svm.createLinearModel(tempInputFileName=inputFile, tempScaling=scaling, tempLogC=costRange, tempProbabilistic=probabilisticOutput)
        if(not outputFiles == False):
            svm.predictFromLinearModel(tempDataFileName=dataPredictionFile, tempModelFileName=outputFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=outputFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=probabilisticOutput)
        result = "."

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(LIBSVMTest)
    return suite
