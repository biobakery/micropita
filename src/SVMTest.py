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
from AbundanceTable import AbundanceTable
from Constants import Constants
from Constants_Testing import Constants_Testing
import numpy as np
import os
from SVM import SVM
import unittest

##
#Tests the Blog object
class SVMTest(unittest.TestCase):

    def ntestConvertAbundanceFileToSVMFileForGoodCase(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"])
        delimiter = Constants.TAB
        labels = [-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]
        firstDataRow = 2
        normalize = False
        skipFirst = True

        #Correct Answer
        answer = "-1 1:1 2:0 3:3 4:0 5:5\n1 1:0 2:10 3:0 4:45 5:0\n-1 1:0 2:0 3:0 4:0 5:0\n1 1:12 2:43 3:29 4:34 5:2\n-1 1:0 2:6 3:0 4:3 5:0\n1 1:6 2:0 3:45 4:0 5:6\n-1 1:0 2:23 3:0 4:0 5:0\n1 1:2 2:0 3:1 4:0 5:1\n-1 1:1 2:1 3:1 4:1 5:1\n1 1:0 2:0 3:0 4:0 5:0\n"

        #Call method
        SVM.convertAbundanceFileToSVMFile(tempInputFile=inputFile, tempOutputSVMFile=outputFile, tempDelimiter=delimiter, tempLabels=labels, sLastMetadataName="STSite", tempSkipFirstColumn=skipFirst, tempNormalize=normalize)
        with open(outputFile,'r') as f:
            result = f.read()
        f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def ntestConvertAbundanceFileToSVMFileForGoodCaseNormalize(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"])
        delimiter = Constants.TAB
        labels = [-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]
        firstDataRow = 2
        normalize = True
        skipFirst = True

        #Correct Answer
        answer = "-1 1:0.111111111111 2:0.0 3:0.333333333333 4:0.0 5:0.555555555556\n1 1:0.0 2:0.181818181818 3:0.0 4:0.818181818182 5:0.0\n-1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n1 1:0.1 2:0.358333333333 3:0.241666666667 4:0.283333333333 5:0.0166666666667\n-1 1:0.0 2:0.666666666667 3:0.0 4:0.333333333333 5:0.0\n1 1:0.105263157895 2:0.0 3:0.789473684211 4:0.0 5:0.105263157895\n-1 1:0.0 2:1.0 3:0.0 4:0.0 5:0.0\n1 1:0.5 2:0.0 3:0.25 4:0.0 5:0.25\n-1 1:0.2 2:0.2 3:0.2 4:0.2 5:0.2\n1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n"

        #Call method
        SVM.convertAbundanceFileToSVMFile(tempInputFile=inputFile, tempOutputSVMFile=outputFile, tempDelimiter=delimiter, tempLabels=labels, sLastMetadataName="STSite", tempSkipFirstColumn=skipFirst, tempNormalize=normalize)
        with open(outputFile,'r') as f:
            result = f.read()
        f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testConvertAbundanceTableToSVMFileForGoodCase(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
#"L_Antecubital_fossa"	"R_Retroauricular_crease"	"L_Retroauricular_crease"	"Subgingival_plaque"	"R_Antecubital_fossa"	"L_Retroauricular_crease"	"R_Retroauricular_crease"	"L_Antecubital_fossa"	"R_Antecubital_fossa"	"Anterior_nares"
        #Correct Answer
        answer = "-1 1:0.111111111111 2:0.0 3:0.333333333333 4:0.0 5:0.555555555556\n1 1:0.0 2:0.181818181818 3:0.0 4:0.818181818182 5:0.0\n-1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n1 1:0.1 2:0.358333333333 3:0.241666666667 4:0.283333333333 5:0.0166666666667\n-1 1:0.0 2:0.666666666667 3:0.0 4:0.333333333333 5:0.0\n1 1:0.105263157895 2:0.0 3:0.789473684211 4:0.0 5:0.105263157895\n-1 1:0.0 2:1.0 3:0.0 4:0.0 5:0.0\n1 1:0.5 2:0.0 3:0.25 4:0.0 5:0.25\n-1 1:0.2 2:0.2 3:0.2 4:0.2 5:0.2\n1 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0\n"

        strOutputFile = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVMFile.txt"])

        abndData = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        SVM.convertAbundanceTableToSVMFile(abndAbundanceTable=abndData, tempOutputSVMFile=strOutputFile, sMetadataLabel=sLastMetadata)

        with open(strOutputFile,'r') as f:
            result = f.read()
        f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def ntestCreateLinearModelForTestDataTrain1(self):

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

    def ntestCreateLinearModelForTestDataSVMGuide1(self):

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
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
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
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVM.txt"
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
    suite = unittest.TestLoader().loadTestsFromTestCase(SVMTest)
    return suite
