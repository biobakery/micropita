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

    def testConvertAbundanceTableToSVMFileForGoodCase(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False

        strOutputFile = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-SVMFile.txt"])

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        lsUniqueLabels = SVM.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abndData, strOutputSVMFile=strOutputFile, sMetadataLabel=sLastMetadata)

        #Correct Answer
        sLAntecubitalFossa = lsUniqueLabels.index("L_Antecubital_fossa")
        sRRetroauricularCrease = lsUniqueLabels.index("R_Retroauricular_crease")
        sLRetroauricularCrease = lsUniqueLabels.index("L_Retroauricular_crease")
        sSubgingivalPlaque = lsUniqueLabels.index("Subgingival_plaque")
        sRAntecubitalFossa = lsUniqueLabels.index("R_Antecubital_fossa")
        sAnteriorNares = lsUniqueLabels.index("Anterior_nares")
        answer = "\n".join([str(sLAntecubitalFossa)+" 1:0.111111111111 2:0.0 3:0.333333333333 4:0.0 5:0.555555555556",
                            str(sRRetroauricularCrease)+" 1:0.0 2:0.181818181818 3:0.0 4:0.818181818182 5:0.0",
                            str(sLRetroauricularCrease)+" 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0",
                            str(sSubgingivalPlaque)+" 1:0.1 2:0.358333333333 3:0.241666666667 4:0.283333333333 5:0.0166666666667",
                            str(sRAntecubitalFossa)+" 1:0.0 2:0.666666666667 3:0.0 4:0.333333333333 5:0.0",
                            str(sLRetroauricularCrease)+" 1:0.105263157895 2:0.0 3:0.789473684211 4:0.0 5:0.105263157895",
                            str(sRRetroauricularCrease)+" 1:0.0 2:1.0 3:0.0 4:0.0 5:0.0",
                            str(sLAntecubitalFossa)+" 1:0.5 2:0.0 3:0.25 4:0.0 5:0.25",
                            str(sRAntecubitalFossa)+" 1:0.2 2:0.2 3:0.2 4:0.2 5:0.2",
                            str(sAnteriorNares)+" 1:0.0 2:0.0 3:0.0 4:0.0 5:0.0"])+"\n"


        with open(strOutputFile,'r') as f:
            result = f.read()
        f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncScaleFeatureForAllZeros(self):

        npaRow = np.array([0,0,0,0,0,0,0,0,0,0])
        #Answer
        answer = npaRow
        #Get answer
        result = SVM.funcScaleFeature(npaRow)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncScaleFeatureForAllOnes(self):

        npaRow = np.array([1,1,1,1,1,1,1,1,1,1])
        #Answer
        answer = npaRow
        #Get answer
        result = SVM.funcScaleFeature(npaRow)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncScaleFeatureFor01(self):

        npaRow = np.array([0,0,0,0,0,0,0,1,0,0])
        #Answer
        answer = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0])
        #Get answer
        result = SVM.funcScaleFeature(npaRow)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncScaleFeatureFor19(self):

        npaRow = np.array([1,2,3,4,5,6,7,8,9])
        #Answer
        answer =np.array([0.0,1.0/8.0,2.0/8.0,3.0/8.0,4.0/8.0,5.0/8.0,6.0/8.0,7.0/8.0,1.0])
        #Get answer
        result = SVM.funcScaleFeature(npaRow)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncScaleFeatureForRandomish(self):

        npaRow = np.array([1,4,3,6,7,4,5,6,6])
        #Answer
        answer = np.array([0.0,3.0/6,2.0/6,5.0/6.0,1.0,3.0/6.0,4.0/6.0,5.0/6.0,5.0/6.0])
        #Get answer
        result = SVM.funcScaleFeature(npaRow)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncWeightLabelsForWeights(self):

        lsLabels =  ["Four","Three","One","One","One","Two","Two","One","Three","One"]

        #Get answer
        result, ldWeightLabels = SVM.funcWeightLabels(lsLabels)
        result = ",".join(sorted([":".join([ldWeightLabels[key],str(result[key])]) for key in result.keys()]))

        #Answer
        answer = {0: 5.0/1, 1: 5.0/2, 2: 5.0/5, 3: 5.0/2}
        answerWeights = ["Four","Three","One","Two"]
        answer = ",".join(sorted([":".join([answerWeights[key],str(answer[key])]) for key in answer.keys()]))

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFuncWeightLabelsForWeights2(self):

        lsLabels =  ["Four","Three","One","One","One","Two","Two","Three","Three","One"]

        #Get answer
        result, ldWeightLabels = SVM.funcWeightLabels(lsLabels)
        result = ",".join(sorted([":".join([ldWeightLabels[key],str(result[key])]) for key in result.keys()]))

        #Answer
        answer = {0: 4.0/1, 1: 4.0/3, 2: 4.0/4, 3: 4.0/2}
        answerWeights = ["Four","Three","One","Two"]
        answer = ",".join(sorted([":".join([answerWeights[key],str(answer[key])]) for key in answer.keys()]))

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFunc10FoldCrossvalidationForGoodCase10(self):

        #SVM
        svm = SVM()

        #Input sample count
        iSampleCount = 10

        #Randomise
        fRandom = False

        #Result
        answer = [[False, True, True, True, True, True, True, True, True, True],
         [True, False, False, False, False, False, False, False, False, False],
         [True, False, True, True, True, True, True, True, True, True],
         [False, True, False, False, False, False, False, False, False, False],
         [True, True, False, True, True, True, True, True, True, True],
         [False, False, True, False, False, False, False, False, False, False],
         [True, True, True, False, True, True, True, True, True, True],
         [False, False, False, True, False, False, False, False, False, False],
         [True, True, True, True, False, True, True, True, True, True],
         [False, False, False, False, True, False, False, False, False, False],
         [True, True, True, True, True, False, True, True, True, True],
         [False, False, False, False, False, True, False, False, False, False],
         [True, True, True, True, True, True, False, True, True, True],
         [False, False, False, False, False, False, True, False, False, False],
         [True, True, True, True, True, True, True, False, True, True],
         [False, False, False, False, False, False, False, True, False, False],
         [True, True, True, True, True, True, True, True, False, True],
         [False, False, False, False, False, False, False, False, True, False],
         [True, True, True, True, True, True, True, True, True, False],
         [False, False, False, False, False, False, False, False, False, True]]

        #Call generator
        result = []
        for liTraining, liValidation in svm.func10FoldCrossvalidation(iTotalSampleCount=iSampleCount, fRandomise=fRandom):
            result.append(liTraining)
            result.append(liValidation)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFunc10FoldCrossvalidationForGoodCase12(self):

        #SVM
        svm = SVM()

        #Input sample count
        iSampleCount = 12

        #Randomise
        fRandom = False

        #Result
        answer = [[False, True, True, True, True, True, True, True, True, True, False, True],
         [True, False, False, False, False, False, False, False, False, False, True, False],
         [True, False, True, True, True, True, True, True, True, True, True, False],
         [False, True, False, False, False, False, False, False, False, False, False, True],
         [True, True, False, True, True, True, True, True, True, True, True, True],
         [False, False, True, False, False, False, False, False, False, False, False, False],
         [True, True, True, False, True, True, True, True, True, True, True, True],
         [False, False, False, True, False, False, False, False, False, False, False, False],
         [True, True, True, True, False, True, True, True, True, True, True, True],
         [False, False, False, False, True, False, False, False, False, False, False, False],
         [True, True, True, True, True, False, True, True, True, True, True, True],
         [False, False, False, False, False, True, False, False, False, False, False, False],
         [True, True, True, True, True, True, False, True, True, True, True, True],
         [False, False, False, False, False, False, True, False, False, False, False, False],
         [True, True, True, True, True, True, True, False, True, True, True, True],
         [False, False, False, False, False, False, False, True, False, False, False, False],
         [True, True, True, True, True, True, True, True, False, True, True, True],
         [False, False, False, False, False, False, False, False, True, False, False, False],
         [True, True, True, True, True, True, True, True, True, False, True, True],
         [False, False, False, False, False, False, False, False, False, True, False, False]]

        #Call generator
        result = []
        for liTraining, liValidation in svm.func10FoldCrossvalidation(iTotalSampleCount=iSampleCount, fRandomise=fRandom):
            result.append(liTraining)
            result.append(liValidation)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),".Received=",str(result),"."]))

    def testFunc10FoldCrossvalidationForGoodCase10ForTrue(self):

        #SVM
        svm = SVM()

        #Input sample count
        iSampleCount = 10

        #Randomise
        fRandom = True

        #Call generator
        resultTraining = []
        resultValidation = []
        error = ""
        for liTraining, liValidation in svm.func10FoldCrossvalidation(iTotalSampleCount=iSampleCount, fRandomise=fRandom):
            resultTraining.append(liTraining.index(False))
            resultValidation.append(liValidation.index(True))
            if (not sum(liTraining) == 9):
                error = error + " liTraining had the sum of "+str(sum(liTraining))
            if (not sum(liValidation) == 1):
                error = error + " liTraining had the sum of "+str(sum(liValidation))

        lTraining = list(set(resultTraining))
        sorted(lTraining)
        lValidation = list(set(resultValidation))
        sorted(lValidation)

        if(not sum(lTraining) == sum([0,1,2,3,4,5,6,7,8,9])):
            error = error + " Training was not the correct sum. Received:"+str(sum(lTraining))
        if(not sum(lValidation) == sum([0,1,2,3,4,5,6,7,8,9])):
            error = error + " Validation was not the correct sum. Received:"+str(sum(lValidation))
        sorted(resultTraining)
        if(str(lTraining) == str(resultTraining)):
            error = error + " Training was not the correct composition. Received:"+str(lTraining)
        sorted(resultValidation)
        if(str(lValidation) == str(resultValidation)):
            error = error + " Validation was not the correct composition. Received:"+str(lValidation)

        #Check result
        self.assertEqual(error,"",error)

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(SVMTest)
    return suite
