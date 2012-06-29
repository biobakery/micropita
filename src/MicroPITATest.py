#######################################################
# Author: Timothy Tickle
# Description: Class to test the RunMicroPITA class
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
from CommandLine import CommandLine
from Constants import Constants
from Constants_Testing import Constants_Testing
from Diversity import Diversity
import mlpy
from MLPYDistanceAdaptor import MLPYDistanceAdaptor
import numpy as np
import os
from MicroPITA import MicroPITA
import re
from SVM import SVM
import unittest
from Utility_Math import Utility_Math
from scikits.learn.cluster import AffinityPropagation
#TODO Get the new import
#from scikits.learn.datasets.samples_generator import make_blobs

##
#Tests the Blog object
class MicroPITATest(unittest.TestCase):

#####Test funcGetTopRankedSamples
    def testfuncGetTopRankedSamplesForGoodCase1(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10]]
        N = 3

        #Correct Answer
        answer = [[9,8,7]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase2(self):
        
        #Inputs
        scores = [[10,9,8,7,6,5,4,3,2,1]]
        N = 3

        #Correct Answer
        answer = [[0,1,2]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase3(self):
        
        #Inputs
        scores = [[10,3,40,56,35,678,3,0,-2366]]
        N = 3

        #Correct Answer
        answer = [[5,3,2]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase4(self):
        
        #Inputs
        scores = [[.1,.4,.2,.5,.6,.7,.46,.9]]
        N = 3

        #Correct Answer
        answer = [[7,5,4]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase5(self):
        
        #Inputs
        scores = [[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]]
        N = 3

        #Correct Answer
        answer = [[9,8,7]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase6(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10],[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1],[.1,.4,.2,.5,.6,.7,.46,.9]]
        N = 3

        #Correct Answer
        answer = [[9,8,7],[9,8,7],[7,5,4]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix = scores, iTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCase6WithSamples(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10],[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1],[.1,.4,.2,.5,.6,.7,.46,.9]]
        names = ["Zero","One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]
        N = 3

        #Correct Answer
        answer = [["Nine","Eight","Seven"],["Nine","Eight","Seven"],["Seven","Five","Four"]]

        #Call method
        result = MicroPITA().funcGetTopRankedSamples(lldMatrix=scores, lsSampleNames=names, iTopAmount=N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData1InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 1

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "['Sample_1']"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData1InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 1

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "['Sample_1']"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData1InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 1

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "['Sample_1']"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData1InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 1

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_1']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData2InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 2

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_1', 'Sample_2']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData3InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 3

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_1', 'Sample_2', 'Sample_3']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData4InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 4

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData5InvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 5

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcNormalize()
        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData1Choa(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 1

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData2Choa(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 2

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50', 'Sample_49']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData3Choa(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 3

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50', 'Sample_49', 'Sample_48']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData4Choa(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 4

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50', 'Sample_49', 'Sample_48', 'Sample_47']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData5Choa(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 5

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50', 'Sample_49', 'Sample_48', 'Sample_47', 'Sample_46']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetTopRankedSamplesForGoodCaseAbridgedData5ChoaInvSimpson(self):

        #Inputs
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_strChao1Diversity,microPITA.c_strInverseSimpsonDiversity]

        #Generate data
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/DiversityTest.pcl"])
        delimiter = Constants.TAB
        sNameRow = "ID"
        sLastMetadata = "Inverse_Simpson"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        iSelectCount = 5

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        abundance = abndData.funcGetAbundanceCopy()
        sampleNames = abndData.funcGetSampleNames()

        #Get results
        metrics = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abundance, lsSampleNames = sampleNames, lsDiversityMetricAlpha = metric)
        result = microPITA.funcGetTopRankedSamples(lldMatrix = metrics, lsSampleNames = sampleNames, iTopAmount = iSelectCount)

        #Correct Answer
        answer = "[['Sample_50', 'Sample_49', 'Sample_48', 'Sample_47', 'Sample_46'], ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


#####Test funcGetBetaMetric Need to test other beta metrics as they come on line
    def testfuncGetBetaMetricForGoodCaseBrayCurtisMetric(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        rawData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        microPITA = MicroPITA()
        metric = microPITA.c_strBrayCurtisDissimilarity

        #Generate data
        abundance = rawData.funcGetAbundanceCopy()
        sampleNames = rawData.funcGetSampleNames()
        abundance = Utility_Math.funcTransposeDataMatrix(abundance)
        abundance = abundance[1:4]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.funcGetBetaMetric(npadAbundancies=abundance, sMetric=metric)

        #Correct Answer
        answer = "[ 1.  1.  1.]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetBetaMetricForGoodCaseBrayCurtisMetric4(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        rawData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        microPITA = MicroPITA()
        metric = microPITA.c_strBrayCurtisDissimilarity

        #Generate data
        abundance = rawData.funcGetAbundanceCopy()
        sampleNames = rawData.funcGetSampleNames()
        abundance = Utility_Math.funcTransposeDataMatrix(abundance)
        abundance = abundance[1:5]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.funcGetBetaMetric(npadAbundancies=abundance, sMetric=metric)

        #Correct Answer
        answer = "[ 1.          1.          0.90697674  1.          0.49714286  1.        ]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetBetaMetricForGoodCaseInvBrayCurtisMetric(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        microPITA = MicroPITA()
        metric = microPITA.c_strInvBrayCurtisDissimilarity

        rawData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        #Generate data
        abundance = rawData.funcGetAbundanceCopy()
        sampleNames = rawData.funcGetSampleNames()
        abundance = Utility_Math.funcTransposeDataMatrix(abundance)
        abundance = abundance[1:4]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.funcGetBetaMetric(npadAbundancies=abundance, sMetric=metric)

        #Correct Answer
        answer = "[ 0.  0.  0.]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncGetBetaMetricForGoodCaseInvBrayCurtisMetric4(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        microPITA = MicroPITA()
        metric = microPITA.c_strInvBrayCurtisDissimilarity

        rawData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        #Generate data
        abundance = rawData.funcGetAbundanceCopy()
        sampleNames = rawData.funcGetSampleNames()
        abundance = Utility_Math.funcTransposeDataMatrix(abundance)
        abundance = abundance[1:5]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.funcGetBetaMetric(npadAbundancies=abundance, sMetric=metric)

        #Correct Answer
        answer = "[ 0.          0.          0.09302326  0.          0.50285714  0.        ]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

#####Test funcGetCentralSamplesByKMedoids
    def testfuncGetCentralSamplesByKMedoidsForGoodCaseBC(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        data = np.array([[0.2,0.8],[1.0,0.0],[0.1,0.9],[.95,0.05],[0.5,0.5],[0.45,0.55],[0.25,0.75],[0.40,0.60],[0.83,0.17]])
        sampleNames = ("One","Two","Three","Four","Five","Six","Seven","Eight","Nine")
        # [0.1,0.9], 3
        #[[0.2,0.8], 1
        # [0.25,0.75], 7
        #
        # [0.5,0.5], 5
        # [0.45,0.55], 6
        # [0.40,0.60], 8
        #
        # [0.83,0.17]] 9
        # [.95,0.05], 4
        # [1.0,0.0], 2

        numberClusters = 3
        numberSamplesReturned = 3

        #Correct Answer
        answer = "['Six', 'One', 'Four']"

        #Call method
        result = microPITA.funcGetCentralSamplesByKMedoids(npaMatrix=data, sMetric=Diversity.c_strBrayCurtisDissimilarity, lsSampleNames = sampleNames, iNumberSamplesReturned = numberSamplesReturned)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

###Test selectExtremeSamplesFromHClust
    def testFuncSelectExtremeSamplesForGoodCase1Different(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six"]
        iSelectSampleCount = 1
        npaAbundanceMatrix = np.array([[1,0,0,0,0],[0,1,1,1,1],[0,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)
        result.sort()

        answer = ["One"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase2Different(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six"]
        iSelectSampleCount = 2
        npaAbundanceMatrix = np.array([[1,0,0,0,0],[0,1,1,1,1],[0,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)
        result.sort()

        answer = ["One","Three"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase1Different1Groups(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four"]
        iSelectSampleCount = 1
        npaAbundanceMatrix = np.array([[1,0,0,0],
                                       [1,1,1,0],
                                       [1,1,0,0],
                                       [0,1,1,0]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)

        result.sort()

        answer = ["One"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase2GradientGroups2Samples(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]
        iSelectSampleCount = 2
        npaAbundanceMatrix = np.array([[1,0,0,0,0],
                                       [1,1,0,0,0],
                                       [1,1,1,0,0],
                                       [1,1,1,1,0],
                                       [1,1,1,1,1],
                                       [0,1,1,1,1],
                                       [0,0,1,1,1],
                                       [0,0,0,1,1],
                                       [0,0,0,0,1]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)

        result.sort()

        answer = ["Four","Nine"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase2GradientGroups4Samples(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]
        iSelectSampleCount = 4
        npaAbundanceMatrix = np.array([[1,0,0,0,0],
                                       [1,1,0,0,0],
                                       [1,1,1,0,0],
                                       [1,1,1,1,0],
                                       [1,1,1,1,1],
                                       [0,1,1,1,1],
                                       [0,0,1,1,1],
                                       [0,0,0,1,1],
                                       [0,0,0,0,1]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)

        result.sort()

        answer = ["Three","Four","Eight","Nine"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase2GradientGroups6Samples(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]
        iSelectSampleCount = 6
        npaAbundanceMatrix = np.array([[1,0,0,0,0],
                                       [1,1,0,0,0],
                                       [1,1,1,0,0],
                                       [1,1,1,1,0],
                                       [1,1,1,1,1],
                                       [0,1,1,1,1],
                                       [0,0,1,1,1],
                                       [0,0,0,1,1],
                                       [0,0,0,0,1]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)

        result.sort()

        answer = ["Two","Three","Four","Seven","Eight","Nine"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectExtremeSamplesForGoodCase3Different3Groups(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strBetaMetric = microPITA.c_strInvBrayCurtisDissimilarity
        lsSampleNames = ["One","Two","Three","Four","Five","Six"]
        iSelectSampleCount = 3
        npaAbundanceMatrix = np.array([[1,0,0,0],
                                       [0,0,0,1],
                                       [1,0,0,1],
                                       [1,1,0,0],
                                       [0,0,1,1],
                                       [1,1,1,1]])

        #Get results
        result = microPITA.funcSelectExtremeSamplesFromHClust(strBetaMetric=strBetaMetric, npaAbundanceMatrix=npaAbundanceMatrix,
                                                 lsSampleNames=lsSampleNames, iSelectSampleCount=iSelectSampleCount)

        result.sort()

        answer = ["One","Four","Five"]
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


###Test funcGetAverageAbundanceSamples
    def testfuncGetAverageAbundanceSamplesForGoodCase1Feature(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]
        fRanked = False

        answer= [["700098980","12.0",1],["700037470","6.0",1],["700098986","1.0",1],["700098988","1.0",1],["700098982","0.0",1]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCase1FeatureRanked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]
        fRanked = True

        answer= [["700037470","2.0",6.0],["700098986","3.0",1.0],["700098980","4.0",12.0],["700098988","5.0",1.0],["700098982","9.0",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCase2Feature(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        fRanked = False

        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361"]
        answer= [["700037470","25.5",1],["700098980","20.5",1],["700098986","2.0",1],["700098988","2.0",1],["700098982","0.0",1]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCase2FeatureRanked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        fRanked = True

        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361"]
        answer= [["700037470","1.5",25.5],["700098986","2.5",2.0],["700098980","3.5",20.5],["700098988","4.0",2.0],["700098982","9.0",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCaseAllFeature(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        fRanked = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                      "Bacteria|unclassified|4904","Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368",
                      "Bacteria|3417"]
        answer= [["700098980","24.0",1],["700037470","11.4",1],["700098988","3.0",1],["700098986","1.8",1],["700098982","0.0",1]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCaseAllFeatureRankedWithTie(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        fRanked = True
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                      "Bacteria|unclassified|4904","Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368",
                      "Bacteria|3417"]
        answer= [["700037470","2.2",11.4],["700098986","2.8",1.8],["700098980","3.0",24.0],["700098988","3.0",3.0],["700098982","9.0",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncGetAverageAbundanceSamplesForGoodCaseAllFeatureRankedWithTies2(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        fRanked = True
        liFeatures = ["Bacteria|unclassified|4904"]
        answer= [["700098980","1.0",43.0],["700098988","4.0",2.0],["700098986","6.0",0.0],["700037470","7.0",0.0],["700098982","9.0",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcGetAverageAbundanceSamples(abndTable=abndData, lsTargetedFeature=liFeatures, fRank=fRanked)
        result = [[resultList[0],"{0:.1f}".format(resultList[1]),resultList[2]] for resultList in result]

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

###Test funcSelectTargetedTaxaSamples
    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect1(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 1
        sMethod = MicroPITA.c_strTargetedAbundance
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700098980"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect2(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 2
        sMethod = MicroPITA.c_strTargetedAbundance
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700098980","700037470"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect3(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 3
        sMethod = MicroPITA.c_strTargetedAbundance
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700098980","700037470","700098986"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect4(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 4
        sMethod = MicroPITA.c_strTargetedAbundance
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700098980","700037470","700098986","700098988"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect5(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 5
        sMethod = MicroPITA.c_strTargetedAbundance
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700098980","700037470","700098986","700098988","700098982"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect1Ranked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 1
        sMethod = MicroPITA.c_strTargetedRanked
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700037470"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect2Ranked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 2
        sMethod = MicroPITA.c_strTargetedRanked
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700037470","700098986"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect3Ranked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 3
        sMethod = MicroPITA.c_strTargetedRanked
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700037470","700098986","700098980"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect4Ranked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 4
        sMethod = MicroPITA.c_strTargetedRanked
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700037470","700098986","700098980","700098988"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testfuncSelectTargetedTaxaSamplesForGoodCase1FeatureSelect5Ranked(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        iSampleCount = 5
        sMethod = MicroPITA.c_strTargetedRanked
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        answer= ["700037470","700098986","700098980","700098988","700098982"]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = microPITA.funcSelectTargetedTaxaSamples(abndMatrix=abndData, lsTargetedTaxa=liFeatures, iSampleSelectionCount=iSampleCount, sMethod=sMethod)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

##### Test funcGetRandomSamples
    def testfuncGetRandomSamplesForGoodCase1of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 1

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().funcGetRandomSamples(lsSamples = samples, iNumberOfSamplesToReturn = N)

        #Check to make sure only N samples are returned
        samplesReturnLength=len(result)

        #Make sure the samples are different (this assumes the samples were initially unique in the test).
        uniqueCounts = len(set(result))

        #Evaluate
        if(samplesReturnLength != N):
            foundError = True
            errorString = errorString + "Should have been "+str(N)+" but received "+str(samplesReturnLength)+"."
        if(uniqueCounts != N):
            foundError = True
            errorString = errorString + "Should have had a unique count of "+str(N)+" but instead received an element count of "+str(uniqueCounts)+". Elements = "+str(result)

        #Check result against answer
        self.assertEqual(str(foundError),str(answer),"".join([str(self),"::",str(errorString),"."]))

    def testfuncGetRandomSamplesForGoodCase3of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 3

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().funcGetRandomSamples(lsSamples = samples, iNumberOfSamplesToReturn = N)

        #Check to make sure only N samples are returned
        samplesReturnLength=len(result)

        #Make sure the samples are different (this assumes the samples were initially unique in the test).
        uniqueCounts = len(set(result))

        #Evaluate
        if(samplesReturnLength != N):
            foundError = True
            errorString = errorString + "Should have been "+str(N)+" but received "+str(samplesReturnLength)+"."
        if(uniqueCounts != N):
            foundError = True
            errorString = errorString + "Should have had a unique count of "+str(N)+" but instead received an element count of "+str(uniqueCounts)+". Elements = "+str(result)

        #Check result against answer
        self.assertEqual(str(foundError),str(answer),"".join([str(self),"::",str(errorString),"."]))

    def testfuncGetRandomSamplesForGoodCase10of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 10

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().funcGetRandomSamples(lsSamples = samples, iNumberOfSamplesToReturn = N)

        #Check to make sure only N samples are returned
        samplesReturnLength=len(result)

        #Make sure the samples are different (this assumes the samples were initially unique in the test).
        uniqueCounts = len(set(result))

        #Evaluate
        if(samplesReturnLength != N):
            foundError = True
            errorString = errorString + "Should have been "+str(N)+" but received "+str(samplesReturnLength)+"."
        if(uniqueCounts != N):
            foundError = True
            errorString = errorString + "Should have had a unique count of "+str(N)+" but instead received an element count of "+str(uniqueCounts)+". Elements = "+str(result)

        #Check result against answer
        self.assertEqual(str(foundError),str(answer),"".join([str(self),"::",str(errorString),"."]))

    ##### funcRunSVM
    def nottestfuncRunSVM(self):

        #Inputs
        #Reading file
        inputFile="./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        outputFile="./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.SVM.txt"
        delimiter=Constants.TAB
        nameRow=0
        sLastMetadataName="STSite"
        skipColumn1=True
        normalize=True
        lowestScaleBound=0
        probabilistic=True
        cRange="-5,-4,-3,-2,-1,0,1,2,3,4,5"
        gRange="-5,-4,-3,-2,-1,0,1,2,3,4,5"

        #Inputs  programmatic
        microPITA = MicroPITA()
        labels = [0,0,0,0,0,1,1,1,1,1]

        #Generate data
        microPITA.funcRunSVM(tempInputFile=inputFile, tempDelimiter=Constants.TAB, tempOutputSVMFile=outputFile, tempMatrixLabels=labels, sLastMetadataName=sLastMetadataName, tempSkipFirstColumn=skipColumn1, tempNormalize=normalize, tempSVMScaleLowestBound=lowestScaleBound, tempSVMLogG=gRange, tempSVMLogC=cRange, tempSVMProbabilistic=probabilistic)

        #Get results
        result = ""

        #Correct Answer
        answer = "[['700037472', '700098984'], ['700037476', '700098980'], ['700037476', '700098980']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncRunSVMForGoodCase(self):

        answer = "."
        result = ""

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### Test runLIBSVM Will be removed
### Test runMLPYSVM
    def testFuncRunMLPYSVMForGoodCase(self):

        answer = "."
        result = ""

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncStoreSVMProbabilityForGoodCaseNoPriorStateProbabilistic(self):

        #Input
        lsValidationSamples = ["Sample1","Sample2","Sample3","Sample4","Sample5"]
        lSVMLabels = [0,1,2]
        ldValidationLabels = [2,1,0,1,2]
        lPredictions = [0,1,2,1,1]
        npaDistances = [[0.0,.7,.3],[.0,.5,.5],[0.0,0.0,1.0],[.4,.1,.5],[.9,.05,.05]]
        dictdProbability = dict()
        dictAllProbabilities = dict()
        dictiPrediction = dict()
        dictAllPredictions = dict()
        fClassifyByProbability = True

        #Make results
        lreturn = MicroPITA()._funcStoreSVMProbability(lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,
                                                       dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)
        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = lreturn

        result = " ".join([str([key,dictdProbability[key]]) for key in sorted(dictdProbability.keys())])
        result = result+", "+" ".join([str([key,dictAllProbabilities[key]]) for key in sorted(dictAllProbabilities.keys())])
        result = result+", "+" ".join([str([key,dictiPrediction[key]]) for key in sorted(dictiPrediction.keys())])
        result = result+", "+" ".join([str([key,dictAllPredictions[key]]) for key in sorted(dictAllPredictions.keys())])

        #Correct answers
        dictdProbabilityAnswer = {"Sample1":.7,"Sample2":.5,"Sample3":1.0,"Sample4":.5,"Sample5":.9}
        dictAllProbabilitiesAnswer = {"Sample1":[1,0.0,.7,.3],"Sample2":[1,.0,.5,.5],"Sample3":[2,0.0,0.0,1.0],"Sample4":[2,.4,.1,.5],"Sample5":[0,.9,.05,.05]}
        dictiPredictionAnswer = {"Sample1":False,"Sample2":True,"Sample3":False,"Sample4":False,"Sample5":False}
        dictAllPredictionsAnswer = {"Sample1":"1","Sample2":"1","Sample3":"2","Sample4":"2","Sample5":"0"}

        answer = " ".join([str([key,dictdProbabilityAnswer[key]]) for key in sorted(dictdProbabilityAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllProbabilitiesAnswer[key]]) for key in sorted(dictAllProbabilitiesAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictiPredictionAnswer[key]]) for key in sorted(dictiPredictionAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllPredictionsAnswer[key]]) for key in sorted(dictAllPredictionsAnswer.keys())])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncStoreSVMProbabilityForGoodCaseNoPriorState(self):

        #Input
        lsValidationSamples = ["Sample1","Sample2","Sample3","Sample4","Sample5"]
        lSVMLabels = [0,1,2]
        ldValidationLabels = [2,1,0,1,2]
        lPredictions = [0,1,2,1,1]
        npaDistances = [[0.0,.7,.3],[.0,.5,.5],[0.0,0.0,1.0],[.4,.1,.5],[.9,.05,.05]]
        dictdProbability = dict()
        dictAllProbabilities = dict()
        dictiPrediction = dict()
        dictAllPredictions = dict()
        fClassifyByProbability = False

        #Make results
        lreturn = MicroPITA()._funcStoreSVMProbability(lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,
                                                       dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)
        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = lreturn

        result = " ".join([str([key,dictdProbability[key]]) for key in sorted(dictdProbability.keys())])
        result = result+", "+" ".join([str([key,dictAllProbabilities[key]]) for key in sorted(dictAllProbabilities.keys())])
        result = result+", "+" ".join([str([key,dictiPrediction[key]]) for key in sorted(dictiPrediction.keys())])
        result = result+", "+" ".join([str([key,dictAllPredictions[key]]) for key in sorted(dictAllPredictions.keys())])

        #Correct answers
        dictdProbabilityAnswer = {"Sample1":.0,"Sample2":.5,"Sample3":1.0,"Sample4":.1,"Sample5":.05}
        dictAllProbabilitiesAnswer = {"Sample1":[0,0.0,.7,.3],"Sample2":[1,.0,.5,.5],"Sample3":[2,0.0,0.0,1.0],"Sample4":[1,.4,.1,.5],"Sample5":[1,.9,.05,.05]}
        dictiPredictionAnswer = {"Sample1":False,"Sample2":True,"Sample3":False,"Sample4":True,"Sample5":False}
        dictAllPredictionsAnswer = {"Sample1":"0","Sample2":"1","Sample3":"2","Sample4":"1","Sample5":"1"}

        answer = " ".join([str([key,dictdProbabilityAnswer[key]]) for key in sorted(dictdProbabilityAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllProbabilitiesAnswer[key]]) for key in sorted(dictAllProbabilitiesAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictiPredictionAnswer[key]]) for key in sorted(dictiPredictionAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllPredictionsAnswer[key]]) for key in sorted(dictAllPredictionsAnswer.keys())])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncStoreSVMProbabilityForGoodCaseWithPriorStateProbabilistic(self):

        #Input
        lsValidationSamples = ["Sample1","Sample2","Sample3","Sample4","Sample5"]
        lSVMLabels = [0,1,2]
        ldValidationLabels = [2,1,0,1,2]
        lPredictions = [0,1,2,1,1]
        npaDistances = [[0.0,.7,.3],[.0,.5,.5],[0.0,0.0,1.0],[.4,.1,.5],[.9,.05,.05]]
        dictdProbability = {"Sample0":.99,"Sample10":.4}
        dictAllProbabilities = {"Sample0":[.99,.1,.0],"Sample10":[.1,.4,.5]}
        dictiPrediction = {"Sample0":True,"Sample10":False}
        dictAllPredictions = {"Sample0":0,"Sample10":1}
        fClassifyByProbability = True

        #Make results
        lreturn = MicroPITA()._funcStoreSVMProbability(lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,
                                                       dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)
        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = lreturn

        result = " ".join([str([key,dictdProbability[key]]) for key in sorted(dictdProbability.keys())])
        result = result+", "+" ".join([str([key,dictAllProbabilities[key]]) for key in sorted(dictAllProbabilities.keys())])
        result = result+", "+" ".join([str([key,dictiPrediction[key]]) for key in sorted(dictiPrediction.keys())])
        result = result+", "+" ".join([str([key,dictAllPredictions[key]]) for key in sorted(dictAllPredictions.keys())])

        #Correct answers
        dictdProbabilityAnswer = {"Sample0":.99,"Sample10":.4,"Sample1":.7,"Sample2":.5,"Sample3":1.0,"Sample4":.5,"Sample5":.9}
        dictAllProbabilitiesAnswer = {"Sample0":[.99,.1,.0],"Sample10":[.1,.4,.5],"Sample1":[1,0.0,.7,.3],"Sample2":[1,.0,.5,.5],"Sample3":[2,0.0,0.0,1.0],"Sample4":[2,.4,.1,.5],"Sample5":[0,.9,.05,.05]}
        dictiPredictionAnswer = {"Sample0":True,"Sample10":False,"Sample1":False,"Sample2":True,"Sample3":False,"Sample4":False,"Sample5":False}
        dictAllPredictionsAnswer = {"Sample0":0,"Sample10":1,"Sample1":"1","Sample2":"1","Sample3":"2","Sample4":"2","Sample5":"0"}

        answer = " ".join([str([key,dictdProbabilityAnswer[key]]) for key in sorted(dictdProbabilityAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllProbabilitiesAnswer[key]]) for key in sorted(dictAllProbabilitiesAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictiPredictionAnswer[key]]) for key in sorted(dictiPredictionAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllPredictionsAnswer[key]]) for key in sorted(dictAllPredictionsAnswer.keys())])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncStoreSVMProbabilityForGoodCaseWithPriorState(self):

        #Input
        lsValidationSamples = ["Sample1","Sample2","Sample3","Sample4","Sample5"]
        lSVMLabels = [0,1,2]
        ldValidationLabels = [2,1,0,1,2]
        lPredictions = [0,1,2,1,1]
        npaDistances = [[0.0,.7,.3],[.0,.5,.5],[0.0,0.0,1.0],[.4,.1,.5],[.9,.05,.05]]
        dictdProbability = {"Sample0":.99,"Sample10":.4}
        dictAllProbabilities = {"Sample0":[.99,.1,.0],"Sample10":[.1,.4,.5]}
        dictiPrediction = {"Sample0":True,"Sample10":False}
        dictAllPredictions = {"Sample0":0,"Sample10":1}
        fClassifyByProbability = False

        #Make results
        lreturn = MicroPITA()._funcStoreSVMProbability(lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,
                                                       dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)
        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = lreturn

        result = " ".join([str([key,dictdProbability[key]]) for key in sorted(dictdProbability.keys())])
        result = result+", "+" ".join([str([key,dictAllProbabilities[key]]) for key in sorted(dictAllProbabilities.keys())])
        result = result+", "+" ".join([str([key,dictiPrediction[key]]) for key in sorted(dictiPrediction.keys())])
        result = result+", "+" ".join([str([key,dictAllPredictions[key]]) for key in sorted(dictAllPredictions.keys())])

        #Correct answers
        dictdProbabilityAnswer = {"Sample0":.99,"Sample10":.4,"Sample1":.0,"Sample2":.5,"Sample3":1.0,"Sample4":.1,"Sample5":.05}
        dictAllProbabilitiesAnswer = {"Sample0":[.99,.1,.0],"Sample10":[.1,.4,.5],"Sample1":[0,0.0,.7,.3],"Sample2":[1,.0,.5,.5],"Sample3":[2,0.0,0.0,1.0],"Sample4":[1,.4,.1,.5],"Sample5":[1,.9,.05,.05]}
        dictiPredictionAnswer = {"Sample0":True,"Sample10":False,"Sample1":False,"Sample2":True,"Sample3":False,"Sample4":True,"Sample5":False}
        dictAllPredictionsAnswer = {"Sample0":0,"Sample10":1,"Sample1":"0","Sample2":"1","Sample3":"2","Sample4":"1","Sample5":"1"}

        answer = " ".join([str([key,dictdProbabilityAnswer[key]]) for key in sorted(dictdProbabilityAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllProbabilitiesAnswer[key]]) for key in sorted(dictAllProbabilitiesAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictiPredictionAnswer[key]]) for key in sorted(dictiPredictionAnswer.keys())])
        answer = answer+", "+" ".join([str([key,dictAllPredictionsAnswer[key]]) for key in sorted(dictAllPredictionsAnswer.keys())])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### test runSupervisedMethods
    def testFuncRunSupervisedMethodsForGoodCase(self):

        answer = "."
        result = ""

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### test _funcSelectSupervisedSamplesFromPredictFile
    def testFuncSelectSupervisedSamplesFromPredictFile2Classes(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 1
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Eight","One"],microPITA.c_strSVMFar:["Three","Seven"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.464303 0.535697",
                                                      "1 0.510597 0.489403",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "0 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "0 1:0.514484",
                                                      "1 1:0.38026 2:0.61974",
                                                      "0 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesDistinct(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = False
        fSelectDistinct = True
        iSelectCount = 1
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMFar:["Three","Seven"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.464303 0.535697",
                                                      "1 0.510597 0.489403",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "0 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "0 1:0.514484",
                                                      "1 1:0.38026 2:0.61974",
                                                      "0 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesDiscriminant(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = False
        iSelectCount = 1
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Eight","One"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.464303 0.535697",
                                                      "1 0.510597 0.489403",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "0 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "0 1:0.514484",
                                                      "1 1:0.38026 2:0.61974",
                                                      "0 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturn8(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 2
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Eight","One","Six","Nine"],microPITA.c_strSVMFar:["Four","Five","Seven","Ten"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.464303 0.535697",
                                                      "1 0.510597 0.489403",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.514484",
                                                      "0 1:0.38026 2:0.61974",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnAll(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 2
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["One","Four","Six","Seven"],microPITA.c_strSVMFar:["Two","Three","Five","Eight"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnMoreThanExists(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 9
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["One","Four","Six","Seven","Two","Three","Five","Eight"],microPITA.c_strSVMFar:["One","Four","Six","Seven","Two","Three","Five","Eight"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnMoreThanExistsMislabeled(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 9
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["One","Four","Six","Seven","Five","Eight"],microPITA.c_strSVMFar:["One","Four","Six","Seven","Five","Eight"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "1 1:0.464303 2:0.535697",
                                                      "0 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnMoreThanExistsClassMislabeled(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 9
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Three","Five","Six","Seven"],microPITA.c_strSVMFar:["Three","Five","Six","Seven"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["1 1:0.482089 2:0.517911",
                                                      "1 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "1 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "1 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnMoreThanExistsClassMislabeled2(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 9
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["One","Two","Four","Eight"],microPITA.c_strSVMFar:["One","Two","Four","Eight"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "0 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "0 1:0.388271 2:0.611729",
                                                      "0 1:0.5 2:0.5",
                                                      "0 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile2ClassesReturnLess0(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 2
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Two","Six","Seven"],microPITA.c_strSVMFar:["Two","Three","Five"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0",
                                                      "0 0.482089 0.517911",
                                                      "0 0.409344 0.590656",
                                                      "1 0.99 0.01",
                                                      "0 0.5 0.5",
                                                      "1 0.988271 0.011729",
                                                      "1 0.5 0.5",
                                                      "1 0.510373 0.489627",
                                                      "0 0.410257 0.589743"])
        strInputFileContents = Constants.ENDLINE.join(["1 1:0.482089 2:0.517911",
                                                      "0 1:0.464303 2:0.535697",
                                                      "1 1:0.510597",
                                                      "1 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "1 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "1 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile3Classes(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 1
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Four","Seven","Six"],microPITA.c_strSVMFar:["One","Five","Three"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0 2",
                                                      "0 0.1 0.5",
                                                      "1 0.2 0.8",
                                                      "2 0.99 0.01",
                                                      "0 0.3 0.3",
                                                      "1 0.988271 0.011729",
                                                      "2 0.3 0.3",
                                                      "1 0.5 0.4",
                                                      "0 0.4 0.5"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "1 1:0.464303 2:0.535697",
                                                      "2 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "2 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile3ClassesSelect2(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 2
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"]

        #Answer
        answer = {microPITA.c_strSVMClose:["Four","Eight","Seven","Two","Six","Three"],microPITA.c_strSVMFar:["One","Eight","Five","Two","Three","Six"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0 2",
                                                      "0 0.1 0.5",
                                                      "1 0.2 0.8",
                                                      "2 0.99 0.01",
                                                      "0 0.3 0.3",
                                                      "1 0.988271 0.011729",
                                                      "2 0.3 0.3",
                                                      "1 0.5 0.4",
                                                      "0 0.4 0.5"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "1 1:0.464303 2:0.535697",
                                                      "2 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "2 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncSelectSupervisedSamplesFromPredictFile3ClassesSelect3(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        strPredictFilePath = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromPredictionFile.txt"])
        strOriginalInputFile = "".join([Constants_Testing.c_strTestingTMP,"TestFuncSelectSupervisedSamplesFromInputFile.txt"])
        fSelectDiscriminant = True
        fSelectDistinct = True
        iSelectCount = 3
        lsSampleNames = ["One","Two","Three","Four","Five","Six","Seven","Eight"]

        #Answer
        answer = {microPITA.c_strSVMClose:["One","Two","Three","Four","Five","Six","Seven","Eight"],microPITA.c_strSVMFar:["One","Two","Three","Four","Five","Six","Seven","Eight"]}

        #Set up. Write rpedict and input files
        strPredFileContents = Constants.ENDLINE.join(["labels 1 0 2",
                                                      "0 0.1 0.5",
                                                      "1 0.2 0.8",
                                                      "2 0.99 0.01",
                                                      "0 0.3 0.3",
                                                      "1 0.988271 0.011729",
                                                      "2 0.3 0.3",
                                                      "1 0.5 0.4",
                                                      "0 0.4 0.5"])
        strInputFileContents = Constants.ENDLINE.join(["0 1:0.482089 2:0.517911",
                                                      "1 1:0.464303 2:0.535697",
                                                      "2 1:0.510597",
                                                      "0 1:0.409344 2:0.590656",
                                                      "1 1:0.388271 2:0.611729",
                                                      "2 1:0.5 2:0.5",
                                                      "1 2:0.489627",
                                                      "0 1:0.410257 2:0.589743"])

        with open(strPredictFilePath, 'w') as f, open(strOriginalInputFile, 'w') as g:
            f.write(strPredFileContents)
            g.write(strInputFileContents)

        #Get answer
        result = microPITA._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=strOriginalInputFile, strPredictFilePath=strPredictFilePath,
                                                                         lsSampleNames=lsSampleNames, iSelectCount=iSelectCount, 
                                                                         fSelectDiscriminant=fSelectDiscriminant, fSelectDistinct=fSelectDistinct)

        #Standardize answer and result
        answer = [[strKey,sorted(answer[strKey])] for strKey in answer]
        result = [[strKey,sorted(result[strKey])] for strKey in result]

        #Delete predict file
        for strFile in [strPredictFilePath,strOriginalInputFile]:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Check result
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### Test run
    def testFuncRunForGoodCase(self):

        answer = "."
        result = ""

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### Test funcWriteSelectionToFile
    def testFuncWriteSelectionToFileForGoodCase(self):

        #Micropita object
        microPITA = MicroPITA()

        dictTest = {"Diversity_C":["Sample_0_D","Sample_1_D","Sample_2_D","Sample_3_D","Sample_4_D","Sample_5_D"],
		"Distinct":["Sample_41_E","Sample_42_E","Sample_43_E","Sample_45_T","Sample_46_T","Sample_47_T"],
		"Extreme_B":["Sample_7_D","Sample_38_E","Sample_8_D","Sample_43_E","Sample_6_D","Sample_39_E"],
		"Discriminant":["Sample_3_D","Sample_5_D","Sample_6_D","Sample_0_D","Sample_1_D","Sample_2_D"],
		"Representative_B":["Sample_38_E","Sample_39_E","Sample_40_E","Sample_43_E","Sample_44_T","Sample_47_T"],
		"Diversity_I":["Sample_45_T","Sample_44_T","Sample_46_T","Sample_13_D","Sample_9_D","Sample_2_D"],
		"Taxa_Defined":["Sample_47_T","Sample_46_T","Sample_44_T","Sample_45_T","Sample_24_R","Sample_19_R"]}
        lsKeys = ["Diversity_C","Distinct","Extreme_B","Discriminant","Representative_B","Diversity_I","Taxa_Defined"]
        sTestFile = "".join([Constants_Testing.c_strTestingTMP,"TempTestSelectFile.txt"])
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"TestSelectFile.txt"])
        answer = ""

        if os.path.exists(sTestFile):
            os.remove(sTestFile)

        microPITA.funcWriteSelectionToFile(dictSelection=dictTest,strOutputFilePath=sTestFile)

        #Read in generated file and answer
        result = ""
        with open(sTestFile) as f, open(sAnswerFile) as g:
            result = f.read()
            answer = g.read()

        if os.path.exists(sTestFile):
            os.remove(sTestFile)

        #Put answer in correct order
        dictresult = dict([(sLine.split(Constants.COLON)[0],sLine.split(Constants.COLON)[1]) for sLine in filter(None,result.split(Constants.ENDLINE))])
        result = Constants.ENDLINE.join([Constants.COLON.join([sKey,dictresult[sKey]]) for sKey in lsKeys])+Constants.ENDLINE

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

### Test funcReadSelectionFileToDictionary
    def testFuncReadSelectionFileToDictionaryForGoodCase(self):

        #Micropita object
        microPITA = MicroPITA()

        dictTest = {"Diversity_C":["Sample_0_D","Sample_1_D","Sample_2_D","Sample_3_D","Sample_4_D","Sample_5_D"],
		"Distinct":["Sample_41_E","Sample_42_E","Sample_43_E","Sample_45_T","Sample_46_T","Sample_47_T"],
		"Extreme_B":["Sample_7_D","Sample_38_E","Sample_8_D","Sample_43_E","Sample_6_D","Sample_39_E"],
		"Discriminant":["Sample_3_D","Sample_5_D","Sample_6_D","Sample_0_D","Sample_1_D","Sample_2_D"],
		"Representative_B":["Sample_38_E","Sample_39_E","Sample_40_E","Sample_43_E","Sample_44_T","Sample_47_T"],
		"Diversity_I":["Sample_45_T","Sample_44_T","Sample_46_T","Sample_13_D","Sample_9_D","Sample_2_D"],
		"Taxa_Defined":["Sample_47_T","Sample_46_T","Sample_44_T","Sample_45_T","Sample_24_R","Sample_19_R"]}
        lsKeys = ["Diversity_C","Distinct","Extreme_B","Discriminant","Representative_B","Diversity_I","Taxa_Defined"]
        sTestFile = "".join([Constants_Testing.c_strTestingTruth,"TestSelectFile.txt"])
        answer = "".join(["".join([sKey,str(dictTest[sKey])]) for sKey in lsKeys])

        dictResults = microPITA.funcReadSelectionFileToDictionary(sTestFile)

        #Put answer in correct order
        result = "".join(["".join([sKey,str(dictResults[sKey])]) for sKey in lsKeys])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncReadSelectionFileToDictionaryAndFuncWriteSelectionToFileForGoodCase(self):

        #Micropita object
        microPITA = MicroPITA()

        dictTest = {"Diversity_C":["Sample_0_D","Sample_1_D","Sample_2_D","Sample_3_D","Sample_4_D","Sample_5_D"],
		"Distinct":["Sample_41_E","Sample_42_E","Sample_43_E","Sample_45_T","Sample_46_T","Sample_47_T"],
		"Extreme_B":["Sample_7_D","Sample_38_E","Sample_8_D","Sample_43_E","Sample_6_D","Sample_39_E"],
		"Discriminant":["Sample_3_D","Sample_5_D","Sample_6_D","Sample_0_D","Sample_1_D","Sample_2_D"],
		"Representative_B":["Sample_38_E","Sample_39_E","Sample_40_E","Sample_43_E","Sample_44_T","Sample_47_T"],
		"Diversity_I":["Sample_45_T","Sample_44_T","Sample_46_T","Sample_13_D","Sample_9_D","Sample_2_D"],
		"Taxa_Defined":["Sample_47_T","Sample_46_T","Sample_44_T","Sample_45_T","Sample_24_R","Sample_19_R"]}
        lsKeys = ["Diversity_C","Distinct","Extreme_B","Discriminant","Representative_B","Diversity_I","Taxa_Defined"]
        sTestFile = "".join([Constants_Testing.c_strTestingTMP,"TempTestSelectFile.txt"])
        answer = "".join(["".join([sKey,str(dictTest[sKey])]) for sKey in lsKeys])

        #Get result
        microPITA.funcWriteSelectionToFile(dictSelection=dictTest,strOutputFilePath=sTestFile)
        dictResults = microPITA.funcReadSelectionFileToDictionary(sTestFile)
        #Put answer in correct order
        result = "".join(["".join([sKey,str(dictResults[sKey])]) for sKey in lsKeys])

        if os.path.exists(sTestFile):
            os.remove(sTestFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroPITATest)
    return suite
