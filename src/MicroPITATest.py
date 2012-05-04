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
from Diversity import Diversity
import mlpy
from MLPYDistanceAdaptor import MLPYDistanceAdaptor
import numpy as np
import os
from MicroPITA import MicroPITA
import re
from SVM import SVM
import unittest
from scikits.learn.cluster import AffinityPropagation
#TODO Get the new import
#from scikits.learn.datasets.samples_generator import make_blobs

##
#Tests the Blog object
class MicroPITATest(unittest.TestCase):

    ##### GetAlphaMetric
    def testGetAlphaMetricForGoodCaseAbridgedData1(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_SIMPSON_A_DIVERSITY

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.getAlphaMetric(tempAbundancies = abundance[sampleNames[0]], tempMetric = metric)

        #Correct Answer
        answer = "0.432098765432"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetAlphaMetricForGoodCaseAbridgedData2(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_INV_SIMPSON_A_DIVERSITY

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.getAlphaMetric(tempAbundancies = abundance[sampleNames[0]], tempMetric = metric)

        #Correct Answer
        answer = "2.31428571429"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetAlphaMetricForGoodCaseAbridgedData3(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_SHANNON_A_DIVERSITY

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.getAlphaMetric(tempAbundancies = abundance[sampleNames[0]], tempMetric = metric)

        #Correct Answer
        answer = "0.936888307539"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetAlphaMetricForGoodCaseAbridgedData4(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        microPITA = MicroPITA()
        metric = microPITA.c_CHAO1_A_DIVERSITY

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.getAlphaMetric(tempAbundancies = abundance[sampleNames[0]], tempMetric = metric)

        #Correct Answer
        answer = "3"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### BuildAlphaMetricsMatrix
    def testBuildAlphaMetricsMatrixForGoodCaseAbridgedData1Metric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_SIMPSON_A_DIVERSITY]

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = metric)

        #Correct Answer
        answer = "[[0.4320987654320988, 0.70247933884297531, 0.0, 0.27736111111111111, 0.55555555555555558, 0.64542936288088648, 1.0, 0.375, 0.20000000000000004, 0.0]]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testBuildAlphaMetricsMatrixForGoodCaseAbridgedDataChaoMetric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        microPITA = MicroPITA()
        metric = [microPITA.c_CHAO1_A_DIVERSITY]

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = metric)

        #Correct Answer
        answer = "[[3, 2, 0, 5, 2, 3, 1, 5.0, 5, 0]]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


    def testBuildAlphaMetricsMatrixForGoodCaseAbridgedData3Metric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_SIMPSON_A_DIVERSITY,microPITA.c_INV_SIMPSON_A_DIVERSITY,microPITA.c_SHANNON_A_DIVERSITY]

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        result = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = metric)

        #Correct Answer
        answer = "[[0.4320987654320988, 0.70247933884297531, 0.0, 0.27736111111111111, 0.55555555555555558, 0.64542936288088648, 1.0, 0.375, 0.20000000000000004, 0.0], [2.3142857142857141, 1.4235294117647057, False, 3.6054081121682522, 1.7999999999999998, 1.5493562231759654, 1.0, 2.6666666666666665, 4.9999999999999991, False], [0.93688830753901586, 0.47413931305783735, 0.0, 1.3667866091157435, 0.63651416829481278, 0.66057888765207562, 0.0, 1.0397207708399179, 1.6094379124341005, 0.0]]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


    ##### GetTopRankedSamples
    def testGetTopRankedSamplesForGoodCase1(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10]]
        N = 3

        #Correct Answer
        answer = [[9,8,7]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase2(self):
        
        #Inputs
        scores = [[10,9,8,7,6,5,4,3,2,1]]
        N = 3

        #Correct Answer
        answer = [[0,1,2]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase3(self):
        
        #Inputs
        scores = [[10,3,40,56,35,678,3,0,-2366]]
        N = 3

        #Correct Answer
        answer = [[5,3,2]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase4(self):
        
        #Inputs
        scores = [[.1,.4,.2,.5,.6,.7,.46,.9]]
        N = 3

        #Correct Answer
        answer = [[7,5,4]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase5(self):
        
        #Inputs
        scores = [[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]]
        N = 3

        #Correct Answer
        answer = [[9,8,7]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase6(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10],[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1],[.1,.4,.2,.5,.6,.7,.46,.9]]
        N = 3

        #Correct Answer
        answer = [[9,8,7],[9,8,7],[7,5,4]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix = scores, tempTopAmount = N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCase6WithSamples(self):
        
        #Inputs
        scores = [[1,2,3,4,5,6,7,8,9,10],[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1],[.1,.4,.2,.5,.6,.7,.46,.9]]
        names = ["Zero","One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]
        N = 3

        #Correct Answer
        answer = [["Nine","Eight","Seven"],["Nine","Eight","Seven"],["Seven","Five","Four"]]

        #Call method
        result = MicroPITA().getTopRankedSamples(tempMatrix=scores, tempSampleNames=names, tempTopAmount=N)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetTopRankedSamplesForGoodCaseAbridgedData3Metric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = [microPITA.c_SIMPSON_A_DIVERSITY,microPITA.c_INV_SIMPSON_A_DIVERSITY,microPITA.c_SHANNON_A_DIVERSITY]

        #Generate data
        abundance = AbundanceTable().textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]

        #Get results
        metrics = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = metric)
        result = microPITA.getTopRankedSamples(tempMatrix = metrics, tempSampleNames = sampleNames, tempTopAmount = 2)

        #Correct Answer
        answer = "[['700037472', '700098984'], ['700037476', '700098980'], ['700037476', '700098980']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### GetBetaMetric Need to test other beta metrics as they come on line
    def testGetBetaMetricForGoodCaseBrayCurtisMetric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_BRAY_CURTIS_B_DIVERSITY

        #Generate data
        rawData = AbundanceTable()
        abundance = rawData.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]
        abundance = rawData.transposeDataMatrix(abundance)
        abundance = abundance[1:4]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.getBetaMetric(tempAbundancies=abundance, tempMetric=metric)

        #Correct Answer
        answer = "[ 1.  1.  1.]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBetaMetricForGoodCaseBrayCurtisMetric4(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_BRAY_CURTIS_B_DIVERSITY

        #Generate data
        rawData = AbundanceTable()
        abundance = rawData.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]
        abundance = rawData.transposeDataMatrix(abundance)
        abundance = abundance[1:5]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.getBetaMetric(tempAbundancies=abundance, tempMetric=metric)

        #Correct Answer
        answer = "[ 1.          1.          0.64166667  1.          0.53484848  1.        ]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBetaMetricForGoodCaseInvBrayCurtisMetric(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY

        #Generate data
        rawData = AbundanceTable()
        abundance = rawData.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]
        abundance = rawData.transposeDataMatrix(abundance)
        abundance = abundance[1:4]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.getBetaMetric(tempAbundancies=abundance, tempMetric=metric)

        #Correct Answer
        answer = "[ 0.  0.  0.]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBetaMetricForGoodCaseInvBrayCurtisMetric4(self):

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = MicroPITA()
        metric = microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY

        #Generate data
        rawData = AbundanceTable()
        abundance = rawData.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        abundance = abundance[0]
        sampleNames = abundance.dtype.names[1:]
        abundance = rawData.transposeDataMatrix(abundance)
        abundance = abundance[1:5]

        #Get results
        #tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
        result = microPITA.getBetaMetric(tempAbundancies=abundance, tempMetric=metric)

        #Correct Answer
        answer = "[ 0.          0.          0.35833333  0.          0.46515152  0.        ]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### GetCentralSamplesByKMedoids
    def testGetCentralSamplesByKMedoidsForGoodCaseBC(self):

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
        result = microPITA.getCentralSamplesByKMedoids(tempMatrix=data, tempMetric=Diversity.c_BRAY_CURTIS_B_DIVERSITY, tempSampleNames = sampleNames, tempNumberClusters = numberClusters, tempNumberSamplesReturned = numberSamplesReturned)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### SelectRepresentiveSamplesFromHClust
    def nottestSelectRepresentativeSamplesFromHClust(self):

        #Micropita object
        microPITA=MicroPITA()

        #Abundance table object to read in and manage data
        rawData = AbundanceTable()

        #Prepare data
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/hq.otu_04-nul-nul-mtd-trn-flt-by-Blood.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/hq.otu_04-nul-nul-mtd-trn-flt-by-Mid_vagina.txt"
#        inputFile = "./testData/microPITA/extremeDissimilarityTest/hq.otu_04-nul-nul-mtd-trn-flt-by-Posterior_fornix.txt"
#        inputFile = "./testData/microPITA/extremeDissimilarityTest/hq.otu_04-nul-nul-mtd-trn-flt-by-Vaginal_introitus.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/lq.phylotype_04-nul-nul-mtd-trn-flt-by-Blood.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/lq.phylotype_04-nul-nul-mtd-trn-flt-by-Mid_vagina.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/lq.phylotype_04-nul-nul-mtd-trn-flt-by-Posterior_fornix.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/lq.phylotype_04-nul-nul-mtd-trn-flt-by-Vaginal_introitus.txt"
##        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
##        inputFile = "./testData/microPITA/extremeDissimilarityTest/lq.phylotype_04-nul-nul-mtd-trn-flt-by-Palatine_Tonsils.txt"
        inputFileMicroPITA = "./input/micropita/src/Testing/Data/microPITA/extremeDissimilarityTest/Insilico.txt"
        inputFile = "./input/micropita/src/Testing/Data/microPITA/extremeDissimilarityTest/InsilicoNorm.txt"

        #Create file names
        prefix = inputFile.splitext()[0]
        figureColorFilePath=prefix+"_ExtremeDissimilarityColor.txt"
        figureLabelFilePath=prefix+"_ExtremeDissimilarityLabel.txt"
        figureDataFile = prefix+"_values.txt"
        outputFileName=prefix+"_ExtremeDissimilarity.pdf"
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        selectCount = 3
        normalize = True
        abundance,metadata = rawData.textToStructuredArray(tempInputFile=inputFileMicroPITA, tempDelimiter=Constants.TAB, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=False)
        sampleNames = abundance.dtype.names[1:]
        abundance = rawData.normalizeColumns(tempStructuredArray=abundance, tempColumns=list(sampleNames))
        tAbundance = rawData.transposeDataMatrix(tempMatrix=abundance, tempRemoveAdornments=True)

        #Beta metric to use
        bMetric=microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY

        #Get samples
        selectedSamples = microPITA.selectRepresentiveSamplesFromHClust(tempBetaMetric=bMetric, tempAbundanceMatrix=tAbundance, tempSampleNames=sampleNames, tempSelectSampleCount=selectCount)
        
        #Generate HClust files
        #If the color file exists, delete
        if(os.path.exists(figureColorFilePath)):
            os.remove(figureColorFilePath)
        #If the label file exists, delete
        if(os.path.exists(figureLabelFilePath)):
            os.remove(figureLabelFilePath)
        #If the data file exists, delete
        if(os.path.exists(figureDataFile)):
            os.remove(figureDataFile)
        #Create file handle to write to files
        with open(figureColorFilePath,'a') as hndlColor:
            colorList = list()
            for selectedSampleName in selectedSamples:
                colorList.append(selectedSampleName+Constants.TAB+"#FF0000")
            hndlColor.write(Constants.ENDLINE.join(colorList))
        hndlColor.close()

        with open(figureLabelFilePath,'a') as hndlLabel:
            labelList = list()
            for selectedSampleName in selectedSamples:
                labelList.append(selectedSampleName+Constants.TAB+selectedSampleName)
            hndlLabel.write(Constants.ENDLINE.join(labelList))
        hndlLabel.close()

        #Create data file
        with open(figureLabelFilePath,'r') as hndlLabel:
            dataContent = hndlLabel.read()
        hndlLabel.close()
        dataContent=dataContent.split(Constants.ENDLINE)
        dataWriteContent = [dataContent[0]]
        if(len(dataContent) > 3):
            dataWriteContent.extend(dataContent[2:(len(dataContent)-1)])

        with open(figureDataFile,'w') as f:
            f.write(Constants.ENDLINE.join(dataWriteContent))
        f.close()

        #Call command
        CommandLine().runCommandLine(["./external/hclust/hclust.py", "--in", figureDataFile, "--out", outputFileName, "--label2cols", figureColorFilePath, "-l", figureLabelFilePath, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "-d", "correlation", "-f", "braycurtis","-y","0.001"])

        #Check result against answer
        result = True
        answer = True
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### reduceToTaxa
    def testReduceToTaxaForGoodCase3(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]
        taxa = ["Bacteria|unclassified|4904","Bacteria|3417","Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0)]"

        #Call method
        result = microPITA.reduceToTaxa(tempMatrix = matrix, tempTaxa = taxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testReduceToTaxaForGoodCase1(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]
        taxa = ["Bacteria|unclassified|4904"]

        #Correct Answer
        answer = "[ ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)]"

        #Call method
        result = microPITA.reduceToTaxa(tempMatrix = matrix, tempTaxa = taxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testReduceToTaxaForGoodCase5(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]
        taxa = ["Bacteria|unclassified|4904","Bacteria|3417","Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72","Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361","Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"]

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0)\n ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)]"

        #Call method
        result = microPITA.reduceToTaxa(tempMatrix = matrix, tempTaxa = taxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testReduceToTaxaForGoodCase0(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]
        taxa = []

        #Correct Answer
        answer = "[]"

        #Call method
        result = microPITA.reduceToTaxa(tempMatrix = matrix, tempTaxa = taxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

#    ##### GetRankAverageSamples
#    def testRankAverageSamplesForGoodCase5(self):
#
#        #Micropita object
#        microPITA = MicroPITA()
#
#        #Inputs
#        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
#        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
#        matrix = matrix[0]
#
#        #Correct Answer
#        answer = "[['700098980', 24.0, 1], ['700037470', 11.4, 2], ['700098984', 11.0, 3], ['700037472', 4.5999999999999996, 4], ['700098986', 1.8, 5], ['700098988', 1.8, 5], ['700037476', 1.0, 6], ['700037474', 0.80000000000000004, 7], ['700098982', 0.0, 8], ['700037478', 0.0, 8]]"
#
#        #Call method
#        result = microPITA.getRankAverageSamples(tempMatrix = matrix)
#
#        #Check result against answer
#        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### GetAverageRankSamples
    def testAverageRankSamplesForGoodCase1(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]

        #Taxa to look for
        lsTaxa = ["Bacteria|3417"]

        #Correct Answer
        answer = "[['700098986', 4.0], ['700037470', 3.0], ['700037474', 3.0], ['700098980', 2.0], ['700098988', 2.0], ['700037472', 2.0], ['700098984', 1.0], ['700098982', 1.0], ['700037476', 1.0], ['700037478', 1.0]]"

        #Call method
        result = microPITA.getAverageRanksSamples(tempMatrix = matrix, tempTargetedTaxa= lsTaxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testAverageRankSamplesForGoodCase2(self):

        #Micropita object
        microPITA = MicroPITA()

        #Inputs
        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
        matrix = matrix[0]

        #Taxa to look for
        lsTaxa = ["Bacteria|unclassified|4904","Bacteria|3417"]

        #Correct Answer
        answer = "[['700098986', 4.0], ['700037470', 3.0], ['700037474', 3.0], ['700098984', 1.5], ['700098980', 1.5], ['700098988', 1.5], ['700037472', 1.5], ['700098982', 1.0], ['700037476', 1.0], ['700037478', 1.0]]"

        #Call method
        result = microPITA.getAverageRanksSamples(tempMatrix = matrix, tempTargetedTaxa= lsTaxa)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

#    def testSelectedTargetedTaxaForGoodCase1(self):

        #Micropita object
#        microPITA = MicroPITA()

        #Inputs
#        inputFile = "./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
#        matrix = AbundanceTable().textToStructuredArray(tempInputFile = inputFile, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 2, tempNormalize = False)
#        matrix = matrix[0]

        #Taxa to look for
#        lsTaxa = ["Bacteria|3417"]

        #Correct Answer
#        answer = "[['700098986', 4.0], ['700037470', 3.0], ['700037474', 3.0], ['700098980', 2.0], ['700098988', 2.0], ['700037472', 2.0], ['700098984', 1.0], ['700098982', 1.0], ['700037476', 1.0], ['700037478', 1.0]]"

        #Call method
#        result = microPITA.selectTargetedTaxaSamples(tempMatrix = matrix, tempTargetedTaxa= lsTaxa)

        #Check result against answer
#        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    ##### GetRandomSamples
    def testGetRandomSamplesForGoodCase1of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 1

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().getRandomSamples(tempSamples = samples, tempNumberOfSamplesToReturn = N)

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

    def testGetRandomSamplesForGoodCase3of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 3

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().getRandomSamples(tempSamples = samples, tempNumberOfSamplesToReturn = N)

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

    def testGetRandomSamplesForGoodCase10of10Samples(self):
        
        #Inputs
        samples = ["one","two","three","four","five","six","seven","eight","nine","ten"]
        N = 10

        #Correct Answer
        answer = False

        #Tracking error information
        foundError = False
        errorString = ""

        #Call method
        result = MicroPITA().getRandomSamples(tempSamples = samples, tempNumberOfSamplesToReturn = N)

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

    def testFuncStratifyDataByMetadataForGoodCase3Metadata(self):
        
        #Inputs
        tdata = [("Name1",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),("Name2",16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),("Name3",31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)]
        dtype=[('id','a5'),('x1',int),('x2',int),('x3',int),('x4',int),('x5',int),('x6',int),('x7',int),('x8',int),('x9',int),('x10',int),('x11',int),('x12',int),('x13',int),('x14',int),('x15',int)]
        rawAbundance = np.array(tdata, dtype=dtype)

        lsMetadata = ["three","three","one","one","two","one","one","three","two","one","one","two","two","three","two"]

        #Correct Answer
        answer = "{'one': array([(3, 4, 6, 7, 10, 11),(18, 19, 21, 22, 25, 26),(33, 34, 36, 37, 40, 41)], "
        answer = answer + "dtype=[('x3', '<i8'), ('x4', '<i8'), ('x6', '<i8'), ('x7', '<i8'), ('x10', '<i8'), ('x11', '<i8')]),"
        answer = answer + "'three': array([(1, 2, 8, 14), (16, 17, 23, 29), (31, 32, 38, 44)], "
        answer = answer + "dtype=[('x1', '<i8'), ('x2', '<i8'), ('x8', '<i8'), ('x14', '<i8')]),"
        answer = answer + "'two': array([(5, 9, 12, 13, 15), (20, 24, 27, 28, 30), (35, 39, 42, 43, 45)], "
        answer = answer + "dtype=[('x5', '<i8'), ('x9', '<i8'), ('x12', '<i8'), ('x13', '<i8'), ('x15', '<i8')])}"
        answer = re.sub(r'\s+',"",answer)

        #Call method
        result = MicroPITA().funcStratifyByMetadata(lsMetadata, rawAbundance)
        result = re.sub(r'\s+',"", str(result))

        #Check result against answer
        self.assertEqual(result,str(answer).strip(' \t\n\r'),"".join([str(self),"::Expected="+str(answer)+".but received="+str(result)+"."]))

    ##### RunSVM
    def nottestRunSVM(self):

        #Inputs
        #Reading file
        inputFile="./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        outputFile="./input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.SVM.txt"
        delimiter=Constants.TAB
        nameRow=0
        firstDataRow=2
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
        microPITA.runSVM(tempInputFile=inputFile, tempDelimiter=Constants.TAB, tempOutputSVMFile=outputFile, tempMatrixLabels=labels, tempFirstDataRow=firstDataRow, tempSkipFirstColumn=skipColumn1, tempNormalize=normalize, tempSVMScaleLowestBound=lowestScaleBound, tempSVMLogG=gRange, tempSVMLogC=cRange, tempSVMProbabilistic=probabilistic)

        #Get results
        result = ""

        #Correct Answer
        answer = "[['700037472', '700098984'], ['700037476', '700098980'], ['700037476', '700098980']]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroPITATest)
    return suite
