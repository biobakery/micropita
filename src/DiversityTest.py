#######################################################
# Author: Timothy Tickle
# Description: Class to test the Diversity class
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
import Diversity
import FileIO
import numpy as np
import unittest

##
#Tests the Diversity object
class DiversityTest(unittest.TestCase):

    def testGetSimpsonsDiversityIndexForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.008,0.0,0.015,0.019,0.000,0.001,0.056,0.031,0.003,0.00,0.010,0.0,0.037,0.04,0.012,0.022,0.051,0.254,0.018,0.0,0.171,0.032,0.031,0.015,0.0,0.037,0.008,0.002,0.001,0.009,0.0,0.0,0.066,0.0,0.031,0.0,0.0,0.0,0.007])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"simpson")
        answer = "0.1144"

        #Call method
        result = Diversity.Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetSimpsonsDiversityIndexForGoodCase2(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.987])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"simpson")
        answer = "0.974338"

        #Call method
        result = Diversity.Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetSimpsonsDiversityIndexForGoodCase3(self):
        
        #Inputs
        sampleAbundancies = np.array([0.0])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"simpson")
        answer = 0.0

        #Call method
        result = Diversity.Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetSimpsonsDiversityIndexForGoodCase4(self):
        
        #Inputs
        sampleAbundancies = np.array([1.0])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"simpson")
        answer = 1.0

        #Call method
        result = Diversity.Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseSimpsonsDiversityIndexForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.008,0.0,0.015,0.019,0.000,0.001,0.056,0.031,0.003,0.00,0.010,0.0,0.037,0.04,0.012,0.022,0.051,0.254,0.018,0.0,0.171,0.032,0.031,0.015,0.0,0.037,0.008,0.002,0.001,0.009,0.0,0.0,0.066,0.0,0.031,0.0,0.0,0.0,0.007])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"inv")
        answer = 1/0.1144

        #Call method
        result = Diversity.Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseSimpsonsDiversityIndexForGoodCase2(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.987])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"inv")
        answer = 1/0.974338

        #Call method
        result = Diversity.Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseSimpsonsDiversityIndexForGoodCase3(self):
        
        #Inputs
        sampleAbundancies = np.array([0.0])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"inv")
        answer = False

        #Call method
        result = Diversity.Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseSimpsonsDiversityIndexForGoodCase4(self):
        
        #Inputs
        sampleAbundancies = np.array([1.0])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"inv")
        answer = 1.0

        #Call method
        result = Diversity.Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetShannonDiversityIndexForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.008,0.0,0.015,0.019,0.000,0.001,0.056,0.031,0.003,0.00,0.010,0.0,0.037,0.04,0.012,0.022,0.051,0.254,0.018,0.0,0.171,0.032,0.031,0.015,0.0,0.037,0.008,0.002,0.001,0.009,0.0,0.0,0.066,0.0,0.031,0.0,0.0,0.0,0.007])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"shannon")
        answer = "2.65980629671"

        #Call method
        result = Diversity.Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetShannonDiversityIndexForGoodCase2(self):
        
        #Inputs
        sampleAbundancies = np.array([0.013,0.987])

        #Correct Answer
        #Checked against the R package Vegan using diversity(data,"shannon")
        answer = "0.0693716084143"

        #Call method
        result = Diversity.Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetShannonDiversityIndexForGoodCase3(self):
        
        #Inputs
        sampleAbundancies = np.array([0.0])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"shannon")
        answer = 0.0

        #Call method
        result = Diversity.Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetShannonDiversityIndexForGoodCase4(self):
        
        #Inputs
        sampleAbundancies = np.array([1.0])

        #Correct Answer
        #Checked against the R package Vegan using 1-diversity(data,"shannon")
        answer = 0.0

        #Call method
        result = Diversity.Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBrayCurtisDissimilarityForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")
        answer = "[ 0.33333333]"

        #Call method
        result = Diversity.Diversity.getBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBrayCurtisDissimilarityForGoodCase2(self):
        
        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0], [2.0,3.0,4.0,5.0,6.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")

        #0-1,0-2,1-2
        answer = "[ 0.33333333  0.14285714  0.2       ]"

        #Call method
        result = Diversity.Diversity.getBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetBrayCurtisDissimilarityForGoodCase3(self):

        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0], [2.0,3.0,4.0,5.0,6.0], [5.0,4.0,3.0,2.0,1.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")
        #answer = [sum(abs(sampleAbundancies[0]-sampleAbundancies[1]))/sum(sampleAbundancies[0]+sampleAbundancies[1]), sum(abs(sampleAbundancies[0]-sampleAbundancies[2]))/sum(sampleAbundancies[0]+sampleAbundancies[2]),sum(abs(sampleAbundancies[0]-sampleAbundancies[3]))/sum(sampleAbundancies[0]+sampleAbundancies[3]),sum(abs(sampleAbundancies[1]-sampleAbundancies[2]))/sum(sampleAbundancies[1]+sampleAbundancies[2]),sum(abs(sampleAbundancies[1]-sampleAbundancies[3]))/sum(sampleAbundancies[1]+sampleAbundancies[3]),sum(abs(sampleAbundancies[2]-sampleAbundancies[3]))/sum(sampleAbundancies[2]+sampleAbundancies[3])]
        answer = "[ 0.33333333  0.14285714  0.4         0.2         0.46666667  0.37142857]"

        #Call method
        result = Diversity.Diversity.getBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseBrayCurtisDissimilarityForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")
        answer = "[ 0.66666667]"

        #Call method
        result = Diversity.Diversity.getInverseBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseBrayCurtisDissimilarityForGoodCase2(self):
        
        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0], [2.0,3.0,4.0,5.0,6.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")

        #0-1,0-2,1-2
        answer = "[ 0.66666667  0.85714286  0.8       ]"

        #Call method
        result = Diversity.Diversity.getInverseBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetInverseBrayCurtisDissimilarityForGoodCase3(self):

        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0], [2.0,3.0,4.0,5.0,6.0], [5.0,4.0,3.0,2.0,1.0]])

        #Correct Answer
        #Checked against the R package Vegan using vegdist(data,method="bray")
        #answer = [sum(abs(sampleAbundancies[0]-sampleAbundancies[1]))/sum(sampleAbundancies[0]+sampleAbundancies[1]), sum(abs(sampleAbundancies[0]-sampleAbundancies[2]))/sum(sampleAbundancies[0]+sampleAbundancies[2]),sum(abs(sampleAbundancies[0]-sampleAbundancies[3]))/sum(sampleAbundancies[0]+sampleAbundancies[3]),sum(abs(sampleAbundancies[1]-sampleAbundancies[2]))/sum(sampleAbundancies[1]+sampleAbundancies[2]),sum(abs(sampleAbundancies[1]-sampleAbundancies[3]))/sum(sampleAbundancies[1]+sampleAbundancies[3]),sum(abs(sampleAbundancies[2]-sampleAbundancies[3]))/sum(sampleAbundancies[2]+sampleAbundancies[3])]
        answer = "[ 0.66666667  0.85714286  0.6         0.8         0.53333333  0.62857143]"

        #Call method
        result = Diversity.Diversity.getInverseBrayCurtisDissimilarity(tempSampleTaxaAbundancies = sampleAbundancies)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetChao1DiversityIndexForGoodCase1NoBias(self):
        
        #Inputs
        sampleAbundancies = np.array([1.0,2.0,3.0,4.0,5.0,0.0,1.0,0.0,1.0,2.5,2.0,4.0,6.0,8.0,10.0,9.0])
        correctBias = False

        #Correct Answer
        #Checked against the R package Fossil using chao1(data)
        answer = 14.0+((3.0**2.0)/(2.0*2.0))
        answer = 16.25

        #Call method
        result = Diversity.Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies, tempCorrectForBias = correctBias)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetChao1DiversityIndexForGoodCase2NoBias(self):
        
        #Inputs
        sampleAbundancies = np.array([1.0])
        correctBias = False

        #Correct Answer
        #Checked against the R package Fossil using chao1(data)
        answer = 1.0

        #Call method
        result = Diversity.Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies, tempCorrectForBias = correctBias)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testGetChao1DiversityIndexForGoodCase3NoBias(self):
        
        #Inputs
        sampleAbundancies = np.array([0.0])
        correctBias = False

        #Correct Answer
        #Checked against the R package Fossil using chao1(data)
        answer = 0.0

        #Call method
        result = Diversity.Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies = sampleAbundancies, tempCorrectForBias = correctBias)

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #This is running the example at http://pycogent.sourceforge.net/examples/unifrac.html 10-21-2011
    def testGetUnifracDistanceForGoodCase1(self):
        
        #Inputs
        sampleAbundancies = np.array([[1.0,2.0,3.0,4.0,5.0], [2.0,4.0,6.0,8.0,10.0]])
        taxTree = "(B:0.2,(C:0.3,D:0.4)E:0.6)F;"
        #Original Abundance matrix
        #Rows = Sample
        #Cols = Sequence
        #[10,11,0]
        #[2,0,9]
        #[2,2,2]
        #Translates to this dictionary
        envs = Constants.Constants.TEST_MICROPITA_DOCUMENTS+"unifracTestData/testQIIMEOTUAbundance.txt"
#        envs = {'B':{'sample1':10, 'sample2':2, 'sample3':2},
#                'C':{'sample1':11,'sample2':0, 'sample3':2},
#                'D':{'sample1':0, 'sample2':9, 'sample3':2}}

        #Correct Answer
        answer = """(array([[ 0.        ,  0.46666667,  0.26666667],
       [ 0.46666667,  0.        ,  0.2       ],
       [ 0.26666667,  0.2       ,  0.        ]]), ['sample1', 'sample2', 'sample3'])"""

        #Call method
        result = Diversity.Diversity.getUnifracDistance(tempSampleTaxaAbundancies=envs, tempTaxonomyTree = taxTree, tempWeighted=False)
        result = str(result['distance_matrix'])

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::\nExpected=\n",str(answer),".\nReceived=\n",str(result),"."]))

    #This is an attempt to run unifrac with HMP data
    def nottestGetUnifracDistanceForGoodCaseHMP(self):
        
        #Newick tree
        taxTree = Constants.Constants.INPUT_DATA_DIRECTORY+"HMPNewickTreeQiimeFormat/rep_set_v35-NoQuote.tre"
        readTree = FileIO.FileIO(taxTree,True,False,False)
        taxTree = readTree.readFullFile()
        readTree.close()

        #Translates to this dictionary
        envs = Constants.Constants.INPUT_DATA_DIRECTORY+"HMPQiimeFormatAbundanceTable/otu_table_psn_v35.red.txt"

        #Call method
        result = Diversity.Diversity.getUnifracDistance(tempSampleTaxaAbundancies=envs, tempTaxonomyTree = taxTree, tempWeighted=False)
        result = "Ran"
        answer = "Ran"

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::\nExpected=\n",str(answer),".\nReceived=\n",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(DiversityTest)
    return suite
