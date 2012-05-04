#######################################################
#
#	Title:		Utility_MathTest
#	Author:		Timothy Tickle 
#	Date:		03/26/2012
#	Purpose:	Test Utility_Math class
#
#######################################################

#Import local code
from AbundanceTable import AbundanceTable
from Constants import Constants
from Constants_Testing import Constants_Testing
import numpy as np
import random
import unittest
from Utility_Math import Utility_Math

class Utility_MathTest(unittest.TestCase):

    ##Set up for tests
    def setUp(self): pass

    def testConvertToBHQValueForGoodCaseR(self):
        methodName = "testConvertToBHFDRForGoodCaseR"

        ldInputPvalues = [0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,
                          0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,
                          0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,
                          0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,
                          0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,
                          0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,
                          0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,
                          0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,
                          0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,
                          0.99,1.00]

        #R was used to produce the following as a result
        #pvalues = seq(0,1,0.01)
        #> pvalues
        #  [1] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14
        # [16] 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29
        # [31] 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44
        # [46] 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59
        # [61] 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74
        # [76] 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89
        # [91] 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00
        #> p.adjust(p=pvalues, method="BH")
        #  [1] 0.0000000 0.5050000 0.6733333 0.7575000 0.8080000 0.8416667 0.8657143
        #  [8] 0.8837500 0.8977778 0.9090000 0.9181818 0.9258333 0.9323077 0.9378571
        # [15] 0.9426667 0.9468750 0.9505882 0.9538889 0.9568421 0.9595000 0.9619048
        # [22] 0.9640909 0.9660870 0.9679167 0.9696000 0.9711538 0.9725926 0.9739286
        # [29] 0.9751724 0.9763333 0.9774194 0.9784375 0.9793939 0.9802941 0.9811429
        # [36] 0.9819444 0.9827027 0.9834211 0.9841026 0.9847500 0.9853659 0.9859524
        # [43] 0.9865116 0.9870455 0.9875556 0.9880435 0.9885106 0.9889583 0.9893878
        # [50] 0.9898000 0.9901961 0.9905769 0.9909434 0.9912963 0.9916364 0.9919643
        # [57] 0.9922807 0.9925862 0.9928814 0.9931667 0.9934426 0.9937097 0.9939683
        # [64] 0.9942188 0.9944615 0.9946970 0.9949254 0.9951471 0.9953623 0.9955714
        # [71] 0.9957746 0.9959722 0.9961644 0.9963514 0.9965333 0.9967105 0.9968831
        # [78] 0.9970513 0.9972152 0.9973750 0.9975309 0.9976829 0.9978313 0.9979762
        # [85] 0.9981176 0.9982558 0.9983908 0.9985227 0.9986517 0.9987778 0.9989011
        # [92] 0.9990217 0.9991398 0.9992553 0.9993684 0.9994792 0.9995876 0.9996939
        # [99] 0.9997980 0.9999000 1.0000000

        strCorrectAnswer = "".join(["[0.0, 0.505, 0.6733333333333333, 0.7575, 0.808, 0.8416666666666668, ",
          "0.8657142857142857, 0.88375, 0.8977777777777778, 0.909, 0.9181818181818183, 0.9258333333333333, ",
          "0.9323076923076923, 0.937857142857143, 0.9426666666666667, 0.9468749999999999, 0.9505882352941176, ",
          "0.953888888888889, 0.9568421052631578, 0.9595, 0.961904761904762, 0.9640909090909091, 0.966086956521739, ",
          "0.9679166666666666, 0.9695999999999999, 0.9711538461538461, 0.9725925925925927, 0.9739285714285716, ",
          "0.9751724137931035, 0.9763333333333333, 0.9774193548387096, 0.9784375, 0.9793939393939394, 0.9802941176470588, ",
          "0.9811428571428572, 0.9819444444444443, 0.9827027027027027, 0.9834210526315789, 0.9841025641025641, 0.98475, ",
          "0.9853658536585367, 0.9859523809523809, 0.9865116279069768, 0.9870454545454546, 0.9875555555555555, ",
          "0.9880434782608696, 0.9885106382978723, 0.9889583333333333, 0.9893877551020408, 0.9898, 0.9901960784313726, ",
          "0.990576923076923, 0.9909433962264151, 0.9912962962962963, 0.9916363636363638, 0.9919642857142857, ",
          "0.992280701754386, 0.9925862068965516, 0.9928813559322034, 0.9931666666666666, 0.9934426229508195, ",
          "0.9937096774193548, 0.993968253968254, 0.99421875, 0.9944615384615385, 0.9946969696969697, 0.9949253731343283, ",
          "0.9951470588235294, 0.9953623188405798, 0.9955714285714286, 0.9957746478873238, 0.9959722222222221, ",
          "0.9961643835616438, 0.9963513513513514, 0.9965333333333333, 0.9967105263157895, 0.9968831168831169, ",
          "0.997051282051282, 0.9972151898734177, 0.9973750000000001, 0.997530864197531, 0.9976829268292683, ",
          "0.9978313253012048, 0.9979761904761905, 0.9981176470588236, 0.9982558139534883, 0.9983908045977011, ",
          "0.9985227272727273, 0.9986516853932583, 0.9987777777777778, 0.998901098901099, 0.9990217391304348, ",
          "0.9991397849462366, 0.9992553191489363, 0.9993684210526316, 0.9994791666666666, 0.9995876288659793, ",
          "0.9996938775510203, 0.9997979797979798, 0.9998999999999999, 1.0]"])

        ldResults = Utility_Math.convertToBHQValue(ldInputPvalues)
        self.assertEquals(strCorrectAnswer, str(ldResults), methodName+" did not give correct result. Expected ."+strCorrectAnswer+". but received ."+str(ldResults)+".")

    def testConvertToBHQValueForGoodCase2(self):
        methodName = "testConvertToBHFDRForGoodCase2"

        ldInputPvalues = [0.0120,0.0330,0.2120,0.9000,0.9800,0.0010,0.9990,0.0003,.00001]

        #R was used to produce the following as a result
        #> pvalues2 = c(0.0120,0.0330,0.2120,0.9000,0.9800,0.0010,0.9990,0.0003,.00001)
        #> p.adjust(p=pvalues2, method="BH")
        #[1] 0.02700 0.05940 0.31800 0.99900 0.99900 0.00300 0.99900 0.00135 0.00009

        strCorrectAnswer = "".join(["[0.027, 0.05940000000000001, 0.318, 1.157142857142857, 1.1025, 0.0030000000000000005, ",
                                    "0.999, 0.0013499999999999999, 9e-05]"])

        ldResults = Utility_Math.convertToBHQValue(ldInputPvalues)
        self.assertEquals(str(strCorrectAnswer), str(ldResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(ldResults)+".")

    def testFuncSampleWithReplacementForGoodCase1(self):
      methodName = "testFuncSampleWithReplacementForGoodCase1"
      iDataLength = 20

      random.seed("NotRandom")
      aData = list(range(iDataLength))
      iSelect = 10

      strCorrectAnswer = str([2, 19, 5, 4, 12, 18, 4, 10, 5, 12])
      lResults = str(list(Utility_Math.funcSampleWithReplacement(aData, iSelect)))
      self.assertEquals(str(strCorrectAnswer), str(lResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(lResults)+".")

    def testFuncSampleWithReplacementForGoodCase2(self):
      methodName = "testFuncSampleWithReplacementForGoodCase2"

      iDataLength = 1

      random.seed("NotRandom")
      aData = list(range(iDataLength))
      iSelect = 8

      strCorrectAnswer = str([0, 0, 0, 0, 0, 0, 0, 0])
      lResults = str(list(Utility_Math.funcSampleWithReplacement(aData, iSelect)))
      self.assertEquals(str(strCorrectAnswer), str(lResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(lResults)+".")

    def testFuncSampleWithReplacementForGoodCase3(self):
      methodName = "testFuncSampleWithReplacementForGoodCase3"
      iDataLength = 1

      random.seed("NotRandom")
      aData = list(range(iDataLength))
      iSelect = 0

      strCorrectAnswer = str([])
      lResults = str(list(Utility_Math.funcSampleWithReplacement(aData, iSelect)))
      self.assertEquals(str(strCorrectAnswer), str(lResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(lResults)+".")

    def testFuncSampleWithReplacementForGoodCase4(self):
      methodName = "testFuncSampleWithReplacementForGoodCase4"
      iDataLength = 10

      random.seed("NotRandom")
      aData = []
      iSelect = 3

      strCorrectAnswer = str([])
      lResults = str(list(Utility_Math.funcSampleWithReplacement(aData, iSelect)))
      self.assertEquals(str(strCorrectAnswer), str(lResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(lResults)+".")

    def testFuncSampleWithReplacementForGoodCase5(self):
      methodName = "testFuncSampleWithReplacementForGoodCase5"
      iDataLength = 10

      random.seed("NotRandom")
      aData = None
      iSelect = 3

      strCorrectAnswer = str([])
      lResults = str(list(Utility_Math.funcSampleWithReplacement(aData, iSelect)))
      self.assertEquals(str(strCorrectAnswer), str(lResults), methodName+" did not give correct result. Expected ."+str(strCorrectAnswer)+". but received ."+str(lResults)+".")

    def testTransposeDataMatrixForGoodCase(self):
        
      #Inputs
      inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
      delimiter = Constants.TAB
      nameRow = 0
      firstDataRow = 2
      normalize = False
      removeAdornment = False

      #Correct Answer
      answer = "[[ 'Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72'\n  'Bacteria|unclassified|4904'\n  'Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361'\n  'Bacteria|3417'\n  'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368']\n ['1.0' '0.0' '3.0' '0.0' '5.0']\n ['0.0' '10.0' '0.0' '45.0' '0.0']\n ['0.0' '0.0' '0.0' '0.0' '0.0']\n ['12.0' '43.0' '29.0' '34.0' '2.0']\n ['0.0' '6.0' '0.0' '3.0' '0.0']\n ['6.0' '0.0' '45.0' '0.0' '6.0']\n ['0.0' '23.0' '0.0' '0.0' '0.0']\n ['2.0' '0.0' '1.0' '0.0' '1.0']\n ['1.0' '1.0' '1.0' '1.0' '1.0']\n ['0.0' '0.0' '0.0' '0.0' '0.0']]"

      #Call method
      result = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)
      result = Utility_Math.transposeDataMatrix(tempMatrix=result.funcGetAbundanceCopy(), tempRemoveAdornments=removeAdornment)

      #Check result against answer
      self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTransposeDataMatrixForGoodCaseRemoveAdornments(self):
        
      #Inputs
      inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
      delimiter = Constants.TAB
      nameRow = 0
      firstDataRow = 2
      removeAdornment = True

      #Correct Answer
      answer = "[[  1.   0.   3.   0.   5.]\n [  0.  10.   0.  45.   0.]\n [  0.   0.   0.   0.   0.]\n [ 12.  43.  29.  34.   2.]\n [  0.   6.   0.   3.   0.]\n [  6.   0.  45.   0.   6.]\n [  0.  23.   0.   0.   0.]\n [  2.   0.   1.   0.   1.]\n [  1.   1.   1.   1.   1.]\n [  0.   0.   0.   0.   0.]]"

      #Call method
      result = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)
      result = Utility_Math.transposeDataMatrix(tempMatrix=result.funcGetAbundanceCopy(), tempRemoveAdornments=removeAdornment)

      #Check result against answer
      self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Create a suite to be called to test
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(Utility_MathTest)
    return suite

#if __name__=='__main__':
#    unittest.main()
