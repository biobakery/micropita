#######################################################
#
#	Title:		ValidateDataTest
#	Author:		Timothy Tickle 
#	Date:		December 8, 2009
#	Purpose:	Test Validate Data containing methods for testing variables
#
#######################################################

#Import local code
import unittest
import numpy as np
from Constants import Constants
from ValidateData import ValidateData

class ValidateDataTest(unittest.TestCase):

    ##Set up for tests
    def setUp(self): pass

    def testIsValidFormatStringForBadCaseNone(self):
        formatString = None
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),False,"Should have received false for a input"+str(formatString)+".")

    def testIsValidFormatStringForBadCaseBlank(self):
        formatString = "  "
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),False,"Should have received false for a input"+str(formatString)+".")

    def testIsValidFormatStringForBadCaseEmpty(self):
        formatString = ""
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),False,"Should have received false for a input"+str(formatString)+".")

    def testIsValidFormatStringForBadCaseWrongType(self):
        formatString = list()
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),False,"Should have received false for a input"+str(formatString)+".")

    def testIsValidFormatStringForGoodCase1(self):
        formatString = "BbcdfHhIiLlPpsx0123456789"
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),True,"Should have received true for a input"+str(formatString)+".")

    def testIsValidFormatStringForBadCaseBadChar(self):
        formatString = "BbcdfHhIiLlPpWQZsx0123456789"
        self.assertEqual(ValidateData.funcIsValidFormatString(formatString),False,"Should have received false for a input"+str(formatString)+".")

    def testIsValidBooleanForGoodCase(self):
        self.assertEqual(ValidateData.funcIsValidBoolean(False),True,"Should have received true for a valid boolean input.")

    def testIsValidBooleanForString(self):
        self.assertEqual(ValidateData.funcIsValidBoolean("False"),False,"Should have received false for a invalid string input.")

    def testIsValidBooleanForNone(self):
        self.assertEqual(ValidateData.funcIsValidBoolean(None),False,"Should have received false for a invalid none input.")

    def testIsValidBooleanForZero(self):
        self.assertEqual(ValidateData.funcIsValidBoolean(0),False,"Should have received false for a invalid boolean input.")

    def testIsValidBooleanForOne(self):
        self.assertEqual(ValidateData.funcIsValidBoolean(1),False,"Should have received false for a invalid boolean input.")

    def testIsTrueForGoodCase(self):
        self.assertEqual(ValidateData.funcIsTrue(True),True,"Should have received true for a valid boolean input.")

    def testIsTrueForFalseCase(self):
        self.assertEqual(ValidateData.funcIsTrue(False),False,"Should have received false for false boolean input.")

    def testIsTrueForString(self):
        self.assertEqual(ValidateData.funcIsTrue("False"),False,"Should have received false for a invalid string input.")

    def testIsTrueForNone(self):
        self.assertEqual(ValidateData.funcIsTrue(None),False,"Should have received false for a invalid none input.")

    def testIsTrueForZero(self):
        self.assertEqual(ValidateData.funcIsTrue(0),False,"Should have received false for a invalid boolean input.")

    def testIsTrueForOne(self):
        self.assertEqual(ValidateData.funcIsTrue(1),False,"Should have received false for a invalid boolean input.")

    def testIsFalseForGoodCase(self):
        self.assertEqual(ValidateData.funcIsFalse(False),True,"Should have received true for a valid boolean input.")

    def testIsFalseForTrueCase(self):
        self.assertEqual(ValidateData.funcIsFalse(True),False,"Should have received false for true boolean input.")

    def testIsFalseForString(self):
        self.assertEqual(ValidateData.funcIsFalse("False"),False,"Should have received false for a invalid string input.")

    def testIsFalseForNone(self):
        self.assertEqual(ValidateData.funcIsFalse(None),False,"Should have received false for a invalid none input.")

    def testIsFalseForZero(self):
        self.assertEqual(ValidateData.funcIsFalse(0),False,"Should have received false for a invalid boolean input.")

    def testIsFalseForOne(self):
        self.assertEqual(ValidateData.funcIsFalse(1),False,"Should have received false for a invalid boolean input.")

    def testIsValidStringForGoodCase(self):
        self.assertEqual(ValidateData.funcIsValidString("Valid"),True,"Should have received true for a valid string input.")

    def testIsValidStringForNoneCase(self):
        self.assertEqual(ValidateData.funcIsValidString(None),False,"Should have received false for a none string input.")

    def testIsValidStringForBlankCase(self):
        self.assertEqual(ValidateData.funcIsValidString("    "),False,"Should have received false for a blank string input.")

    def testIsValidStringForEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidString(""),False,"Should have received false for a empty string input.")

    def testIsValidStringForWrongCase(self):
        self.assertEqual(ValidateData.funcIsValidString(1),False,"Should have received false for a numeric input.")

    def testIsValidStringTypeForGoodCase(self):
        self.assertEqual(ValidateData.funcIsValidStringType("Valid"),True,"Should have received true for a valid string input.")

    def testIsValidStringTypeForNoneCase(self):
        self.assertEqual(ValidateData.funcIsValidStringType(None),False,"Should have received false for a none string input.")

    def testIsValidStringTypeForBlankCase(self):
        self.assertEqual(ValidateData.funcIsValidStringType("    "),True,"Should have received true for a blank string input.")

    def testIsValidStringTypeForEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidStringType(""),True,"Should have received true for a empty string input.")

    def testIsValidStringTypeForWrongCase(self):
        self.assertEqual(ValidateData.funcIsValidStringType(1),False,"Should have received false for a numeric input.")

    def testIsValidCharForGoodCase(self):
        self.assertEqual(ValidateData.funcIsValidChar('g'),True,"Should have received true for a valid char input.")

    def testIsValidCharForNoneCase(self):
        self.assertEqual(ValidateData.funcIsValidChar(None),False,"Should have received false for a none char input.")

    def testIsValidCharForBlankCase(self):
        self.assertEqual(ValidateData.funcIsValidChar("    "),False,"Should have received false for a blank char input.")

    def testIsValidCharForEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidChar(""),False,"Should have received false for a empty char input.")

    def testIsValidCharForWrongTypeCase(self):
        self.assertEqual(ValidateData.funcIsValidChar(46764735),False,"Should have received false for a wrong type input.")

    def testIsValidPositiveIntegerForGoodCase(self):
        self.assertEqual(ValidateData.funcIsValidPositiveInteger(89),True,"Should have received true for a valid positive integer input.")

    def testIsValidPositiveIntegerForNegativeCase(self):
        self.assertEqual(ValidateData.funcIsValidPositiveInteger(-45),False,"Should have received false for a invalid negative integer input.")

    def testIsValidPositiveIntegerForLargeNegativeCase(self):
        self.assertEqual(ValidateData.funcIsValidPositiveInteger(-999999999999999999999999999999999),False,"Should have received false for a invalid negative integer input.")

    def testIsValidPositiveIntegerForStringCase(self):
        self.assertEqual(ValidateData.funcIsValidPositiveInteger("45"),False,"Should have received false for a invalid string integer input.")

    def testIsValidPositiveIntegerForNoneCase(self):
        self.assertEqual(ValidateData.funcIsValidPositiveInteger(None),False,"Should have received false for a invalid none integer input.")

    def testIsValidListForGoodEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidList(list()),True,"Should have received true for a valid list input.")

    def testIsValidListForGoodCase(self):
        testingList = ["1","2","3","4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),True,"Should have received true for a valid list input.")

    def testIsValidListForGoodCaseListWithinList(self):
        testingList = ["1","2","3",[1,2,3,4],"4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),True,"Should have received true for a valid list within list input.")

    def testIsValidListForGoodCaseListWithinListinList(self):
        testingList = ["1","2","3",[1,2,["a","r","w"],3,4],"4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),True,"Should have received true for a valid list within list within list input.")

    def testIsValidListForNoneListCase(self):
        self.assertEqual(ValidateData.funcIsValidList(None),False,"Should have received false for a invalid none input.")

    def testIsValidListForListWithNoneCase(self):
        testingList = ["1","2","3",None,"4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidListForWrongTypeCase(self):
        self.assertEqual(ValidateData.funcIsValidList(8),False,"Should have received false for a invalid type input.")

    def testIsValidListForListWithListOfNoneCase(self):
        testingList = ["1","2","3",[None,None,None],"4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidListForListWithListOfNoneCase2(self):
        testingList = ["1","2","3",["sge",[None,None,None]],"4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidStringListForEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidStringList(list()),True,"Should have received true for a valid list input.")

    def testIsValidStringListForGoodCase(self):
        testingList = ["1","2","3","4","5","6","7"]
        self.assertEqual(ValidateData.funcIsValidStringList(testingList),True,"Should have received true for a valid list input.")

    def testIsValidStringListForBadCaseListWithinList(self):
        testingList = ["1","2","3",["1","2","3","4"],"4","5","6","7"]
        self.assertEqual(ValidateData.funcIsValidStringList(testingList),False,"Should have received true for a valid list within list input.")

    def testIsValidStringListForNoneListCase(self):
        self.assertEqual(ValidateData.funcIsValidStringList(None),False,"Should have received false for a invalid none input.")

    def testIsValidStringListForListWithNoneCase(self):
        testingList = ["1","2","3",None,"4","5","6","7"]
        self.assertEqual(ValidateData.funcIsValidStringList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidStringListForWrongTypeCase(self):
        self.assertEqual(ValidateData.funcIsValidStringList(8),False,"Should have received false for a invalid type input.")

    def testIsValidStringListForListWithWrongTypeCase(self):
        testingList = ["1","2","3","4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidStringList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidNumericListForEmptyCase(self):
        self.assertEqual(ValidateData.funcIsValidNumericList(list()),True,"Should have received true for a valid list input.")

    def testIsValidNumericListForGoodCase(self):
        testingList = [1,2,3,4,5,6,7]
        self.assertEqual(ValidateData.funcIsValidNumericList(testingList),True,"Should have received true for a valid list input.")

    def testIsValidNumericListForBadCaseListWithinList(self):
        testingList = [1,2,3,[1,2,3,4],4,5,6,7]
        self.assertEqual(ValidateData.funcIsValidNumericList(testingList),False,"Should have received true for a valid list within list input.")

    def testIsValidNumericListForNoneListCase(self):
        self.assertEqual(ValidateData.funcIsValidNumericList(None),False,"Should have received false for a invalid none input.")

    def testIsValidNumericListForListWithNoneCase(self):
        testingList = [1,2,3,None,4,5,6,7]
        self.assertEqual(ValidateData.funcIsValidNumericList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidNumericListForWrongTypeCase(self):
        self.assertEqual(ValidateData.funcIsValidNumericList(8),False,"Should have received false for a invalid type input.")

    def testIsValidNumericListForListWithWrongTypeCase(self):
        testingList = ["1","2","3","4",5,6,7]
        self.assertEqual(ValidateData.funcIsValidNumericList(testingList),False,"Should have received false for a invalid list with none input.")

    def testIsValidNPArrayForGoodCaseEmpty(self):
        testingArray = np.array([1,2,3,4,5,6])
        self.assertEquals(ValidateData.funcIsValidNPArray(testingArray),True,"Should have received True for a np array")

    def testIsValidNPArrayForBadCaseList(self):
        testingArray = [1,2,3,4,5,6]
        self.assertEquals(ValidateData.funcIsValidNPArray(testingArray),False,"Should have received True for a list")

    def testIsValidNPArrayForBadCaseNone(self):
        testingArray = None
        self.assertEquals(ValidateData.funcIsValidNPArray(testingArray),False,"Should have received True for a None")

    def testIsValidNPArrayForGoodCaseWrongType(self):
        testingArray = True
        self.assertEquals(ValidateData.funcIsValidNPArray(testingArray),False,"Should have received True for a boolean")

    def testIsValidDictionaryForGoodCaseEmpty(self):
        testingDictionary = dict()
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),True,"Should have received True for an empty dictionary")

    def testIsValidDictionaryForGoodCase(self):
        testingDictionary = dict(x=56, y=23, z=45)
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),True,"Should have received True for a valid dictionary")

    def testIsValidDictionaryForGoodCase2(self):
        testingDictionary = dict([[1,56], [2,23], [3,45]])
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),True,"Should have received True for a valid dictionary")

    def testIsValidDictionaryForGoodCase3(self):
        testingDictionary = dict([[1,[2,56]], [2,23], [3,[34,87,45]]])
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),True,"Should have received True for a valid dictionary")

    def testIsValidDictionaryForNoneKey(self):
        testingDictionary = dict([[None,[2,56]], [2,23], [3,[34,87,45]]])
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),False,"Should have received false for a none key")

    def testIsValidDictionaryForNoneItem(self):
        testingDictionary = dict([[1,[2,56]], [2,None], [3,[34,87,45]]])
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),False,"Should have received false for a none item")

    def testIsValidDictionaryForNoneItemInList(self):
        testingDictionary = dict([[1,[2,56]], [2,23], [3,[34,None,45]]])
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),False,"Should have received false for a none item")

    def testIsValidDictionaryForNoneDictionary(self):
        testingDictionary = None
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),False,"Should have received false for a none dictionary")

    def testIsValidDictionaryForWrongTypeDictionary(self):
        testingDictionary = list()
        self.assertEquals(ValidateData.funcIsValidDictionary(testingDictionary),False,"Should have received false for a wrong type")

    def testIsValidNucleotideBaseForAGoodCase(self):
        testingSequence = "A"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for an A")

    def testIsValidNucleotideBaseForaGoodCase(self):
        testingSequence = "a"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for an a")

    def testIsValidNucleotideBaseForUGoodCase(self):
        testingSequence = "U"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for an U")

    def testIsValidNucleotideBaseForuGoodCase(self):
        testingSequence = "u"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for an u")

    def testIsValidNucleotideBaseForTGoodCase(self):
        testingSequence = "T"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a T")

    def testIsValidNucleotideBaseFortGoodCase(self):
        testingSequence = "t"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a t")

    def testIsValidNucleotideBaseForCGoodCase(self):
        testingSequence = "C"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a C")

    def testIsValidNucleotideBaseForcGoodCase(self):
        testingSequence = "c"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a c")

    def testIsValidNucleotideBaseForGGoodCase(self):
        testingSequence = "G"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a G")

    def testIsValidNucleotideBaseForgGoodCase(self):
        testingSequence = "g"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),True,"Should have received True for a g")

    def testIsValidNucleotideBaseForNoneCase(self):
        testingSequence = None
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),False,"Should have received True for a None")

    def testIsValidNucleotideBaseForEmptyCase(self):
        testingSequence = ""
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),False,"Should have received True for an empty")

    def testIsValidNucleotideBaseForBlankCase(self):
        testingSequence = "     "
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),False,"Should have received True for a blank")

    def testIsValidNucleotideBaseForBadStringCase(self):
        testingSequence = "$"
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),False,"Should have received True for a $")

    def testIsValidNucleotideBaseForNumberCase(self):
        testingSequence = 1
        self.assertEquals(ValidateData.funcIsValidNucleotideBase(testingSequence),False,"Should have received True for a 1")

    def testIsValidDNASequenceForGoodCase(self):
        testingSequence = "AGCTGCATCGATCGATCGATTAGTACAGTTACCAGAGATCAGTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good lowercases sequence")

    def testIsValidDNASequenceForGoodLowerCase(self):
        testingSequence = "actgaatcgacgacagcaggtatcatgcgacgtacgatcgtacgtat"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good sequence")

    def testIsValidDNASequenceForGoodCaseA(self):
        testingSequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good A sequence")

    def testIsValidDNASequenceForGoodCaseT(self):
        testingSequence = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good T sequence")

    def testIsValidDNASequenceForGoodCaseC(self):
        testingSequence = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good C sequence")

    def testIsValidDNASequenceForGoodCaseG(self):
        testingSequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good G sequence")

    def testIsValidDNASequenceForGoodCaseCG(self):
        testingSequence = "CGCGCGCG"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good CG sequence")

    def testIsValidDNASequenceForGoodCaseAT(self):
        testingSequence = "ATATATATATAT"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),True,"Should have received True for a good AT sequence")

    def testIsValidDNASequenceForNumberCase(self):
        testingSequence = "AGCTGCATCGATCGA1TCGATTAGTACAGTTACCAGAGATCAGTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a bad sequence with a number")

    def testIsValidDNASequenceForYCase(self):
        testingSequence = "AGCTGCATCGATCGATCGATTYAGTACAGTTACCAGAGATCAGTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a bad sequence with a Y")

    def testIsValidDNASequenceForSpaceCase(self):
        testingSequence = "AGCTGCATCGATCGATCGATTAGTACAGTTACCAGAGATCA   GTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a bad sequence with a space")

    def testIsValidDNASequenceForLeadingSpaceCase(self):
        testingSequence = "    AGCTGCATCGATCGATCGATTAGTACAGTTACCAGAGATCAGTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a bad sequence with a leading space")

    def testIsValidDNASequenceForEndingSpaceCase(self):
        testingSequence = "AGCTGCATCGATCGATCGATTAGTACAGTTACCAGAGATCAGTCTCGA     "
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a bad sequence with an ending space")

    def testIsValidDNASequenceForBadLetterCase(self):
        testingSequence = "AGCTGCATCGATCGATCGAWTTAGTACAGTTACCAGAGATCAGTCTCGA"
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for sequence with an inccorect letter with an ending space")

    def testIsValidDNASequenceForNoneCase(self):
        testingSequence = None
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a none sequence.")

    def testIsValidDNASequenceForWrongCase(self):
        testingSequence = 1
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a wrong case sequence.")

    def testIsValidDNASequenceForBlankCase(self):
        testingSequence = ""
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a blank case sequence.")

    def testIsValidDNASequenceForEmptyCase(self):
        testingSequence = "    "
        self.assertEquals(ValidateData.funcIsValidDNASequence(testingSequence),False,"Should have received False for a empty case sequence.")

    def testIsValidPositiveNumberCharForGoodCase(self):
        testingChar = '5.0'
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),True,"Should have received True for a valid char representing a positive number.")

    def testIsValidPositiveNumberCharForGoodCase0(self):
        testingChar = '0.0'
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),True,"Should have received True for a valid char representing a 0.")

    def testIsValidPositiveNumberCharForNegativeNumberCharCase(self):
        testingChar = '-5'
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received false for a negative number char.")

    def testIsValidPositiveNumberCharForGoodStringCase(self):
        testingChar = "5"
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),True,"Should have received True for a valid string representing a positive number.")

    def testIsValidPositiveNumberCharForGoodCase0String(self):
        testingChar = "0"
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),True,"Should have received True for a valid string representing a 0.")

    def testIsValidPositiveNumberCharForNegativeNumberStringCase(self):
        testingChar = "-5"
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received False for a negative number string.")

    def testIsValidPositiveNumberCharForNoneCase(self):
        testingChar = None
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received False for a none char.")

    def testIsValidPositiveNumberCharForPosNumberCase(self):
        testingChar = 5
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received False for a positive number.")

    def testIsValidPositiveNumberCharForZeroNumberCase(self):
        testingChar = 0
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received False for a zero number.")

    def testIsValidPositiveNumberCharForNegativeNumberCase(self):
        testingChar = -5
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received false for a negative number.")

    def testIsValidPositiveNumberCharForEmptyCase(self):
        testingChar = ""
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received false for an empty case.")

    def testIsValidPositiveNumberCharForBlankCase(self):
        testingChar = "    "
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received false for a blank case.")

    def testIsValidPositiveNumberCharForWrongTypeCase(self):
        testingChar = list()
        self.assertEquals(ValidateData.funcIsValidPositiveNumberChar(testingChar),False,"Should have received false for a wrong type case.")

    def testIsValidFlagCharFor0Case(self):
        testingChar = "0"
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),True,"Should have received true for a string 0 case.")

    def testIsValidFlagCharFor0CharCase(self):
        testingChar = '0'
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),True,"Should have received true for a char 0 case.")

    def testIsValidFlagCharFor1Case(self):
        testingChar = "1"
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),True,"Should have received true for a string 1 case.")

    def testIsValidFlagCharFor1CharCase(self):
        testingChar = '1'
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),True,"Should have received true for a char 1 case.")

    def testIsValidFlagCharForNoneCase(self):
        testingChar = None
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),False,"Should have received false for a none case.")

    def testIsValidFlagCharForEmptyCase(self):
        testingChar = ""
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),False,"Should have received false for an empty case.")

    def testIsValidFlagCharForBlankCase(self):
        testingChar = "    "
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),False,"Should have received false for a blank case.")

    def testIsValidFlagCharForWrongCase(self):
        testingChar = True
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),False,"Should have received boolean for a none case.")

    def testIsValidFlagCharForBadStringCase(self):
        testingChar = "True"
        self.assertEquals(ValidateData.funcIsValidFlagChar(testingChar),False,"Should have received false for a bad string case.")

    def testIsValidIntegerForGoodCase(self):
        testingInt = 568
        self.assertEquals(ValidateData.funcIsValidInteger(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidIntegerForGoodNegCase(self):
        testingInt = -9568
        self.assertEquals(ValidateData.funcIsValidInteger(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidIntegerForGoodZeroCase(self):
        testingInt = 0
        self.assertEquals(ValidateData.funcIsValidInteger(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidIntegerForWrongTypeCase(self):
        testingInt = list()
        self.assertEquals(ValidateData.funcIsValidInteger(testingInt),False,"Should have received False for a "+str(testingInt)+".")

    def testIsValidIntegerForNoneCase(self):
        testingInt = None
        self.assertEquals(ValidateData.funcIsValidInteger(testingInt),False,"Should have received False for a "+str(testingInt)+".")

    def testIsValidBoundedNumberCharForGoodCase1(self):
        testingChar = "568"
        boundingIntOne = 0
        boundingIntTwo = 6000
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForGoodCase2(self):
        testingChar = "568"
        boundingIntOne = 6000
        boundingIntTwo = 0
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForGoodCase3(self):
        testingChar = "-578"
        boundingIntOne = 45
        boundingIntTwo = -6754
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForGoodCase4(self):
        testingChar = "-578"
        boundingIntOne = -6754
        boundingIntTwo = 45
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForGoodCase5(self):
        testingChar = "0"
        boundingIntOne = 0
        boundingIntTwo = 0
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForGoodCase6(self):
        testingChar = "-34"
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),True,"Should have received True for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForWrongTypeCase1(self):
        testingChar = set()
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+str(testingChar)+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForWrongTypeCase2(self):
        testingChar = "-34"
        boundingIntOne = list()
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForWrongTypeCase3(self):
        testingChar = "-34"
        boundingIntOne = -34
        boundingIntTwo = dict()
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForNoneCase1(self):
        testingChar = None
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+str(testingChar)+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForNoneCase2(self):
        testingChar = "-34"
        boundingIntOne = None
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForNoneCase3(self):
        testingChar = "-34"
        boundingIntOne = -34
        boundingIntTwo = None
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForEmptyCase(self):
        testingChar = ""
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForBlankCase(self):
        testingChar = "   "
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidBoundedNumberCharForStringCase(self):
        testingChar = "Negative thirty four"
        boundingIntOne = -34
        boundingIntTwo = -34
        self.assertEquals(ValidateData.funcIsValidBoundedIntegerChar(testingChar,boundingIntOne,boundingIntTwo),False,"Should have received False for a "+testingChar+" bounded by: "+str(boundingIntOne)+" and "+str(boundingIntTwo)+".")

    def testIsValidNumericForGoodCase(self):
        testingInt = 568
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForGoodNegCase(self):
        testingInt = -9568
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForGoodZeroCase(self):
        testingInt = 0
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForFloatGoodCase(self):
        testingInt = 568.0
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForFloatGoodNegCase(self):
        testingInt = -9568.1
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForFloatGoodZeroCase(self):
        testingInt = 0.98765435678987654346789876543
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForComplexGoodCase(self):
        testingInt = 1+0j
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForComplexGoodNegCase(self):
        testingInt = -1+0j
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForComplexGoodZeroCase(self):
        testingInt = 0j
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForLongGoodCase(self):
        testingInt = 999999999999999999L
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForLongGoodNegCase(self):
        testingInt = -777777777777777777777L
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForLongGoodZeroCase(self):
        testingInt = 0L
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),True,"Should have received True for a "+str(testingInt)+".")

    def testIsValidNumericForWrongTypeCase(self):
        testingInt = list()
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),False,"Should have received False for a "+str(testingInt)+".")

    def testIsValidNumericForNoneCase(self):
        testingInt = None
        self.assertEquals(ValidateData.funcIsValidNumeric(testingInt),False,"Should have received False for a "+str(testingInt)+".")

    def testIsValidClassForNoneCase(self):
        testingClass = None
        self.assertEquals(ValidateData.funcIsValidClass(parameterValue = testingClass, strCorrectName = "None"),False,"Should have received False for a "+str(testingClass)+".")

    def testIsValidClassForGoodCase(self):
        testingClass = ValidateData()
        self.assertEquals(ValidateData.funcIsValidClass(parameterValue = testingClass, strCorrectName = "ValidateData"),True,"Should have received False for a "+str(testingClass)+".")

    def testIsValidClassForBadCaseWrongName(self):
        testingClass = []
        self.assertEquals(ValidateData.funcIsValidClass(parameterValue = testingClass, strCorrectName = "List"),False,"Should have received False for a "+str(testingClass)+".")

    def testIsValidClassForBadCaseNoneName(self):
        testingClass = []
        self.assertEquals(ValidateData.funcIsValidClass(parameterValue = testingClass, strCorrectName = None),False,"Should have received False for a "+str(testingClass)+".")

    def testIsValidClassForBadCaseWrongTypeName(self):
        testingClass = []
        self.assertEquals(ValidateData.funcIsValidClass(parameterValue = testingClass, strCorrectName = 1),False,"Should have received False for a "+str(testingClass)+".")

##
#Create a suite to be called to test
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(ValidateDataTest)
    return suite

#if __name__=='__main__':
#    unittest.main()
