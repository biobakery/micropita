#######################################################
#
#	Title:		ValidateData
#	Author:		Timothy Tickle 
#	Date:		August 12, 2009
#	Purpose:	Validate Data containing methods for testing variables
#
#######################################################

#Import local code
from types import *
import decimal
import os
import re
import string

##
#Includes methods that help validating parameters
#Created 8/12/09
class ValidateData:

    ####
    #Validate Primitives
    ####

    ##
    #Validates boolean parameters
    #Created 11/9/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a boolean
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isValidBoolean(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is BooleanType:
            return False
        return True

    ##
    #Validates a boolean parameter as true
    #@param ParameterValue value to be evaluated as a True
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isTrue(parameterValue):
        if(ValidateData.isValidBoolean(parameterValue)):
            if(parameterValue == True):
                return True
        return False

    ##
    #Validates a boolean parameter as false
    #@param ParameterValue value to be evaluated as a False
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isFalse(parameterValue):
        if(ValidateData.isValidBoolean(parameterValue)):
            if(parameterValue == False):
                return True
        return False

    ##
    #Validates interger parameters
    #Created 4/13/10
    #Tested 4/13/10
    #@param ParameterValue value to be evaluated as an integer
    #@return bool Indicator of integer parameter validity
    @staticmethod
    def isValidInteger(parameterValue):

        #Check to make sure it is not null
        if (parameterValue == None):
            return False

        #Check to make sure it is an integer
        if not type(parameterValue) is IntType:
            return False

        return True

    ##
    #Validates integer parameters making sure they are greater than 0
    #Created 11/9/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a positive integer
    #@param tempZero Value to return for Zero
    #@return bool Indicator of integer parameter validity
    @staticmethod
    def isValidPositiveInteger(parameterValue, tempZero = False):

        #Check to make sure it is not null
        if not ValidateData.isValidInteger(parameterValue):
            return False

        #Check to see it is positive
        if (parameterValue < 0):
            return False

        #Check for zero value
        if(parameterValue == 0):
            return tempZero
        return True

    ##
    #Validates numeric parameters
    #Created 6-24-2012
    #Tested 6-24-2012
    #@param ParameterValue value to be evaluated as a numeric
    #@return bool Indicator of integer parameter validity
    @staticmethod
    def isValidNumeric(parameterValue):
        #Check to make sure it is not null
        if (parameterValue == None):
            return False
        #Check to make sure it is an integer
        if((type(parameterValue) == IntType)or(type(parameterValue) == LongType)or(type(parameterValue) == FloatType)or(type(parameterValue) == ComplexType)or(str(type(parameterValue)) == "<type 'numpy.float64'>")):
            if(not type(parameterValue) == BooleanType):
                return True
        return False

    ##
    #Validates string parameters, allows string to be blank or empty
    #Created 10/4/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a string
    #@return bool Indicator of string parameter validity
    @staticmethod
    def isValidStringType(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is StringType:
            return False

        return True

    ##
    #Validates string parameters
    #Created 10/4/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a string
    #@return bool Indicator of string parameter validity
    @staticmethod
    def isValidString(parameterValue):

        #Type check
        if not ValidateData.isValidStringType(parameterValue):
            return False

        #Check to see it is not blank
        if parameterValue.strip() == "":
            return False
        return True

    ##
    #Validates format strings used for binary data I/O
    #Created 5/24/2010
    #Tested 6-24-2012
    #@param ParameterValue value to be evaluated as a format string
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidFormatString(parameterValue):
        lettersValid = False
        if ValidateData.isValidString(parameterValue):
            validChars = "BbcdfHhIiLlPpsx0123456789"
            for letter in parameterValue:
                lettersValid = letter in validChars
                if(not lettersValid):
                    break
        return lettersValid

    ##
    #Validates char parameters
    #Created 11/15/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a character
    #@return bool Indicator of char parameter validity
    @staticmethod
    def isValidChar(parameterValue):
        return ValidateData.isValidString(parameterValue)

    ##
    #Validates char parameters that must be a number of 0 or more
    #Created 4/11/10
    #Tested 4/12/10
    #@param ParameterValue value to be evaluated as a char representing a decimal equal to or greater than 0
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidPositiveNumberChar(parameterValue):
        #Check to make sure is a valid string
        if not ValidateData.isValidString(parameterValue):
            return False

        #Try to convert to decimal
        try:
            decimalConversion = decimal.Decimal(parameterValue)
            if decimalConversion < 0:
                return False
        except:
            return False
        return True

    ##
    #Validates char parameters that must be a 0 or 1
    #Created 4/11/10
    #Tested 4/12/10
    #@param ParameterValue value to be evaluated as a boolean
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidFlagChar(parameterValue):
        if parameterValue == '0' or parameterValue == "0" or parameterValue == '1' or parameterValue == "1":
            return True
        return False

    ##
    #Validates char parameters that must represent a integer between and including the two given integers.
    #The lower integer bounds the number at the lower end and the higher bounds at the higher.
    #The integers can be given in highest or lowest first. The bounding includes the integer.
    #Created 4/13/10
    #Tested 4/13/10
    #@param ParameterValue value to be evaluated as an integer
    #@param tempValueOne integer representing the one value that bounds the paramterValue
    #@param tempValueTwo integer representing the a second value that bounds the paramterValue
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidBoundedIntegerChar(parameterValue, tempValueOne, tempValueTwo):
        #Check to make sure is a valid string
        if not ValidateData.isValidString(parameterValue):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.isValidInteger(tempValueOne):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.isValidInteger(tempValueTwo):
            return False

        #Try to convert to decimal
        try:
            intConversion = int(parameterValue)
            if(tempValueOne < tempValueTwo):
                if ((intConversion >= tempValueOne) and (intConversion <= tempValueTwo)):
                    return True
                return False
            if(tempValueTwo < tempValueOne):
                if ((intConversion >= tempValueTwo) and (intConversion <= tempValueOne)):
                    return True
                return False
            if(tempValueOne == tempValueTwo):
                if (intConversion == tempValueOne):
                    return True
                return False
        except:
            return False
    ####
    #Collections
    ####

    ##
    #Validates list parameters
    #Created 10/4/09
    #Tested 12/12/09
    #@param ParameterValue value to be evaluated as a list
    #@return bool Indicator of list parameter validity
    @staticmethod
    def isValidList(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is ListType:
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in range(0,listSize):
            if parameterValue[i] == None:
                return False
            if type(parameterValue[i]) is ListType:
                if ValidateData.isValidList(parameterValue[i]) == False:
                    return False
        return True

    ##
    #Validates tuple parameters
    #Tested 11/25/2011
    #@param ParameterValue value to be evaluated as a tuple
    #@return bool Indicator of list parameter validity
    @staticmethod
    def isValidTuple(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is TupleType:
            return False

        #Check elements
        tupleSize = len(parameterValue)
        for i in range(0,tupleSize):
            if parameterValue[i] == None:
                return False
        return True

    ##
    #Validates list parameters with only numeric types
    #Tests created 6-24-2012
    #Tested 6-24-2012
    #@param ParameterValue value to be evaluated as a numeric list
    #@return bool Indicator of list parameter validity
    @staticmethod
    def isValidNumericList(parameterValue):
        #Check is valid list
        if(not ValidateData.isValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in xrange(0,listSize):
            if(not ValidateData.isValidNumeric(parameterValue[i])):
                return False
        return True

    ##
    #Validates list parameters with only string types
    #Created tests 6-24-2012
    #Tested 6-24-2012
    #@param ParameterValue value to be evaluated as a string list
    #@return bool Indicator of list parameter validity
    @staticmethod
    def isValidStringList(parameterValue):
        #Check is valid list
        if(not ValidateData.isValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in xrange(0,listSize):
            if(not ValidateData.isValidString(parameterValue[i])):
                return False
        return True

    ##
    #Makes sure the object is not none and has a type of ndarray
    #Created 2/25/2011
    @staticmethod
    def isValidStructuredArray(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a structure array
        if not str(type(parameterValue)) == "<type 'numpy.ndarray'>":
            return False

        return True

    ##
    #Validates list parameters
    #Created 2/8/2010
    #Tested 2/8/2010
    #@param ParameterValue value to be evaluated as a dictionary
    #@return bool Indicator of dictionary parameter validity
    @staticmethod
    def isValidDictionary(parameterValue):

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is DictType:
            return False

        #Check key elements
        keyList = parameterValue.keys()
        keyListSize = len(keyList)
        for i in range(0,keyListSize):
            if keyList[i] == None:
                return False
            if type(keyList[i]) is ListType:
                if validateData.isValidList(keyList[i]) == False:
                    return False

        #Check key elements
        itemList = parameterValue.values()
        itemListSize = len(itemList)

        for i in range(0,itemListSize):
            if itemList[i] == None:
                return False
            if type(itemList[i]) is ListType:
                if ValidateData.isValidList(itemList[i]) == False:
                    return False
        return True

    ####
    #Bioinformatics Primitives
    ####

    ##
    #Validates DNA Sequences parameters
    #Created 2/8/2010
    #Tested 2/8/2010
    #@param ParameterValue value to be evaluated as a DNA sequence
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidDNASequence(parameterValue):
        if ValidateData.isValidString(parameterValue):
            expression = re.compile(r'[^atcgATCG]')
            if not None == expression.search(parameterValue):
                return False
            return True
        return False

    ##
    #Validates Nucleotide bases
    #Created 5/10/2010
    #Tested 5/10/2010
    #@param ParameterValue value to be evaluated as a nucleotide base
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidNucleotideBase(parameterValue):
        if (ValidateData.isValidDNASequence(parameterValue) or (parameterValue == 'u') or (parameterValue == "U")):
            if (len(parameterValue) == 1):
                return True
        return False

    ####
    #File Structure
    ####
    

    ####
    #Classes, Functions, and Instances
    ####
    ##
    #Validates string file or directory names that exist
    #Created 10/04/11
    #@param ParameterValue value to be evaluated as a file name
    #@return bool Indicator of char parameter validity
    @staticmethod
    def isValidFileName(parameterValue):
        if(ValidateData.isValidString(parameterValue)):
            return os.path.exists(parameterValue)
        return False

    ##
    #Validates classes to be of the give type
    #Created 6-24-2012
    #Tested 6-24-2012
    #@param tempClassInstance Instance of the class to be evaluated
    #@param tempCorrectName TestName of the class to evaluate the instance against
    #@return bool Indicator of parameter validity
    @staticmethod
    def isValidClass(tempClassInstance, tempCorrectName):
        if(tempClassInstance==None):
            return False
        if not ValidateData.isValidString(tempCorrectName):
            return False
        classType = type(tempClassInstance).__name__
        if(classType == tempCorrectName):
            return True
        if(classType == 'instance'):
            if(tempClassInstance.__class__.__name__==tempCorrectName):
                return True
            else:
                return False
        return False

    ##
    #Validates a parameter as a function
    #@param ParameterValue value to be evaluated as a function
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isValidFunctionType(parameterValue):
        if parameterValue == None:
            return False
        if (not (type(parameterValue) is FunctionType)):
            return False
        return True

    ##
    #Validates a parameter as an instance
    #@param ParameterValue value to be evaluated as an instance
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isValidInstance(parameterValue):
        if parameterValue == None:
            return False
        if (not (type(parameterValue) is InstanceType)):
            return False
        return True

    ##
    #Validates a parameter as a method
    #@param ParameterValue value to be evaluated as a method
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isValidMethod(parameterValue):
        if parameterValue == None:
            return False
        if (not (type(parameterValue) is MethodType)):
            return False
        return True

    ##
    #Validates a parameter as an instance of a method
    #@param ParameterValue value to be evaluated as an instance of a method
    #@return bool Indicator of boolean parameter validity
    @staticmethod
    def isValidInstanceMethod(parameterValue):
        if parameterValue == None:
            return False
        if (not (str(type(parameterValue)) == "<type 'instancemethod'>")):
            return False
        return True
