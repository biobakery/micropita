#######################################################
#
#	Title:		ValidateData
#	Author:		Timothy Tickle 
#	Date:		August 12, 2009
#	Purpose:	Validate Data containing methods for testing variables
#
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

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

    @staticmethod
    def funcIsValidBoolean(parameterValue):
        """
        Validates a boolean parameter as a valid boolean.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid boolean.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is BooleanType:
            return False
        return True

    @staticmethod
    def funcIsTrue(parameterValue):
        """
        Validates a boolean parameter as true.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is True.
        :type	Boolean
        """

        if(ValidateData.funcIsValidBoolean(parameterValue)):
            if(parameterValue == True):
                return True
        return False

    @staticmethod
    def funcIsFalse(parameterValue):
        """
        Validates a boolean parameter as false.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is False.
        :type	Boolean
        """

        if(ValidateData.isValidBoolean(parameterValue)):
            if(parameterValue == False):
                return True
        return False

    @staticmethod
    def funcIsValidInteger(parameterValue):
        """
        Validates a boolean parameter as an integer.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is an integer.
        :type	Boolean
        """

        #Check to make sure it is not null
        if (parameterValue == None):
            return False

        #Check to make sure it is an integer
        if not type(parameterValue) is IntType:
            return False

        return True

    @staticmethod
    def funcIsValidPositiveInteger(parameterValue, tempZero = False):
        """
        Validates a boolean parameter as false.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :param	tempZero:	Allows one to set what the value for zero should return.
        :type	Boolean	The return value for zero.
        :return	Boolean:	True indicates the parameter is a positive integer.
        :type	Boolean
        """

        #Check to make sure it is not null
        if not ValidateData.funcIsValidInteger(parameterValue):
            return False

        #Check to see it is positive
        if (parameterValue < 0):
            return False

        #Check for zero value
        if(parameterValue == 0):
            return tempZero
        return True

    @staticmethod
    def funcIsValidNumeric(parameterValue):
        """
        Validates a boolean parameter as an integer.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a numeric.
        :type	Boolean
        """

        #Check to make sure it is not null
        if (parameterValue == None):
            return False
        #Check to make sure it is an integer
        if((type(parameterValue) == IntType)or(type(parameterValue) == LongType)or(type(parameterValue) == FloatType)or(type(parameterValue) == ComplexType)or(str(type(parameterValue)) == "<type 'numpy.float64'>")):
            if(not type(parameterValue) == BooleanType):
                return True
        return False

    @staticmethod
    def funcIsValidStringType(parameterValue):
        """
        Validates a boolean parameter as a string. This allows the string to be blank or empty.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a string type.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is StringType:
            return False

        return True

    @staticmethod
    def funcIsValidString(parameterValue):
        """
        Validates a boolean parameter as a string. Does NOT allow string to be blank or empty.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a string.
        :type	Boolean
        """

        #Type check
        if not ValidateData.funcIsValidStringType(parameterValue):
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
    def funcIsValidFormatString(parameterValue):

        lettersValid = False
        if ValidateData.funcIsValidString(parameterValue):
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
    def funcIsValidChar(parameterValue):
        return ValidateData.funcIsValidString(parameterValue)

    ##
    #Validates char parameters that must be a number of 0 or more
    #Created 4/11/10
    #Tested 4/12/10
    #@param ParameterValue value to be evaluated as a char representing a decimal equal to or greater than 0
    #@return bool Indicator of parameter validity
    @staticmethod
    def funcIsValidPositiveNumberChar(parameterValue):
        #Check to make sure is a valid string
        if not ValidateData.funcIsValidString(parameterValue):
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
    def funcIsValidFlagChar(parameterValue):
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
    def funcIsValidBoundedIntegerChar(parameterValue, tempValueOne, tempValueTwo):
        #Check to make sure is a valid string
        if not ValidateData.funcIsValidString(parameterValue):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.funcIsValidInteger(tempValueOne):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.funcIsValidInteger(tempValueTwo):
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

    @staticmethod
    def funcIsValidList(parameterValue):
        """
        Validates a boolean parameter as a list.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a list
        :type	Boolean
        """

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
                if ValidateData.funcIsValidList(parameterValue[i]) == False:
                    return False
        return True

    ##
    #Validates tuple parameters
    #Tested 11/25/2011
    #@param ParameterValue value to be evaluated as a tuple
    #@return bool Indicator of list parameter validity
    @staticmethod
    def funcIsValidTuple(parameterValue):

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
    def funcIsValidNumericList(parameterValue):
        #Check is valid list
        if(not ValidateData.funcIsValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in xrange(0,listSize):
            if(not ValidateData.funcIsValidNumeric(parameterValue[i])):
                return False
        return True

    ##
    #Validates list parameters with only string types
    #Created tests 6-24-2012
    #Tested 6-24-2012
    #@param ParameterValue value to be evaluated as a string list
    #@return bool Indicator of list parameter validity
    @staticmethod
    def funcIsValidStringList(parameterValue):
        #Check is valid list
        if(not ValidateData.funcIsValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in xrange(0,listSize):
            if(not ValidateData.funcIsValidString(parameterValue[i])):
                return False
        return True

    ##
    #Makes sure the object is not none and has a type of ndarray
    #Created 2/25/2011
    @staticmethod
    def funcIsValidNPArray(parameterValue):

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
    def funcIsValidDictionary(parameterValue):

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
                if validateData.funcIsValidList(keyList[i]) == False:
                    return False

        #Check key elements
        itemList = parameterValue.values()
        itemListSize = len(itemList)

        for i in range(0,itemListSize):
            if itemList[i] == None:
                return False
            if type(itemList[i]) is ListType:
                if ValidateData.funcIsValidList(itemList[i]) == False:
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
    def funcIsValidDNASequence(parameterValue):
        if ValidateData.funcIsValidString(parameterValue):
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
    def funcIsValidNucleotideBase(parameterValue):
        if (ValidateData.funcIsValidDNASequence(parameterValue) or (parameterValue == 'u') or (parameterValue == "U")):
            if (len(parameterValue) == 1):
                return True
        return False

    ####
    #File Structure
    ####
    

    ####
    #Classes, Functions, and Instances
    ####

    @staticmethod
    def funcIsValidFileName(parameterValue):
        """
        Validates a boolean parameter as a valid file name.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid file path.
        :type	Boolean
        """

        if(ValidateData.funcIsValidString(parameterValue)):
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
    def funcIsValidClass(tempClassInstance, tempCorrectName):
        if(tempClassInstance==None):
            return False
        if not ValidateData.funcIsValidString(tempCorrectName):
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
