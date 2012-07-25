"""
Author: Timothy Tickle
Description: Utility class for generic math functions.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libaries
import itertools
import numpy as np
import operator
import random
from ValidateData import ValidateData

##
#Utility function for data generation
class Utility_Math():
    """
    Class to perform misc math methods.
    """

    ##
    #Happy path test 2
    @staticmethod
    def funcConvertToBHQValue(ldPValues, iNumberOfTests=None):
        """
        Convert a list of p-value to a list of q-values.

        :param	ldPValues:	List of doubles (p-values) to convert.
        :type	List
        :param	iNumberOfTests:	Number of (multiple) tests if different than the ldValue length. If not set the length of ldPValues is used.
        :type	Integer
        :return	List:	List of Q-values made with a BH modification.
        """

        #If the number of tests is not specified, use the number of pvalues
        if(iNumberOfTests == None):
            iNumberOfTests = len(ldPValues)
        #Used to hold the pvalues as they are being manipulated
        lsConvertToQValues = list()
        #Is used to set the ordr of the pvalues as they are placed in the lsConvertToQValues
        dOrder = 1
        for dValue in ldPValues:
            lsConvertToQValues.append([dValue,dOrder,None])
            dOrder = dOrder + 1

        #Sort by pvalue
        lsConvertToQValues.sort(key=lambda x: x[0])

        #Used to keep track of the current test number
        iTest = 1
        for dConvValue in lsConvertToQValues:
            dConvValue[2] = dConvValue[0] * iNumberOfTests / iTest
            iTest = iTest + 1

        #Sort by original order
        lsConvertToQValues.sort(key=lambda x: x[1])

        #return just 1 dimension (the qvalue)
        return [ldValues[2] for ldValues in lsConvertToQValues]

    #Happy path tested 5
    @staticmethod
    def funcSampleWithReplacement(aData, iSelect):
        """
        Sample from a vector of data (aData) with replacement iSelect many objects.

        :param	aData:	Data to sample from with replacement.
        :type	List
        :param	iSelect:	Amount of data to select from the original data population.
        :type	Integer.
        :return	List:	List of sampled data.
                        Returns an empty list on error.
        """

        if iSelect and aData:
            iDataSize = len(aData)
            funcRandom, funcInt = random.random, int
            lsSampling =  operator.itemgetter(*[funcInt(funcRandom() * iDataSize) for selected in itertools.repeat(None, iSelect)])(aData)
            if isinstance(lsSampling, basestring):
                lsSampling = [lsSampling]
            return lsSampling
        return []

    #Happy Path Tested 2
    @staticmethod
    def funcSumRowsOfColumns(npaAbundance, lsSampleNames):
        """
        Takes the column names of a npArray and sums the rows into one column.

        :param	npaAbundance:	Array of data to sum.
        :type	Numpy Array
        :param	lsSampleNames:	List of sample names.
        :type	List	List of strings.
        :return	List	List of data summed at each row.
        """

        #Compress by data name
        npPooledSample = npaAbundance[lsSampleNames[0]]
        for strSampleName in lsSampleNames[1:]:
            #When combining, combine counts by summing
            npPooledSample = npPooledSample + npaAbundance[strSampleName]
        return list(npPooledSample)

    #Testing Status: Light happy path testing 2
    @staticmethod
    def funcTransposeDataMatrix(npaMatrix, fRemoveAdornments=False):
        """
        Transposes a numpy array.

        :param	npaMatrix:	Data matrix to transpose.
        :type	Numpy Array	
        :param	fRemoveAdornments:	Remove the first column before transposing.
        :type	Boolean	True indicates removing the column.
        :return	Boolean or Numpy Array:	Transposed array or a boolean indicating error.
                                   Boolean	False is returned on error.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidNPArray(npaMatrix)):
            print "".join(["Utility_Math:transposeDataMatrix::Error, transposeDataMatrix was an invalid structured array. Value =",str(npaMatrix)])
            return False
        if(not ValidateData.funcIsValidBoolean(fRemoveAdornments)):
            print "".join(["Utility_Math:transposeDataMatrix::Error, fRemoveAdornments was an invalid boolean. Value =",str(fRemoveAdornments)])
            return False

        #Change to samples x taxa as is needed for the compute method below
        #Also remove the first row which is taxa identification
        startColumnElement = 0
        if(fRemoveAdornments == True):
            startColumnElement = 1
        conversionMatrix = list()
        for row in npaMatrix:
            conversionMatrix.append(list(row)[startColumnElement:])
        npaMatrix = None
        conversionMatrix = np.array(conversionMatrix)
        conversionMatrix = conversionMatrix.transpose()
        return conversionMatrix

