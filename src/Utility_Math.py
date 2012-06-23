#######################################################
#
#	Title:		Utility_Math
#	Author:		Timothy Tickle
#	Date:		03/26/2012
#	Purpose:	Utility class for genric math functions.
#
#######################################################

#Import libaries
import itertools
import numpy as np
import operator
import random
from ValidateData import ValidateData

##
#Utility function for data generation
class Utility_Math():

  ##
  #Contructor
  def __init__(): pass

  ##
  #Happy path test 2
  #Generate matrix for microPITA
  #@param tempOutPutFile
  @staticmethod
  def convertToBHQValue(ldPValues, iNumberOfTests=None):
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
  #Sample from a vector of data (aData) with replacement iSelect many objects
  @staticmethod
  def funcSampleWithReplacement(aData, iSelect):
    if iSelect and aData:
      iDataSize = len(aData)
      funcRandom, funcInt = random.random, int
      lsSampling =  operator.itemgetter(*[funcInt(funcRandom() * iDataSize) for selected in itertools.repeat(None, iSelect)])(aData)
      if isinstance(lsSampling, basestring):
        lsSampling = [lsSampling]
      return lsSampling
    return []

  #Happy Path Tested 2
  #Takes the column indices of a npArray and sums the rows into one column
  #Returns a list which is the row sums of the column
  @staticmethod
  def funcSumRowsOfColumns(npaAbundance, lsSampleNames):
    #Compress by data name
    npPooledSample = npaAbundance[lsSampleNames[0]]
    for strSampleName in lsSampleNames[1:]:
      #When combining, combine counts by summing
      npPooledSample = npPooledSample + npaAbundance[strSampleName]
    return list(npPooledSample)

  #Testing Status: Light happy path testing 2
  #Transposes a matrix.
  #Removes the first column before transposing if tempRemoveAdornments = True
  #@params tempMatrix Structured array to transpose
  #@params tempRemoveAdornments Remove first column before transposing
  #@return Returns Transposed structured array (list of lists, list=previous columns) with potentially the first column removed.
  @staticmethod
  def transposeDataMatrix(tempMatrix, tempRemoveAdornments=False):
    #Validate parameters
    if(not ValidateData.isValidNPArray(tempMatrix)):
      print "".join(["Utility_Math:transposeDataMatrix::Error, transposeDataMatrix was an invalid structured array. Value =",str(tempMatrix)])
      return False
    if(not ValidateData.isValidBoolean(tempRemoveAdornments)):
      print "".join(["Utility_Math:transposeDataMatrix::Error, tempRemoveAdornments was an invalid boolean. Value =",str(tempRemoveAdornments)])
      return False

    #Change to samples x taxa as is needed for the compute method below
    #Also remove the first row which is taxa identification
    startColumnElement = 0
    if(tempRemoveAdornments == True):
      startColumnElement = 1
    conversionMatrix = list()
    for row in tempMatrix:
      conversionMatrix.append(list(row)[startColumnElement:])
    tempMatrix = None
    conversionMatrix = np.array(conversionMatrix)
    conversionMatrix = conversionMatrix.transpose()
    return conversionMatrix

