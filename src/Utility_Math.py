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
import operator
import random

##
#Utility function for data generation
class Utility_Math():

  ##
  #Contructor
  def __init__(): pass

  ##
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

  ##
  #Bootstrap and perform the given method on a population of data
  #Bootstrap, sample with replacement (each iteration of sampling may not have unique selection).
#  @staticmethod
#  def funcBootstrapRowsFromMatrix(aMatrix, iSelectRowCount, iIterationCount, funcFunction, funcSummaryFunction):
#    iColumnCount = dim(aMatrix)
#    iRowCount = 0
#    if iColumnCount:
#      iRowCount = len(aMatrix[0])
#    return [funcSummaryFunction([funcFunction(funcSampleWithReplacement(row,iSelectRowCount)) for iIter in xrange(iIteractionCount)]) for row in aMatrix]

  #Sample from a vector of data (aData) with replacement iSelect many objects
  @staticmethod
  def funcSampleWithReplacement(aData, iSelect):
    if iSelect and aData:
      iDataSize = len(aData)
      funcRandom, funcInt = random.random, int
      return operator.itemgetter(*[funcInt(funcRandom() * iDataSize) for selected in itertools.repeat(None, iSelect)])(aData)
    return []

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

