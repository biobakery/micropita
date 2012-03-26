#######################################################
#
#	Title:		Utility_Math
#	Author:		Timothy Tickle
#	Date:		03/26/2012
#	Purpose:	Utility class for genric math functions.
#
#######################################################

#Import libaries

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
