#######################################################
# Author: Timothy Tickle
# Description: Class to Allow Support Vector Machine 
# analysis and to contain associated scripts.
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Libraries
from AbundanceTable import AbundanceTable
from Constants import Constants
from CommandLine import CommandLine
import math
import operator
import os
from random import shuffle
from ValidateData import ValidateData

class SVM:
    """
    Class which holds generic methods for SVM use.
    """

    #1 Happy Path tested
    @staticmethod
    def funcConvertAbundanceTableToSVMFile(abndAbundanceTable, strOutputSVMFile, sMetadataLabel):
        """
        Converts abundance files to input SVM files.

        :param strOutputSVMFile: File to save SVM data to when converted from the abundance table.
        :type	String
        :param	sMetadataLabel: The name of the last row in the abundance table representing metadata.
        :type	String
        :return	lsUniqueLabels:	List of unique labels.
        :type	List	List of strings
        """

        #Validate parameters
        if abndAbundanceTable == None:
            print "Error, invalid Abundance table."
            return False
        if(not ValidateData.funcIsValidString(strOutputSVMFile)):
            print "Error, file not valid. File:"+str(strOutputSVMFile)
            return False

        #If output file exists, delete
        if(os.path.exists(strOutputSVMFile)):
            os.remove(strOutputSVMFile)

        #Create data matrix
        dataMatrix = zip(*abndAbundanceTable.funcGetAbundanceCopy())

        #Add labels
        llData = []
        lsLabels = abndAbundanceTable.funcGetMetadata(sMetadataLabel)
        lsUniqueLabels = list(set(lsLabels))
        dictLabels = dict([[str(lenuLabels[1]),str(lenuLabels[0])] for lenuLabels in enumerate(lsUniqueLabels)])
        lsLabels = [dictLabels[sLabel] for sLabel in lsLabels]

        iRowIndex = 0
        for dataRow in dataMatrix[1:]:
            llData.append(" ".join([lsLabels[iRowIndex]]+[Constants.COLON.join([str(enuSamples[0]+1),str(enuSamples[1])])
                            for enuSamples in enumerate(dataRow)])+Constants.ENDLINE)
            iRowIndex = iRowIndex + 1

        #Output file
        with open(strOutputSVMFile,'a') as f:
            (f.write("".join(llData)))

        return lsUniqueLabels

    #Tested
    @staticmethod
    def funcScaleFeature(npdData):
        """
        Scale a feature between 0 and 1. Using 01 and not 01,1 because it keeps te sparsity of the data and may save time.

        :param	npdData:	Feature data to scale.
        :type	Numpy Array	Scaled feature data.
        """

        if sum(npdData) == 0 or len(set(npdData))==1:
            return npdData
        dMin = min(npdData)
        return (npdData-dMin)/float(max(npdData-dMin))

    #Tested
    @staticmethod
    def funcWeightLabels(lLabels):
        """
        Returns weights for labels based on how balanced the labels are. Weights try to balance unbalanced results.

        :params	lLabels:	List of labels to use for measure how balanced the comparison is.
        :type	List
        :return	List:		[dictWeights ({"label":weight}),lUniqueLabels (unique occurences of original labels)]
        :type	List
        """

        #Convert to dict
        lUniqueLabels = list(set(lLabels))
        dictLabels = dict(zip(lUniqueLabels, range(len(lUniqueLabels))))

        #Build a dict of weights per label {label:weight, label:weight}
        #Get the occurence of each label
        dictWeights = dict()
        for sLabelKey in dictLabels:
            sCurLabel = dictLabels[sLabelKey]
            dictWeights[sCurLabel] = lLabels.count(sLabelKey)

        #Divide the highest occurence each occurence
        iMaxOccurence = max(dictWeights.values())
        for sWeightKey in dictWeights:
            dictWeights[sWeightKey]=iMaxOccurence/float(dictWeights[sWeightKey])

        return [dictWeights,lUniqueLabels]

    #Tested 3/4 cases could add in test 12 with randomize True
    def func10FoldCrossvalidation(self, iTotalSampleCount, fRandomise = False):
        """
        Generator.
        Generates the indexes for a 10 fold crossvalidation given a sample count.
        If there are less than 10 samples, it uses the sample count as the K-fold crossvalidation
        as a leave one out method.

	:param	iTotalSampleCount:	Total Sample Count
	:type	Integer	Sample Count
	:param	fRandomise:	Random sample indices
	:type	Boolean	True indicates randomise (Default False)
        """

        #Make indices and shuffle if needed
        liindices = range(iTotalSampleCount)
        if fRandomise:
            shuffle(liindices)

        #For 10 times
        iKFold = 10
        if iTotalSampleCount < iKFold:
            iKFold = iTotalSampleCount
        for iiteration in xrange(iKFold):
            lfTraining = [iindex % iKFold != iiteration for iindex in liindices]
            lfValidation = [not iindex for iindex in lfTraining]
            yield lfTraining, lfValidation
