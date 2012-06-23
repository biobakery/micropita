#######################################################
# Author: Timothy Tickle
# Description: Class to Allow Support Vector Machine 
# analysis and to contain associated scripts.
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
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

    #1 Happy Path tested
    #Converts abundance files to input SVM files.
    #@tempInputFile Abundance file to read (should be a standard Qiime output abundance table)
    #@tempOutputSVMFile File to save SVM data to when converted from teh abundance table
    #@tempDelimiter Delimiter of the Abundance table
    #@tempLabels Ordered labels to use to classify the samples in the abundance table
    #@sLastMetadataName The name of the last row in the abundance table representing metadata
    #@tempSkipFirstColumn Boolean Indicates to skip the first column (true) (for instance if it contains taxonomy identifiers)
    #@tempNormalize Boolean to indicate if the abundance data should be normalized (true) before creating the file (normalized by total sample abundance)
    @staticmethod
    def convertAbundanceTableToSVMFile(abndAbundanceTable, tempOutputSVMFile, sMetadataLabel):
        #Validate parameters
        if abndAbundanceTable == None:
            print "Error, invalid Abundance table."
            return False
        if(not ValidateData.isValidString(tempOutputSVMFile)):
            print "Error, file not valid. File:"+str(tempOutputSVMFile)
            return False

        #If output file exists, delete
        if(os.path.exists(tempOutputSVMFile)):
            os.remove(tempOutputSVMFile)

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
        with open(tempOutputSVMFile,'a') as f:
            (f.write("".join(llData)))
        f.close()
        return lsUniqueLabels

    #Tested
    @staticmethod
    def funcScaleFeature(npdData):
        if sum(npdData) == 0 or len(set(npdData))==1:
            return npdData
        dMin = min(npdData)
        return (npdData-dMin)/float(max(npdData-dMin))

    #Tested
    @staticmethod
    def funcWeightLabels(lLabels):
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
        Generates the indexes for a 10 fold crossvalidation given a sample count.
        If there are less than 10 samples, it uses the sample count as the K-fold crossvalidation
        as a leave one out method.

	:param	iTotalSampleCount:	Total Sample Count
	:type	int:	Sample Count
	:param	fRandomise:	Random sample indices
	:type	boolean:	True indicates randomise (Default False)
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
