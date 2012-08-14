#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Run analysis for the microPITA paper
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

import sys
import argparse
from src.breadcrumbs.AbundanceTable import AbundanceTable
from src.breadcrumbs.ConstantsBreadCrumbs import ConstantsBreadCrumbs
from src.breadcrumbs.Metric import Metric
from src.breadcrumbs.KMedoids import Kmedoids
from src.breadcrumbs.MLPYDistanceAdaptor import MLPYDistanceAdaptor
from src.breadcrumbs.SVM import SVM
from src.breadcrumbs.UtilityMath import UtilityMath

from src.ConstantsMicropita import ConstantsMicropita
import csv
#import itertools
import logging
import math
import mlpy
import numpy as np
import operator
import os
import random
#import re
import scipy.cluster.hierarchy as hcluster
#import sys
from types import *

class MicroPITA:
    """
    Selects samples from a first tier of a multi-tiered study to be used in a second tier.
    Different methods can be used for selection.
    The expected input is an abundance table (and potentially a text file of targeted features,
    if using the targeted features option). Output is a list of samples exhibiting the
    characteristics of interest.
    """

####Group 1## Diversity
    #Testing: Happy path Testing (8)
    def funcGetTopRankedSamples(self, lldMatrix = None, lsSampleNames = None, iTopAmount = None):
	"""
	Given a list of lists of measurements, for each list the indices of the highest values are returned. If lsSamplesNames is given
        it is treated as a list of string names that is in the order of the measurements in each list. Indices are returned or the sample
        names associated with the indices.
	
	:param	lldMatrix:	List of lists [[value,value,value,value],[value,value,value,value]].
	:type	List of lists	List of measurements. Each list is a different measurement. Each measurement in possionally related to a sample.
	:param	lsSampleNames:	List of sample names positionally related (the same) to each list (Optional).
	:type	List of strings	List of strings.
	:param	iTopAmount:	The amount of top measured samples (assumes the higher measurements are better).
	:type	integer	Integer amount of sample names/ indices to return.
	:return	List:	List of samples to be selected.
	"""

        topRankList = []
        #If none are given return and empty list.
        if(len(lldMatrix)<1):
            return topRankList
        for rowMetrics in lldMatrix:
            #Create 2 d array to hold value and index and sort
            rowsMetricsLength = len(rowMetrics)
            indexX = [rowMetrics,range(rowsMetricsLength)]
            indexX[1].sort(key = indexX[0].__getitem__,reverse = True)
            if(lsSampleNames == None):
                topRankList.append(indexX[1][:iTopAmount])
            else:
                sampleIndexes = indexX[1][:iTopAmount]
                addSamplesToRank = []
                for index in sampleIndexes:
                    addSamplesToRank.append(lsSampleNames[index])
                topRankList.append(addSamplesToRank)
        return topRankList

####Group 2## Representative Dissimilarity

    #Testing: Happy Path Tested for BrayCurtis and Inverse BrayCurtis
    def funcGetBetaMetric(self, npadAbundancies=None, sMetric=None):
	"""
	Takes a matrix of values and returns a beta metric matrix. The metric returned is indicated by name (sMetric).
	
	:param	npadAbundancies:	Numpy array of sample abundances to measure against.
	:type	Numpy Array	Numpy array where row=samples and columns = features.
	:param	sMetric:	String name of beta metric. Possibilities are listed in microPITA.
	:type	String	String name of beta metric. Possibilities are listed in microPITA.
	:return	Double:	Measurement indicated by metric for given abundance list
	"""

        if sMetric == ConstantsMicropita.c_strBrayCurtisDissimilarity:
            return Metric.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        elif(sMetric == ConstantsMicropita.c_strInvBrayCurtisDissimilarity):
            return Metric.funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        else:
            return False

    #Testing: Happy path tested 1
    def funcGetCentralSamplesByKMedoids(self, npaMatrix=None, sMetric=None, lsSampleNames=None, iNumberSamplesReturned=0):
	"""
	Gets centroid samples by k-medoids clustering of a given matrix.
	
	:param	npaMatrix:	Numpy array where row=features and columns=samples
	:type	Numpy array	Abundance Data.
	:param	sMetric:	String name of beta metric used as the distance metric.
	:type	String	String name of beta metric. Possibilities are listed in microPITA.
	:param	lsSampleNames:	The names of the sample
	:type	List	List of strings
	:param	iNumberSamplesReturned:	Number of samples to return, each will be a centroid of a sample.
	:type	Integer	Number of samples to return
	:return	List:	List of selected samples.
	"""

        #Validate parameters
        if(iNumberSamplesReturned < 0):
            logging.error("MicroPITA.funcGetCentralSamplesByKMedoids. Number of samples to return must be atleast 1.")
            return False

        #Count of how many rows
        sampleCount = npaMatrix.shape[0]
        if(iNumberSamplesReturned > sampleCount):
            logging.error("".join(["MicroPITA.funcGetCentralSamplesByKMedoids. There are not enough samples to return the amount of samples specified. Return sample count = ",str(iNumberSamplesReturned),". Sample number = ",str(sampleCount),"."]))
            return False

        #Samples to return
        returningSamples = list()

        #If the cluster count is equal to the sample count return all samples
        if(sampleCount == iNumberSamplesReturned):
            return list(lsSampleNames)

        #Get distance matrix
        distanceMatrix=self.funcGetBetaMetric(npadAbundancies=npaMatrix, sMetric=sMetric)
        if type(distanceMatrix) is BooleanType:
            if distanceMatrix == False:
                logging.error("MicroPITA.funcGetCentralSamplesByKMedoids. Received false for betaMetrix matrix generation, returning false.")
            return False

        #Log distance matrix
        logging.debug("".join(["Distance matrix for representative selection using metric=",str(sMetric)]))

        distance = MLPYDistanceAdaptor(npaDistanceMatrix=distanceMatrix, fIsCondensedMatrix=True)

        #Create object to determine clusters/medoids
        medoidsMaker = Kmedoids(k=iNumberSamplesReturned, dist=distance)
        #medoidsData includes(1d numpy array, medoids indexes; 
        #              1d numpy array, non-medoids indexes;
        #              1d numpy array, cluster membership for non-medoids;
        #              double, cost of configuration)
        #npaMatrix is samples x rows
        #Build a matrix of lists of indicies to pass to the distance matrix
        indicesMatrix = []
        for indexPosition in xrange(0,len(npaMatrix)):
            indicesMatrix.append([indexPosition])
        medoidsData = medoidsMaker.compute(np.array(indicesMatrix))
        logging.debug("Results from the kmedoid method in representative selection:")
        logging.debug(str(medoidsData))

        #If returning the same amount of clusters and samples
        #Return centroids
        selectedIndexes = medoidsData[0]
        for index in xrange(0,iNumberSamplesReturned):
            returningSamples.append(lsSampleNames[selectedIndexes[index]])
        return returningSamples

####Group 3## Highest Dissimilarity
    #Testing: Happy path tested
    def funcSelectExtremeSamplesFromHClust(self, strBetaMetric, npaAbundanceMatrix, lsSampleNames, iSelectSampleCount):
	"""
	Select extreme samples from HClustering.
	
	:param	strBetaMetric:	The beta metric to use for distance matrix generation.
	:type	String	The name of the beta metric to use.
	:param	npaAbundanceMatrix:	Numpy array where row=samples and columns=features.
	:type	Numpy Array	Abundance data.
	:param	lsSampleNames:	The names of the sample.
	:type	List	List of strings.
	:param	iSelectSampleCount:	Number of samples to select (return).
	:type	Integer	Integer number of samples returned.
	:return	Samples:	List of samples.
	"""

        #If they want all the sample count, return all sample names
        sampleCount=len(npaAbundanceMatrix[:,0])
        if(iSelectSampleCount==sampleCount):
          return(lsSampleNames)

        #Holds the samples to be returned
        returnSamples = []

        #Holds the roots to each cluster
        clusterRoots = []
        #List of lists. Each internal list is a cluster containing all of the membership's indices
        clusterIndices = []

        #Generate beta matrix
        #Returns condensed matrix
        tempDistanceMatrix = self.funcGetBetaMetric(npadAbundancies=npaAbundanceMatrix, sMetric=strBetaMetric)

        #Feed beta matrix to linkage to cluster
        #Send condensed matrix
        linkageMatrix = hcluster.linkage(tempDistanceMatrix, method=ConstantsMicropita.c_strHierarchicalClusterMethod)

        #Extract cluster information from dendrogram
        #The linakge matrix is of the form
        #[[int1 int2 doube int3],...]
        #int1 and int1 are the paired samples indexed at 0 and up.
        #each list is an entry for a branch that is number starting with the first
        #list being sample count index + 1
        #each list is then named by an increment as they appear
        #this means that if a number is in the list and is = sample count or greater it is not
        #terminal and is instead a branch.
        #This method just takes the lowest metric measurement (highest distance pairs/clusters)
        #Works much better than the original technique
        #get total number of samples
        iSampleCount = len(lsSampleNames)

        iCurrentSelectCount = 0
        for row in linkageMatrix:
            #Get nodes ofthe lowest pairing (so the furthest apart pair)
            iNode1 = int(row[0])
            iNode2 = int(row[1])
            #Make sure the nodes are a terminal node (sample) and not a branch in the dendrogram
            #The branching in the dendrogram will start at the number of samples and increment higher.
            #Add each of the pair one at a time breaking when enough samples are selected.
            if(iNode1<iSampleCount):
                returnSamples.append(lsSampleNames[iNode1])
                iCurrentSelectCount = iCurrentSelectCount + 1
            if iCurrentSelectCount == iSelectSampleCount:
                break
            if(iNode2<iSampleCount):
                returnSamples.append(lsSampleNames[iNode2])
                iCurrentSelectCount = iCurrentSelectCount + 1
            if iCurrentSelectCount == iSelectSampleCount:
                break

        #Return selected samples
        return returnSamples

####Group 4## Rank Average of user Defined Taxa
    #Testing: Happy Path Tested
    def funcGetAverageAbundanceSamples(self, abndTable, lsTargetedFeature, fRank=False):
	"""
	Averages feature abundance or ranked abundance. Expects a column 0 of taxa id that is skipped.
	
	:param	abndTable:	Abundance Table to analyse
	:type	AbundanceTable	Abundance Table
	:param	lsTargetedFeature:	String names
	:type	list	list of string names of features (bugs) which are measured after ranking against the full sample
	:param  fRank:	Indicates to rank the abundance before getting the average abundance of the features (default false)
	:type   boolean	Flag indicating ranking abundance before calculating average feature measurement (false= no ranking)
	:return	List of lists or boolean:	List of lists or False on error. One internal list per sample indicating the sample,
            feature average abundance or ranked abundance. Lists will already be sorted.
            For not Ranked [[sample,average abundance of selected feature,1]]
        	For Ranked [[sample,average ranked abundance, average abundance of selected feature]]
			Error Returns false
	"""
        llAbundance = abndTable.funcGetAverageAbundancePerSample(lsTargetedFeature)
        if not llAbundance:
            logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not get average abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table."]))
            return False
        #Add a space for ranking if needed
        #Not ranked will be [[sSample,average abundance,1]]
        #(where 1 will not discriminant ties if used in later functions, so this generalizes)
        #Ranked will be [[sSample, average rank, average abundance]]
        #Here average abundance can be used incase of tie breaking.
        llRetAbundance = [llist+[1] for llist in llAbundance]

        #Rank if needed
        if fRank:
            abndRanked = abndTable.funcRankAbundance()
            if abndRanked == None:
                logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not rank the abundance table, returned false."]))
                return False
            dictRankedAbundance = dict([tuple(lRankedItem) for lRankedItem in abndRanked.funcGetAverageAbundancePerSample(lsTargetedFeature)])
            if not dictRankedAbundance:
                logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not get average ranked abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table."]))
                return False

            #Check for the case where all features are 0
            #Rank could be 1 and these could be selected
            #Give any sample with an average feature abundance of 0 the worst rank
            iWorstRank = float(abndRanked.funcGetFeatureCount())
            llRetAbundance = [[llist[0],iWorstRank] if llist[1] == 0.0 
                             else [llist[0],dictRankedAbundance[llist[0]],llist[1]]
                             for llist in llRetAbundance]
            
        #Sort first for ties and then for the main feature
        if ConstantsMicropita.c_fBreakRankTiesByDiversity:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[2], reverse = True)
        if fRank:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[1], reverse = False)
        else:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[1], reverse = True)

        #return
        return llRetAbundance

    #Testing: Happy Path Tested
    def funcSelectTargetedTaxaSamples(self, abndMatrix, lsTargetedTaxa, iSampleSelectionCount, sMethod):
      """
      Selects samples with the highest ranks or abundance of targeted features.
      If ranked, select the highest abundance for tie breaking

      :param	abndMatrix:	Abundance table to analyse 
      :type	AbundanceTable	Abundance table
      :param	lsTargetedTaxa:	List of features
      :type	list	list of strings
      :param	iSampleSelectionCount:	Number of samples to select
      :type	integer	integer
      :param	sMethod:	Method to select targeted features
      :type	string	String (Can be values found in microPITA)
      :return	List of strings:	List of sample names which were selected
      List of strings	Empty list is returned on an error.
      """

      #Check data
      if(len(lsTargetedTaxa) < 1):
        logging.error("MicroPITA.funcSelectTargetedTaxaSamples. Taxa defined selection was requested but no features were given.")
        return []
      if not sMethod:
        logging.error("MicroPITA.funcSelectTargetedTaxaSamples. Taxa defined selection was requested but no Method were given.")
        return []

      #Call function
      lsTargetedSamples = False
      if sMethod.lower() == ConstantsMicropita.c_strTargetedRanked.lower():
          lsTargetedSamples = self.funcGetAverageAbundanceSamples(abndTable=abndMatrix, lsTargetedFeature=lsTargetedTaxa, fRank=True)
      elif sMethod.lower() == ConstantsMicropita.c_strTargetedAbundance.lower():
          lsTargetedSamples = self.funcGetAverageAbundanceSamples(abndTable=abndMatrix, lsTargetedFeature=lsTargetedTaxa, fRank=False)

      #If an error occured or the key word for the method was not recognized
      if lsTargetedSamples == False: 
          logging.error("MicroPITA::funcSelectTargetedTaxaSamples: Was not able to select for the features given. So targeted feature selection was performed. Check to make sure the features are spelled correctly and exist in the abundance file.")
          return []

      #Select from results
      return [sSample[0] for sSample in lsTargetedSamples[:iSampleSelectionCount]]


####Group 5## Random

    #Testing: Happy path Tested
    def funcGetRandomSamples(self, lsSamples=None, iNumberOfSamplesToReturn=0):
	"""
	Returns random sample names of the number given. No replacement.
	
	:param	lsSamples:	List of sample names 
	:type	list	list of strings
	:param	iNumberOfSamplesToReturn:	Number of samples to select
	:type	integer	integer.
	:return	List:	List of selected samples (strings).
	"""

        #Input matrix sample count
        sampleCount = len(lsSamples)

        #Validate Number of samples is not negative and not greater than the matrix
        if((iNumberOfSamplesToReturn < 0) or (iNumberOfSamplesToReturn > sampleCount)):
            logging.error("".join(["MicroPITA.funcGetRandomSamples. TempNumberOfSamples was not useful for deriving subset. ",str(iNumberOfSamplesToReturn)]))
            return False

        #Return the full matrix if they ask for a return matrix where length == original
        if(iNumberOfSamplesToReturn == sampleCount):
            return lsSamples

        #Get the random indices for the sample (without replacement)
        randomIndices = random.sample(range(sampleCount), iNumberOfSamplesToReturn)

        #Create a boolean array of if indexes are to be included in the reduced array
        randomIndicesBoolean = list()
        for rIndex in range(sampleCount):
            randomIndicesBoolean.append(rIndex in randomIndices)

        #Reduce array to just what is randomly sampled
        return np.compress(condition=randomIndicesBoolean, a=lsSamples)

####Group 6## Supervised
    def funcRunMLPYSVM(self, abndAbundanceTable, sMetadataForLabel, strInputSVMFile, strPredictionFile):
	"""
	Runs a linear SVM using MLPY.
	
	:param	abndAbundanceTable:	Abundance table of data.
	:type	Abundance table	Abundance table of data.
	:param	sMetadataForLabel:	Metadata for label used to supervised learning.
	:type	string	string
	:param	strInputSVMFile:	File name for the file generated to mock the input SVM file.
	:type	string	string
	:param	strPredictionFile	LIBSVM style output file of probabilistic predictions.
	:type	string	File path
	:return	Analysis Results:	Dictionary of files generated (input and prediction files)
         Return false on error.
	"""

        #SVM general class
        svm = SVM()

        #Convert abundancies file to SVM file
        lsUniqueLabelOrder = SVM.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abndAbundanceTable, strOutputSVMFile=strInputSVMFile, sMetadataLabel=sMetadataForLabel)
        if not lsUniqueLabelOrder:
            logging.error("MicroPITA.funcRunMLPYSVM: Received an error when creating the input SVM file in the MLPY LIBSVM analysis pipeline.")
            return False

        #Will hold the SVM related return data including files generated during the process.
        #Although these files are not necessary for MLPY Libsvm, they are made to mirror the LIBSVM process and
        #and can be used to troubleshoot.
        dictSVMReturn = {}
        dictSVMReturn[ConstantsMicropita.c_strKeywordInputFile] = strInputSVMFile

        #Cost values to optimize to
        liCost = ConstantsMicropita.c_lCostRange

        #Get labels
        lsLabels = abndAbundanceTable.funcGetMetadata(sMetadataForLabel)

        #Get weights for labels
        dictWeights, ldWeightLabels = SVM.funcWeightLabels(lsLabels)
        #Create a new label list but coded as integers
        ldLabels = np.array([ldWeightLabels.index(sLabel) for sLabel in lsLabels])
        logging.debug("".join(["ldLabels ", str(ldLabels)]))

        #Get 2D array of data (has feature names)
        npData = abndAbundanceTable.funcGetAbundanceCopy()

        #Scale features (has feature names)
        for iRowIdex, npadRow in enumerate(npData):
            sRowName = npadRow[0]
            npadRow = np.array(list(npadRow)[1:])
            npData[iRowIdex] = tuple([sRowName]+list(SVM.funcScaleFeature(npadRow)))

        #Get each sample data (this drops the feature names)
        lsSampleNames = np.array(abndAbundanceTable.funcGetSampleNames())

        #SampleName labels
        dictLabels = dict(zip(lsSampleNames,ldLabels))

        #Best Accuracy
        dBestAccuracy = 0
        iBestCost = None
        dictIBestPredictions = None
        dictDBestProbabilities = None
        dictAllProbabilities = None
        lSVMLabels = None

        #Flag used to progress ties up the validation half the time
        #This allows the cost value selected in a ties situation to be in the
        #middle and not an extreme of the cost value plateau created by the ties
        fTiesProgress = True

        #For each cost
        for iCost in liCost:

            dictiPrediction = dict()
            dictdProbability = dict()
            dictAllProbabilities = dict()
            dictAllPredictions = dict()

            #For each cross validation fold for the cost
            for lfTraining, lfValidation in svm.func10FoldCrossvalidation(iTotalSampleCount = len(lsSampleNames), fRandomise = True):

                #Create SVM object
                svmMLPY = mlpy.LibSvm(C=math.pow(2,iCost), probability=True, weight=dictWeights)

                #Training set
                lsTrainingSamples = lsSampleNames.compress(lfTraining)
                ldTrainingLabels = ldLabels.compress(lfTraining)
                llTrainingData = [list(npData[sSamples]) for sSamples in lsTrainingSamples]

                #Validation set
                lsValidationSamples = lsSampleNames.compress(lfValidation)
                ldValidationLabels = ldLabels.compress(lfValidation)
                llValidationData = [list(npData[sSamples]) for sSamples in lsValidationSamples]

                #Create classifier
                #x = [sample of the same number as labels, features]
                #y = [labels]
                svmMLPY.learn(x=llTrainingData,y=ldTrainingLabels)

                #Get Classes
                #[samples,features]
                lPredictions = svmMLPY.pred(llValidationData)
                lSVMLabels = list(svmMLPY.labels())

                #Get probabilities
                #svm functions require [samples,features]
                npaDistances = svmMLPY.pred_probability(llValidationData)

                #Store the fold update
                ldictReturn = self._funcStoreSVMProbability(lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,
                                                            dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions)
                dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = ldictReturn

            #Get Accuracy and if is the best accuracy store
            dAccuracy = sum(dictiPrediction.values())/float(len(lsSampleNames))
            if dAccuracy > dBestAccuracy:
                iBestCost = int(iCost)
                dBestAccuracy = dAccuracy
                dictIBestPredictions = dictiPrediction
                dictDBestProbabilities = dictdProbability
            elif dAccuracy == dBestAccuracy:
                if fTiesProgress:
                    iBestCost = int(iCost)
                    dBestAccuracy = dAccuracy
                    dictIBestPredictions = dictiPrediction
                    dictDBestProbabilities = dictdProbability
                fTiesProgress = not fTiesProgress

        #Create SVM object
        svmMLPY = mlpy.LibSvm(C=math.pow(2,iBestCost), probability=True, weight=dictWeights)

        #Get full data
        llFullData = [list(npData[sSamples]) for sSamples in lsSampleNames]

        #Create classifier
        #x = [sample of the same number as labels, features]
        #y = [labels]
        svmMLPY.learn(x=llFullData,y=ldLabels)

        #Get Classes
        #[samples,features]
        lPredictions = svmMLPY.pred(llFullData)
        lSVMLabels = list(svmMLPY.labels())

        #Get probabilities
        #svm functions require [samples,features]
        npaDistances = svmMLPY.pred_probability(llFullData)

        #Store all sample data
        ldictReturn = self._funcStoreSVMProbability(lsSampleNames,ldLabels,lSVMLabels,npaDistances,lPredictions,
                                                        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions)

        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = ldictReturn

        #Log best
        logging.debug("".join(["funcRunMLPYSVM::Best Accuracy=", str(dBestAccuracy)]))
        logging.debug("".join(["funcRunMLPYSVM::Best Cost=", str(iBestCost)]))

        #Create output prediction file
        strPredictionOutput = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace.join(["labels"]+[str(iLabel) for iLabel in lSVMLabels])+ConstantsMicropita.c_strEndline
        strPredictionOutput = strPredictionOutput + ConstantsMicropita.c_strEndline.join([ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace.join([str(dProb) for dProb in dictAllProbabilities[sSampleName]])
            for sSampleName in lsSampleNames])

        #Write prediction file to file
        with open(strPredictionFile, 'w') as f:
            f.write(strPredictionOutput)
        dictSVMReturn[ConstantsMicropita.c_strKeywordPredFile] = strPredictionFile

        #Return
        return dictSVMReturn

    #Testing: Happy path tested
    def _funcStoreSVMProbability(self,lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions):
        """
        Takes the validation probabilities and stores in a collection of dictionaries

        :param	lsValidationSamples:	List of samples (strings) being evaluated
        :type	List of strings	Order is important, should match npaDistances
        :param	ldValidationLabels:	Original labels of each sample
        :type	List of int	Order is important, should lsValidationSamples
        :param	lSVMLabels:	Unique labels in order of the probabilities/distances in (npaDistances)
        :type	List of int	Labeled ordering
        :param	npaDistances:	Distances from the MLPY pred and probablity functions
        :type	Numpy Array
        :param	lPredictions:	List of labels which are the prediction from MLPY from the pred()
        :type	List of int	Order matching lsValidationSamples
        :param	dictdProbability:	Dictionary to add probabilities {"SampleName":probability}
        :type	Dictionary	Order is important, should match npaDistances
        :param	dictAllProbabilities:	Dictionary to add all returned probabilities {"SampleName":[probability label 1,probablitity label 2]...}
        :type	Dictionary	Order is important, should math npaDistances
        :param	dictiPrediction:	Dictionary to add label prdictions {"SampleName":label}
        :type	Dictionary	Order is important, should math npaDistances
        :param	dictAllPredictions:	Dictionary to add all returned predictions {"SampleName":[label 1,label 2]...}
        :type	Dictionary	Order is important, should math npaDistances
        :return	Analysis Results:	[dictIBestPredictions (dictionary of predictions),
                                         dictDBestProbabilities (dictionary of probabilities),
                                         dictLabels (list of original labels),
                                         lsLabels (list of the labels in order of the MLPY output, these are unique values]
        """

        #Store fold results
        for indexSamples, sSampleName in enumerate(lsValidationSamples):
            #Determine label by highest probability
            dMaxProbability = max(npaDistances[indexSamples])
            dictdProbability[sSampleName] = dMaxProbability
            iLabel = lSVMLabels[list(npaDistances[indexSamples]).index(dMaxProbability)]
            dictiPrediction[sSampleName] = ldValidationLabels[indexSamples] == iLabel
            dictAllProbabilities[sSampleName] = [int(iLabel)]+list(npaDistances[indexSamples])
            dictAllPredictions[sSampleName] = str(int(iLabel))

        return [dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions]

    #Happy path tested (case 3)
    def funcGetAveragePopulation(self, abndTable, lfCompress):
        """
        Get the average row per column in the abndtable.

        :param abndTable: AbundanceTable of data to be averaged
        :type AbudanceTable
        :param lfCompress: List of boolean flags (false means to remove sample before averaging
        :type List of floats
        :return List of doubles: 
        """
        #Get the average populations
        lAverage = []

        for sFeature in abndTable.funcGetAbundanceCopy():
            sFeature = list(sFeature)[1:]
            sFeature=np.compress(lfCompress,sFeature,axis=0)
            #If there are no samples then return empty list.
            if len(sFeature) == 0:
                lAverage.append(0.0)
            #Return average of feature
            else:
                lAverage.append(sum(sFeature)/float(len(sFeature)))
        return lAverage

    #Happy path tested (2 cases)
    def funcGetDistanceFromAverage(self, abndTable,ldAverage,lsSamples,lfSelected):
        """
        Given an abundance table and an average sample, this returns the distance of each sample
        (measured using brays-curtis dissimilarity) from the average.
        The distances are reduced by needing to be in the lsSamples and being a true in the lfSelected
        (which is associated with the samples in the order of the samples in the abundance table;
        use abundancetable.funcGetSampleNames() to see the order if needed).

        :param abndTable: Abundance table holding the data to be analyzed.
        :type AbundanceTable
        :param ldAverage: Average population (Average features of the abundance table of samples)
        :type List of doubles which represent the average population
        :param lsSamples: These are the only samples used in the analysis
        :type List of strings (sample ids)
        :param lfSelected: Samples to be included in the analysis
        :type List of boolean (true means include)
        :return: List of distances (doubles)
        """
        #Get the distance from label 1 of all samples in label0 splitting into selected and not selected lists
        ldSelectedDistances = []
        ldNotSelectedDistances = []

        for sSampleName in [sSample for iindex, sSample in enumerate(lsSamples) if lfSelected[iindex]]:
            #Get the sample measurements
            ldSelectedDistances.append(Metric.funcGetBrayCurtisDissimilarity(np.array([abndTable.funcGetSample(sSampleName),ldAverage]))[0])
        return ldSelectedDistances

    #Happy path tested (1 case)
    def funcMeasureDistanceFromLabelToAverageOtherLabel(self, abndTable, lfGroup, lfGroupOther):
        """
        Get the distance of samples from one label from the average sample of not the label.
        Note: This assumes 2 classes.  

        :param abndTable: Table of data to work out of.
        :type: Abundace Table
        :param lfGroup: Boolean indicator of the sample being in the frist group
        :type: It is assumed that False is the other label and there is not 3rd class (like misclassfied).
        :return List of List of doubles: [list of tuples (string sample name,double distance) for the selected population, list of tuples for the not selected population]
        """
        #Hold data for combined boxplots
        llBoxplotData = []
        lsBoxplotLabels = []

        #Get all sample names
        lsAllSamples = abndTable.funcGetSampleNames()

        #Get average populations
        lAverage = self.funcGetAveragePopulation(abndTable=abndTable, lfCompress=lfGroup)
        lAverageOther = self.funcGetAveragePopulation(abndTable=abndTable, lfCompress=lfGroupOther)

        #Get the distance from the average of the other label (label 1)
        ldSelectedDistances = self.funcGetDistanceFromAverage(abndTable = abndTable, ldAverage = lAverageOther,
                                                    lsSamples = lsAllSamples, lfSelected = lfGroup)
        ldNotSelectedDistances = self.funcGetDistanceFromAverage(abndTable = abndTable,ldAverage = lAverage,
                                                    lsSamples = lsAllSamples, lfSelected = lfGroupOther)

        ldSelectedDistances = zip([lsAllSamples[iindex] for iindex, fGroup in enumerate(lfGroup) if fGroup],ldSelectedDistances)
        ldNotSelectedDistances = zip([lsAllSamples[iindex] for iindex, fGroupOther in enumerate(lfGroupOther) if fGroupOther],ldNotSelectedDistances)
        return [ldSelectedDistances,ldNotSelectedDistances]

    #Happy path tested (1 test case)
    def funcPerformDistanceSelection(self, abndTable, iSelectionCount, sLabel):
        """
        Given labels, labels from one label are measured from the average (centroid) value of another group.
        Currently runs on 2 labels.

        :params  abndTable: Abundance of measurements
        :type AbundanceTable: 
        :params iSelectionCount: The number of samples selected per sample.
        :type Integer: Integer greater than 0
        :params sLabel: Label used by supervised methods.
        :type String: 
        """
        llBoxplotData = []

        #Get labels
        liLabels = SVM.funcMakeLabels(abndTable.funcGetMetadata(sLabel))
        #Check labels
        if not len(set(liLabels)) == 2:
           print "Did not get two labels, errors"
           return False
        liUniqueLabels = list(set(liLabels))
  
        #Get boolean indicator of labels
        lfLabels1 = [liUniqueLabels[0] == iLabel for iLabel in liLabels]
        lfLabels2 = [not fLabel for fLabel in lfLabels1]

        #Selection in 16S
        #Get the distances for 16S data from two different groups to the average of the other
        ldLabel1Distances, ldLabel2Distances = self.funcMeasureDistanceFromLabelToAverageOtherLabel(abndTable, lfLabels1, lfLabels2)
        ldLabel1Distances = sorted(ldLabel1Distances,key=operator.itemgetter(1))
        ldLabel2Distances = sorted(ldLabel2Distances,key=operator.itemgetter(1))

        #Get the closest and farthest distances
        ltupleDiscriminantSamples0 = ldLabel1Distances[:iSelectionCount]
        ltupleDiscriminantSamples1 = ldLabel2Distances[:iSelectionCount]
        ltupleDistinctSamples0 = ldLabel1Distances[iSelectionCount*-1:]
        ltupleDistinctSamples1 = ldLabel2Distances[iSelectionCount*-1:]

        #Remove the selected samples from the larger population of distances (better visualization)
        ldSelected = [tpleSelected[0] for tpleSelected in ltupleDiscriminantSamples0+ltupleDiscriminantSamples1+ltupleDistinctSamples0+ltupleDistinctSamples1]

        #Set up data and labels (label 0)
        return [ltupleDiscriminantSamples0, ltupleDistinctSamples0,
               [tplData for tplData in ldLabel1Distances if tplData[0] not in ldSelected],
               ltupleDiscriminantSamples1, ltupleDistinctSamples1,
               [tplData for tplData in ldLabel2Distances if tplData[0] not in ldSelected]]

###
    #Run the supervised method surrounding distance from centroids
    #Happy path tested (3 test cases)
    def funcRunSupervisedDistancesFromCentroids(self, abundanceTable, fRunDistinct, fRunDiscriminant,
                                       strOuputSVMFile, strPredictSVMFile, strSupervisedMetadata,
                                       iSampleSVMSelectionCount):
        """
	Runs supervised methods based on measuring distances of one label from the centroid of another.
        This does make output files like the other supervised methods excluding one difference,
        the output "prediction" file which is where one would find the list of differences for each sample
        in the order of the sample names (after trimming away samples with invalid metadata values), has
        only one label listed in the header "other". This is because the method only generates one distance
        for a sample, it's distance from the other label centroid.
	
	:param	abundanceTable:	AbundanceTable
	:type	AbudanceTable	Data to analyze
	:param	fRunDistinct:	Run distinct selection method
	:type	Boolean	boolean (true runs method)
	:param	fRunDiscriminant:	Run discriminant method
	:type	Boolean	boolean (true runs method)
	:param	strOutputSVMFile:	File output from  SVM (scaled input file in the style of LIBSVM)
	:type	String	String
	:param	strPredictSVMFile:	File label prediction from  SVM
	:type	String	String
	:param	iSampleSVMSelectionCount:	Number of samples to select
	:type	Integer	int sample selection count
	:return	Selected Samples:	A dictionary of selected samples by selection ID
        Dictionary	{"Selection Method":["SampleID","SampleID"...]}
	"""

        #Run supervised blocks
        #Select supervised (using SVM)
        #Will contain the samples selected to return
        dictSelectedSamples = dict()

        #For now perform the selection here, should eventually be performed in micropta
        ltupleDisc0, ltupleDist0, ltpleNotSelected0, ltupleDisc1, ltupleDist1, ltpleNotSelected1 = self.funcPerformDistanceSelection(abndTable=abundanceTable,
                                                                                                      iSelectionCount=iSampleSVMSelectionCount,
                                                                                                      sLabel=strSupervisedMetadata)

        #Make expected output files for supervised methods
        #1. Output file which is similar to an input file for SVMs
        #2. Output file that is similar to the probabilitic output of a SVM (LibSVM)
        #Manly for making output of supervised methods (SVM and Distince from Centroid) similar
        #MicropitaVis needs some of these files
        lsAllSampleNames = abundanceTable.funcGetSampleNames()
        lsUniqueLabels = SVM.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abundanceTable, strOutputSVMFile=strOuputSVMFile, sMetadataLabel=strSupervisedMetadata)
        lsLabels = SVM.funcReadLabelsFromFile(sSVMFile=strOuputSVMFile, lsAllSampleNames=lsAllSampleNames, isPredictFile=False)
        lsLabels = [(sSample,sLabel) for sLabel in lsLabels.keys() for sSample in lsLabels[sLabel]]
        dictDistances = dict(ltupleDisc0+ltupleDist0+ltpleNotSelected0+ltupleDisc1+ltupleDist1+ltpleNotSelected1)
        dictLabels = dict(lsLabels)

        strPredictionOutput = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace.join(["labels"]+["Other"])+ConstantsMicropita.c_strEndline
        strPredictionOutput = strPredictionOutput + ConstantsMicropita.c_strEndline.join([ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace.join([str(dictLabels[sSample]),str(dictDistances[sSample])]) for sSample in lsAllSampleNames])

        #Write prediction file to file
        with open(strPredictSVMFile, 'w') as f:
            f.write(strPredictionOutput)

        if fRunDistinct:
            dictSelectedSamples[ConstantsMicropita.c_strDistinct] = [ltple[0] for ltple in ltupleDist0 + ltupleDist1]
        if fRunDiscriminant:
            dictSelectedSamples[ConstantsMicropita.c_strDiscriminant] = [ltple[0] for ltple in ltupleDisc0 + ltupleDisc1]
        return dictSelectedSamples

    #Run the supervised methods
    def funcRunSupervisedMethods(self, abundanceTable, fRunDistinct, fRunDiscriminant,
                                       strOuputSVMFile, strPredictSVMFile, strSupervisedMetadata,
                                       iSampleSVMSelectionCount):
        """
	Runs supervised methods.
	
	:param	abundanceTable:	AbundanceTable
	:type	AbudanceTable	Data to analyze
	:param	fRunDistinct:	Run distinct selection method
	:type	Boolean	boolean (true runs method)
	:param	fRunDiscriminant:	Run discriminant method
	:type	Boolean	boolean (true runs method)
	:param	strOutputSVMFile:	File output from  SVM (scaled input file in the style of LIBSVM)
	:type	String	String
	:param	strPredictSVMFile:	File label prediction from  SVM
	:type	String	String
        :param  strSupervisedMetadata:	Label for supervised method
        :type   String:	String
	:param	iSampleSVMSelectionCount:	Number of samples to select
	:type	Integer	int sample selection count
	:return	Selected Samples:	A dictionary of selected samples by selection ID
        Dictionary	{"Selection Method":["SampleID","SampleID"...]}
	"""

        #Run supervised blocks
        #Will contain the samples selected to return
        dictSelectedSamples = dict()
        #SVM related output from the SVM call
        svmRelatedData = None

        #Remove all files associated with supervised methods
        for f in [strOuputSVMFile,strPredictSVMFile]:
            if os.path.exists(f):
                os.remove(f)

        #True = Select supervised (using SVM)
        #False = Select supervised (using distances)
        if ConstantsMicropita.fRunSVM:
            #Run MLPY SVM
            svmRelatedData = self.funcRunMLPYSVM(abndAbundanceTable=abundanceTable, sMetadataForLabel=strSupervisedMetadata,
                                             strInputSVMFile=strOuputSVMFile, strPredictionFile=strPredictSVMFile)

            #Read in prediction file and select samples
            if svmRelatedData:
                dictSelectedSamples = self._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=svmRelatedData[ConstantsMicropita.c_strKeywordInputFile],
                                                      strPredictFilePath=svmRelatedData[ConstantsMicropita.c_strKeywordPredFile], lsSampleNames=abundanceTable.funcGetSampleNames(),
                                                      iSelectCount=iSampleSVMSelectionCount, fSelectDiscriminant = fRunDiscriminant, fSelectDistinct = fRunDistinct)
            return dictSelectedSamples

        else:
            #Run supervised blocks
            #Select supervised (using distances from the centroid other label)
            return self.funcRunSupervisedDistancesFromCentroids(abundanceTable=abundanceTable, fRunDistinct=fRunDistinct, fRunDiscriminant=fRunDiscriminant,
                                       strOuputSVMFile=strOuputSVMFile, strPredictSVMFile=strPredictSVMFile, strSupervisedMetadata=strSupervisedMetadata,
                                       iSampleSVMSelectionCount=iSampleSVMSelectionCount)

    #Testing: Happy path tested
    def _funcSelectSupervisedSamplesFromPredictFile(self, strOriginalInputFile, strPredictFilePath, lsSampleNames, iSelectCount, fSelectDiscriminant, fSelectDistinct):
        """
        Selects sample for the supervised distinct and discriminant methods from an SVM predict file.

        :param	strOriginalInputFile:	The original input file to parse the original class labels.
        :type	String
        :param	strPredictFilePath:	A file path to the predict file generated by the SVM methods.
        :type	String
        :param	lsSampleNames:	List of string sample ids.
        :type	String
        :param	iSelectCount:	Amount of samples to select per label.
        :type	Integer
        :param	fSelectDiscriminant:	Indicates samples should be selected for the disriminant method.
        :type	Boolean			
        :param	fSelectDistinct:	Boolean. Indicates samples should be selected for the distinct method.
        :type	Boolean
        :return    Dictionary    Dictionary of selected samples {"technique":[sampleName1,sampleName2,sampleName3...]}
        """

        #Create and array to hold difference of the samples probabilities from the central probability
        dictDistanceFromHyperplane = dict()

        #Holds selected ids of selected samples
        dictSelectedSamples = dict()

        #Holds labels to compare to the predictions
        lsOriginalLabels = None
        #Open prediction file and input file and get labels to compare to the predictions
        with open(strOriginalInputFile,'r') as f, open(strPredictFilePath,'r') as g:
            reader = csv.reader(f, delimiter=ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace, quoting=csv.QUOTE_NONE)
            lsOriginalLabels = [row[0] for row in reader]
            predictionLists = g.read()
            predictionLists = [filter(None,strPredictionList) for strPredictionList in predictionLists.split(ConstantsMicropita.c_strEndline)]

        #Get label count (meaning the number of label categories)(-1 to not count the predicted first entry)
        labelCount = len(predictionLists[0].split(ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace))-1

        #Central probability
        centralProbability = 1.0 / float(labelCount)
        logging.debug("centralProbability")
        logging.debug(centralProbability)

        #For each line in the file subtract all but the first item from the central value
        #The first item being the header line not a row of probabilities
        #Selected_label prob_for_label1 prob_ForLabel2....
        #Remove prediction label header
        iCurrectlyClassifiedCount = 0
        predictionLists = predictionLists[1:]
        for lineIndex in xrange(0,len(predictionLists)):
            logging.debug("lineIndex")
            logging.debug(lineIndex)

            #Split line into elements by whitespace and remove first element (the label)
            lineElements = predictionLists[lineIndex]
            #Skip blank entries
            if lineElements.strip() == '':
                continue

            #Get label and probabilities
            lineElements = lineElements.split(ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)
            iCurLabel = str(lineElements[0])
            lineElements = lineElements[1:]

            #Only work with samples that are correctly predicted
            if iCurLabel == lsOriginalLabels[lineIndex]:
                iCurrectlyClassifiedCount = iCurrectlyClassifiedCount + 1
                logging.debug("Correctly predicted")
                #Sum the absolute values of the differences
                #Store the index and the deviation from the central probability
                deviation = sum([math.pow((float(strPrediction)-centralProbability),2) for strPrediction in lineElements])
                curSVMData = dictDistanceFromHyperplane.get(iCurLabel, [])
                curSVMData.append([deviation,lineIndex])
                dictDistanceFromHyperplane[iCurLabel] = curSVMData
                logging.debug("curSVMData")
                logging.debug(curSVMData)

            logging.debug("".join(["iCurrectlyClassifiedCount=",str(iCurrectlyClassifiedCount)]))

        #Sort sample by distance from the center and take the top N indexes
        for scurKey in dictDistanceFromHyperplane:
            lcurLabelSamples = dictDistanceFromHyperplane[scurKey]
            dictDistanceFromHyperplane[scurKey] = sorted(lcurLabelSamples, key=operator.itemgetter(0))
            logging.debug("dictDistanceFromHyperplane[scurKey]")
            logging.debug(dictDistanceFromHyperplane[scurKey])

        selectedSamplesIndicesClose = list()
        selectedSamplesIndicesFar = list()
        #If the amount of samples needed are greater than what was analyzed with the SVM,
        #Return them all for both far and near.
        #Get the samples closeRuns to hyperplane
        #Get samples farthest from hyperplane
        #Make this balanced to the labels so for each sample label select the iSelectCount count of samples
        for scurKey in dictDistanceFromHyperplane:
            lcurDeviations = dictDistanceFromHyperplane[scurKey]
            iLengthCurDeviations = len(lcurDeviations)
            if(iSelectCount > iLengthCurDeviations):
                licurIndices = [measurement[1] for measurement in lcurDeviations]
                selectedSamplesIndicesClose.extend(licurIndices)
                selectedSamplesIndicesFar.extend(licurIndices)
            else:
                selectedSamplesIndicesClose.extend([measurement[1] for measurement in lcurDeviations[0:iSelectCount]])
                selectedSamplesIndicesFar.extend([measurement[1] for measurement in lcurDeviations[(iLengthCurDeviations-iSelectCount):iLengthCurDeviations]])
        logging.debug("selectedSamplesIndicesClose")
        logging.debug(selectedSamplesIndicesClose)
        logging.debug("selectedSamplesIndicesFar")
        logging.debug(selectedSamplesIndicesFar)
        logging.debug("sampleNames")
        logging.debug(lsSampleNames)

        if fSelectDiscriminant:
            SVMSamples = list()
            for selectedSampleIndex in selectedSamplesIndicesClose:
                SVMSamples.append(lsSampleNames[selectedSampleIndex])
            dictSelectedSamples[ConstantsMicropita.c_strSVMClose]=SVMSamples
            logging.debug("SVMSamples Close")
            logging.debug(SVMSamples)
        if fSelectDistinct:
            SVMSamples = list()
            for selectedSampleIndex in selectedSamplesIndicesFar:
                SVMSamples.append(lsSampleNames[selectedSampleIndex])
            dictSelectedSamples[ConstantsMicropita.c_strSVMFar]=SVMSamples
            logging.debug("SVMSamples Far")
            logging.debug(SVMSamples)

        return dictSelectedSamples

    def _funcRunNormalizeSensitiveMethods(self, abndData, iSampleSelectionCount, dictSelectedSamples, lsAlphaMetrics, lsBetaMetrics, lsInverseBetaMetrics,
                                                fRunDiversity, fRunRepresentative, fRunExtreme):
        """
        Manages running methods that are sensitive to normalization. This is called twice, once for the set of methods which should not be normalized and the other
        for the set that should be normalized.

        :param	abndData:	Abundance table object holding the samples to be measured.
        :type	AbundanceTable
        :param	iSampleSelectionCount	The number of samples to select per method.
        :type	Integer
        :param	dictSelectedSamples	Will be added to as samples are selected {"Method:["strSelectedSampleID","strSelectedSampleID"...]}.
        :type	Dictionary
        :param	lsAlphaMetrics:	List of alpha metrics to use on alpha metric dependent assays (like highest diversity).
        :type	List of strings
        :param	lsBetaMetrics:	List of beta metrics to use on beta metric dependent assays (like most representative).
        :type	List of strings
        :param	lsInverseBetaMetrics:	List of inverse beta metrics to use on inverse beta metric dependent assays (like most dissimilar).
        :type	List of strings
        :return	Dictionary:	Returns dictSelectedSamples with the addition of any newly measured samples.
        :param	fRunDiversity:	Run Diversity based methods (true indicates run).
        :type	Boolean	
        :param	fRunRepresentative:	Run Representative based methods (true indicates run).
        :type	Boolean	
        :param	fRunExtreme:	Run Extreme based methods (true indicates run).
        :type	Boolean		
        """

        #Sample ids/names
        lsSampleNames = abndData.funcGetSampleNames()

        #Generate alpha metrics and get most diverse
        if(fRunDiversity):
            #If the table is summed get just the terminal taxa
            internalAlphaMatrix = None
            if abndData.funcIsSummed():
                lsFileElements = os.path.splitext(abndData.funcGetName())
                npaTerminalAbundance = abndData.funcGetFeatureAbundanceTable(abndData.funcGetTerminalNodes()).funcGetAbundanceCopy()
                #Get Alpha metrics matrix
                #Expects Observations (Taxa (row) x sample (column))
                #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                internalAlphaMatrix = Metric.funcBuildAlphaMetricsMatrix(npaSampleAbundance = npaTerminalAbundance, lsSampleNames = lsSampleNames, lsDiversityMetricAlpha = lsAlphaMetrics)
            else:
                #Get Alpha metrics matrix
                #Expects Observations (Taxa (row) x sample (column))
                #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                internalAlphaMatrix = Metric.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abndData.funcGetAbundanceCopy(), lsSampleNames = lsSampleNames, lsDiversityMetricAlpha = lsAlphaMetrics)

            #Get top ranked alpha diversity by most diverse
            #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
            #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
            mostDiverseAlphaSamplesIndexes = self.funcGetTopRankedSamples(lldMatrix=internalAlphaMatrix, lsSampleNames=lsSampleNames, iTopAmount=iSampleSelectionCount)

            #Add to results
            for index in xrange(0,len(lsAlphaMetrics)):
                astrSelectionMethod = ConstantsMicropita.dictConvertAMetricDiversity[lsAlphaMetrics[index]]
                if not astrSelectionMethod in dictSelectedSamples:
                    dictSelectedSamples[astrSelectionMethod]=list()
                dictSelectedSamples[astrSelectionMethod].extend(mostDiverseAlphaSamplesIndexes[index])

        logging.info("Selected Samples 1b")
        logging.info(dictSelectedSamples)

        #Generate beta metrics and 
        if((fRunRepresentative)or(fRunExtreme)):

            #Abundance matrix transposed
            npaTransposedAbundance = UtilityMath.funcTransposeDataMatrix(abndData.funcGetAbundanceCopy(), fRemoveAdornments=True)

            #Get center selection using clusters/tiling
            #This will be for beta metrics in normalized space
            if(fRunRepresentative):
                logging.info("Performing representative selection on normalized data.")
                for bMetric in lsBetaMetrics:

                    #Get representative dissimilarity samples
                    medoidSamples=self.funcGetCentralSamplesByKMedoids(npaMatrix=npaTransposedAbundance, sMetric=bMetric, lsSampleNames=lsSampleNames, iNumberSamplesReturned=iSampleSelectionCount)

                    if(not medoidSamples == False):
                        astrSelectionMethod = ConstantsMicropita.dictConvertBMetricRepresentative[bMetric]
                        if not astrSelectionMethod in dictSelectedSamples:
                            dictSelectedSamples[astrSelectionMethod]=list()
                        dictSelectedSamples[astrSelectionMethod].extend(medoidSamples)

            #Get extreme selection using clusters, tiling
            if(fRunExtreme):
                logging.info("Performing extreme selection on normalized data.")
                #Run KMedoids with inverse custom distance metric in normalized space
                for bMetric in lsInverseBetaMetrics:

                    #Samples for repersentative dissimilarity
                    #This involves inverting the distance metric,
                    #Taking the dendrogram level of where the number cluster == the number of samples to select
                    #Returning a repersentative sample from each cluster
                    extremeSamples = self.funcSelectExtremeSamplesFromHClust(strBetaMetric=bMetric, npaAbundanceMatrix=npaTransposedAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount=iSampleSelectionCount)

                    #Add selected samples
                    if(not extremeSamples == False):
                        astrSelectionMethod = ConstantsMicropita.dictConvertBMetricExtreme[bMetric]
                        if not astrSelectionMethod in dictSelectedSamples:
                            dictSelectedSamples[astrSelectionMethod]=list()
                        dictSelectedSamples[astrSelectionMethod].extend(extremeSamples)

            logging.info("Selected Samples 2,3b")
            logging.info(dictSelectedSamples)
            return dictSelectedSamples

    #Start micropita selection
    def funcRun(self, fIsAlreadyNormalized, fCladesAreSummed, sMetadataID, sLastMetadataName, strInputAbundanceFile,
                      strInputPredictFile, strPredictPredictFile, strCheckedAbndFile, strOutputFile="MicroPITAOutput.txt",
                      cDelimiter = ConstantsMicropita.c_cTab, cFeatureNameDelimiter = "|",
                      strUserDefinedTaxaFile=None, iSampleSelectionCount=0, iSupervisedSampleCount=0,
                      strSelectionTechnique=None, strLabel=None, strStratify=None, fSumData=True, sFeatureSelectionMethod=None):
	"""
	Writes the selection of samples by method to an output file.
	
	:param	fIsAlreadyNormalized:	Indicates if the abundance table is normalized.
	:type	boolean	boolean indicator if the table is normalized (true= normalized).
	:param	fCladesAreSummed:	Indicates if the abundance table is summed.
	:type	boolean	boolean indicator if the table is summed (true= summed).
	:param	strOutputFile:	File to store selection data.
	:type	String	String file path.
	:param	cDelimiter:	Delimiter of abundance table.
	:type	Character	Char (default TAB).
	:param	cFeatureNameDelimiter:	Delimiter of the name of features (for instance if they contain consensus lineages indicating clades).
	:type	Character	Char (default |).
	:param	strInputAbundanceFile:	Abundance table data file.
	:type	String	String path to abundance table file.
	:param	strUserDefinedTaxaFile:	File containing features to select for.
	:type	String	String path to existing file.
	:param	strCheckedAbndFile:	After the input file is checked it will be saved as this file name.
	:type	String String file path.
	:param	iSampleSelectionCount:	Number of samples to select with unsupervised methods.
	:type	Integer	integer.
	:param	iSupervisedSampleCount:	Number of samples to select with supervised methods.
	:type	Integer	integer.
	:param	strSelectionTechnique:	List of strings indicating selection techniques.
	:type	String	List of strings each a selection technique.
	:param	strLabel:	The metadata used for supervised labels.
	:type	String	String (metadata id).
	:param	strStratify:	The metadata used to stratify unsupervised data.
	:type	String	String (metadata id).
	:param	sMetadataID:	The id of the metadata used as an id for each sample.
	:type	String	String metadata id.
	:param	sLastMetadataName:	The id of the metadata positioned last in the abundance table.
	:type	String	String metadata id.
	:param	sFeatureSelectionMethod:	Which method to use to select features in a targeted manner (Using average ranked abundance or abundance).
	:type	String	String (specific values indicated in microPITA).
	:return	Selected Samples:	Samples selected by methods.
            Dictionary	{"Selection Method":["SampleID","SampleID","SampleID",...]}
	"""

        #Holds the top ranked samples from different metrics
        #dict[metric name] = [samplename,samplename...]
        selectedSamples = dict()

        #Check parameters
        if not sMetadataID or not sLastMetadataName or not cDelimiter or not cFeatureNameDelimiter:
          if not sMetadataID:
            logging.error("MicroPITA.funcRun. Please specify Metadata ID. Stopped.")
          if not sLastMetadataName:
            logging.error("MicroPITA.funcRun. Please specify Last meta data name. Stopped.")
          if not cDelimiter:
            logging.error("MicroPITA.funcRun. Please specify delimiter. Stopped.")
          if not cFeatureNameDelimiter:
            logging.error("MicroPITA.funcRun. Please specify feature name delimiter. Stopped.")
          return False

        #If a target feature file is given make sure that targeted feature is in the selection methods, if not add
        if strUserDefinedTaxaFile or ConstantsMicropita.c_strTaxa in strSelectionTechnique:
          if not (strUserDefinedTaxaFile and sFeatureSelectionMethod):
            logging.error("MicroPITA::Did not receive both the Targeted feature file and the feature selection method. MicroPITA did not run.")
            return False
          elif ConstantsMicropita.c_strTaxa not in strSelectionTechnique:
            strSelectionTechnique.append(ConstantsMicropita.c_strTaxa)

        #Check supervised related arguments
        #If there is a superversied count or a label or a supervised method, all three previously mentioned variables must be valid.
        #Otherwise stop
        iActualSupervisedMethodCount = len(set(ConstantsMicropita.c_lsAllSupervisedMethods)&set(strSelectionTechnique))
        if ((iActualSupervisedMethodCount > 0 or strLabel or iSupervisedSampleCount)):
          if not strLabel:
            logging.error("MicroPITA::Did not received a value for sLabel. MiroPITA did not run. Received="+str(strLabel))
            return False
          if iSupervisedSampleCount < 1:
            logging.error("".join(["MicroPITA::Did not receive a selection count for supervised selection above 0, received=",
                          str(iSupervisedSampleCount),". Did not continue analysis."]))
            return False
          if iActualSupervisedMethodCount < 1:
            logging.error("MicroPITA::Did not receive a supervised method even though a supervised count or label was given.")
            return False

        #Check unsupervised methods
        #If an unsupervised count is specified so should a method.
        #Otherwise stop
        iActualUnsupervisedMethodCount = len(set(ConstantsMicropita.c_lsAllUnsupervisedMethods)&set(strSelectionTechnique))
        if ((iActualUnsupervisedMethodCount > 0) or (iSampleSelectionCount)):
          if iSampleSelectionCount < 1:
            logging.error("".join(["MicroPITA::Did not receive a selection count for unsupervised selection above 0, received=",
                                   str(iSampleSelectionCount),". Did not continue analysis."]))
            return False
          if iActualUnsupervisedMethodCount < 1:
            logging.error("MicroPITA::Did not receive a selection count but did receive selection methods. MicroPITA did not run.")
            return False

        #If unsupervised stratify is given, make sure an unsupervised method is given
        if strStratify:
           if iActualUnsupervisedMethodCount < 1:
              logging.error("MicroPITA::Did not receive an unsupervised method when aske dto stratify an unsupervised method. MicroPITA did not run.")
              return False

        #After checking parameters, if still no selection counts or methods are indicated return
        if iSampleSelectionCount+iSupervisedSampleCount < 1:
          logging.error("MicroPITA.funcRun. Please specify a selection amount. Stopped.")
          return False
        if len(strSelectionTechnique) < 1:
          logging.error("MicroPITA.funcRun. Please specify a selection technique. Stopped.")
          return False

        #Create file paths if not already given
        lsFilesToMake = []
        if not strCheckedAbndFile:
          strCheckedAbndFile = os.path.splitext(strInputAbundanceFile)[0]+"-checked.pcl"
          lsFilesToMake.append(strCheckedAbndFile)

        if strLabel and iSupervisedSampleCount:
          if not strInputPredictFile:
            strInputPredictFile = "".join([os.path.splitext(strCheckedAbndFile)[0],"-select-",str(iSupervisedSampleCount),"-",str(strLabel),"-SVM.txt"])
            lsFilesToMake.append(strInputPredictFile)
          if not strPredictPredictFile:
            strPredictPredictFile = "".join([os.path.splitext(strCheckedAbndFile)[0],"-select-",str(iSupervisedSampleCount),"-",str(strLabel),"-SVM.predict"])
            lsFilesToMake.append(strPredictPredictFile)

        #If the directories do not already exist, create
        for f in [strInputAbundanceFile,strOutputFile] + lsFilesToMake:
            strDir = os.path.dirname(f)
            if not os.path.exists(strDir):
                os.makedirs(strDir)

        #microPITA object
        microPITA = MicroPITA()

        #SVM parameters
        #Constants associated with the abundance to SVM input file conversion
        c_ABUNDANCE_DELIMITER=ConstantsMicropita.c_cTab
        c_NORMALIZE_RELATIVE_ABUNDANCY=True
        c_SKIP_FIRST_COLUMN=True

        #Diversity metrics to run
        diversityMetricsAlpha = [ConstantsMicropita.c_strInverseSimpsonDiversity]
        diversityMetricsBeta = [ConstantsMicropita.c_strBrayCurtisDissimilarity]
        inverseDiversityMetricsBeta = [ConstantsMicropita.c_strInvBrayCurtisDissimilarity]
        diversityMetricsAlphaNoNormalize = []#ConstantsMicropita.c_strChao1Diversity]
        diversityMetricsBetaNoNormalize = []
        inverseDiversityMetricsBetaNoNormalize = []

        #Targeted taxa
        userDefinedTaxa = []

        #Perform different flows flags
        c_RUN_MAX_DIVERSITY_1 = False
        if ConstantsMicropita.c_strDiversity in strSelectionTechnique and iSampleSelectionCount:
            c_RUN_MAX_DIVERSITY_1 = True
        c_RUN_REPRESENTIVE_DISSIMILARITY_2 = False
        if ConstantsMicropita.c_strRepresentativeDissimilarity in strSelectionTechnique and iSampleSelectionCount:
            c_RUN_REPRESENTIVE_DISSIMILARITY_2 = True
        c_RUN_MAX_DISSIMILARITY_3 = False
        if ConstantsMicropita.c_strExtremeDissimilarity in strSelectionTechnique and iSampleSelectionCount:
            c_RUN_MAX_DISSIMILARITY_3 = True
        c_RUN_RANK_AVERAGE_USER_4 = False
        if ConstantsMicropita.c_strTaxa in strSelectionTechnique and iSampleSelectionCount:
            c_RUN_RANK_AVERAGE_USER_4 = True
            if not strUserDefinedTaxaFile:
                logging.error("MicroPITA.funcRun. No taxa file was given for taxa selection.") 
                return False
            if not os.path.exists(strUserDefinedTaxaFile):
                logging.error("MicroPITA.funcRun. The taxa file given for selection does not exist.") 
                return False
            if not sFeatureSelectionMethod:
                logging.error("MicroPITA.funcRun. No feature selection method was given for taxa selection.") 
                return False
            if c_RUN_RANK_AVERAGE_USER_4:
                #Read in taxa list, break down to lines and filter out empty strings
                with open(strUserDefinedTaxaFile,'r') as fhndlTaxaInput:
                    userDefinedTaxa = filter(None,fhndlTaxaInput.read().split(ConstantsMicropita.c_strEndline))
                if not sFeatureSelectionMethod:
                    sFeatureSelectionMethod = ConstantsMicropita.c_strTargetedRanked

        c_RUN_RANDOM_5 = False
        if ConstantsMicropita.c_strRandom in strSelectionTechnique:
            c_RUN_RANDOM_5 = True
        c_RUN_DISTINCT = False
        if (ConstantsMicropita.c_strDistinct in strSelectionTechnique) and strLabel and iSupervisedSampleCount:
            c_RUN_DISTINCT = True
        c_RUN_DISCRIMINANT = False
        if (ConstantsMicropita.c_strDiscriminant in strSelectionTechnique) and strLabel and iSupervisedSampleCount:
            c_RUN_DISCRIMINANT = True

        #Input file path components
        inputFileComponents = os.path.splitext(strInputAbundanceFile)
        inputFilePrefix = inputFileComponents[0]

        sampleSVMSelectionCount = iSupervisedSampleCount

        #Holds the alpha diversity metrics for samples
        internalAlphaMatrix = []
        #Holds the beta diversity metrics for samples
        internalBetaMatrix = dict()
        userRankedSamples = None
        randomlySelectedSamples = None

        #Check/reduce raw abundance data
        #If already normalized you cant run the occurence filter so make None to turn off.
        liOccurenceFilter = ConstantsMicropita.c_liOccurenceFilter
        if fIsAlreadyNormalized:
            liOccurenceFilter = None

        strInputAbundanceFile = AbundanceTable.funcCheckRawDataFile(strReadDataFileName=strInputAbundanceFile, sLastMetadataName=sLastMetadataName, lOccurenceFilter = liOccurenceFilter, strOutputFileName=strCheckedAbndFile)

        #Read in abundance data
        #Abundance is a structured array. Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
        #Abundance table object to read in and manage data
        totalAbundanceTable = AbundanceTable.funcMakeFromFile(strInputFile=strInputAbundanceFile, fIsNormalized=fIsAlreadyNormalized, fIsSummed=fCladesAreSummed,
                                   cDelimiter=cDelimiter, sMetadataID=sMetadataID, sLastMetadata=sLastMetadataName, cFeatureNameDelimiter=cFeatureNameDelimiter)

        if not totalAbundanceTable:
            logging.error("MicroPITA.funcRun. Could not read abundance table. Stopped.")
            return False

        if fSumData:
            totalAbundanceTable.funcSumClades()

###TODO REMOVE
#        totalAbundanceTable.funcReduceFeaturesToCladeLevel(6)   

        dictTotalMetadata = totalAbundanceTable.funcGetMetadataCopy()

        #Log metadata keys
        logging.debug(" ".join(["Micropita:funcRun.","Received metadata keys=",str(dictTotalMetadata.keys())]))

        #If there is only 1 unique value for the labels, do not run the Supervised methods
        if len(set(dictTotalMetadata.get(strLabel,[]))) < 2:
            c_RUN_DISCRIMINANT = False
            c_RUN_DISTINCT = False
            logging.debug("".join(["The label ",str(strLabel)," did not have 2 or more values. Labels found="]+dictTotalMetadata.get(strLabel,[])))

        logging.debug(" ".join(["Micropita:funcRun.","Received metadata=",str(dictTotalMetadata)]))

        #Run unsupervised methods###
        #Stratify the data if need be and drop the old data
        lStratifiedAbundanceTables = None
        if (not strStratify == None) and (not strStratify == "None"):
            lStratifiedAbundanceTables = totalAbundanceTable.funcStratifyByMetadata(strStratify,fWriteToFile="".join([os.path.split(strOutputFile)[0],"/"]))
        else:
            lStratifiedAbundanceTables = [totalAbundanceTable]

        #For each stratified abundance block or for the unstratfified abundance
        #Run the unsupervised blocks
        for stratAbundanceTable in lStratifiedAbundanceTables:

            #Check to make sure the stratification is not just NA values which can happen
            #If so skip them
            if (not strStratify == None) and (not strStratify == "None"):
                if len(set(stratAbundanceTable.funcGetMetadata(strStratify)) & set(ConstantsMicropita.lNAs)) > 0:
                    continue

            logging.info("Running abundance block:"+stratAbundanceTable.funcGetName())

            #Only perform if the data is not yet normalized
            if not stratAbundanceTable.funcIsNormalized():

                #Need to first work with unnormalized data
                if((c_RUN_MAX_DIVERSITY_1)or(c_RUN_REPRESENTIVE_DISSIMILARITY_2) or (c_RUN_MAX_DISSIMILARITY_3)):

                    self._funcRunNormalizeSensitiveMethods(abndData=stratAbundanceTable, iSampleSelectionCount=iSampleSelectionCount,
                                                     dictSelectedSamples=selectedSamples, lsAlphaMetrics=diversityMetricsAlphaNoNormalize,
                                                     lsBetaMetrics=diversityMetricsBetaNoNormalize,
                                                     lsInverseBetaMetrics=inverseDiversityMetricsBetaNoNormalize,
                                                     fRunDiversity=c_RUN_MAX_DIVERSITY_1,fRunRepresentative=c_RUN_REPRESENTIVE_DISSIMILARITY_2,
                                                     fRunExtreme=c_RUN_MAX_DISSIMILARITY_3)

                #Normalize data at this point
                fNormalizeSuccess = stratAbundanceTable.funcNormalize()
                if not fNormalizeSuccess:
                    logging.error("MicroPITA.funcRun. Error occured during normalizing data. Stopped.")
                    return False

            #Need to first work with unnormalized data
            if((c_RUN_MAX_DIVERSITY_1)or(c_RUN_REPRESENTIVE_DISSIMILARITY_2)or(c_RUN_MAX_DISSIMILARITY_3)):

                self._funcRunNormalizeSensitiveMethods(abndData=stratAbundanceTable, iSampleSelectionCount=iSampleSelectionCount,
                                                 dictSelectedSamples=selectedSamples, lsAlphaMetrics=diversityMetricsAlpha,
                                                 lsBetaMetrics=diversityMetricsBeta,
                                                 lsInverseBetaMetrics=inverseDiversityMetricsBeta,
                                                 fRunDiversity=c_RUN_MAX_DIVERSITY_1,fRunRepresentative=c_RUN_REPRESENTIVE_DISSIMILARITY_2,
                                                 fRunExtreme=c_RUN_MAX_DISSIMILARITY_3)

            #Write file For debugging
#TODO REMOVE            stratAbundanceTable.funcWriteToFile(os.path.split(strOutputFile)[0]+"-sumnorm.pcl")

            #Generate selection by the rank average of user defined taxa
            #Expects (Taxa (row) by Samples (column))
            #Expects a column 0 of taxa id that is skipped
            #Returns [(sample name,average,rank)]
            if(c_RUN_RANK_AVERAGE_USER_4):
              if not ConstantsMicropita.c_strUserRanked in selectedSamples:
                  selectedSamples[ConstantsMicropita.c_strUserRanked]=list()
              selectedSamples[ConstantsMicropita.c_strUserRanked].extend(microPITA.funcSelectTargetedTaxaSamples(abndMatrix=stratAbundanceTable, lsTargetedTaxa=userDefinedTaxa, iSampleSelectionCount=iSampleSelectionCount, sMethod=sFeatureSelectionMethod))
            logging.info("Selected Samples 4")
            logging.info(selectedSamples)

            #5::Select randomly
            #Expects sampleNames = List of sample names [name, name, name...]
            if(c_RUN_RANDOM_5):

                #Select randomly from sample names
                randomlySelectedSamples = microPITA.funcGetRandomSamples(lsSamples=stratAbundanceTable.funcGetSampleNames(), iNumberOfSamplesToReturn=iSampleSelectionCount)
                if not ConstantsMicropita.c_strRandom in selectedSamples:
                    selectedSamples[ConstantsMicropita.c_strRandom]=list()
                selectedSamples[ConstantsMicropita.c_strRandom].extend(list(randomlySelectedSamples))

            logging.info("Selected Samples 5")
            logging.info(selectedSamples)

        #Run supervised methods#
        lStratifiedAbundanceTables = None
        totalAbundanceTable.funcNormalize()

        ##Remove NA entries from the abundance table for the metadata label
        ##Warning this modifies the table itself and does NOT return a copy
        ##Any analysis after the point may be working with a subset of samples
        ##ALL UNSUPERVISED SELECTION SHOULD HAPPEN BEFORE HERE
        #This is valid for running supervised methods on 1 label but not multiple labels
        #Or adding additional downstream analysis after this point (unless it is contengent on the supervised label).
        if not strLabel is None:
            fRemoveSuccess = totalAbundanceTable.funcRemoveSamplesByMetadata(strLabel,ConstantsMicropita.lNAs)

            if fRemoveSuccess:
                if(c_RUN_DISTINCT or c_RUN_DISCRIMINANT):
                    selectedSamples.update(self.funcRunSupervisedMethods(abundanceTable=totalAbundanceTable,fRunDistinct=c_RUN_DISTINCT, fRunDiscriminant=c_RUN_DISCRIMINANT,
                                       strOuputSVMFile=strInputPredictFile,strPredictSVMFile=strPredictPredictFile,
                                       strSupervisedMetadata=strLabel, iSampleSVMSelectionCount=sampleSVMSelectionCount))
                    logging.info("Selected Samples Unsupervised")
                    logging.info(selectedSamples)
            else:
                logging.error("MicroPITA.funcRun. Error occured when cleaning up table for metadata label values. Supervised methods not ran.")
        elif(c_RUN_DISTINCT or c_RUN_DISCRIMINANT):
            logging.error("MicroPITA.funcRun. strLabel was not given so supervised methods were not ran.")
        return selectedSamples

    #Testing: Happy path tested
    @staticmethod
    def funcWriteSelectionToFile(dictSelection,strOutputFilePath):
	"""
	Writes the selection of samples by method to an output file.
	
	:param	dictSelection:	The dictionary of selections by method to be written to a file.
	:type	Dictionary	The dictionary of selections by method {"method":["sample selected","sample selected"...]}
	:param	strOutputFilePath:	String path to file to output dictionary.
	:type	String	String path to file to write to
	"""

        #Holds the output content
        strOutputContent = ""
        #Create output content from dictionary
        for sKey in dictSelection:
            strOutputContent = "".join([strOutputContent,sKey,ConstantsMicropita.c_cTab,ConstantsMicropita.c_cTab.join(dictSelection[sKey]),ConstantsMicropita.c_strEndline])

        #Write to file
        if(not strOutputContent == ""):
            with open(strOutputFilePath,'w') as fHndlOutput:
                fHndlOutput.write(str(strOutputContent))
        logging.debug("".join(["Selected samples output to file:",strOutputContent]))

    #Testing: Happy Path tested
    @staticmethod
    def funcReadSelectionFileToDictionary(strInputFile):
	"""
	Reads in an output selection file from micropita and formats it into a dictionary.
	
	:param	strInputFile:	String path to file to read and translate into a dictionary.
                                {"method":["sample selected","sample selected"...]}
	:type	String	String path to file to read and translate.
	:return    Dictionary:    Samples selected by methods.
                Dictionary    {"Selection Method":["SampleID","SampleID","SampleID",...]}
	"""

        #Read in selection file
        strSelection = ""

        #Check for file
        if not os.path.exists(strInputFile):
            return False

        with open(strInputFile,'r') as fHndlInput:
            strSelection = fHndlInput.read()

        #Dictionary to hold selection data
        dictSelection = dict()
        for strSelectionLine in filter(None,strSelection.split(ConstantsMicropita.c_strEndline)):
            astrSelectionMethod = strSelectionLine.split(ConstantsMicropita.c_cTab)
            dictSelection[astrSelectionMethod[0]] = astrSelectionMethod[1:]

        #Return dictionary
        return dictSelection

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicroPITA.py", 
    description = """Selects samples from abundance tables based on various selection schemes.""" )

#Arguments
#Optional parameters
#Logging
argp.add_argument(ConstantsMicropita.c_strLoggingArgument, dest="strLogLevel", metavar= "Log level", default="INFO", 
                  choices=ConstantsMicropita.c_lsLoggingChoices, help= ConstantsMicropita.c_strLoggingHelp)

#Abundance associated
argp.add_argument(ConstantsMicropita.c_strIDNameArgument, dest="sIDName", metavar= "Sample ID Metadata Name", default="ID", help= ConstantsMicropita.c_strIDNameHelp)
argp.add_argument(ConstantsMicropita.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "Last Metadata Name",
                  help= ConstantsMicropita.c_strLastMetadataNameHelp)
argp.add_argument(ConstantsMicropita.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store_true", default=False,
                  help= ConstantsMicropita.c_strIsNormalizedHelp)
argp.add_argument(ConstantsMicropita.c_strIsSummedArgument, dest="fIsSummed", action = "store_true", default= False, help= ConstantsMicropita.c_strIsSummedHelp)
argp.add_argument(ConstantsMicropita.c_strSumDataArgument, dest="fSumData", action = "store_false", default = True, help= ConstantsMicropita.c_strSumDataHelp)
argp.add_argument(ConstantsMicropita.c_strTargetedFeatureMethodArgument, dest="sFeatureSelection", metavar= "Feature Selection Method", default=ConstantsMicropita.lsTargetedFeatureMethodValues[0], 
                  choices=ConstantsMicropita.lsTargetedFeatureMethodValues, help= ConstantsMicropita.c_strTargetedFeatureMethodHelp)
argp.add_argument(ConstantsMicropita.c_strFileDelimiterArgument, dest= "cFileDelimiter", action= "store", metavar="File Delimiter", default=ConstantsMicropita.c_cTab, help=ConstantsMicropita.c_strFileDelimiterHelp) 
argp.add_argument(ConstantsMicropita.c_strFeatureNameDelimiterArgument, dest= "cFeatureNameDelimiter", action= "store", metavar="Feature Name Delimiter", default=ConstantsMicropita.c_cPipe, help=ConstantsMicropita.c_strFeatureNameDelimiterHelp) 

#Select count
argp.add_argument(ConstantsMicropita.c_strUnsupervisedCountArgument, dest="iUnsupervisedSelectionCount", metavar = "Number Samples To Select (Unsupervised)", default=0, type = int, help = ConstantsMicropita.c_strUnsupevisedCountHelp)
argp.add_argument(ConstantsMicropita.c_strTargetedSelectionFileArgument, dest="strFileTaxa", metavar = "Targeted Feature File", default=None, help = ConstantsMicropita.c_strTargetedSelectionFileHelp)
argp.add_argument(ConstantsMicropita.c_strUnsupervisedStratifyMetadataArgument, dest="sUnsupervisedStratify", metavar= "Metadata to Stratify Unsupervised Selection", default=None, 
                  help= ConstantsMicropita.c_strUnsupervisedStratifyMetadataHelp)

#SVM label
#Label parameter to be used with SVM
argp.add_argument(ConstantsMicropita.c_strSupervisedLabelArgument, dest="sLabel", metavar= "Supervised Label Metadata Name", default=None, help= ConstantsMicropita.c_strSupervisedLabelCountHelp)
argp.add_argument(ConstantsMicropita.c_strSupervisedLabelCountArgument, dest="iSupervisedCount", metavar= "Supervised Sample Selection Count", default=0, type=int,
                  help= ConstantsMicropita.c_strSupervisedLabelCountHelp)

#Files
argp.add_argument(ConstantsMicropita.c_strCheckedAbundanceFileArgument, dest="strCheckedFile", metavar = "Checked Abundance File Path", default="", help = ConstantsMicropita.c_strCheckedAbundanceFileHelp)
argp.add_argument(ConstantsMicropita.c_strLoggingFileArgument, dest="strLoggingFile", metavar = "Logging File Path.", default="", help = ConstantsMicropita.c_strLoggingFileHelp)
argp.add_argument(ConstantsMicropita.c_strSupervisedInputFile, dest="strInputPredictFile", metavar = "Output File Path of the Scaled Values for Supervised Predictions.", default="", help = ConstantsMicropita.c_strSupervisedInputFileHelp)
argp.add_argument(ConstantsMicropita.c_strSupervisedPredictedFile, dest="strPredictFile", metavar = "Output File Path of the Supervised Predictions.", default="", help = ConstantsMicropita.c_strSupervisedPredictedFileHelp)

#Required
#Input
#Abundance file
argp.add_argument("strFileAbund", metavar = "Abundance file", help = ConstantsMicropita.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument("strOutFile", metavar = "Selection Output File", help = ConstantsMicropita.c_strGenericOutputDataFileHelp)
#Selection parameter
argp.add_argument("strSelection", metavar = "Selection Methods", help = ConstantsMicropita.c_strSelectionTechniquesHelp, nargs="*")


__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    if not args.strLoggingFile:
        args.strLoggingFile = "".join([os.path.splitext(args.strOutFile)[0],".log"])

    #If the logging directory does not already exist, create
    for f in [args.strLoggingFile]:
        strDir = os.path.dirname(f)
        if not os.path.exists(strDir):
            os.makedirs(strDir)
    logging.basicConfig(filename = args.strLoggingFile, filemode = 'w', level=iLogLevel)

    #Run micropita
    logging.info("Start microPITA")
    microPITA = MicroPITA()

    dictSelectedSamples = microPITA.funcRun(fIsAlreadyNormalized=args.fIsNormalized,
                                        fCladesAreSummed=args.fIsSummed,
                                        sMetadataID=args.sIDName,
                                        sLastMetadataName=args.sLastMetadataName,
                                        strInputAbundanceFile=args.strFileAbund,
                                        strInputPredictFile=args.strInputPredictFile,
                                        strPredictPredictFile=args.strPredictFile,
                                        strCheckedAbndFile=args.strCheckedFile,
                                        strOutputFile=args.strOutFile,
                                        cDelimiter=args.cFileDelimiter,
                                        cFeatureNameDelimiter = args.cFeatureNameDelimiter,
                                        strUserDefinedTaxaFile=args.strFileTaxa,
                                        iSampleSelectionCount=args.iUnsupervisedSelectionCount,
                                        iSupervisedSampleCount=args.iSupervisedCount,
                                        strSelectionTechnique=args.strSelection,
                                        strLabel=args.sLabel,
                                        strStratify=args.sUnsupervisedStratify,
                                        fSumData=args.fSumData,
                                        sFeatureSelectionMethod=args.sFeatureSelection)

    if dictSelectedSamples == False:
        logging.error("".join(["MicroPITA::Error, did not get a result from analysis."]))
        return -1
    logging.info("End microPITA")

    #Log output for debugging
    logging.debug("".join(["Returned the following samples:",str(dictSelectedSamples)]))

    #Write selection to file
    microPITA.funcWriteSelectionToFile(dictSelection=dictSelectedSamples,strOutputFilePath=args.strOutFile)

if __name__ == "__main__":
    _main( )

