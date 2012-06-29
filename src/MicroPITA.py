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

#External libraries
from AbundanceTable import AbundanceTable
import argparse
from Constants import Constants
from Constants_Arguments import Constants_Arguments
import csv
from Diversity import Diversity
import itertools
from KMedoids import Kmedoids
from LIBSVM import LIBSVM
import logging
import math
import mlpy
from MLPYDistanceAdaptor import MLPYDistanceAdaptor
import numpy as np
import operator
import os
from Pycluster import *
import random
import re
import scipy.cluster.hierarchy as hcluster
from SVM import SVM
import sys
from Utility_Math import Utility_Math
from ValidateData import ValidateData

class MicroPITA:
    """
    Performs analysis associated with sample selection.
    """

    #Constants
    #Diversity metrics Alpha
    c_strInverseSimpsonDiversity = Diversity.c_strInvSimpsonDiversity
    c_strChao1Diversity = Diversity.c_strChao1Diversity

    #Diversity metrics Beta
    c_strBrayCurtisDissimilarity = Diversity.c_strBrayCurtisDissimilarity

    #Additive inverses of diversity metrics beta
    c_strInvBrayCurtisDissimilarity = Diversity.c_strInvBrayCurtisDissimilarity

    #Selelection methods
    c_strDiversity = "Diversity"
    c_strExtremeDissimilarity = "Extreme"
    c_strDiscriminant = "Discriminant"
    c_strDistinct = "Distinct"
    c_strRandom = "Random"
    c_strRepresentativeDissimilarity = "Representative"
    c_strTaxa = "Taxa_Defined"
    c_lsAllSupervisedMethods = [c_strDiversity,c_strExtremeDissimilarity,c_strRandom,c_strRepresentativeDissimilarity,c_strTaxa]
    c_lsAllUnsupervisedMethods = [c_strDiscriminant,c_strDistinct]

    #Technique Names
    c_strDiversity1 = "".join([c_strDiversity,"_I"])
    c_strDiversity2 = "".join([c_strDiversity,"_C"])
    c_strExtremeDissimiarity1 = "".join([c_strExtremeDissimilarity,"_B"])
    c_strRepresentativeDissimilarity1 = "".join([c_strRepresentativeDissimilarity,"_B"])
    c_strUserRanked = c_strTaxa
    c_strSVMClose = c_strDiscriminant
    c_strSVMFar = c_strDistinct

    #Targeted feature settings
    c_strTargetedRanked = Constants_Arguments.c_strTargetedRanked
    c_strTargetedAbundance = Constants_Arguments.c_strTargetedAbundance

    #Technique groupings
    c_lsDiversityMethods = [c_strDiversity1,c_strDiversity2]

    #Converts ecology metrics into standardized method selection names
    dictConvertAMetricDiversity = {c_strInverseSimpsonDiversity:c_strDiversity1, c_strChao1Diversity:c_strDiversity2}
    dictConvertMicroPITAToAMetric = {c_strDiversity1:c_strInverseSimpsonDiversity, c_strDiversity2:c_strChao1Diversity}
    dictConvertBMetricRepresentative = {c_strBrayCurtisDissimilarity:c_strRepresentativeDissimilarity1}
    dictConvertBMetricExtreme = {c_strInvBrayCurtisDissimilarity:c_strExtremeDissimiarity1}

    #Linkage used in the Hierarchical clustering
    c_strHierarchicalClusterMethod = 'average'

####Group 1## Diversity

    #Testing: Happy path Testing (8)
    def funcGetTopRankedSamples(self, lldMatrix = None, lsSampleNames = None, iTopAmount = None):
	"""
	Given a list of lists of measurements, for each list the indices of the highest values are returned. If lsSamplesNames is given
        it is treated as a list of string names that is in the order of the measurements in each list. Indices are returned or the sample
        names positionally associated with the indices.
	
	:param	lldMatrix:	List of lists [[value,value,value,value],[value,value,value,value]].
	:type	List of lists	List of measurements. Each list is a different measurement. Each measurement in possionally related to a sample.
	:param	lsSampleNames:	List of sample names positionally related (the same) to each list (Optional).
	:type	List of strings	List of strings.
	:param	iTopAmount:	The amount of top measured samples (assumes the higher measurements are better).
	:type	integer	Integer amount of sample names/ indices to return.
        :return	List:	List of samples to be selected.
        :type	List	List of strings.
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
	
	:param	npadAbundancies:	Numpy arrayof sample abundances to measure against.
	:type	Numpy Array	Numpy array where row=samples and columns = features.
	:param	sMetric:	String name of beta metric. Possibilities are listed in microPITA.
	:type	String	String name of beta metric. Possibilities are listed in microPITA.
        :return	Double:	Get Measurement indicated by metric for given abundancies
        :type	Double	Measurement
	"""

        if(sMetric == self.c_strBrayCurtisDissimilarity):
            return Diversity.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        elif(sMetric == self.c_strInvBrayCurtisDissimilarity):
            return Diversity.funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        else:
            return False

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
        :type	List	List of strings, returns false on error.
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
        if(ValidateData.funcIsFalse(distanceMatrix)):
          logging.error("MicroPITA.funcGetCentralSamplesByKMedoids. Received false for betaMetrix matrix generation, returning false.")
          return False

        #Log distance matrix
        logging.debug("".join(["Distance matrix for representative selection using metric=",str(sMetric)]))

#        #Run MLPY implementation
#        if not fMLPY:
#            #Perform Kmedoid clustering in pycluster
#            lclusterid, error, nfound = kmedoids (distance=distanceMatrix, nclusters=iNumberSamplesReturned, npass=1000, initialid=None)
#            #TODO make sure you are getting condensed for unifrac and braycurtis
#            logging.debug("Results from the kmedoid method in representative selection:")
#            logging.debug("clusterid:"+str(list(lclusterid)))
#            logging.debug("error:"+str(error))
#            logging.debug("nfound:"+str(nfound))
#
#            #Convert centroid indexes to names and return
#            #Return centroids
#            return [lsSampleNames[iindex] for iindex in set(lclusterid)]

#        #Run Pycluster implementation
#        else:
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
        :type	List of strings	Selected samples.
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
        linkageMatrix = hcluster.linkage(tempDistanceMatrix, method=self.c_strHierarchicalClusterMethod)

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
            if(iNode1<iSampleCount):
                returnSamples.append(lsSampleNames[iNode2])
                iCurrentSelectCount = iCurrentSelectCount + 1
            if iCurrentSelectCount == iSelectSampleCount:
                break

        #Return selected samples
        return returnSamples

####Group 4## Rank Average of user Defined Taxa
    #Happy Path Tested
    def funcGetAverageAbundanceSamples(self, abndTable, lsTargetedFeature, fRank=False):
	"""
	Averages feature abundance or ranked abundance. Expects a column 0 of taxa id that is skipped.
	
	:param	abndTable:	Abundance Table to analyse
	:type	AbundanceTable	Abundance Table
	:param	lsTargetedFeature:	String names
	:type	list	list of string names of features (bugs) which are measured after ranking against the full sample
        :param  fRank:	Indicates to rank the abundance before getting the average abundance of the features (default false)
        :type   boolean	Flag indicating ranking abundance before calculating average feature measurement (false= no ranking)
        :return	List of lists or boolean:	List of lists or False on error. One internal list per sample indicating the sample, feature average abundance or ranked abundance and somethign to break ties
						Lists will already be sorted.
        :type	list	For not Ranked [[sample,average abundance of selected feature,1]]
        		For Ranked [[sample,average ranked abundance, average abundance of selected feature]]
			Error Returns false
	"""

        llAbundance = abndTable.funcGetAverageAbundancePerSample(lsTargetedFeature)
        if not llAbundance:
            logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not get average abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table."]))

        #Add a space for ranking if needed
        #Not ranked will be [[sSample,average abundance,1]]
        #(where 1 will not discriminant ties if used in later functions, so this generalizes)
        #Ranked will be [[sSample, average rank, average abundance]]
        #Here average abundance can be used incase of tie breaking.
        llRetAbundance = [llist+[1] for llist in llAbundance]

        #Rank if needed
        if fRank:
            abndRanked = abndTable.funcRankAbundance()
            if abndTable == None:
                logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not rank the abundance table, returned false."]))
                return False
            dictRankedAbundance = dict([tuple(lRankedItem) for lRankedItem in abndRanked.funcGetAverageAbundancePerSample(lsTargetedFeature)])
            if not dictRankedAbundance:
                logging.error("".join(["MicroPITA.funcGetAverageAbundanceSamples. Could not get average ranked abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table."]))
                return False

            #Check for the case where all features are 0
            #Rank could be 1 and these could be selected
            #Make them the worst ranked by giving them feature count + feature count-rank
            #The ideas that if the sample is choosen then the samples with the most occurences of other features
            #Should be the ones with the 0 features having the highest ranks (assuming minimal tieing)
            #So this amount is subtracted from the total featurecount and added to featurecount
            #To make the 0 feature samples at the end of the list but ranked by how much abundance are in the
            #OTHER features.....This could be done by ranking 0s by total sample abundance
            #Also add in the rank information as the second element and move the abundance information to the 3rd element
            iWorstRank = float(abndRanked.funcGetFeatureCount())
            llRetAbundance = [[llist[0],iWorstRank+iWorstRank-dictRankedAbundance[llist[0]],llist[1]] if llist[1] == 0.0 
                             else [llist[0],dictRankedAbundance[llist[0]],llist[1]]
                             for llist in llRetAbundance]

        #Sort first for ties and then for the main feature
        if Constants.fBreakRankTiesByDiversity:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[2], reverse = True)
        if fRank:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[1], reverse = False)
        else:
            llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[1], reverse = True)

        #return
        return llRetAbundance

    #Happy Path Tested
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
      :type	List of strings	Empty list is returned on an error.
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
      if sMethod.lower() == self.c_strTargetedRanked.lower():
          lsTargetedSamples = self.funcGetAverageAbundanceSamples(abndTable=abndMatrix, lsTargetedFeature=lsTargetedTaxa, fRank=True)
      elif sMethod.lower() == self.c_strTargetedAbundance.lower():
          lsTargetedSamples = self.funcGetAverageAbundanceSamples(abndTable=abndMatrix, lsTargetedFeature=lsTargetedTaxa, fRank=False)

      #If an error occured or the key word for the method was not recognized
      if lsTargetedSamples == False: 
          logging.error("MicroPITA::funcSelectTargetedTaxaSamples: Was not able to select for the features given. So targeted feature selection was performed. Check to make sure the features are spelled correctly and exist in the abundance file.")
          return []

      #Select from results
      return [sSample[0] for sSample in lsTargetedSamples[:iSampleSelectionCount]]


####Group 5## Random

    #Happy path Tested
    def funcGetRandomSamples(self, lsSamples=None, iNumberOfSamplesToReturn=0):
	"""
	Returns random sample names of the number given. No replacement.
	
	:param	lsSamples:	List of sample names 
	:type	list	list of strings
	:param	iNumberOfSamplesToReturn:	Number of samples to select
	:type	integer	integer.
        :return	List:	List of selected samples.
        :type	List	List of strings.
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
    #@return Dictionary of file paths which were generated by the model and prediction steps
    def funcRunLIBSVM(self, abndTable, strOutputSVMFile, sLabelID, strTempDir):
	"""
	Runs a linear SVM (Adapted from the easy.py script included in the standard libsvm install).

	:param	abndTable:	AbundanceTable of data.
	:type	AbundanceTable	Data.
	:param	strOutputSVMFile:	Out path used for basing ouput files generated by LIBSVM.
	:type	String	List of strings
	:param	sLabelID:	Metadata data ID for the labels.
	:type	String	Metadta ID
        :param	strTempDir:	String directroy path for temporary files.
        :type	String
        :return	Files:	Dictionary of generated files, double cost value used, and double best accuracy.
	"""

        #If the output should be probabilistic
        strSVMProbabilistic = Constants.c_fProbabilitistic

        #The cost range
        lSVMLogC = Constants.c_lCostRange

        #The lower bound on scaling
        iSVMScaleLowestBound = Constants.c_intScaleLowerBound

        #Create SVM object
        svm = LIBSVM()

        #Holds files generated by SVM code
        #Files generated by modeling
        modelFiles = None
        #Files generated by prediction
        predictionFiles = None

        #Convert abundancies file to SVM file
        noError = svm.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abndTable, strOutputSVMFile=strOutputSVMFile, sMetadataLabel=sLabelID)
        modelFiles = False
        predictionFiles = False
        #Run SVM
        if(not noError==False):
            modelFiles = svm.createLinearModel(tempInputFileName=strOutputSVMFile, tempTMPDirectory=strTempDir, tempScaling=iSVMScaleLowestBound, tempLogC=lSVMLogC, tempProbabilistic=strSVMProbabilistic)
            logging.debug("modelFiles")
            logging.debug(str(modelFiles))
            if(not modelFiles == False):
                predictionFiles = svm.predictFromLinearModel(tempDataFileName=strOutputSVMFile, tempModelFileName=modelFiles[Constants.c_strKeywordModelFile], tempRangeFileName=modelFiles[Constants.c_strKeywordRangeFile], tempProbabilistic=strSVMProbabilistic)
            else:
                logging.error("".join(["MicroPITA.funcRunLIBSVM: Could not run prediction, model file was false."]))
        #Combine output dictionarys and return
        if(not predictionFiles == False):
            for pKey in predictionFiles:
                modelFiles[pKey]=(predictionFiles[pKey])
        else:
            logging.error("MicroPITA.funcRunLIBSVM: Prediction files were false.")
        return modelFiles

    #@return Dictionary of file paths which were generated by the model and prediction steps
    def funcRunMLPYSVM(self, abndAbundanceTable, sMetadataForLabel, strInputSVMFile, strPredictionFile, fPredictWithFullData = True, fClassifyByProbability = True):
	"""
	Runs a linear SVM using MLPY.
	
	:param	abndAbundanceTable:	Abundance table of data.
	:type	Abundance table	Abundance table of data.
	:param	sMetadataForLabel:	Metadata for label used to supervised learning.
	:type	string	string
	:param	strOutputSVMFile:	File name for the file generated to mock the input SVM file.
	:type	string	string
        :param	strPredictionFile	LIBSVM style output file of probabilistic predictions.
        :type	string	File path
        :return	Analysis Results:	Dictionary of files generated (input and prediction files)
                                        Return false on error.
	"""

        #SVM general class
        svm = SVM()

        #Convert abundancies file to SVM file
        noError = svm.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abndAbundanceTable, strOutputSVMFile=strInputSVMFile, sMetadataLabel=sMetadataForLabel)
        if not noError:
            logging.error("MicroPITA.funcRunMLPYSVM: Received an error when creating the input SVM file in the MLPY LIBSVM analysis pipeline.")
            return False

        #Will hold the SVM related return data including files generated during the process.
        #Although these files are not necessary for MLPY Libsvm, they are made to mirror the LIBSVM process and
        #and can be used to troubleshoot.
        dictSVMReturn = {}
        dictSVMReturn[Constants.c_strKeywordInputFile] = strInputSVMFile

        #Cost values to optimize to
        liCost = Constants.c_lCostRange

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
        #middle and not an extreme of the cost value plateau created by the tieing
        fTiesProgress = True

        #For each cost
        for iCost in liCost:

#            print "iCost:", iCost
#            print "Cost Value used: ", math.pow(2,int(iCost))

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
                                                            dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)
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

        if fPredictWithFullData:
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

            print "lPredictions", lPredictions
            print "lSVMLabels", lSVMLabels
            print "npaDistances", npaDistances
            #Store all sample data
            ldictReturn = self._funcStoreSVMProbability(lsSampleNames,ldLabels,lSVMLabels,npaDistances,lPredictions,
                                                        dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability)

            dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions = ldictReturn

        #Log best
        logging.debug("".join(["funcRunMLPYSVM::Best Accuracy=", str(dBestAccuracy)]))
        logging.debug("".join(["funcRunMLPYSVM::Best Cost=", str(iBestCost)]))

        #Create output prediction file
        strPredictionOutput = " ".join(["labels"]+[str(iLabel) for iLabel in lSVMLabels])+Constants.ENDLINE
        strPredictionOutput = strPredictionOutput + Constants.ENDLINE.join([" ".join([str(dProb) for dProb in dictAllProbabilities[sSampleName]])
            for sSampleName in lsSampleNames])

        #Write prediction file to file
        with open(strPredictionFile, 'w') as f:
            f.write(strPredictionOutput)
        dictSVMReturn[Constants.c_strKeywordPredFile] = strPredictionFile

        #Return
        return dictSVMReturn

    #Happy path tested
    def _funcStoreSVMProbability(self,lsValidationSamples,ldValidationLabels,lSVMLabels,npaDistances,lPredictions,dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions,fClassifyByProbability):
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
        :param	fClassifyByProbability: Indicates if the distance from the hyperplane is measured from the probabilities or the predictions
        :type	Boolean	True = Probabilities, False = Predictions
        :return	Analysis Results:	[dictIBestPredictions (dictionary of predictions),
                                         dictDBestProbabilities (dictionary of probabilities),
                                         dictLabels (list of original labels),
                                         lsLabels (list of the labels in order of the MLPY output, these are unique values]
        :type	List	[dictIBestPredictions,dictDBestProbabilities,dictLabels,lsLabels]
        """

        #Store fold results
        for indexSamples, sSampleName in enumerate(lsValidationSamples):
            #Determine label by highest probability
            if fClassifyByProbability:
                dMaxProbability = max(npaDistances[indexSamples])
                dictdProbability[sSampleName] = dMaxProbability
                iLabel = lSVMLabels[list(npaDistances[indexSamples]).index(dMaxProbability)]
                dictiPrediction[sSampleName] = ldValidationLabels[indexSamples] == iLabel
                dictAllProbabilities[sSampleName] = [int(iLabel)]+list(npaDistances[indexSamples])
                dictAllPredictions[sSampleName] = str(int(iLabel))
            #Determine label by pred method
            else:
                dictiPrediction[sSampleName] = ldValidationLabels[indexSamples] == lPredictions[indexSamples]
                dictAllProbabilities[sSampleName] = [int(lPredictions[indexSamples])]+list(npaDistances[indexSamples])
                dictAllPredictions[sSampleName] = str(int(lPredictions[indexSamples]))
                dictdProbability[sSampleName] = npaDistances[indexSamples][lSVMLabels.index(lPredictions[indexSamples])]

        return [dictdProbability,dictAllProbabilities,dictiPrediction,dictAllPredictions]

#TODO Normalize? When is this done?
    #Run the supervised methods
    def funcRunSupervisedMethods(self, abundanceTable, fRunDistinct, fRunDiscriminant,
                                   strOuputSVMFile, strTMPDirectory, strSupervisedMetadata, iSampleSVMSelectionCount, fLibSVM=False):
        """
	Runs supervised methods.
	
	:param	abundanceTable:	AbundanceTable
	:type	AbudanceTable	Data to analyze
	:param	fRunDistinct:	Run distinct selection method
	:type	Boolean	boolean (true runs method)
	:param	fRunDiscriminant:	Run discriminant method
	:type	Boolean	boolean (true runs method)
	:param	strOutputSVMFile:	File output from  SVM
	:type	String	String
	:param	strTMPDirectory:	Directory to place temporary files
	:type	String	String
	:param	strSupervisedMetadata:	Label for supervised selection
	:type	String	String
	:param	iSampleSVMSelectionCount:	Number of samples to select
	:type	Integer	int sample selection count
        :return	Selected Samples:	A dictionary of selected samples by selection ID
        :type	Dictionary	{"Selection Method":["SampleID","SampleID"...]}
	"""

        #Run supervised blocks
        #Select supervised (using SVM)
        #Will contain the samples selected to return
        dictSelectedSamples = dict()
        #SVM related output from the SVM call
        svmRelatedData = None

        #Remove all files associated with supervised methods
        #Create file names and delete if they exist
        fileNamePrefix = os.path.splitext(strOuputSVMFile)[0]
        fileDirectory = os.path.split(strOuputSVMFile)[0]
        scaledFile,rangeFile,cvOutFile,modelFile,scaledPredictFile,predictFile = ["".join([fileNamePrefix,strExtension]) for strExtension in
                   [Constants.c_PREDICT_FILE_EXT,Constants.c_SCALING_PARAMETERS,Constants.c_CV_FILE_EXT,Constants.c_MODEL_FILE_EXT,
                    Constants.c_SCALED_FOR_PREDICTION_FILE_EXT,Constants.c_PREDICT_FILE_EXT]]
        for f in [scaledFile,rangeFile,cvOutFile,modelFile,scaledPredictFile,predictFile]:
            if os.path.exists(f):
                os.remove(f)

        #Run either the MLPY or our wrappers for LIBSVM
        if fLibSVM:
            logging.debug("".join(["Using LibSVM Support Vector machine package"]))

            #Run linear SVM
            svmRelatedData = self.funcRunLIBSVM(abndTable=abundanceTable, strOutputSVMFile=strOuputSVMFile, sLabelID=strSupervisedMetadata, strTempDir=strTMPDirectory)

        else:
            logging.debug("".join(["Using MLPY Support Vector machine package"]))

            #Run MLPY SVM
            svmRelatedData = self.funcRunMLPYSVM(abndAbundanceTable=abundanceTable, sMetadataForLabel=strSupervisedMetadata,
                                             strInputSVMFile=strOuputSVMFile, strPredictionFile=predictFile)

        #Read in prediction file and select samples
        if svmRelatedData:
            dictSelectedSamples = self._funcSelectSupervisedSamplesFromPredictFile(strOriginalInputFile=svmRelatedData[Constants.c_strKeywordInputFile],
                                                      strPredictFilePath=svmRelatedData[Constants.c_strKeywordPredFile], lsSampleNames=abundanceTable.funcGetSampleNames(),
                                                      iSelectCount=iSampleSVMSelectionCount,
                                                      fSelectDiscriminant = fRunDiscriminant, fSelectDistinct = fRunDistinct)
        return dictSelectedSamples

    def _funcSelectSupervisedSamplesFromPredictFile(self, strOriginalInputFile, strPredictFilePath, lsSampleNames, iSelectCount, fSelectDiscriminant, fSelectDistinct):
        """
        Selects sample for the supervised distinct and discriminant methods from an SVM predict file.

        :param	strOriginalInputFile:	The original input file to parse the original class labels.
        :type	String
        :param	strPredictFilePath:	A file path to the predict file generated by the SVM methods.
        :type	String
        :param	lsSampleNames:	List of string sample ids.
        :type	String
        :param	iSelectCount:	Amount os samples to select per label.
        :type	Integer
        :param	fSelectDiscriminant:	Indicates samples should be selected for the disriminant method.
        :type	Boolean			
        :param	fSelectDistinct:	Indicates samples should be selected for the distinct method.
        :type	Boolean			
        """

        #Create and array to hold difference of the samples probabilities from the central probability
        dictDistanceFromHyperplane = dict()

        #Holds selected ids of selected samples
        dictSelectedSamples = dict()

        #Holds labels to compare to the predictions
        lsOriginalLabels = None
        #Open perdiction file and input file and get labels to compare to the predictions
        with open(strOriginalInputFile,'r') as f, open(strPredictFilePath,'r') as g:
            reader = csv.reader(f, delimiter=Constants.WHITE_SPACE, quoting=csv.QUOTE_NONE)
            lsOriginalLabels = [row[0] for row in reader]
            predictionLists = g.read()
            predictionLists = [filter(None,strPredictionList) for strPredictionList in predictionLists.split(Constants.ENDLINE)]

        #Get label count (meaning the number of label categories)(-1 to not count the predicted first entry)
        labelCount = len(predictionLists[0].split(Constants.WHITE_SPACE))-1

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
            lineElements = lineElements.split(Constants.WHITE_SPACE)
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
            dictSelectedSamples[self.c_strSVMClose]=SVMSamples
            logging.debug("SVMSamples Close")
            logging.debug(SVMSamples)
        if fSelectDistinct:
            SVMSamples = list()
            for selectedSampleIndex in selectedSamplesIndicesFar:
                SVMSamples.append(lsSampleNames[selectedSampleIndex])
            dictSelectedSamples[self.c_strSVMFar]=SVMSamples
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
            #Get Alpha metrics matrix
            #Expects Observations (Taxa (row) x sample (column))
            #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
            internalAlphaMatrix = Diversity.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abndData.funcGetAbundanceCopy(), lsSampleNames = lsSampleNames, lsDiversityMetricAlpha = lsAlphaMetrics)
            #Get top ranked alpha diversity by most diverse
            #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
            #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
            mostDiverseAlphaSamplesIndexes = self.funcGetTopRankedSamples(lldMatrix=internalAlphaMatrix, lsSampleNames=lsSampleNames, iTopAmount=iSampleSelectionCount)

            #Add to results
            for index in xrange(0,len(lsAlphaMetrics)):
                astrSelectionMethod = self.dictConvertAMetricDiversity[lsAlphaMetrics[index]]
                if not astrSelectionMethod in dictSelectedSamples:
                    dictSelectedSamples[astrSelectionMethod]=list()
                dictSelectedSamples[astrSelectionMethod].extend(mostDiverseAlphaSamplesIndexes[index])

        logging.info("Selected Samples 1b")
        logging.info(dictSelectedSamples)

        #Generate beta metrics and 
        if((fRunRepresentative)or(fRunExtreme)):

            #Abundance matrix transposed
            npaTransposedAbundance = Utility_Math.funcTransposeDataMatrix(abndData.funcGetAbundanceCopy(), fRemoveAdornments=True)

            #Get center selection using clusters/tiling
            #This will be for beta metrics in normalized space
            if(fRunRepresentative):
                logging.info("Performing representative selection on normalized data.")
                for bMetric in lsBetaMetrics:

                    #Get representative dissimilarity samples
                    medoidSamples=self.funcGetCentralSamplesByKMedoids(npaMatrix=npaTransposedAbundance, sMetric=bMetric, lsSampleNames=lsSampleNames, iNumberSamplesReturned=iSampleSelectionCount)

                    if(not medoidSamples == False):
                        astrSelectionMethod = self.dictConvertBMetricRepresentative[bMetric]
                        if not astrSelectionMethod in dictSelectedSamples:
                            dictSelectedSamples[astrSelectionMethod]=list()
                        dictSelectedSamples[astrSelectionMethod].extend(medoidSamples)

            #Get extreme selection using clusters, tiling
            if(fRunExtreme):
                logging.info("Performing extreme selection on normalized data.")
                #Run KMedoids with inverse custom distance metric in normalized space
                for bMetric in lsInverseBetaMetrics:

                    #Samples for repersentative dissimilaritystrUserDefinedTaxaFile
                    #This involves inverting the distance metric,
                    #Taking the dendrogram level of where the number cluster == the number of samples to select
                    #Returning a repersentative sample from each cluster
                    extremeSamples = self.funcSelectExtremeSamplesFromHClust(strBetaMetric=bMetric, npaAbundanceMatrix=npaTransposedAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount=iSampleSelectionCount)

                    #Add selected samples
                    if(not extremeSamples == False):
                        astrSelectionMethod = self.dictConvertBMetricExtreme[bMetric]
                        if not astrSelectionMethod in dictSelectedSamples:
                            dictSelectedSamples[astrSelectionMethod]=list()
                        dictSelectedSamples[astrSelectionMethod].extend(extremeSamples)

            logging.info("Selected Samples 2,3b")
            logging.info(dictSelectedSamples)
            return dictSelectedSamples

    #Start micropita selection
    def funcRun(self, fIsAlreadyNormalized, fCladesAreSummed, strOutputFile="MicroPITAOutput.txt", cDelimiter = Constants.TAB, cFeatureNameDelimiter = "|", strInputAbundanceFile=None,
            strUserDefinedTaxaFile=None, strTemporaryDirectory="./TMP", strCheckedAbndFile = "", iSampleSelectionCount=0, iSupervisedSampleCount=1,
            strSelectionTechnique=None, strLabel=None, strStratify=None, sMetadataID=None, sLastMetadataName=None, fSumData=True, sFeatureSelectionMethod=None, sCostRange="0"):
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
	:param	strTemporaryDirectory:	Directory that will be used to store secondary files important to analysis but not the direct deliverable.
	:type	String	String directory path.
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
        :param	sCostRange:	A comma delimited list of intergers (positive and negative) to be used as the cost range in the supervised methods (in MLPY only the highest is used).
        :type	String	String (comma delimited integers).
        :return	Selected Samples:	Samples selected by methods.
        :type	Dictionary	{"Selection Method":["SampleID","SampleID","SampleID",...]}
	"""

        #microPITA object
        microPITA = MicroPITA()

        #SVM parameters
        #Constants associated with the abundance to SVM input file conversion
        c_ABUNDANCE_DELIMITER=Constants.TAB
        c_NORMALIZE_RELATIVE_ABUNDANCY=True
        c_SKIP_FIRST_COLUMN=True

        #Diversity metrics to run
        diversityMetricsAlpha = [microPITA.c_strInverseSimpsonDiversity]
        diversityMetricsBeta = []#[microPITA.c_strBrayCurtisDissimilarity]
        inverseDiversityMetricsBeta = [microPITA.c_strInvBrayCurtisDissimilarity]
        diversityMetricsAlphaNoNormalize = [microPITA.c_strChao1Diversity]
        diversityMetricsBetaNoNormalize = [microPITA.c_strBrayCurtisDissimilarity]#[]
        inverseDiversityMetricsBetaNoNormalize = [microPITA.c_strInvBrayCurtisDissimilarity]#[]

        #Targeted taxa
        userDefinedTaxa = []

        #Perform different flows flags
        c_RUN_MAX_DIVERSITY_1 = False
        if(microPITA.c_strDiversity in strSelectionTechnique):
            c_RUN_MAX_DIVERSITY_1 = True
        c_RUN_REPRESENTIVE_DISSIMILARITY_2 = False
        if(microPITA.c_strRepresentativeDissimilarity in strSelectionTechnique):
            c_RUN_REPRESENTIVE_DISSIMILARITY_2 = True
        c_RUN_MAX_DISSIMILARITY_3 = False
        if(microPITA.c_strExtremeDissimilarity in strSelectionTechnique):
            c_RUN_MAX_DISSIMILARITY_3 = True
        c_RUN_RANK_AVERAGE_USER_4 = False
        if(microPITA.c_strTaxa in strSelectionTechnique):
            c_RUN_RANK_AVERAGE_USER_4 = True
            if(strUserDefinedTaxaFile == None):
                c_RUN_RANK_AVERAGE_USER_4 = False
                logging.error("MicroPITA.funcRun. No taxa file was given for taxa selection.")
            #Read in taxa list, break down to lines and filter out empty strings
            with open(strUserDefinedTaxaFile,'r') as fhndlTaxaInput:
                userDefinedTaxa = filter(None,fhndlTaxaInput.read().split(Constants.ENDLINE))

        c_RUN_RANDOM_5 = False
        if(microPITA.c_strRandom in strSelectionTechnique):
            c_RUN_RANDOM_5 = True
        c_RUN_DISTINCT = False
        if((microPITA.c_strDistinct in strSelectionTechnique) and (not strLabel == None)):
            c_RUN_DISTINCT = True
        c_RUN_DISCRIMINANT = False
        if((microPITA.c_strDiscriminant in strSelectionTechnique) and (not strLabel == None)):
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

        #Holds the top ranked samples from different metrics
        #dict[metric name] = [samplename,samplename...]
        selectedSamples = dict()

        #Check/reduce raw abundance data
        strInputAbundanceFile = AbundanceTable.funcCheckRawDataFile(strReadDataFileName=strInputAbundanceFile, sLastMetadataName=sLastMetadataName, strOutputFileName=strCheckedAbndFile)

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
        dictTotalMetadata = totalAbundanceTable.funcGetMetadataCopy()

        #Log metadata keys
        logging.debug(" ".join(["Micropita:funcRun.","Received metadata keys=",str(dictTotalMetadata.keys())]))

        #If there is only 1 unique value for the labels, do not run the Supervised methods
        if len(set(dictTotalMetadata.get(strLabel,[]))) < 2:
            c_RUN_DISCRIMINANT = False
            c_RUN_DISTINCT = False
            logging.error("".join(["The label ",strLabel," did not have 2 or more values. Labels found="]+dictTotalMetadata.get(strLabel,[])))

        logging.debug(" ".join(["Micropita:funcRun.","Received metadata=",str(dictTotalMetadata)]))

        #Run unsupervised methods###
        #Stratify the data if need be and drop the old data
        lStratifiedAbundanceTables = None
        if (not strStratify == None) and (not strStratify == "None"):
            lStratifiedAbundanceTables = totalAbundanceTable.funcStratifyByMetadata(strStratify,xWriteToFile="".join([os.path.split(strOutputFile)[0],"/"]))
        else:
            lStratifiedAbundanceTables = [totalAbundanceTable]

        #For each stratified abundance block or for the unstratfified abundance
        #Run the unsupervised blocks
        for stratAbundanceTable in lStratifiedAbundanceTables:
            logging.info("Running abundance block:"+stratAbundanceTable.funcGetName())

            #Only perform if the data is not yet normalized
            if not stratAbundanceTable.funcIsNormalized():

                #Need to first work with unnormalized data
                if((c_RUN_MAX_DIVERSITY_1)or(c_RUN_REPRESENTIVE_DISSIMILARITY_2)or(c_RUN_MAX_DISSIMILARITY_3)):

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

            #Generate selection by the rank average of user defined taxa
            #Expects (Taxa (row) by Samples (column))
            #Expects a column 0 of taxa id that is skipped
            #Returns [(sample name,average,rank)]
            if(c_RUN_RANK_AVERAGE_USER_4):
              if not microPITA.c_strUserRanked in selectedSamples:
                  selectedSamples[microPITA.c_strUserRanked]=list()
              selectedSamples[microPITA.c_strUserRanked].extend(microPITA.funcSelectTargetedTaxaSamples(abndMatrix=stratAbundanceTable, lsTargetedTaxa=userDefinedTaxa, iSampleSelectionCount=iSampleSelectionCount, sMethod=sFeatureSelectionMethod))
            logging.info("Selected Samples 4")
            logging.info(selectedSamples)

            #5::Select randomly
            #Expects sampleNames = List of sample names [name, name, name...]
            if(c_RUN_RANDOM_5):

                #Select randomly from sample names
                randomlySelectedSamples = microPITA.funcGetRandomSamples(lsSamples=stratAbundanceTable.funcGetSampleNames(), iNumberOfSamplesToReturn=iSampleSelectionCount)
                if not microPITA.c_strRandom in selectedSamples:
                    selectedSamples[microPITA.c_strRandom]=list()
                selectedSamples[microPITA.c_strRandom].extend(list(randomlySelectedSamples))

            logging.info("Selected Samples 5")
            logging.info(selectedSamples)

        #Run supervised methods#
        lStratifiedAbundanceTables = None
        totalAbundanceTable.funcNormalize()
        if(c_RUN_DISTINCT or c_RUN_DISCRIMINANT):
            selectedSamples.update(self.funcRunSupervisedMethods(abundanceTable=totalAbundanceTable,fRunDistinct=c_RUN_DISTINCT, fRunDiscriminant=c_RUN_DISCRIMINANT,
                                   strOuputSVMFile="".join([strTemporaryDirectory,"/",os.path.splitext(os.path.basename(totalAbundanceTable.funcGetName()))[0],"-SVM.txt"]),
                                   strTMPDirectory=strTemporaryDirectory,
                                   strSupervisedMetadata=strLabel, iSampleSVMSelectionCount=sampleSVMSelectionCount))
            logging.info("Selected Samples Unsupervised")
            logging.info(selectedSamples)

        return selectedSamples

    #Happy path tested
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
            strOutputContent = "".join([strOutputContent,sKey,Constants.COLON,", ".join(dictSelection[sKey]),Constants.ENDLINE])

        #Write to file
        if(not strOutputContent == ""):
            with open(strOutputFilePath,'w') as fHndlOutput:
                fHndlOutput.write(str(strOutputContent))
        logging.debug("".join(["Selected samples output to file:",strOutputContent]))

    #Happy Path tested
    #Reads in an output selection file from micropita and formats it into a dictionary
    #Dictionary formatted as- {"method":["sample selected","sample selected"...]}
    @staticmethod
    def funcReadSelectionFileToDictionary(strInputFile):
	"""
	Reads in an output selection file from micropita and formats it into a dictionary.
	
	:param	strInputFile:	String path to file to read and translate into a dictionary.
                                {"method":["sample selected","sample selected"...]}
	:type	String	String path to file to read and translate
	"""

        #Read in selection file
        strSelection = ""
        with open(strInputFile,'r') as fHndlInput:
            strSelection = fHndlInput.read()

        #Dictionary to hold selection data
        dictSelection = dict()
        for strSelectionLine in filter(None,strSelection.split(Constants.ENDLINE)):
            astrSelectionMethod = strSelectionLine.split(Constants.COLON)
            dictSelection[astrSelectionMethod[0].split()[0]] = [strSample.split()[0] for strSample in filter(None,astrSelectionMethod[1].split(Constants.COMMA))]

        #Return dictionary
        return dictSelection

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicroPITA.py", 
    description = """Selects samples from abundance tables based on various selection schemes.""" )

#Arguments
#Optional parameters
#Logging
argp.add_argument(Constants_Arguments.c_strLoggingArgument, dest="strLogLevel", metavar= "Log level", default="INFO", 
                  choices=Constants_Arguments.c_lsLoggingChoices, help= Constants_Arguments.c_strLoggingHelp)

#Abundance associated
argp.add_argument(Constants_Arguments.c_strIDNameArgument, dest="sIDName", metavar= "Sample ID Metadata Name", help= Constants_Arguments.c_strIDNameHelp)
argp.add_argument(Constants_Arguments.c_strLastMetadataNameArgument, dest="sLastMetadataName", metavar= "Last Metadata Name",
                  help= Constants_Arguments.c_strLastMetadataNameHelp)
argp.add_argument(Constants_Arguments.c_strIsNormalizedArgument, dest="fIsNormalized", action = "store_true", default=False,
                  help= Constants_Arguments.c_strIsNormalizedHelp)
argp.add_argument(Constants_Arguments.c_strIsSummedArgument, dest="fIsSummed", action = "store_true", default= False, help= Constants_Arguments.c_strIsSummedHelp)
argp.add_argument(Constants_Arguments.c_strSumDataArgument, dest="fSumData", action = "store_false", default = True, help= Constants_Arguments.c_strSumDataHelp)
argp.add_argument(Constants_Arguments.c_strTargetedFeatureMethodArgument, dest="sFeatureSelection", metavar= "Feature Selection Method", default=Constants_Arguments.lsTargetedFeatureMethodValues[0], 
                  choices=Constants_Arguments.lsTargetedFeatureMethodValues, help= Constants_Arguments.c_strTargetedFeatureMethodHelp)
argp.add_argument(Constants_Arguments.c_strFileDelimiterArgument, dest= "cFileDelimiter", action= "store", metavar="File Delimiter", default=Constants.TAB, help=Constants_Arguments.c_strFileDelimiterHelp) 
argp.add_argument(Constants_Arguments.c_strFeatureNameDelimiterArgument, dest= "cFeatureNameDelimiter", action= "store", metavar="Feature Name Delimiter", default=Constants.PIPE, help=Constants_Arguments.c_strFeatureNameDelimiterHelp) 

#Select count
argp.add_argument(Constants_Arguments.c_strUnsupervisedCountArgument, dest="iUnsupervisedSelectionCount", metavar = "Number Samples To Select (Unsupervised)", default=0, type = int, help = Constants_Arguments.c_strUnsupevisedCountHelp)
argp.add_argument(Constants_Arguments.c_strTargetedSelectionFileArgument, dest="strFileTaxa", metavar = "Targeted Feature File", default=None, help = Constants_Arguments.c_strTargetedSelectionFileHelp)
argp.add_argument(Constants_Arguments.c_strUnsupervisedStratifyMetadataArgument, dest="sUnsupervisedStratify", metavar= "Metadata to Stratify Unsupervised Selection", default=None, 
                  help= Constants_Arguments.c_strUnsupervisedStratifyMetadataHelp)

#SVM label
#Label parameter to be used with SVM
argp.add_argument(Constants_Arguments.c_strSupervisedLabelArgument, dest="sLabel", metavar= "Supervised Label Metadata Name", default=None, help= Constants_Arguments.c_strSupervisedLabelCountHelp)
argp.add_argument(Constants_Arguments.c_strSupervisedLabelCountArgument, dest="iSupervisedCount", metavar= "Supervised Sample Selection Count", default=0, type=int,
                  help= Constants_Arguments.c_strSupervisedLabelCountHelp)
argp.add_argument(Constants_Arguments.c_strCostArgument, dest="sCostRange", metavar= "Cost for Supervised Selection Actual cost = 2^sCostRange.", default="1",
                  help= Constants_Arguments.c_strCostHelp)

#Output
argp.add_argument(Constants_Arguments.c_strTemporaryDirectoryArgument, dest="strTMPDir", metavar = "Temporary Directory", default=None, help = Constants_Arguments.c_genericTMPDirLocationHelp)
argp.add_argument(Constants_Arguments.c_strCheckedAbundanceFileArgument, dest="strCheckedFile", metavar = "Checked Abundance File Name", default="", help = Constants_Arguments.c_strCheckedAbundanceFileHelp)

#Required
#Input
#Abundance file
argp.add_argument("strFileAbund", metavar = "Abundance file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument("strOutFile", metavar = "Selection Output File", help = Constants_Arguments.c_strGenericOutputDataFileHelp)
#Selection parameter
argp.add_argument("strSelection", metavar = "Selection Methods", help = Constants_Arguments.c_strSelectionTechniquesHelp, nargs="*")


__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    #TODO does this stop the full analysis process? Not if the selection file already exists...
    if not args.sIDName:
        logging.error("MicroPITA::Did not received a value for sIDName. MiroPITA did not run.")
        return False
    if not args.sLastMetadataName:
        logging.error("MicroPITA::Did not received a value for sIDName. MiroPITA did not run.")
        return False

    #Run micropita
    logging.info("Start microPITA")
    microPITA = MicroPITA()

    #Check inputs
    #If a target feature file is given make sure that targeted feature is in the selection methods, if not add
    if args.strFileTaxa:
        if microPITA.c_strTaxa not in args.strSelection:
            args.selection.append(microPITA.c_strTaxa)

    #If a supervised selection method is indicated make sure supervised selection count is indicated and above 0
    #Otherwise stop
    if len(set(microPITA.c_lsAllSupervisedMethods)&set(args.strSelection))>0:
        if args.iSupervisedCount < 1:
            logging.error("".join(["MicroPITA::Did not receive a selection count for supervised selection above 0, received=",
                                   str(args.iSupervisedCount),". Did not continue analysis."]))
            return -1

    #If an unsupervised selection method is indicated make sure the unsupervised selection count is above 0
    #Otherwise stop
    if len(set(microPITA.c_lsAllUnsupervisedMethods)&set(args.strSelection))>0:
        if args.iUnsupervisedSelectionCount < 1:
            logging.error("".join(["MicroPITA::Did not receive a selection count for unsupervised selection above 0, received=",
                                   str(args.iUnsupervisedSelectionCount),". Did not continue analysis."]))
            return -2

    #If the tmp directory is not made, make
    if not args.strTMPDir:
        args.strTMPDir = "."


    dictSelectedSamples = microPITA.funcRun(fIsAlreadyNormalized=args.fIsNormalized, fCladesAreSummed=args.fIsSummed, strOutputFile=args.strOutFile,
                                        cDelimiter=args.cFileDelimiter, cFeatureNameDelimiter = args.cFeatureNameDelimiter, strInputAbundanceFile=args.strFileAbund,
                                        strUserDefinedTaxaFile=args.strFileTaxa, strTemporaryDirectory=args.strTMPDir, strCheckedAbndFile = args.strCheckedFile,
                                        iSampleSelectionCount=args.iUnsupervisedSelectionCount,
                                        iSupervisedSampleCount=args.iSupervisedCount, strLabel=args.sLabel,
                                        strStratify=args.sUnsupervisedStratify, strSelectionTechnique=args.strSelection,
                                        sMetadataID=args.sIDName, sLastMetadataName=args.sLastMetadataName, fSumData=args.fSumData,
                                        sFeatureSelectionMethod=args.sFeatureSelection,sCostRange=args.sCostRange)
    logging.info("End microPITA")

    #Log output for debugging
    logging.debug("".join(["Returned the following samples:",str(dictSelectedSamples)]))

    #Write selection to file
    microPITA.funcWriteSelectionToFile(dictSelection=dictSelectedSamples,strOutputFilePath=args.strOutFile)

if __name__ == "__main__":
    _main( )

