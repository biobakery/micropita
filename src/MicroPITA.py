#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Run analysis for the microPITA paper
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = "1.0"
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
import logging
import math
import mlpy
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
    c_SHANNON_A_DIVERSITY = Diversity.c_SHANNON_A_DIVERSITY
    c_SIMPSON_A_DIVERSITY = Diversity.c_SIMPSON_A_DIVERSITY
    c_INV_SIMPSON_A_DIVERSITY = Diversity.c_INV_SIMPSON_A_DIVERSITY
    c_CHAO1_A_DIVERSITY = Diversity.c_CHAO1_A_DIVERSITY

    #Diversity metrics Beta
    c_UNIFRAC_B_DIVERSITY = Diversity.c_UNIFRAC_B_DIVERSITY
    c_WEIGHTED_UNIFRAC_B_DIVERSITY = Diversity.c_WEIGHTED_UNIFRAC_B_DIVERSITY
    c_BRAY_CURTIS_B_DIVERSITY = Diversity.c_BRAY_CURTIS_B_DIVERSITY

    #Additive inverses of diversity metrics beta
    c_INVERSE_BRAY_CURTIS_B_DIVERSITY = Diversity.c_INVERSE_BRAY_CURTIS_B_DIVERSITY
    c_INVERSE_UNIFRAC_B_DIVERSITY = Diversity.c_INVERSE_UNIFRAC_B_DIVERSITY
    c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY = Diversity.c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY

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
    c_DIVERSITY_1 = "".join([c_strDiversity,"_I"])
    c_DIVERSITY_2 = "".join([c_strDiversity,"_C"])
    c_EXTREME_DISSIMILARITY_1 = "".join([c_strExtremeDissimilarity,"_B"])
    c_EXTREME_DISSIMILARITY_2 = "".join([c_strExtremeDissimilarity,"_U"])
    c_EXTREME_DISSIMILARITY_3 = "".join([c_strExtremeDissimilarity,"_W"])
    c_REPRESENTATIVE_DISSIMILARITY_1 = "".join([c_strRepresentativeDissimilarity,"_B"])
    c_REPRESENTATIVE_DISSIMILARITY_2 = "".join([c_strRepresentativeDissimilarity,"_U"])
    c_REPRESENTATIVE_DISSIMILARITY_3 = "".join([c_strRepresentativeDissimilarity,"_W"])
    c_RANDOM = c_strRandom
    c_USER_RANKED = c_strTaxa
    c_SVM_CLOSE = c_strDiscriminant
    c_SVM_FAR = c_strDistinct

    #Targeted feature settings
    c_TARGETED_METHOD_RANKED = Constants_Arguments.c_TARGETED_METHOD_RANKED
    c_TARGETED_METHOD_ABUNDANCE = Constants_Arguments.c_TARGETED_METHOD_ABUNDANCE

    #Technique groupings
    c_lsDiversityMethods = [c_DIVERSITY_1,c_DIVERSITY_2]

    #Converts ecology metrics into standardized method selection names
    convertAMetricDiversity = {c_INV_SIMPSON_A_DIVERSITY:c_DIVERSITY_1, c_CHAO1_A_DIVERSITY:c_DIVERSITY_2}
    convertMicroPITAToAMetric = {c_DIVERSITY_1:c_INV_SIMPSON_A_DIVERSITY, c_DIVERSITY_2:c_CHAO1_A_DIVERSITY}
    convertBMetricRepresentative = {c_BRAY_CURTIS_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_1, c_UNIFRAC_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_2, c_WEIGHTED_UNIFRAC_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_3}
    convertBMetricExtreme = {c_INVERSE_BRAY_CURTIS_B_DIVERSITY:c_EXTREME_DISSIMILARITY_1, c_INVERSE_UNIFRAC_B_DIVERSITY:c_EXTREME_DISSIMILARITY_2, c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY:c_EXTREME_DISSIMILARITY_3}

    #Used to store scores and pvalues
    c_SCORE_INDEX = 0
    c_PVALUE_INDEX = 1

    #Linkage used in the Hierarchical clustering
    c_HIERARCHICAL_CLUSTER_METHOD = 'average'

####Group 1## Diversity

    #Testing: Happy path Testing (8)
    def getTopRankedSamples(self, tempMatrix = None, tempSampleNames = None, tempTopAmount = None):
	"""
	Given a list of lists of measurements, for each list the indices of the highest values are returned. If tempSamplesNames is given
        it is treated as a list of string names that is in the order of the measurements in each list. Indices are returned or the sample
        names positionally associated with the indices.
	
	:param	tempMatrix:	List of lists [[value,value,value,value],[value,value,value,value]].
	:type	List of lists:	List of measurements. Each list is a different measurement. Each measurement in possionally related to a sample.
	:param	tempSampleNames:	List of sample names positionally related (the same) to each list (Optional).
	:type	List of strings:	List of strings
	:param	tempTopAmount:		The amount of top measured samples (assumes the higher measurements are better)
	:type	integer:	Integer amount of sample names/ indices to return.
	"""

        topRankList = []
        if(len(tempMatrix)<1):
            return topRankList
        for rowMetrics in tempMatrix:
            #Create 2 d array to hold value and index and sort
            rowsMetricsLength = len(rowMetrics)
            indexX = [rowMetrics,range(rowsMetricsLength)]
            indexX[1].sort(key = indexX[0].__getitem__,reverse = True)
            if(tempSampleNames == None):
                topRankList.append(indexX[1][:tempTopAmount])
            else:
                sampleIndexes = indexX[1][:tempTopAmount]
                addSamplesToRank = []
                for index in sampleIndexes:
                    addSamplesToRank.append(tempSampleNames[index])
                topRankList.append(addSamplesToRank)
        return topRankList

####Group 2## Representative Dissimilarity

    #Testing: Happy Path Tested for BrayCurtis and Inverse BrayCurtis
    def getBetaMetric(self, tempAbundancies=None, tempMetric=None):
	"""
	Takes a matrix of values and returns a beta metric matrix. The metric returned is indicated by name (tempMetric).
	
	:param	tempAbundancies:	Numpy array where row=samples and columns = features
	:type	Matrix:	Numpy array
	:param	tempMetric:	String name of beta metric. Possibilities are listed in microPITA.
	:type	String:	String name of beta metric. Possibilities are listed in microPITA.
	"""

        if(not ValidateData.isValidString(tempMetric)):
            return False
        elif(tempMetric == self.c_BRAY_CURTIS_B_DIVERSITY):
            return Diversity.getBrayCurtisDissimilarity(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == self.c_INVERSE_BRAY_CURTIS_B_DIVERSITY):
            return Diversity.getInverseBrayCurtisDissimilarity(tempSampleTaxaAbundancies=tempAbundancies)
        #Needs NOT normalized abundances
        elif(tempMetric == self.c_WEIGHTED_UNIFRAC_B_DIVERSITY):
            return Diversity.getUnifracDistance(tempSampleTaxaAbundancies=tempAbundancies, tempTaxonomyTree=None, tempWeighted=True)
        elif(tempMetric == self.c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY):
            return False
#            return Diversity.getUnifracDistance(tempSampleTaxaAbundancies=tempAbundancies, tempTaxonomyTree=None, tempWeighted=True)
        #Needs NOT normalized abundances
        elif(tempMetric == self.c_UNIFRAC_B_DIVERSITY):
            return Diversity.getUnifracDistance(tempSampleTaxaAbundancies=tempAbundancies, tempTaxonomyTree=None, tempWeighted=False)
        elif(tempMetric == self.c_INVERSE_UNIFRAC_B_DIVERSITY):
            return False
#            return Diversity.getUnifracDistance(tempSampleTaxaAbundancies=tempAbundancies, tempTaxonomyTree=None, tempWeighted=False)
        else:
            return False

    #Testing: Needs Testing
    def getCentralSamplesByKMedoids(self, tempMatrix=None, tempMetric=None, tempSampleNames=None, tempNumberSamplesReturned=0):
	"""
	Gets centroid samples by k-medoids clustering of a given matrix.
	
	:param	tempMatrix:	Numpy array where row=features and columns=samples
	:type	Matrix:	Numpy array
	:param	tempMetric:	String name of beta metric used as the distance metric.
	:type	String:	String name of beta metric. Possibilities are listed in microPITA.
	:param	tempSampleNames:	The names of the sample
	:type	List:	List of strings
	:param	tempNumberSamplesReturned:	Number of samples to return, each will be a centroid of a sample.
	:type	Integer:	Number of samples to return
	"""

        #Validate parameters
        if(tempNumberSamplesReturned < 0):
            logging.error("MicroPITA.getCentralSamplesByKMedoids. Number of samples to return must be atleast 1.")
            return False

        #Count of how many rows
        sampleCount = tempMatrix.shape[0]
        if(tempNumberSamplesReturned > sampleCount):
            logging.error("".join(["MicroPITA.getCentralSamplesByKMedoids. There are not enough samples to return the amount of samples specified. Return sample count = ",str(tempNumberSamplesReturned),". Sample number = ",str(sampleCount),"."]))
            return False

        #Samples to return
        returningSamples = list()

        #If the cluster count is equal to the sample count return all samples
        if(sampleCount == tempNumberSamplesReturned):
            return list(tempSampleNames)

        #Get distance matrix
        distanceMatrix=self.getBetaMetric(tempAbundancies=tempMatrix, tempMetric=tempMetric)
        if(ValidateData.isFalse(distanceMatrix)):
          logging.error("MicroPITA.getCentralSamplesByKMedoids. Received false for betaMetrix matrix generation, returning false.")
          return False

        #Log distance matrix
        logging.debug("".join(["Distance matrix for representative selection using metric=",str(tempMetric)]))

        if(( tempMetric==Diversity.c_UNIFRAC_B_DIVERSITY ) or ( tempMetric==Diversity.c_WEIGHTED_UNIFRAC_B_DIVERSITY )):
          distanceMatrix = distanceMatrix['distance_matrix'][0]

        #Perform Kmedoid clustering in pycluster
        lclusterid, error, nfound = kmedoids (distance=distanceMatrix, nclusters=tempNumberSamplesReturned, npass=1000, initialid=None)
        #TODO make sure you are getting condensed for unifrac and braycurtis
        #TODO figure out a way to select optimal k medoid iteration 
        logging.debug("Results from the kmedoid method in representative selection:")
        logging.debug("clusterid:"+str(list(lclusterid)))
        logging.debug("error:"+str(error))
        logging.debug("nfound:"+str(nfound))

        #Convert centroid indexes to names and return
        #Return centroids
        return [tempSampleNames[iindex] for iindex in set(lclusterid)]

####Group 3## Highest Dissimilarity
    #Testing: Needs testing
    def selectExtremeSamplesFromHClust(self, tempBetaMetric, tempAbundanceMatrix, tempSampleNames, tempSelectSampleCount):
	"""
	Select extreme samples from HClustering.
	
	:param	tempBetaMetric:	The beta metric to use for distance matrix generation
	:type	String:	The name of the beta metric to use
	:param	tempAbundanceMatrix:	Numpy array where row=samples and columns=features
	:type	Matrix:	Numpy Array
	:param	tempSampleNames:	The names of the sample
	:type	List:	List of strings
	:param	tempSelectSampleCount:	Number of samples to select (return)
	:type	Integer:	Integer number of samples returned
	"""

        #If they want all the sample count, return all sample names
        sampleCount=len(tempAbundanceMatrix[:,0])
        if(tempSelectSampleCount==sampleCount):
          return(tempSampleNames)

        #Holds the samples to be returned
        returnSamples = []

        #Holds the roots to each cluster
        clusterRoots = []
        #List of lists. Each internal list is a cluster containing all of the membership's indices
        clusterIndices = []

        #Generate beta matrix
        #Returns condensed matrix
        tempDistanceMatrix = self.getBetaMetric(tempAbundancies=tempAbundanceMatrix, tempMetric=tempBetaMetric)

        #Feed beta matrix to linkage to cluster
        #Send condensed matrix
        linkageMatrix = hcluster.linkage(tempDistanceMatrix, method=self.c_HIERARCHICAL_CLUSTER_METHOD)

        #Extract cluster information from cluster
        #Sort by distance (column index = 2)
#TODO is this always sorted by distance, it looks like it but I need to know for sure
#The rest of the code assumes so
#       linkageMatrix = linkageMatrix[linkageMatrix[:,2].argsort()]

        #This method just takes the lowest metric measurement (highest distance pairs/clusters)
        #Works much better than the original technique
        iSampleCount = len(tempSampleNames)
        iSelectedSampleCount = 0
        for row in linkageMatrix:
            if(iSelectedSampleCount < tempSelectSampleCount):
                iNode = int(row[0])
                if(iNode<iSampleCount):
                    returnSamples.append(tempSampleNames[iNode])
                    iSelectedSampleCount = iSelectedSampleCount + 1
            if(iSelectedSampleCount < tempSelectSampleCount):
                iNode = int(row[1])
                if(iNode<iSampleCount):
                    returnSamples.append(tempSampleNames[iNode])
                    iSelectedSampleCount = iSelectedSampleCount + 1

        #Return selected samples
        return returnSamples

####Group 4## Rank Average of user Defined Taxa

    def getAverageAbundanceSamples(self, tempMatrix, tempTargetedTaxa, fRank=False):
	"""
	Averages feature abundance. Expects a column 0 of taxa id that is skipped.
	
	:param	tempMatrix:	Abundance Table to analyse
	:type	AbundanceTable:	Abundance Table
	:param	tempTargetedTaxa:	String names
	:type	list:	list of string names of taxa which are measured after ranking against the full sample
        :param  fRank:	Indicates to rank the abundance before getting the average abundance of the features (default false)
        :type   boolean:	Flag indicating ranking abundance before calculating average feature measurement (false= no ranking)
	"""

        #Sample rank averages [[sample,average abundance of selected taxa]]
        #Returned
        sampleAbundanceAverages = []
        
        #Get sample names
        sampleNames = tempMatrix.funcGetSampleNames()
        #Get taxa names
        allTaxaNames = tempMatrix.funcGetFeatureNames()

        #Rank if needed
        if fRank:
            tempMatrix = tempMatrix.funcRankAbundance()
            if tempMatrix == None:
                logging.error("".join(["MicroPITA.getAverageAbundanceSamples. Could not rank the abundance table, returned false."]))
                return False  

        #Get an abundance table compressed to features of interest
        abndReducedTable = tempMatrix.funcGetFeatureAbundanceTable(tempTargetedTaxa)
        if abndReducedTable == None:
            logging.error("".join(["MicroPITA.getAverageAbundanceSamples. Could not reduce the abundance table to the given taxa, returned false."]))
            return False  

        #If the taxa to be selected are not in the list
        #Return nothing and log
        lsMissing = []
        for tempTaxa in tempTargetedTaxa:
            if not tempTaxa in allTaxaNames:
                lsMissing.append(tempTaxa)
            else:
                #Check to make sure the taxa of interest is not average abundance of 0
                iFeatureSum = tempMatrix.funcGetFeatureSumAcrossSamples(tempTaxa)
                if (iFeatureSum == None) or (iFeatureSum == 0):
                    lsMissing.append(tempTaxa)

        if len(lsMissing) > 0:
            logging.error("".join(["MicroPITA.getAverageAbundanceSamples. The following feature/s is not in the data set or has the abundance of 0. ",",".join(lsMissing)]))
            return False        

        #For each sample name get average abundance
        for sName in sampleNames:
            npaFeaturesSample = abndReducedTable.funcGetSample(sName)
            sampleAbundanceAverages.append([sName,sum(npaFeaturesSample)/float(len(npaFeaturesSample))])

        #Sort based on average
        sampleAbundanceAverages = sorted(sampleAbundanceAverages, key = lambda sampleData: sampleData[1], reverse = True)

        #return
        return sampleAbundanceAverages

    def selectTargetedTaxaSamples(self, tempMatrix, tempTargetedTaxa, sampleSelectionCount, sSampleIDName, sMethod):
      """
      Selects samples with the highest ranks or abundance of targeted features.

      :param	tempMatrix:	Abundance table to analyse 
      :type	AbundanceTable:	Abundance table
      :param	tempTargetedTaxa:	List of features
      :type	list:	list of strings
      :param	sampleSelectionCount:	Number of samples to select
      :type	integer:	integer
      :param	sSampleIDName:	ID for the feature name
      :type	string:	String of the metadata indicating the feature name
      :param	sMethod:	Method to select targeted features
      :type	string:	String (Can be values found in microPITA)
      """

      if(len(tempTargetedTaxa) < 1):
        logging.error("MicroPITA.selectTargetedTaxaSamples. Taxa defined selection was requested but no features were given.")
        return []
      if not sMethod:
        logging.error("MicroPITA.selectTargetedTaxaSamples. Taxa defined selection was requested but no Method were given.")
        return []

      lsTargetedSamples = None
      #Lower is better
      if sMethod.lower() == self.c_TARGETED_METHOD_RANKED.lower():
          lsTargetedSamples = self.getAverageAbundanceSamples(tempMatrix=tempMatrix, tempTargetedTaxa=tempTargetedTaxa, fRank=True)
          return [diversity[0] for diversity in lsTargetedSamples[(sampleSelectionCount*-1):]]
      #Higher is better
      elif sMethod.lower() == self.c_TARGETED_METHOD_ABUNDANCE.lower():
          lsTargetedSamples = self.getAverageAbundanceSamples(tempMatrix=tempMatrix, tempTargetedTaxa=tempTargetedTaxa, fRank=False)
          return [diversity[0] for diversity in lsTargetedSamples[:sampleSelectionCount]]
      
      return []

####Group 5## Random

    #Returns random sample names of the number given
    #@return A list of sample name lists with each list having a randomly selected 
    #amount of sample names as defined by tempNumberOfSamplesToReturn [[randomName1, randomName2,...],[randomName1, randomName2,...]]
    def getRandomSamples(self, tempSamples=None, tempNumberOfSamplesToReturn=0):
	"""
	Returns random sample names of the number given. No replacement.
	
	:param	tempSamples:	List of sample names 
	:type	list:	list of strings
	:param	tempNumberOfSamplesToReturn:	Number of samples to select
	:type	integer:	integer
	"""
        #Input matrix sample count
        sampleCount = len(tempSamples)

        #Validate Number of samples is not negative and not greater than the matrix
        if((tempNumberOfSamplesToReturn < 0) or (tempNumberOfSamplesToReturn > sampleCount)):
            logging.error("".join(["MicroPITA.getRandomSamples. TempNumberOfSamples was not useful for deriving subset. ",str(tempNumberOfSamplesToReturn)]))
            return False

        #Return the full matrix if they ask for a return matrix where length == original
        if(tempNumberOfSamplesToReturn == sampleCount):
            return tempSamples

        #Get the random indices for the sample (without replacement)
        randomIndices = random.sample(range(sampleCount), tempNumberOfSamplesToReturn)

        #Create a boolean array of if indexes are to be included in the reduced array
        randomIndicesBoolean = list()
        for rIndex in range(sampleCount):
            randomIndicesBoolean.append(rIndex in randomIndices)

        #Reduce array to just what is randomly sampled
        return np.compress(condition=randomIndicesBoolean, a=tempSamples)

####Group 6## Supervised
        
    #@return Dictionary of file paths which were generated by the model and prediction steps
    def runLIBSVM(self, abndTable, tempOutputSVMFile, sLabelID, tempSVMScaleLowestBound = 0, tempSVMLogG="-5,-4,-3,-2,-1,0,1,2,3,4,5", tempSVMLogC="-5,-4,-3,-2,-1,0,1,2,3,4,5", tempSVMProbabilistic=True):
	"""
	Runs a linear SVM (Adapted from the easy.py script included in the standard libsvm intall).

	:param	tempOutputSVMFile:	String File path and name to output in SVM format based on the input file
	:type	string:	String file path
	:param	tempMatrixLabels:	List of labels (does not have to be strings, just have to be appropriate for labels when casted to string type
	:type	List:	List of strings
	:param	sLastMetadataName:	String Found immediately before the first row to read (skips header rows)
	:type	string:	string
	:param	tempSkipFirstColumn:	Boolean True indicates the first column will be skipped (due to it containing row identifying information like OTU names).
	:type	boolean:	boolean
	:param	tempNormalize:	Boolean True indicates normalizes to relative abundancy per sample (column)
	:type	boolean:	boolean
	:param	tempSVMScaleLowestBound:	Lowest value (0 or -1) the values were scaled to (max = 1)
	:type	integer:	integer
	:param	tempSVMLogC:	Comma delimited string giving values for C cross validation
	:type	string:	Comma delimited string
	:param	tempSVMProbabilistic:	Get probablistic outcome for SVM
	:type	boolean:	boolean
	"""

        #Create SVM object
        svm = SVM()

        #Holds files generated by SVM code
        #Files generated by modeling
        modelFiles = None
        #Files generated by prediction
        predictionFiles = None

        #Convert abundancies file to SVM file
        noError = svm.convertAbundanceTableToSVMFile(abndAbundanceTable=abndTable, tempOutputSVMFile=tempOutputSVMFile, sMetadataLabel=sLabelID)

        modelFiles = False
        predictionFiles = False
        #Run SVM
        if(noError==True):
            modelFiles = svm.createLinearModel(tempInputFileName=tempOutputSVMFile, tempScaling=tempSVMScaleLowestBound, tempLogC=tempSVMLogC, tempProbabilistic=tempSVMProbabilistic)
            if(not modelFiles == False):
                predictionFiles = svm.predictFromLinearModel(tempDataFileName=tempOutputSVMFile, tempModelFileName=modelFiles[svm.c_KEYWORD_MODEL_FILE], tempRangeFileName=modelFiles[svm.c_KEYWORD_RANGE_FILE], tempProbabilistic=tempSVMProbabilistic)
            else:
                logging.error("".join(["MicroPITA.runSVM: Could not run prediction, model file was false."]))
        #Combine output dictionarys and return
        if(not predictionFiles == False):
            for pKey in predictionFiles:
                modelFiles[pKey]=(predictionFiles[pKey])
        else:
            logging.error("MicroPITA.runSVM: Prediction files were false.")
        return modelFiles

    #@return Dictionary of file paths which were generated by the model and prediction steps
    def runMLPYSVM(self, abndAbundanceTable, sMetadataForLabel, strOutputModelFile, strPredictionFile, iCost = 1, fSVMProbabilistic=True):
	"""
	Runs a linear SVM using MLPY.
	
	:param	abndAbundanceTable:	Abundance table of data
	:type	Abundance table:	Abundance table of data
	:param	sMetadataForLabel:	Metadata for label used to supervised learning
	:type	string:	string
	:param	strOutputModelFile:	String File path and name to output the SVM model
	:type	string:	String file path
        :param	iCost	Integer which will be used to define the cost as 2^iCost
        :type	integer	Integer
	:param	fSVMProbabilistic:	Get probablistic outcome for SVM
	:type	boolean:	boolean
	"""

        #Get labels
        lsLabels = abndAbundanceTable.funcGetMetadata(sMetadataForLabel)

        #Create a new label list but coded as integers
        setsLabels = set(lsLabels)
        dictLabels = dict(zip(setsLabels, range(len(setsLabels))))
        ldLabels = [dictLabels[sLabel] for sLabel in lsLabels]

        logging.debug("".join(["dictLabels ", str(dictLabels)]))
        logging.debug("".join(["ldLabels ", str(ldLabels)]))

        #Build a dict of weights per label {label:weight, label:weight}
        #Get the occurence of each label
        dictWeights = dict()
        for sLabelKey in dictLabels:
            sCurLabel = dictLabels[sLabelKey]
            dictWeights[sCurLabel] = ldLabels.count(sCurLabel)

        #Divide the highest occurence each occurence
        iMaxOccurence = max(dictWeights.values())
        for sWeightKey in dictWeights:
            dictWeights[sWeightKey]=iMaxOccurence/float(dictWeights[sWeightKey])

        logging.debug("".join(["dictWeights ", str(dictWeights)]))

        #Create SVM object
        svm = mlpy.LibSvm(C=math.pow(2,int(iCost)), probability=fSVMProbabilistic, weight=dictWeights)

        #Get 2D array of data
        npData = abndAbundanceTable.funcGetAbundanceCopy()
        lsSampleNames = abndAbundanceTable.funcGetSampleNames()
        llData = [list(npData[sSamples]) for sSamples in lsSampleNames]

        #Create classifier
        #x = [sample of the same number as labels, features]
        #y = [labels]
        svm.learn(x=llData,y=ldLabels)

        #Get Classes
        #[samples,features]
        retPred = svm.pred(llData)

        #Get probabilities
        #svm functions require [samples,features]
        retDistance = None
        if fSVMProbabilistic:
            retDistance = svm.pred_probability(llData)
        else:
            retDistance = svm.pred_values(llData)

        #Create output prediction file
        liSortedLabels = sorted(list(set(ldLabels)))
        strPredictionOutput = " ".join(["labels"]+[str(iLabel) for iLabel in liSortedLabels])+Constants.ENDLINE
        strPredictionOutput = strPredictionOutput + Constants.ENDLINE.join([" ".join([str(retPred[iIndexPredictions])]+[str(dPred) for dPred in list(retDistance[iIndexPredictions])])
                                for iIndexPredictions in xrange(len(retPred))])
        print "strPredictionFile", strPredictionFile 
        with open(strPredictionFile, 'w') as f:
            f.write(strPredictionOutput)
        f.close()

        #Save model to file
        svm.save_model(strOutputModelFile)

        return [retPred,retDistance,dictLabels,lsLabels]

    #Run the supervised methods
    def runSupervisedMethods(self, abundanceTable, fRunDistinct, fRunDiscriminant,
                                   strOuputSVMFile, strSupervisedMetadata, sampleSVMSelectionCount,
                                   sLastMetadataName, fSkipFirstColumn, fNormalize,
                                   iScaleLowestBound, strCostRange, fProbabilitic, fLibSVM=False):
        """
	Runs supervised methods.
	
	:param	abundanceTable:	AbundanceTable
	:type	AbudanceTable:	Data to analyze
	:param	fRunDistinct:	Run distinct selection method
	:type	boolean:	boolean (true runs method)
	:param	fRunDiscriminant:	Run discriminant method
	:type	boolean:	boolean (true runs method)
	:param	strOutputSVMFile:	File output from  SVM
	:type	string:	String
	:param	strSupervisedMetadata:	Label for supervised selection
	:type	string:	String
	:param	sampleSVMSelectionCount:	Number of samples to select
	:type	int:	int sample selection count
	:param	sLastMetadataName:	The metadata id for the last metadata in the abundance table
	:type	string:	String metadata id
	:param	fSkipFirstColumn:	Indicates skipping the first row
	:type	boolean:	boolean (true=skip)
	:param	fNormalize:	Indicates the need to normalize
	:type	boolean:	
	:param	iScaleLowestBound:	Integer min value data is scaled to (-1 0r 0)
	:type	integer: -1 or 0
	:param	strCostRange:	String of comma delimited values to cross validate cost values with
	:type	string:
	:param	fProbabilistic:	Indicator or getting probabilistic ouput from the SVM.
	:type	boolean:	
	"""

        #Run supervised blocks
        #Select supervised (using SVM)
        #Expects input file's matrix to be Taxa (row) by Sample (col) with a taxa id column (index=0)
        #Get sample names
        sampleNames = abundanceTable.funcGetSampleNames()
        #Will contain the samples selected to return
        dictSelectedSamples = dict()
        #Create and array to hold difference of the samples probabilities from the central probability
        centralDeviation = dict()

        #Remove all files associated with supervised methods
        #Create file names and delete if they exist
        fileNamePrefix = os.path.splitext(strOuputSVMFile)[0]
        fileDirectory = os.path.split(strOuputSVMFile)[0]
        scaledFile = "".join([fileNamePrefix,Constants.c_SCALED_FILE_EXT])
        rangeFile = "".join([fileNamePrefix,Constants.c_SCALING_PARAMETERS])
        cvOutFile = "".join([fileNamePrefix,Constants.c_CV_FILE_EXT])
        modelFile = "".join([fileNamePrefix,Constants.c_MODEL_FILE_EXT])
        scaledPredictFile = "".join([fileNamePrefix,Constants.c_SCALED_FOR_PREDICTION_FILE_EXT])
        predictFile = "".join([fileNamePrefix,Constants.c_PREDICT_FILE_EXT])
        if os.path.exists(scaledFile):
            os.remove(scaledFile)
        if os.path.exists(rangeFile):
            os.remove(rangeFile)
        if os.path.exists(cvOutFile):
            os.remove(cvOutFile)
        if os.path.exists(modelFile):
            os.remove(modelFile)
        if os.path.exists(scaledPredictFile):
            os.remove(scaledPredictFile)
        if os.path.exists(predictFile):
            os.remove(predictFile)
        if not os.path.exists(fileDirectory):
            os.makedirs(fileDirectory)

        if fLibSVM:
            #Run linear SVM
            svmRelatedData = self.runLIBSVM(abndTable=abundanceTable, tempOutputSVMFile=strOuputSVMFile, sLabelID=strSupervisedMetadata, tempSVMScaleLowestBound=iScaleLowestBound,
                                     tempSVMLogC=strCostRange, tempSVMProbabilistic=fProbabilitic)
            #Read in prediction file and select samples
            if svmRelatedData:
                if(SVM.c_KEYWORD_PREDICTION_FILE in svmRelatedData):
                    #Get prediction file path
                    predictionFile = svmRelatedData[SVM.c_KEYWORD_PREDICTION_FILE]
                    #Get the input file for the SVMS which has the original labels before classification
                    strSVMInputFile = svmRelatedData[SVM.c_KEYWORD_INPUT_FILE]
                    #Holds labels to compare to the predictions
                    lsOriginalLabels = None
                    #Open input file and get labels to compare to the predictions
                    with open(strSVMInputFile,'r') as f:
                        reader = csv.reader(f, delimiter=Constants.WHITE_SPACE, quoting=csv.QUOTE_NONE)
                        lsOriginalLabels = [row[0] for row in reader]
                    f.close()
                    #TODO test selection of samples
                    #Read in file
                    with open(predictionFile,'r') as f:
                        predictionLists = f.read()
                    f.close()
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
                        logging.debug("iCurLabel")
                        logging.debug(iCurLabel)
                        logging.debug("lsOriginalLabels[lineIndex]")
                        logging.debug(lsOriginalLabels[lineIndex])
                        if iCurLabel == lsOriginalLabels[lineIndex]:
                            logging.debug("Correctly predicted")
                            #Sum the absolute values of the differences
                            #Store the index and the deviation from the central probability
                            deviation = 0
                            for prediction in lineElements:
                                deviation += math.pow((float(prediction)-centralProbability),2)
                            if(not iCurLabel in centralDeviation):
                                centralDeviation[iCurLabel] = []
                            curSVMData = centralDeviation[iCurLabel]
                            curSVMData.append([deviation,lineIndex,iCurLabel])
                            centralDeviation[iCurLabel] = curSVMData
                            logging.debug("curSVMData")
                            logging.debug(curSVMData)

        #if not fLibSVM:
        else:
            #Can only run 1 cost so in MLPY mode run the highest cost
            iCost = max([int(strCost) for strCost in strCostRange.split(",")])

            #Run MLPY SVM
            svmRelatedData = self.runMLPYSVM(abndAbundanceTable=abundanceTable, sMetadataForLabel=strSupervisedMetadata,
                                             strOutputModelFile=strOuputSVMFile, strPredictionFile=predictFile,
                                             iCost=iCost, fSVMProbabilistic=True)
            #Returned data
            ldLabels = svmRelatedData[0]
            lldDistances = svmRelatedData[1]
            dictCovertToLabels = svmRelatedData[2]
            lsOriginalLabels = svmRelatedData[3]

            #Flip keys and values in dictionary
            dictConvertFromLabels = dict((iValue, strKey) for strKey, iValue in dictCovertToLabels.iteritems())

            #Get central probability
            centralProbability = 1.0 / float(len(lsOriginalLabels))

        #Go through predictions 
        #Check to make sure the predictions are correct
        iCorrectlyClassified = 0
        for iDistance, lDistances in enumerate(lldDistances):
            strCurLabel = dictConvertFromLabels[ldLabels[iDistance]]
            if strCurLabel == lsOriginalLabels[iDistance]:
                logging.debug("Correctly predicted")
                iCorrectlyClassified = iCorrectlyClassified + 1
                #Sum the absolute values of the differences
                #Store the index and the deviation from the central probability
                deviation = 0
                for prediction in lDistances:
                    deviation += math.pow((float(prediction)-centralProbability),2)
                if(not strCurLabel in centralDeviation):
                    centralDeviation[strCurLabel] = []
                curSVMData = centralDeviation[strCurLabel]
                curSVMData.append([deviation,iDistance,strCurLabel])
                centralDeviation[strCurLabel] = curSVMData
                logging.debug("curSVMData")
                logging.debug(curSVMData)

        logging.debug("".join(["iCorrectlyClassified",str(iCorrectlyClassified)]))

        #Sort sample by summed absolute deviations from the center and take the top N indexes
        for scurKey in centralDeviation:
            lcurLabeSamples = centralDeviation[scurKey]
            centralDeviation[scurKey] = sorted(lcurLabeSamples, key=operator.itemgetter(0))
            logging.debug("centralDeviation[scurKey]")
            logging.debug(centralDeviation[scurKey])

        selectedSamplesIndicesClose = list()
        selectedSamplesIndicesFar = list()
        #If the amount of samples needed are greater than what was analyzed with the SVM,
        #Return them all for both far and near.
        #Get the samples closeRuns to hyperplane
        #Get samples farthest from hyperplane
        #Make this balanced to the labels so for each sample label select the sampleSVMSelectionCount count of samples
        for scurKey in centralDeviation:
            lcurDeviations = centralDeviation[scurKey]
            iLengthCurDeviations = len(lcurDeviations)
            if(sampleSVMSelectionCount > iLengthCurDeviations):
                licurIndices = [measurement[1] for measurement in lcurDeviations]
                selectedSamplesIndicesClose.extend(licurIndices)
                selectedSamplesIndicesFar.extend(licurIndices)
            else:
                selectedSamplesIndicesClose.extend([measurement[1] for measurement in lcurDeviations[0:sampleSVMSelectionCount]])
                selectedSamplesIndicesFar.extend([measurement[1] for measurement in lcurDeviations[(iLengthCurDeviations-sampleSVMSelectionCount):iLengthCurDeviations]])
        logging.debug("selectedSamplesIndicesClose")
        logging.debug(selectedSamplesIndicesClose)
        logging.debug("selectedSamplesIndicesFar")
        logging.debug(selectedSamplesIndicesFar)
        logging.debug("sampleNames")
        logging.debug(sampleNames)
        #Take indicies and translate to sample names
        #Select close to hyperplane
        if fRunDiscriminant:
            SVMSamples = list()
            for selectedSampleIndex in selectedSamplesIndicesClose:
                SVMSamples.append(sampleNames[selectedSampleIndex])
            dictSelectedSamples[self.c_SVM_CLOSE]=SVMSamples
            logging.debug("SVMSamples Close")
            logging.debug(SVMSamples)
        #Select far from hyperplane
        if fRunDistinct:
            SVMSamples = list()
            for selectedSampleIndex in selectedSamplesIndicesFar:
                SVMSamples.append(sampleNames[selectedSampleIndex])
            dictSelectedSamples[self.c_SVM_FAR]=SVMSamples
            logging.debug("SVMSamples Far")
            logging.debug(SVMSamples)
        return dictSelectedSamples

    #Start micropita selection
    def run(self, fIsAlreadyNormalized, fCladesAreSummed, strOutputFile="MicroPITAOutput.txt", cDelimiter = Constants.TAB, cFeatureNameDelimiter = "|", strInputAbundanceFile=None,
            strUserDefinedTaxaFile=None, strTemporaryDirectory="./TMP", iSampleSelectionCount=0, iSupervisedSampleCount=1,
            strSelectionTechnique=None, strLabel=None, strStratify=None, sMetadataID=None, sLastMetadataName=None, fSumData=True, sFeatureSelectionMethod=None, sCostRange="0"):
	"""
	Writes the selection of samples by method to an output file.
	
	:param	fIsAlreadyNormalized:	Indicates if the abundance table is normalized
	:type	boolean:	boolean indicator if the table is normalized (true= normalized)
	:param	fCladesAreSummed:	Indicates if the abundance table is summed
	:type	boolean:	boolean indicator if the table is summed (true= summed)
	:param	strOutputFile:	File to store selection data
	:type	string:	String file path
	:param	cDelimiter:	Delimiter of abundance table
	:type	char:	Char (default TAB)
	:param	cFeatureNameDelimiter:	Delimiter of the name of features (for instance if they contain consensus lineages indicating clades).
	:type	char:	Char (default |)
	:param	strInputAbundanceFile:	Abundance table data file
	:type	string:	String path to abundance table file
	:param	strUserDefinedTaxaFile:	File containing features to select for
	:type	string:	String path to existing file
	:param	strTemporaryDirectory:	Directory that will be used to store secondary files important to analysis but not the direct deliverable
	:type	string:	String directory path
	:param	iSampleSelectionCount:	Number of samples to select with unsupervised methods
	:type	integer:	integer
	:param	iSupervisedSampleCount:	Number of samples to select with supervised methods
	:type	integer:	integer
	:param	strSelectionTechnique:	List of strings indicating selection techniques
	:type	string:	List of strings each a selection technique
	:param	strLabel:	The metadata used for supervised labels
	:type	string:	String (metadata id)
	:param	strStratify:	The metadata used to stratify unsupervised data
	:type	string:	String (metadata id)
	:param	sMetadataID:	The id of the metadata used as an id for each sample
	:type	string:	String metadata id
	:param	sLastMetadataName:	The id of the metadata positioned last in the abundance table.
	:type	string:	String metadata id
	:param	sFeatureSelectionMethod:	Which method to use to select features in a targeted manner (Using average ranked abundance or abundance)
	:type	string:	String (specific values indicated in microPITA)
        :param	sCostRange:	A comma delimited list of intergers (positive and negative) to be used as the cost range in the supervised methods (in MLPY only the highest is used)
        :type	string:	String (comma delimited integers)	
	"""

        #microPITA object
        microPITA = MicroPITA()

        #SVM parameters
        #Constants associated with the abundance to SVM input file conversion
        c_ABUNDANCE_DELIMITER=Constants.TAB
        c_NORMALIZE_RELATIVE_ABUNDANCY=True
        c_SKIP_FIRST_COLUMN=True
        #Constants associated with the running of the linear SVM
        c_SVM_COST_RANGE = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10"
        c_SVM_PROBABILISTIC = True
        c_SVM_SCALING_LOWER_BOUND = 0

        #Diversity metrics to run
        #TODO make diversity metrics separable
        diversityMetricsAlpha = [microPITA.c_INV_SIMPSON_A_DIVERSITY]
        diversityMetricsAlphaNoNormalize = [microPITA.c_CHAO1_A_DIVERSITY]
        diversityMetricsBeta = [microPITA.c_BRAY_CURTIS_B_DIVERSITY]
        inverseDiversityMetricsBeta = [microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY]
        diversityMetricsBetaNoNormalize = []#microPITA.c_UNIFRAC_B_DIVERSITY]#, microPITA.c_WEIGHTED_UNIFRAC_B_DIVERSITY]
        inverseDiversityMetricsBetaNoNormalize = []

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
                logging.error("MicroPITA.run. No taxa file was given for taxa selection.")
            #Read in taxa list, break down to lines and filter out empty strings
            fhndlTaxaInput = open(strUserDefinedTaxaFile,'r')
            userDefinedTaxa = filter(None,fhndlTaxaInput.read().split(Constants.ENDLINE))
            fhndlTaxaInput.close()

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

        sampleSelectionCount = iSampleSelectionCount
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
        if(not os.path.exists("".join([inputFilePrefix,"-checked",inputFileComponents[1]]))):
            strInputAbundanceFile = AbundanceTable.funcCheckRawDataFile(strReadDataFileName=strInputAbundanceFile, sLastMetadataName=sLastMetadataName, strOutputFileName="".join([inputFilePrefix,"-checked",inputFileComponents[1]]))
        else:
            strInputAbundanceFile = "".join([inputFilePrefix,"-checked",inputFileComponents[1]])

        #Read in abundance data
        #Abundance is a structured array. Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
        #Abundance table object to read in and manage data
        totalAbundanceTable = AbundanceTable.makeFromFile(strInputFile=strInputAbundanceFile, fIsNormalized=fIsAlreadyNormalized, fIsSummed=fCladesAreSummed,
                                   cDelimiter=cDelimiter, sMetadataID=sMetadataID, sLastMetadata=sLastMetadataName, cFeatureNameDelimiter=cFeatureNameDelimiter)
        
        if not totalAbundanceTable:
            logging.error("MicroPITA.run. Could not read abundance table. Stopped.")
            return False

        if fSumData:
            totalAbundanceTable.funcSumClades()
        dictTotalMetadata = totalAbundanceTable.funcGetMetadataCopy()

        #Log metadata keys
        logging.debug(" ".join(["Micropita:run.","Received metadata keys=",str(dictTotalMetadata.keys())]))

        #If there is only 1 unique value for the labels, do not run the Supervised methods
        if len(set(dictTotalMetadata.get(strLabel,[]))) < 2:
            c_RUN_DISCRIMINANT = False
            c_RUN_DISTINCT = False
            logging.error("".join(["The label ",strLabel," did not have 2 or more values. Labels found="]+dictTotalMetadata.get(strLabel,[])))

        logging.debug(" ".join(["Micropita:run.","Received metadata=",str(dictTotalMetadata)]))

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

            #Sample names
            lsSampleNames = stratAbundanceTable.funcGetSampleNames()

            #Only perform if the data is not yet normalized
            if not stratAbundanceTable.funcIsNormalized():

                #Need to first work with unnormalized data
                if((c_RUN_MAX_DIVERSITY_1)or(c_RUN_REPRESENTIVE_DISSIMILARITY_2)or(c_RUN_MAX_DISSIMILARITY_3)):
                    if(c_RUN_MAX_DIVERSITY_1):
                        #Must first generate metrics that do not want normalization before normalization occurs
                        #Expects Observations (Taxa (row) x sample (column))
                        #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                        internalAlphaMatrix = Diversity.buildAlphaMetricsMatrix(tempSampleAbundance = stratAbundanceTable.funcGetAbundanceCopy(), tempSampleNames = lsSampleNames, tempDiversityMetricAlpha = diversityMetricsAlphaNoNormalize)
                        #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
                        #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
                        mostDiverseAlphaSamplesIndexesNoNorm = microPITA.getTopRankedSamples(tempMatrix=internalAlphaMatrix, tempSampleNames=lsSampleNames, tempTopAmount=sampleSelectionCount)
                        #Add to results
                        for index in xrange(0,len(diversityMetricsAlphaNoNormalize)):
                            astrSelectionMethod = microPITA.convertAMetricDiversity[diversityMetricsAlphaNoNormalize[index]]
                            if not astrSelectionMethod in selectedSamples:
                                selectedSamples[astrSelectionMethod]=list()
                            selectedSamples[astrSelectionMethod].extend(mostDiverseAlphaSamplesIndexesNoNorm[index])

                    #Abundance matrix transposed
                    npaTransposedUnnormalizedAbundance = Utility_Math.transposeDataMatrix(stratAbundanceTable.funcGetAbundanceCopy(), tempRemoveAdornments=True)

                    if(c_RUN_REPRESENTIVE_DISSIMILARITY_2):
                        logging.info("Performing representative selection on unnormalized data.")
                        #Run KMedoids with custom distance metric in unnormalized space
                        for bMetric in diversityMetricsBetaNoNormalize:

                            #Get representative dissimilarity samples
                            medoidSamples=microPITA.getCentralSamplesByKMedoids(tempMatrix=npaTransposedUnnormalizedAbundance, tempMetric=bMetric, tempSampleNames=lsSampleNames, tempNumberSamplesReturned=sampleSelectionCount)

                            if(not medoidSamples == False):
                                astrSelectionMethod = microPITA.convertBMetricRepresentative[bMetric]
                                if not astrSelectionMethod in selectedSamples:
                                    selectedSamples[astrSelectionMethod]=list()
                                selectedSamples[astrSelectionMethod].extend(medoidSamples)

                    if(c_RUN_MAX_DISSIMILARITY_3):
                        logging.info("Performing extreme selection on unnormalized data.")
                        #Run HClust with inverse custom distance metric in unnormalized space

                        for bMetric in inverseDiversityMetricsBetaNoNormalize:
                            #Samples for repersentative dissimilarity
                            #This involves inverting the distance metric,
                            #Taking the dendrogram level of where the number cluster == the number of samples to select
                            #Returning a representative sample from each cluster
                            extremeSamples = microPITA.selectExtremeSamplesFromHClust(tempBetaMetric=bMetric, tempAbundanceMatrix=npaTransposedUnnormalizedAbundance, tempSampleNames=lsSampleNames, tempSelectSampleCount=sampleSelectionCount)

                            #Add selected samples
                            if not extremeSamples == False:
                                astrSelectionMethod = microPITA.convertBMetricExtreme[bMetric]
                                if not astrSelectionMethod in selectedSamples:
                                    selectedSamples[astrSelectionMethod]=list()
                                selectedSamples[astrSelectionMethod].extend(extremeSamples)

                logging.info("Selected Samples 1,2,3a")
                logging.info(selectedSamples)

                #Normalize data at this point
                fNormalizeSuccess = stratAbundanceTable.funcNormalize()
                if not fNormalizeSuccess:
                    logging.error("MicroPITA.run. Error occured during normalizing data. Stopped.")
                    return False

            #Generate alpha metrics and get most diverse
            if(c_RUN_MAX_DIVERSITY_1):
                #Get Alpha metrics matrix
                #Expects Observations (Taxa (row) x sample (column))
                #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                internalAlphaMatrix = Diversity.buildAlphaMetricsMatrix(tempSampleAbundance = stratAbundanceTable.funcGetAbundanceCopy(), tempSampleNames = lsSampleNames, tempDiversityMetricAlpha = diversityMetricsAlpha)
                #Get top ranked alpha diversity by most diverse
                #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
                #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
                mostDiverseAlphaSamplesIndexes = microPITA.getTopRankedSamples(tempMatrix=internalAlphaMatrix, tempSampleNames=lsSampleNames, tempTopAmount=sampleSelectionCount)

                #Add to results
                for index in xrange(0,len(diversityMetricsAlpha)):
                    astrSelectionMethod = microPITA.convertAMetricDiversity[diversityMetricsAlpha[index]]
                    if not astrSelectionMethod in selectedSamples:
                        selectedSamples[astrSelectionMethod]=list()
                    selectedSamples[astrSelectionMethod].extend(mostDiverseAlphaSamplesIndexes[index])

            logging.info("Selected Samples 1b")
            logging.info(selectedSamples)

            #Generate beta metrics and 
            if((c_RUN_REPRESENTIVE_DISSIMILARITY_2)or(c_RUN_MAX_DISSIMILARITY_3)):

                #Abundance matrix transposed
                npaTransposedUnnormalizedAbundance = Utility_Math.transposeDataMatrix(stratAbundanceTable.funcGetAbundanceCopy(), tempRemoveAdornments=True)

                #Get center selection using clusters/tiling
                #This will be for beta metrics in normalized space
                if(c_RUN_REPRESENTIVE_DISSIMILARITY_2):
                    logging.info("Performing representative selection on normalized data.")
                    for bMetric in diversityMetricsBeta:

                        #Get representative dissimilarity samples
                        medoidSamples=microPITA.getCentralSamplesByKMedoids(tempMatrix=npaTransposedUnnormalizedAbundance, tempMetric=bMetric, tempSampleNames=lsSampleNames, tempNumberSamplesReturned=sampleSelectionCount)

                        if(not medoidSamples == False):
                            astrSelectionMethod = microPITA.convertBMetricRepresentative[bMetric]
                            if not astrSelectionMethod in selectedSamples:
                                selectedSamples[astrSelectionMethod]=list()
                            selectedSamples[astrSelectionMethod].extend(medoidSamples)

                #Get extreme selection using clusters, tiling
                if(c_RUN_MAX_DISSIMILARITY_3):
                    logging.info("Performing extreme selection on normalized data.")
                    #Run KMedoids with inverse custom distance metric in normalized space
                    for bMetric in inverseDiversityMetricsBeta:

                        #Samples for repersentative dissimilaritystrUserDefinedTaxaFile
                        #This involves inverting the distance metric,
                        #Taking the dendrogram level of where the number cluster == the number of samples to select
                        #Returning a repersentative sample from each cluster
                        extremeSamples = microPITA.selectExtremeSamplesFromHClust(tempBetaMetric=bMetric, tempAbundanceMatrix=npaTransposedUnnormalizedAbundance, tempSampleNames=lsSampleNames, tempSelectSampleCount=sampleSelectionCount)

                        #Add selected samples
                        if(not extremeSamples == False):
                            astrSelectionMethod = microPITA.convertBMetricExtreme[bMetric]
                            if not astrSelectionMethod in selectedSamples:
                                selectedSamples[astrSelectionMethod]=list()
                            selectedSamples[astrSelectionMethod].extend(extremeSamples)

            logging.info("Selected Samples 2,3b")
            logging.info(selectedSamples)

            #Generate selection by the rank average of user defined taxa
            #Expects (Taxa (row) by Samples (column))
            #Expects a column 0 of taxa id that is skipped
            #Returns [(sample name,average,rank)]
            if(c_RUN_RANK_AVERAGE_USER_4):
              if not microPITA.c_USER_RANKED in selectedSamples:
                  selectedSamples[microPITA.c_USER_RANKED]=list()
              selectedSamples[microPITA.c_USER_RANKED].extend(microPITA.selectTargetedTaxaSamples(tempMatrix=stratAbundanceTable, tempTargetedTaxa=userDefinedTaxa, sampleSelectionCount=sampleSelectionCount, sSampleIDName=stratAbundanceTable.funcGetIDMetadataName(),sMethod=sFeatureSelectionMethod))
            logging.info("Selected Samples 4")
            logging.info(selectedSamples)

            #5::Select randomly
            #Expects sampleNames = List of sample names [name, name, name...]
            if(c_RUN_RANDOM_5):

                #Select randomly from sample names
                randomlySelectedSamples = microPITA.getRandomSamples(tempSamples=lsSampleNames, tempNumberOfSamplesToReturn=sampleSelectionCount)
                if not microPITA.c_RANDOM in selectedSamples:
                    selectedSamples[microPITA.c_RANDOM]=list()
                selectedSamples[microPITA.c_RANDOM].extend(list(randomlySelectedSamples))

            logging.info("Selected Samples 5")
            logging.info(selectedSamples)

        #Run supervised methods#
        lStratifiedAbundanceTables = None
        totalAbundanceTable.funcNormalize()
        if(c_RUN_DISTINCT or c_RUN_DISCRIMINANT):
            selectedSamples.update(self.runSupervisedMethods(abundanceTable=totalAbundanceTable,fRunDistinct=c_RUN_DISTINCT, fRunDiscriminant=c_RUN_DISCRIMINANT,
                                   strOuputSVMFile="".join([strTemporaryDirectory,"/",os.path.splitext(os.path.basename(totalAbundanceTable.funcGetName()))[0],"-SVM.txt"]),
                                   strSupervisedMetadata=strLabel, sampleSVMSelectionCount=sampleSVMSelectionCount, sLastMetadataName=sLastMetadataName,
                                   fSkipFirstColumn=c_SKIP_FIRST_COLUMN, fNormalize=c_NORMALIZE_RELATIVE_ABUNDANCY,
                                   iScaleLowestBound=c_SVM_SCALING_LOWER_BOUND, strCostRange=sCostRange, fProbabilitic=c_SVM_PROBABILISTIC))
            logging.info("Selected Samples Unsupervised")
            logging.info(selectedSamples)

        return selectedSamples

    @staticmethod
    def funcWriteSelectionToFile(dictSelection,strOutputFilePath):
	"""
	Writes the selection of samples by method to an output file.
	
	:param	dictSelection:	The dictionary of selections by method to be written to a file.
	:type	Dictionary:	The dictionary of selections by method {"method":["sample selected","sample selected"...]}
	:param	strOutputFilePath:	String path to file to output dictionary.
	:type	String:	String path to file to write to
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
            fHndlOutput.close()
        logging.debug("".join(["Selected samples output to file:",strOutputContent]))

    #Reads in an output selection file from micropita and formats it into a dictionary
    #Dictionary formatted as- {"method":["sample selected","sample selected"...]}
    @staticmethod
    def funcReadSelectionFileToDictionary(strInputFile):
	"""
	Reads in an output selection file from micropita and formats it into a dictionary.
	
	:param	strInputFile:	String path to file to read and translate into a dictionary.
                                {"method":["sample selected","sample selected"...]}
	:type	String:	String path to file to read and translate
	"""

        #Read in selection file
        strSelection = ""
        with open(strInputFile,'r') as fHndlInput:
            strSelection = fHndlInput.read()
        fHndlInput.close()

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
argp.add_argument(Constants_Arguments.c_strCostArgument, dest="sCostRange", metavar= "Cost for Supervised Selection (actual cost = 2^sCostRange).", default="0",
                  help= Constants_Arguments.c_strCostHelp)

#Output
argp.add_argument(Constants_Arguments.c_strTemporaryDirectoryArgument, dest="strTMPDir", metavar = "Temporary Directory", default=None, help = Constants_Arguments.c_genericTMPDirLocationHelp)

#Required
#Input
#Abundance file
argp.add_argument("strFileAbund", metavar = "Abundance file", help = Constants_Arguments.c_strAbundanceFileHelp)
#Outputfile
argp.add_argument("strOutFile", metavar = "Selection Output File", help = Constants_Arguments.c_genericOutputDataFileHelp)
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


    dictSelectedSamples = microPITA.run(fIsAlreadyNormalized=args.fIsNormalized, fCladesAreSummed=args.fIsSummed, strOutputFile=args.strOutFile,
                                        cDelimiter=args.cFileDelimiter, cFeatureNameDelimiter = args.cFeatureNameDelimiter, strInputAbundanceFile=args.strFileAbund,
                                        strUserDefinedTaxaFile=args.strFileTaxa, strTemporaryDirectory=args.strTMPDir, iSampleSelectionCount=args.iUnsupervisedSelectionCount,
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



