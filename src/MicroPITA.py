#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Run analysis for the microPITA paper
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from AbundanceTable import AbundanceTable
import argparse
from Constants import Constants
from Diversity import Diversity
from FileIO import FileIO
import logging
import math
import mlpy
from MLPYDistanceAdaptor import MLPYDistanceAdaptor
import numpy as np
import operator
import os
import random
import re
import scipy.cluster.hierarchy as hcluster
from SVM import SVM
import sys
from ValidateData import ValidateData


class MicroPITA:
    """
    Micropita class
    """
    #Constants
    #Controls if the bias correcting Chao1 is used (false = not used, uses intead original non bias correcting formula)
    c_CHAO_CORRECT_BIAS = False

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

    #Converts ecology metrics into standardized method selection names
    convertAMetricDiversity = {c_INV_SIMPSON_A_DIVERSITY:c_DIVERSITY_1, c_CHAO1_A_DIVERSITY:c_DIVERSITY_2}
    convertBMetricRepresentative = {c_BRAY_CURTIS_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_1, c_UNIFRAC_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_2, c_WEIGHTED_UNIFRAC_B_DIVERSITY:c_REPRESENTATIVE_DISSIMILARITY_3}
    convertBMetricExtreme = {c_INVERSE_BRAY_CURTIS_B_DIVERSITY:c_EXTREME_DISSIMILARITY_1, c_INVERSE_UNIFRAC_B_DIVERSITY:c_EXTREME_DISSIMILARITY_2, c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY:c_EXTREME_DISSIMILARITY_3}

    #Used to store scores and pvalues
    c_SCORE_INDEX = 0
    c_PVALUE_INDEX = 1

    #Linkage used in the Hierarchical clustering
    c_HIERARCHICAL_CLUSTER_METHOD = 'average'

####General
    #Happy path tested
    #Expectes a structured array for npData (rows = Taxa/OTU)
    #Expectes the first entry of every row to be an id that is ignored but carried forward
    #Metadata is used to collapse by
    def funcStratifyDataByMetadata(self,lsMetadata, npData):
        dictAbundanceBlocks = dict()
        setValues = set(lsMetadata)
        lsNames = npData.dtype.names
        #Get index of values to break up
        for value in setValues:
            fDataIndex = [sData==value for sData in lsMetadata]
            #The true is added to keep the first column which should be the feature id
            dictAbundanceBlocks[value] = npData[np.compress([True]+fDataIndex,lsNames)]
        return dictAbundanceBlocks

    def funcStratifyMetadataByMetadata(self, lsMetadata, dictMetadataToStratify):
        dictMetadataBlocks = dict()
        setValues = set(lsMetadata)
        #Get index of values to break up
        for value in setValues:
            fDataIndex = [sData==value for sData in lsMetadata]
            dictBrokenMetadata = dict()
            for metadataType in dictMetadataToStratify:
                dictValues = dictMetadataToStratify[metadataType]
                dictBrokenMetadata[metadataType] = np.compress(fDataIndex,dictValues).tolist()
            #The true is added to keep the first column which should be the feature id
            dictMetadataBlocks[value] = dictBrokenMetadata
        return dictMetadataBlocks

####Group 1## Diversity

    #Testing: Happy Path Tested (4)
    #TODO Need to figure out how to combine the non normalized and normalized metric values going in and going out of metric creation
    #Get alpha abundance of the metric for the vector
    #@params tempAbundancies List of values to compute diversity
    #@params tempMetric Alpha metric to use to define diversity
    #@return float
    def getAlphaMetric(self, tempAbundancies=None, tempMetric=None):
        if(not ValidateData.isValidString(tempMetric)):
            return False
        elif(tempMetric == self.c_SHANNON_A_DIVERSITY):
            return Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == self.c_SIMPSON_A_DIVERSITY):
            return Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == self.c_INV_SIMPSON_A_DIVERSITY):
            return Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        #Needs NOT Normalized Abundance
        elif(tempMetric == self.c_CHAO1_A_DIVERSITY):
            return Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies=tempAbundancies, tempCorrectForBias=self.c_CHAO_CORRECT_BIAS)
        else:
            return False

    #Testing: Happy path Testing (3)
    #Build a matrix of alpha diversity metrics for each sample
    #Row = metric, column = sample
    #@params tempSampleAbundance Observations (Taxa (row) x sample (column))
    #@params tempSampleNames List of sample names of samples to measure (do not include the taxa id  column name or other column names which should not be read)
    #@params tempDiversityMetricAlpha List of diversity metrics to use in measuring
    #@return A lists of lists. Each internal list is a list of (floats) indicating a specific metric measurement method measuring multiple samples
    #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
    def buildAlphaMetricsMatrix(self, tempSampleAbundance = None, tempSampleNames = None, tempDiversityMetricAlpha = None):
        #Create return
        returnMetricsMatrix = []
        [returnMetricsMatrix.append(list()) for index in tempDiversityMetricAlpha]

        #Get amount of metrics
        metricsCount = len(tempDiversityMetricAlpha)

        #For each sample get all metrics
        #Place in list of lists
        #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
        for sample in tempSampleNames:
            sampleAbundance = tempSampleAbundance[sample]
            for metricIndex in xrange(0,metricsCount):
                returnMetricsMatrix[metricIndex].append(self.getAlphaMetric(tempAbundancies = sampleAbundance, tempMetric = tempDiversityMetricAlpha[metricIndex]))
        return returnMetricsMatrix


    #Testing: Happy path Testing (8)
    #Returns the indices of the top numbers
    #@params tempMatrix List of lists
    #[[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
    #@params tempSampleNames Sample Names to select from, if no sample names are given, indices are returned
    #@params tempTopAmount The amount of top ranked samples to return
    #@return list of lists which contain indices or samples names of the top ranked samples
    #[[indexOfSample1, indexOfSample2, indexOfSampleN],[indexOfSample1, indexOfSample2, indexOfSampleN]]
    #[[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
    def getTopRankedSamples(self, tempMatrix = None, tempSampleNames = None, tempTopAmount = None):
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
    #Generate a beta metric matrix between abundancies for the given beta metrics
    #@params tempAbundancies Abundancies to be measured. Matrix should be where row = samples and columns = organisms
    #@params tempMetric String metric to be returned
    #@returns
    def getBetaMetric(self, tempAbundancies=None, tempMetric=None):
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
    #Samples are rows, taxa are columns
    #@params tempMatrix Data matrix taxa (rows) by samples (columns).
    #Note: Not a structured array, a matrix stripped of adornment.
    #@params tempSampleNames The names of the samples
    #@params tempNumberClusters The K number of clusters to divide the space
    #@params tempNumberSamplesReturned The number of samples needing to be returned (Almost completely ignored currently
    #TODO stop ignoring tempNumberSamplesReturned
    #otherwise a distance object is expected and will be used.
    def getCentralSamplesByKMedoids(self, tempMatrix=None, tempMetric=None, tempSampleNames=None, tempNumberClusters=0, tempNumberSamplesReturned=0):
        #Validate parameters
        if(tempNumberClusters > tempNumberSamplesReturned):
            logging.error("MicroPITA.getCentralSamplesByKMedoids. Number of clusters should be equal to or less than the number samples returned. We will not represent a cluster otherwise.")
            return False

        #Count of how many rows
        sampleCount = tempMatrix.shape[0]
        if(tempNumberClusters > sampleCount):
            logging.error("".join(["MicroPITA.getCentralSamplesByKMedoids. There are not enough samples to make that many clusters. Cluster number = ",str(tempNumberClusters),". Sample number = ",str(sampleCount),"."]))
            return False

        #Samples to return
        returningSamples = list()

        #If the cluster count is equal to the sample count return all samples
        if(sampleCount == tempNumberClusters):
            returningSamples = dict()
            for index in xrange(1,(tempNumberClusters+1)):
                returningSamples[str(index)] = [tempSampleNames[(index-1)]]
            return returningSamples

        #Get distance matrix
        distanceMatrix=self.getBetaMetric(tempAbundancies=tempMatrix, tempMetric=tempMetric)
        if(ValidateData.isFalse(distanceMatrix)):
          logging.error("MicroPITA.getCentralSamplesByKMedoids. Received false for betaMetrix matrix generation, returning false.")
          return False

        #Log distance matrix
        logging.debug("".join(["Distance matrix for representative selection using metric=",str(tempMetric)]))

        if(( tempMetric==Diversity.c_UNIFRAC_B_DIVERSITY ) or ( tempMetric==Diversity.c_WEIGHTED_UNIFRAC_B_DIVERSITY )):
          distanceMatrix = distanceMatrix['distance_matrix'][0]
        #TODO make sure you are getting condensed for unifrac and braycurtis
        #Mane the adaptor for the MLPY Kmediods method to use custom distanc matrices
        distance = MLPYDistanceAdaptor(tempDistanceMatrix=distanceMatrix, tempIsCondensedMatrix=True)

        #Create object to determine clusters/medoids
        medoidsMaker = mlpy.Kmedoids(k=tempNumberClusters, dist=distance)

        #medoidsData includes(1d numpy array, medoids indexes; 
        #              1d numpy array, non-medoids indexes;
        #              1d numpy array, cluster membership for non-medoids;
        #              double, cost of configuration)
        #tempMatrix is samples x rows
        #Build a matrix of lists of indicies to pass to the distance matrix
        indicesMatrix = []
        for indexPosition in xrange(0,len(tempMatrix)):
            indicesMatrix.append([indexPosition])
        medoidsData = medoidsMaker.compute(np.array(indicesMatrix))
        logging.debug("Results from the kmedoid method in representative selection:")
        logging.debug(str(medoidsData))

        #If returning the same amount of clusters and samples
        #Return centroids
        if(tempNumberClusters == tempNumberSamplesReturned):
            selectedIndexes = medoidsData[0]
            for index in xrange(0,tempNumberClusters):
                returningSamples.append(tempSampleNames[selectedIndexes[index]])
            return returningSamples
        return False

####Group 3## Highest Dissimilarity
    #Testing: Needs testing
    #Select a given amount of representative samples, one per cluster
    #from a hierarchical cluster given a beta metric.
    #@params tempBetaMetric Betric Metric to use for distance matrix generation
    #@params tempAbundanceMatrix Abundance matrix (sample (row) x taxa (column))
    #@params tempSelectSampleCount Select
    #@params tempNormalze Should be None or a number to normalize by
    def selectExtremeSamplesFromHClust(self, tempBetaMetric, tempAbundanceMatrix, tempSampleNames, tempSelectSampleCount):

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

    #Testing: Happy Path Testing
    #Reduce by taxa
    #@params tempMatrix Structured array [("taxa label",numeric,numeric,numeric),...]
    #@params tempTaxa List of string which are the taxa labels
    #@return Structured array of just the given taxa
    def reduceToTaxa(self, tempMatrix = None, tempTaxa = None):
        #Validate parameters

        #To indicate how to reduce rows
        reduceIndex = []
        for row in tempMatrix:
            taxa = row[0]
            reduceIndex.append(taxa in tempTaxa)

        #Reduce
        return np.compress(reduceIndex, tempMatrix, axis = 0)

    #Returns Rank Average of samples given abundancy in input taxa
    #Expects (Taxa (row) by Samples (column))
    #Expects a column 0 of taxa id that is skipped
    #Allows ties
    #@param tempMatrix [taxa (row) x sample(column)] Matrix with first column as a taxa id column which is ignored
    #@return [(sample name,average,rank)]
#    def getRankAverageSamples(self, tempMatrix = None):
#        #Validate Matrix
#
#        #List to return with [(sample name,average,rank)]
#        rankAverageSamplesReturn = []
#        
#        #Get sample names
#        sampleNames = tempMatrix.dtype.names[1:]
#
#        #For each sample name get average and add to return list
#        for name in sampleNames:
#            rankAverageSamplesReturn.append([name,np.average(tempMatrix[name]),-1])
#
#        #Sort based on average
#        rankAverageSamplesReturn = sorted(rankAverageSamplesReturn, key = lambda sampleData: sampleData[1], reverse = True)
#
#        #Add ranks
#        rank = 1
#        currentValue = rankAverageSamplesReturn[0][1]
#        for sampleData in rankAverageSamplesReturn:
#            sampleAverage = sampleData[1]
#            #Error samples are out of order
#            if(sampleAverage > currentValue):
#                return False
#            #Allow ties
#            if(sampleAverage == currentValue):
#                sampleData[2] = rank
#            else:
#                #Update/set rank and value
#                currentValue = sampleAverage
#                rank = rank + 1
#                sampleData[2] = rank
#
#        #return
#        return rankAverageSamplesReturn

#TODO TEST AND CHECK
    #Ranks the taxa by abundance and then averages thier ranks in the samples
    #Expects (Taxa (row) by Samples (column))
    #Expects a column 0 of taxa id that is skipped
    #Allows ties
    #@param tempMatrix [taxa (row) x sample(column)] Matrix with first column as a taxa id column which is ignored
    #@param tempTargetedTaxa list of string names of taxa which are measured after ranking against the full sample.
    #@return [(sample name,average rank)]
    def getAverageRanksSamples(self, tempMatrix, tempTargetedTaxa, sSampleIDName):

        #Sample rank averages [[sample,average rank of selected taxa]]
        #Returned
        sampleRankAverages = []
        
        #Get sample names
        sampleNames = tempMatrix.dtype.names[1:]
        #Get taxa names
        allTaxaNames = tempMatrix[tempMatrix.dtype.names[0]]

        #For each sample name get the ranks
        for name in sampleNames:
            #Lists of taxa with the following information [[taxa name,value,rank]]
            ranks = []

            #For each taxa
            for taxaIDIndex in xrange(0,len(allTaxaNames)):
                currentAbundance = tempMatrix[name][taxaIDIndex]
                ranks.append([tempMatrix[sSampleIDName][taxaIDIndex],currentAbundance,-1])

            #Sort based on abundance
            ranks = sorted(ranks, key = lambda sampleData: sampleData[1], reverse = True)

            #Add ranks
            rank = 1
            currentValue = ranks[0][1]
            for sampleData in ranks:
                sampleAbundance = sampleData[1]
                #Error samples are out of order
                if(sampleAbundance > currentValue):
                    return False
                #Allow ties
                if(sampleAbundance == currentValue):
                    sampleData[2] = rank
                else:
                    #Update/set rank and value
                    currentValue = sampleAbundance
                    rank = rank + 1
                    sampleData[2] = rank

            #Get average ranks
            sumRank = 0
            countRank = 0.0
            for rankData in ranks:
               if(rankData[0] in tempTargetedTaxa):
                   sumRank = sumRank + rankData[2] 
                   countRank = countRank + 1.0

            #Save average rank by sample name
            if(not countRank == 0):
                sampleRankAverages.append([name,sumRank/countRank])
                print([name,sumRank/countRank])
            else:
                logging.error("".join(["MicroPITA.getAverageRanksSamples. Found no taxa for sample=",str(name)]))

        #Sort based on average
        sampleRankAverages = sorted(sampleRankAverages, key = lambda sampleData: sampleData[1], reverse = True)
            
        #return
        return sampleRankAverages

    def selectTargetedTaxaSamples(self, tempMatrix, tempTargetedTaxa, sampleSelectionCount, sSampleIDName):
      if(len(tempTargetedTaxa) < 1):
        logging.error("MicroPITA.getAverageRanksSamples. Taxa defined selection was requested but no taxa were given.")
      #Rank the samples
      userRankedSamples = self.getAverageRanksSamples(tempMatrix=tempMatrix, tempTargetedTaxa=tempTargetedTaxa, sSampleIDName=sSampleIDName)

      #Select the top samples
#      topRankedSamples = userRankedSamples[0:(sampleSelectionCount-1):]
      topRankedSamples = userRankedSamples[(sampleSelectionCount*-1):]
      topRankedSamplesNames = np.compress([True,False],topRankedSamples,axis=1)
      return [item for sublist in topRankedSamplesNames for item in sublist]

####Group 5## Random

    #Returns random sample names of the number given
    #No replacement
    #@params tempSamples list of sample names [name, name, name...]
    #@params tempNumberOfSamplesToReturn = Number of samples to return 1 to count(samples)
    #@return A list of sample name lists with each list having a randomly selected 
    #amount of sample names as defined by tempNumberOfSamplesToReturn [[randomName1, randomName2,...],[randomName1, randomName2,...]]
    def getRandomSamples(self, tempSamples=None, tempNumberOfSamplesToReturn=0):
        #Validate Matrix

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
        
    #Adapted from the easy.py script included in the standard Libsvm install
    #Runs a SVC linear model.
    #@params tempInputFile String File path to Qiime-like output of abundance data to be converted to SVM format
    #@params tempDelimiter Delimiter to use to parse the input file
    #@params tempOutputSVMFile String File path and name to output in SVM format based on the input file
    #@params tempMatrixLabels List of labels (does not have to be strings, just have to be appropriate for labels when casted to string type
    #@params tempFirstDataRow Integer First row to read (skips header rows) (0-based)
    #@params tempSkipFirstColumn Boolean True indicates the first column will be skipped (due to it containing row identifying information like OTU names).
    #@params tempNormalize Boolean True indicates normalizes to relative abundancy per sample (column)
    #@return Dictionary of file paths which were generated by the model and prediction steps
    def runSVM(self, tempInputFile=None, tempDelimiter=Constants.TAB, tempOutputSVMFile=None, tempMatrixLabels=None, tempFirstDataRow=2, tempSkipFirstColumn=True, tempNormalize=True, tempSVMScaleLowestBound = 0, tempSVMLogG="-5,-4,-3,-2,-1,0,1,2,3,4,5", tempSVMLogC="-5,-4,-3,-2,-1,0,1,2,3,4,5", tempSVMProbabilistic=True):
        #Validate data

        #Create SVM object
        svm = SVM()

        #Holds files generated by SVM code
        #Files generated by modeling
        modelFiles = None
        #Files generated by prediction
        predictionFiles = None

        #Convert metadata to labels
        #Create a dictionary converting labels to indexes
        metadataLabelToIndex = dict()
        indexCount = 0
        for metadataLabel in tempMatrixLabels:
            if(not metadataLabel in metadataLabelToIndex):
                metadataLabelToIndex[metadataLabel] = indexCount
                indexCount += 1
        #Create a new label list but coded as integers
        metadataLabelsAsIntegerCodes = list()
        for data in tempMatrixLabels:
            metadataLabelsAsIntegerCodes.append(metadataLabelToIndex[data])

        #Convert abundancies file to SVM file
        noError = svm.convertAbundanceFileToSVMFile(tempInputFile=tempInputFile, tempOutputSVMFile=tempOutputSVMFile, tempDelimiter=tempDelimiter, tempLabels=metadataLabelsAsIntegerCodes, tempFirstDataRow=tempFirstDataRow, tempSkipFirstColumn=tempSkipFirstColumn, tempNormalize=tempNormalize)

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

    #Start micropita selection
    def run(self, strOutputFile="MicroPITAOutput.txt", strInputAbundanceFile=None, strUserDefinedTaxaFile=None, strTemporaryDirectory="./TMP", iSampleSelectionCount=0, iSupervisedSampleCount=1, strSelectionTechnique=None, strLabel=None, strStratify=None, iSampleNameRow=0, iFirstDataRow=1):
        #microPITA object0
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

        #Abundance table object to read in and manage data
        totalData = AbundanceTable()

        #Input file path components
        inputFileComponents = os.path.split(strInputAbundanceFile)
        inputFilePrefix = inputFileComponents[0]

        #Read in abundance data
        #Abundance is a structured array. Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
        rawAbundance,metadata = totalData.textToStructuredArray(tempInputFile=strInputAbundanceFile, tempDelimiter=Constants.TAB, tempNameRow=iSampleNameRow, tempFirstDataRow=iFirstDataRow, tempNormalize=False)

        logging.debug(" ".join(["Micropita:run.","Received metadata=",str(metadata)]))
        print(" ".join(["Micropita:run.","Received metadata=",str(metadata)]))

        sampleSelectionCount = iSampleSelectionCount
        sampleSVMSelectionCount = iSupervisedSampleCount
        #TODO
        clusterCount = iSampleSelectionCount

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
        #TODO remove this check
        if(not os.path.exists(inputFilePrefix+"-checked.txt")):
            inputFile = totalData.checkRawDataFile(strInputAbundanceFile)
        else:
            inputFile = inputFilePrefix+"-checked.txt"

        #Stratify the data if need be
        dictAbundanceBlocks = dict()
        if (strStratify == None) or (strStratify.lower() == "none"):
          dictAbundanceBlocks["Total"]=rawAbundance
        elif strStratify in metadata:
          lsMetadata = metadata[strStratify]
          dictAbundanceBlocks = microPITA.funcStratifyDataByMetadata(lsMetadata=lsMetadata, npData=rawAbundance)
          dictMetadataBlocks = microPITA.funcStratifyMetadataByMetadata(lsMetadata=lsMetadata, dictMetadataToStratify=metadata)
          #Write to file to record
          lsOutputPathElements = os.path.split(strOutputFile)
          for strStratifiedGroup in dictAbundanceBlocks:
              strMetadataHeader = ""
              strOutputBase = filter(None,[lsSample for lsSample in lsOutputPathElements[0].split("/")])[-1]
              strOutputFileName = "".join([lsOutputPathElements[0],"/",strOutputBase,"-StratBy-",strStratifiedGroup,".",lsOutputPathElements[1]])
              #Write metadata
              dictCurMetadata = dictMetadataBlocks[strStratifiedGroup]
              for sMetaKey in dictCurMetadata:
                  strMetadataHeader = strMetadataHeader + Constants.TAB.join([sMetaKey]+dictCurMetadata[sMetaKey])+Constants.ENDLINE
                  curAbundance = dictAbundanceBlocks[strStratifiedGroup]
                  lsCurNames = curAbundance.dtype.names
                  curAbundance = curAbundance.tolist()
              with open(strOutputFileName, 'w') as f:
                  f.write(Constants.TAB.join(lsCurNames)+Constants.ENDLINE)
                  f.write(strMetadataHeader)
                  lsOutput = list()
                  for curAbundanceRow in curAbundance:
                      lsOutput.append(Constants.TAB.join([str(curAbundanceElement) for curAbundanceElement in curAbundanceRow]))
                  f.write(Constants.ENDLINE.join(lsOutput))
                  f.close()

        #For each stratified abundance block or for the unstratfified abundance
        #Run the unsupervised blocks
        for abundanceBlockValue in dictAbundanceBlocks:
            print("Running abundance block:"+abundanceBlockValue)
            abundance = dictAbundanceBlocks[abundanceBlockValue]
            #Get sample names excluding the taxa id column name
            sampleNames = abundance.dtype.names[1:]
            sampleID = abundance.dtype.names[0]

            #Need to first work with unnormalized data
            if((c_RUN_MAX_DIVERSITY_1)or(c_RUN_REPRESENTIVE_DISSIMILARITY_2)or(c_RUN_MAX_DISSIMILARITY_3)):
                if(c_RUN_MAX_DIVERSITY_1):
                    #Must first generate metrics that do not want normalization before normalization occurs
                    #Expects Observations (Taxa (row) x sample (column))
                    #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                    internalAlphaMatrix = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = diversityMetricsAlphaNoNormalize)
                    #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
                    #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
                    mostDiverseAlphaSamplesIndexesNoNorm = microPITA.getTopRankedSamples(tempMatrix=internalAlphaMatrix, tempSampleNames=sampleNames, tempTopAmount=sampleSelectionCount)
                    #Add to results
                    for index in xrange(0,len(diversityMetricsAlphaNoNormalize)):
                        astrSelectionMethod = microPITA.convertAMetricDiversity[diversityMetricsAlphaNoNormalize[index]]
                        if not astrSelectionMethod in selectedSamples:
                            selectedSamples[astrSelectionMethod]=list()
                        selectedSamples[astrSelectionMethod].extend(mostDiverseAlphaSamplesIndexesNoNorm[index])

                if(c_RUN_REPRESENTIVE_DISSIMILARITY_2):
                    logging.info("Performing representative selection on unnormalized data.")
                    #Run KMedoids with custom distance metric in unnormalized space
                    for bMetric in diversityMetricsBetaNoNormalize:

                        #Get representative dissimilarity samples
                        transposedAbundance = totalData.transposeDataMatrix(abundance, tempRemoveAdornments=True)
                        medoidSamples=microPITA.getCentralSamplesByKMedoids(tempMatrix=transposedAbundance, tempMetric=bMetric, tempSampleNames=sampleNames, tempNumberClusters=clusterCount, tempNumberSamplesReturned=sampleSelectionCount)

                        if(not medoidSamples == False):
                            astrSelectionMethod = microPITA.convertBMetricRepresentative[bMetric]
                            if not astrSelectionMethod in selectedSamples:
                                selectedSamples[astrSelectionMethod]=list()
                            selectedSamples[astrSelectionMethod].extend(medoidSamples)

                if(c_RUN_MAX_DISSIMILARITY_3):
                    logging.info("Performing extreme selection on unnormalized data.")
                    #Run HClust with inverse custom distance metric in unnormalized space
                    #TODO centralize all transpose needs
                    #TODO transpose is being performed
                    tAbundance = totalData.transposeDataMatrix(tempMatrix=abundance, tempRemoveAdornments=True)

                    for bMetric in inverseDiversityMetricsBetaNoNormalize:
                        #Samples for repersentative dissimilarity
                        #This involves inverting the distance metric,
                        #Taking the dendrogram level of where the number cluster == the number of samples to select
                        #Returning a representative sample from each cluster
                        extremeSamples = microPITA.selectExtremeSamplesFromHClust(tempBetaMetric=bMetric, tempAbundanceMatrix=tAbundance, tempSampleNames=sampleNames, tempSelectSampleCount=sampleSelectionCount)

                        #Add selected samples
                        if not extremeSamples == False:
                            astrSelectionMethod = microPITA.convertBMetricExtreme[bMetric]
                            if not astrSelectionMethod in selectedSamples:
                                selectedSamples[astrSelectionMethod]=list()
                            selectedSamples[astrSelectionMethod].extend(extremeSamples)

            logging.info("Selected Samples 1,2,3a")
            logging.info(selectedSamples)

            #Normalize data at this point
            abundance = totalData.normalizeColumns(tempStructuredArray=abundance, tempColumns=list(sampleNames))
            if(abundance == False):
                logging.error("MicroPITA.run. Error occured during normalizing data. Stopped.")
                return False

            #Generate alpha metrics and get most diverse
            if(c_RUN_MAX_DIVERSITY_1):
                #Get Alpha metrics matrix
                #Expects Observations (Taxa (row) x sample (column))
                #Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
                internalAlphaMatrix = microPITA.buildAlphaMetricsMatrix(tempSampleAbundance = abundance, tempSampleNames = sampleNames, tempDiversityMetricAlpha = diversityMetricsAlpha)
                #Get top ranked alpha diversity by most diverse
                #Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
                #Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
                mostDiverseAlphaSamplesIndexes = microPITA.getTopRankedSamples(tempMatrix=internalAlphaMatrix, tempSampleNames=sampleNames, tempTopAmount=sampleSelectionCount)

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

                #Transpose data
                transposedAbundance = totalData.transposeDataMatrix(abundance, tempRemoveAdornments=True)

                #Get center selection using clusters/tiling
                #This will be for beta metrics in normalized space
                if(c_RUN_REPRESENTIVE_DISSIMILARITY_2):
                    logging.info("Performing representative selection on normalized data.")
                    for bMetric in diversityMetricsBeta:

                        #Get representative dissimilarity samples
                        medoidSamples=microPITA.getCentralSamplesByKMedoids(tempMatrix=transposedAbundance, tempMetric=bMetric, tempSampleNames=sampleNames, tempNumberClusters=clusterCount, tempNumberSamplesReturned=sampleSelectionCount)

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
                        extremeSamples = microPITA.selectExtremeSamplesFromHClust(tempBetaMetric=bMetric, tempAbundanceMatrix=transposedAbundance, tempSampleNames=sampleNames, tempSelectSampleCount=sampleSelectionCount)

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
              selectedSamples[microPITA.c_USER_RANKED].extend(microPITA.selectTargetedTaxaSamples(tempMatrix=abundance, tempTargetedTaxa=userDefinedTaxa, sampleSelectionCount=sampleSelectionCount, sSampleIDName=sampleID))
            logging.info("Selected Samples 4")
            logging.info(selectedSamples)

            #5::Select randomly
            #Expects sampleNames = List of sample names [name, name, name...]
            if(c_RUN_RANDOM_5):

                #Select randomly from sample names
                randomlySelectedSamples = microPITA.getRandomSamples(tempSamples=sampleNames, tempNumberOfSamplesToReturn=sampleSelectionCount)
                if not microPITA.c_RANDOM in selectedSamples:
                    selectedSamples[microPITA.c_RANDOM]=list()
                selectedSamples[microPITA.c_RANDOM].extend(list(randomlySelectedSamples))

            logging.info("Selected Samples 5")
            logging.info(selectedSamples)

        #Run supervised blocks
        #Select supervised (using SVM)
        #Expects input file's matrix to be Taxa (row) by Sample (col) with a taxa id column (index=0)
        if(c_RUN_DISTINCT or c_RUN_DISCRIMINANT):
            #Get file name without extention
            strTail = os.path.split(strInputAbundanceFile)[1]
            #Get sample names
            sampleNames = rawAbundance.dtype.names[1:]
            #Run linear SVM
            svmRelatedData = microPITA.runSVM(tempInputFile=strInputAbundanceFile, tempDelimiter=c_ABUNDANCE_DELIMITER, tempOutputSVMFile="".join([strTemporaryDirectory,"/",os.path.splitext(strTail)[0],"-SVM.txt"]), tempMatrixLabels=metadata[strLabel], tempFirstDataRow=iFirstDataRow, tempSkipFirstColumn=c_SKIP_FIRST_COLUMN, tempNormalize=c_NORMALIZE_RELATIVE_ABUNDANCY, tempSVMScaleLowestBound=c_SVM_SCALING_LOWER_BOUND, tempSVMLogC=c_SVM_COST_RANGE, tempSVMProbabilistic=c_SVM_PROBABILISTIC)
            #Read in prediction file and select samples
            #Selecting samples most ambiguous to all hyperplanes
            if(not svmRelatedData == False):
                #TODO make the constants static so this doesnt happen
                if(SVM().c_KEYWORD_PREDICTION_FILE in svmRelatedData):
                    predictionFile = svmRelatedData[SVM().c_KEYWORD_PREDICTION_FILE]
                    
                    #TODO test selection of samples
                    #Read in file
                    readIn = FileIO(predictionFile,True,False,False)
                    predictionLists = readIn.readFullFile()
                    readIn.close()
                    predictionLists = predictionLists.split(Constants.ENDLINE)
                    #Get label count (meaning the number of label categories)
                    labelCount = len(predictionLists[0].split(Constants.WHITE_SPACE))-1
                    #Central probability
                    centralProbability = 1.0 / float(labelCount)
                    #Create and array to hold difference of the samples probabilities from the central probability
                    centralDeviation = dict()

                    #For each line in the file subtract all but the first item from the central value
                    #The first item being the header line not a row of probabilities
                    #Selected_label prob_for_label1 prob_ForLabel2....
                    for lineIndex in xrange(1,len(predictionLists)-1):
                        #Split line into elements by whitespace and remove first element (the label)
                        lineElements = predictionLists[lineIndex]
                        lineElements = lineElements.split(Constants.WHITE_SPACE)
                        iCurLabel = str(lineElements[0])
                        lineElements = lineElements[1:]

                        #Sum the absolute values of the differences
                        #Store the index and the summed absolute differences
                        deviation = 0
                        for prediction in lineElements:
                            deviation += math.pow((float(prediction)-centralProbability),2)
                        if(not iCurLabel in centralDeviation):
                            centralDeviation[iCurLabel] = []
                        curSVMData = centralDeviation[iCurLabel]
                        curSVMData.append([deviation,lineIndex,iCurLabel])
                        centralDeviation[iCurLabel] = curSVMData

                    #Sort sample by summed absolute deviations from the center and take the top N indexes
                    for scurKey in centralDeviation:
                        lcurLabeSamples = centralDeviation[scurKey]
                        centralDeviation[scurKey] = sorted(lcurLabeSamples, key=operator.itemgetter(0))

                    selectedSamplesIndicesClose = list()
                    selectedSamplesIndicesFar = list()
                    #If the amount of samples needed are greater than what was analyzed with the SVM,
                    #Return them all for both far and near.
                    #Get the samples closes to hyperplanestrUserDefinedTaxaFile
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
                    #Take indicies and translate to sample names
                    #TODO A little scary, requires no one mess with the sampleNames...maybe a different way to do this?
                    #Select close to hyperplane
                    if(c_RUN_DISCRIMINANT):
                        SVMSamples = list()
                        for selectedSampleIndex in selectedSamplesIndicesClose:
                            SVMSamples.append(sampleNames[selectedSampleIndex-1])
                        selectedSamples[microPITA.c_SVM_CLOSE]=SVMSamples
                    #Select far from hyperplane
                    if(c_RUN_DISTINCT):
                        SVMSamples = list()
                        for selectedSampleIndex in selectedSamplesIndicesFar:
                            SVMSamples.append(sampleNames[selectedSampleIndex-1])
                        selectedSamples[microPITA.c_SVM_FAR]=SVMSamples

        logging.info("Selected Samples 6")
        logging.info(selectedSamples)
        return selectedSamples

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicroPITA.py", 
    description = """Selects samples from abundance tables based on various selection schemes.""" )

#Arguments
#Optional parameters
#Logging
argp.add_argument("-l", dest="strLogLevel", metavar= "Loglevel", default="INFO", 
                  choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], 
                  help= "Logging level which will be logged to a .log file with the same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
#Abundance associated
argp.add_argument("-n", dest="iSampleNameRow", metavar= "SampleNameRow", default=0, 
                  help= "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering.")
argp.add_argument("-d", dest="iFirstDataRow", metavar= "FirstDataRow", default=1, 
                  help= "The row in the abundance file that is the first row to contain abundance data. This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata.")
argp.add_argument("-s", dest="iSupervisedCount", metavar= "CountSupervisedSamplesSelected", default=1, 
                  help= "The count of labeled data to select per label (default =1)")
#SVM label
#Label parameter to be used with SVM
argp.add_argument("-p", dest="strLabel", metavar= "Label", default="", 
                  help= "The name of the phenotype data row on which to perform supervised methods.")
argp.add_argument("-u", dest="strUnsupervisedStratify", metavar= "UnsupervisedStratify", default=None, 
                  help= "The metatdata to stratify unsupervised analysis.")
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = "An optional output file" )
#Abundance file
argp.add_argument( "strFileAbund", metavar = "Abundance_file", help = "An abundance table." )
#Taxa file
argp.add_argument( "strFileTaxa", metavar = "Taxa_file",
    help = "A file containing taxa to be used in taxa-directed selection." )
#Temporary folder
argp.add_argument( "strTMPDir", metavar = "Temporary_Directory", help = "Directory to place temporary and intermediate files.")
#Select count
argp.add_argument( "icount", metavar = "number", type = int, help = "The number of samples to select (An integer greater than 0.)." )
#Selection parameter
argp.add_argument("strSelection", metavar = "Selection_Methods", help = "Select techniques listed one after another.", nargs="*")

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )

    #Massage input
    if args.strFileTaxa == "None":
      args.strFileTaxa = None

    #Set up logger
    iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
    if not isinstance(iLogLevel, int):
        raise ValueError('Invalid log level: %s. Try DEBUG, INFO, WARNING, ERROR, or CRITICAL.' % strLogLevel)
    logging.basicConfig(filename="".join([os.path.splitext(args.strOutFile)[0],".log"]), filemode = 'w', level=iLogLevel)

    logging.info("Start microPITA")
    dictSelectedSamples = MicroPITA().run(strOutputFile=args.strOutFile, strInputAbundanceFile=args.strFileAbund, strUserDefinedTaxaFile=args.strFileTaxa, strTemporaryDirectory=args.strTMPDir, iSampleSelectionCount=int(args.icount), iSupervisedSampleCount=int(args.iSupervisedCount), strLabel=args.strLabel, strStratify=args.strUnsupervisedStratify, strSelectionTechnique=args.strSelection, iSampleNameRow=int(args.iSampleNameRow), iFirstDataRow=int(args.iFirstDataRow))
    logging.info("End microPITA")

    logging.debug("".join(["Returned the following samples:",str(dictSelectedSamples)]))

    strOutputContent = ""
    for sKey in dictSelectedSamples:
        strOutputContent = "".join([strOutputContent,sKey,Constants.COLON,", ".join(dictSelectedSamples[sKey]),Constants.ENDLINE])

    #Write to file
    if(not strOutputContent == ""):
        fHndlOutput = open(args.strOutFile,'w')
        fHndlOutput.write(str(strOutputContent))
        fHndlOutput.close()
    logging.debug("".join(["Selected samples output to file:",args.strOutFile]))

if __name__ == "__main__":
    _main( )



