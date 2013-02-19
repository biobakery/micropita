#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Run analysis for the microPITA paper
"""

#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "MIT"
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
import logging
import math
import mlpy
import numpy as np
import operator
import os
import random
import scipy.cluster.hierarchy as hcluster
import scipy.spatial.distance
from types import *

class MicroPITA:
	"""
	Selects samples from a first tier of a multi-tiered study to be used in a second tier.
	Different methods can be used for selection.
	The expected input is an abundance table (and potentially a text file of targeted features,
	if using the targeted features option). Output is a list of samples exhibiting the
	characteristics of interest.
	"""

	#Constants
	#Diversity metrics Alpha
	c_strInverseSimpsonDiversity = Metric.c_strInvSimpsonDiversity
	c_strChao1Diversity = Metric.c_strChao1Diversity

	#Diversity metrics Beta
	c_strBrayCurtisDissimilarity = Metric.c_strBrayCurtisDissimilarity

	#Additive inverses of diversity metrics beta
	c_strInvBrayCurtisDissimilarity = Metric.c_strInvBrayCurtisDissimilarity

	#Technique Names
	ConstantsMicropita.c_strDiversity2 = ConstantsMicropita.c_strDiversity+"_C"

	#Targeted feature settings
	c_strTargetedRanked = ConstantsMicropita.c_strTargetedRanked
	c_strTargetedAbundance = ConstantsMicropita.c_strTargetedAbundance

	#Technique groupings
#	c_lsDiversityMethods = [ConstantsMicropita.c_strDiversity,ConstantsMicropita.c_strDiversity2]

	#Converts ecology metrics into standardized method selection names
	dictConvertAMetricDiversity = {c_strInverseSimpsonDiversity:ConstantsMicropita.c_strDiversity, c_strChao1Diversity:ConstantsMicropita.c_strDiversity2}
#	dictConvertMicroPITAToAMetric = {ConstantsMicropita.c_strDiversity:c_strInverseSimpsonDiversity, ConstantsMicropita.c_strDiversity2:c_strChao1Diversity}
	dictConvertBMetricToMethod = {c_strBrayCurtisDissimilarity:ConstantsMicropita.c_strRepresentative}
	dictConvertInvBMetricToMethod = {c_strBrayCurtisDissimilarity:ConstantsMicropita.c_strExtreme}

	#Linkage used in the Hierarchical clustering
	c_strHierarchicalClusterMethod = 'average'

####Group 1## Diversity
	#Testing: Happy path Testing (8)
	def funcGetTopRankedSamples(self, lldMatrix = None, lsSampleNames = None, iTopAmount = None):
		"""
		Given a list of lists of measurements, for each list the indices of the highest values are returned. If lsSamplesNames is given
			it is treated as a list of string names that is in the order of the measurements in each list. Indices are returned or the sample
			names associated with the indices.
		
		:param	lldMatrix:	List of lists [[value,value,value,value],[value,value,value,value]].
		:type:	List of lists	List of measurements. Each list is a different measurement. Each measurement in positionally related to a sample.
		:param	lsSampleNames:	List of sample names positionally related (the same) to each list (Optional).
		:type:	List of strings	List of strings.
		:param	iTopAmount:	The amount of top measured samples (assumes the higher measurements are better).
		:type:	integer	Integer amount of sample names/ indices to return.
		:return	List:	List of samples to be selected.
		"""
		topRankListRet = []
		for rowMetrics in lldMatrix:
			#Create 2 d array to hold value and index and sort
			liIndexX = [rowMetrics,range(len(rowMetrics))]
			liIndexX[1].sort(key = liIndexX[0].__getitem__,reverse = True)

			if lsSampleNames:
				topRankListRet.append([lsSampleNames[iIndex] for iIndex in liIndexX[1][:iTopAmount]])
			else:
				topRankListRet.append(liIndexX[1][:iTopAmount])

		return topRankListRet
	
	####Group 2## Representative Dissimilarity
	#Testing: Happy path tested 1
	def funcGetCentralSamplesByKMedoids(self, npaMatrix=None, sMetric=None, lsSampleNames=None, iNumberSamplesReturned=0, istmBetaMatrix=None, istrmTree=None, istrmEnvr=None):
		"""
		Gets centroid samples by k-medoids clustering of a given matrix.
		
		:param	npaMatrix:	Numpy array where row=features and columns=samples
		:type:	Numpy array	Abundance Data.
		:param	sMetric:	String name of beta metric used as the distance metric.
		:type:	String	String name of beta metric.
		:param	lsSampleNames:	The names of the sample
		:type:	List	List of strings
		:param	iNumberSamplesReturned:	Number of samples to return, each will be a centroid of a sample.
		:type:	Integer	Number of samples to return
		:return	List:	List of selected samples.
		:param	istmBetaMatrix: File with beta-diversity matrix
		:type:	File stream or file path string
		"""

		#Count of how many rows
		sampleCount = npaMatrix.shape[0]
		if iNumberSamplesReturned > sampleCount:
			logging.error("MicroPITA.funcGetCentralSamplesByKMedoids:: There are not enough samples to return the amount of samples specified. Return sample count = "+str(iNumberSamplesReturned)+". Sample number = "+str(sampleCount)+".")
			return False

		#If the cluster count is equal to the sample count return all samples
		if sampleCount == iNumberSamplesReturned:
			return list(lsSampleNames)

		#Get distance matrix
		distanceMatrix=scipy.spatial.distance.squareform(Metric.funcReadMatrixFile(istmMatrixFile=istmBetaMatrix,lsSampleOrder=lsSampleNames)[0]) if istmBetaMatrix else Metric.funcGetBetaMetric(npadAbundancies=npaMatrix, sMetric=sMetric, istrmTree=istrmTree, istrmEnvr=istrmEnvr, lsSampleOrder=lsSampleNames)
		if type(distanceMatrix) is BooleanType:
			logging.error("MicroPITA.funcGetCentralSamplesByKMedoids:: Could not read in the supplied distance matrix, returning false.")
			return False
	
		#Log distance matrix
		logging.debug("MicroPITA.funcGetCentralSamplesByKMedoids:: Distance matrix for representative selection using metric="+str(sMetric))
	
		distance = MLPYDistanceAdaptor(npaDistanceMatrix=distanceMatrix, fIsCondensedMatrix=True)
	
		#Create object to determine clusters/medoids
		medoidsMaker = Kmedoids(k=iNumberSamplesReturned, dist=distance)
		#medoidsData includes(1d numpy array, medoids indexes; 
		#			  1d numpy array, non-medoids indexes;
		#			  1d numpy array, cluster membership for non-medoids;
		#			  double, cost of configuration)
		#npaMatrix is samples x rows
		#Build a matrix of lists of indicies to pass to the distance matrix
		lliIndicesMatrix = [[iIndexPosition] for iIndexPosition in xrange(0,len(npaMatrix))]
		medoidsData = medoidsMaker.compute(np.array(lliIndicesMatrix))
		logging.debug("MicroPITA.funcGetCentralSamplesByKMedoids:: Results from the kmedoid method in representative selection:")
		logging.debug(str(medoidsData))
	
		#If returning the same amount of clusters and samples
		#Return centroids
		selectedIndexes = medoidsData[0]
		return [lsSampleNames[selectedIndexes[index]] for index in xrange(0,iNumberSamplesReturned)]
	
	####Group 3## Highest Dissimilarity
	#Testing: Happy path tested
	def funcSelectExtremeSamplesFromHClust(self, strBetaMetric, npaAbundanceMatrix, lsSampleNames, iSelectSampleCount, istmBetaMatrix=None, istrmTree=None, istrmEnvr=None):
		"""
		Select extreme samples from HClustering.
		
		:param	strBetaMetric:	The beta metric to use for distance matrix generation.
		:type:	String	The name of the beta metric to use.
		:param	npaAbundanceMatrix:	Numpy array where row=samples and columns=features.
		:type:	Numpy Array	Abundance data.
		:param	lsSampleNames:	The names of the sample.
		:type:	List	List of strings.
		:param	iSelectSampleCount:	Number of samples to select (return).
		:type:	Integer	Integer number of samples returned.
		:return	Samples:	List of samples.
		:param	istmBetaMatrix: File with beta-diversity matrix
		:type:	File stream or file path string
		"""
	
		#If they want all the sample count, return all sample names
		iSampleCount=len(npaAbundanceMatrix[:,0])
		if iSelectSampleCount==iSampleCount:
		  return lsSampleNames
	
		#Holds the samples to be returned
		lsReturnSamplesRet = []
	
		#Generate beta matrix
		#Returns condensed matrix
		tempDistanceMatrix = scipy.spatial.distance.squareform(Metric.funcReadMatrixFile(istmMatrixFile=istmBetaMatrix,lsSampleOrder=lsSampleNames)[0]) if istmBetaMatrix else Metric.funcGetBetaMetric(npadAbundancies=npaAbundanceMatrix, sMetric=strBetaMetric, istrmTree=istrmTree, istrmEnvr=istrmEnvr, lsSampleOrder=lsSampleNames, fAdditiveInverse = True)

		if type(tempDistanceMatrix) is BooleanType:
			logging.error("MicroPITA.funcSelectExtremeSamplesFromHClust:: Could not read in the supplied distance matrix, returning false.")
			return False
		if istmBetaMatrix:
			tempDistanceMatrix = 1-tempDistanceMatrix

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
	
		iCurrentSelectCount = 0
		for row in linkageMatrix:
			#Get nodes ofthe lowest pairing (so the furthest apart pair)
			iNode1 = int(row[0])
			iNode2 = int(row[1])
			#Make sure the nodes are a terminal node (sample) and not a branch in the dendrogram
			#The branching in the dendrogram will start at the number of samples and increment higher.
			#Add each of the pair one at a time breaking when enough samples are selected.
			if iNode1<iSampleCount:
				lsReturnSamplesRet.append(lsSampleNames[iNode1])
				iCurrentSelectCount = iCurrentSelectCount + 1
			if iCurrentSelectCount == iSelectSampleCount:
				break
			if iNode2<iSampleCount:
				lsReturnSamplesRet.append(lsSampleNames[iNode2])
				iCurrentSelectCount = iCurrentSelectCount + 1
			if iCurrentSelectCount == iSelectSampleCount:
				break
	
		#Return selected samples
		return lsReturnSamplesRet
	
	####Group 4## Rank Average of user Defined Taxa
		#Testing: Happy Path Tested
	def funcGetAverageAbundanceSamples(self, abndTable, lsTargetedFeature, fRank=False):
		"""
		Averages feature abundance or ranked abundance. Expects a column 0 of taxa id that is skipped.
		
		:param	abndTable:	Abundance Table to analyse
		:type:	AbundanceTable	Abundance Table
		:param	lsTargetedFeature:	String names
		:type:	list	list of string names of features (bugs) which are measured after ranking against the full sample
		:param  fRank:	Indicates to rank the abundance before getting the average abundance of the features (default false)
		:type:   boolean	Flag indicating ranking abundance before calculating average feature measurement (false= no ranking)
		:return	List of lists or boolean:	List of lists or False on error. One internal list per sample indicating the sample,
				feature average abundance or ranked abundance. Lists will already be sorted.
				For not Ranked [[sample,average abundance of selected feature,1]]
				For Ranked [[sample,average ranked abundance, average abundance of selected feature]]
				Error Returns false
		"""
		
		llAbundance = abndTable.funcGetAverageAbundancePerSample(lsTargetedFeature)
		if not llAbundance:
			logging.error("MicroPITA.funcGetAverageAbundanceSamples:: Could not get average abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table.")
			return False
		#Add a space for ranking if needed
		#Not ranked will be [[sSample,average abundance,1]]
		#(where 1 will not discriminant ties if used in later functions, so this generalizes)
		#Ranked will be [[sSample, average rank, average abundance]]
		llRetAbundance = [[llist[0],-1,llist[1]] for llist in llAbundance]
		#Rank if needed
		if fRank:
			abndRanked = abndTable.funcRankAbundance()
			if abndRanked == None:
				logging.error("MicroPITA.funcGetAverageAbundanceSamples:: Could not rank the abundance table, returned false.")
				return False
			llRetRank = abndRanked.funcGetAverageAbundancePerSample(lsTargetedFeature)
			if not llRetRank:
				logging.error("MicroPITA.funcGetAverageAbundanceSamples:: Could not get average ranked abundance, returned false. Make sure the features (bugs) are spelled correctly and in the abundance table.")
				return False
			dictRanks = dict(llRetRank)
			llRetAbundance = [[a[0],dictRanks[a[0]],a[2]] for a in llRetAbundance]
			
		#Sort first for ties and then for the main feature
 		if not fRank or ConstantsMicropita.c_fBreakRankTiesByDiversity:
			llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[2], reverse = not fRank)
		if fRank:
			llRetAbundance = sorted(llRetAbundance, key = lambda sampleData: sampleData[1], reverse = not fRank)
		return llRetAbundance
	
	#Testing: Happy Path Tested
	def funcSelectTargetedTaxaSamples(self, abndMatrix, lsTargetedTaxa, iSampleSelectionCount, sMethod = ConstantsMicropita.lsTargetedFeatureMethodValues[0]):
	  """
	  Selects samples with the highest ranks or abundance of targeted features.
	  If ranked, select the highest abundance for tie breaking
	
	  :param	abndMatrix:	Abundance table to analyse 
	  :type:	AbundanceTable	Abundance table
	  :param	lsTargetedTaxa:	List of features
	  :type:	list	list of strings
	  :param	iSampleSelectionCount:	Number of samples to select
	  :type:	integer	integer
	  :param	sMethod:	Method to select targeted features
	  :type:	string	String (Can be values found in ConstantsMicropita.lsTargetedFeatureMethodValues)
	  :return	List of strings:	List of sample names which were selected
	  List of strings	Empty list is returned on an error.
	  """
	
	  #Check data
	  if(len(lsTargetedTaxa) < 1):
		logging.error("MicroPITA.funcSelectTargetedTaxaSamples. Taxa defined selection was requested but no features were given.")
		return []

	  lsTargetedSamples = self.funcGetAverageAbundanceSamples(abndTable=abndMatrix, lsTargetedFeature=lsTargetedTaxa,
	  	fRank=sMethod.lower() == self.c_strTargetedRanked.lower())
	  #If an error occured or the key word for the method was not recognized
	  if lsTargetedSamples == False: 
		  logging.error("MicroPITA.funcSelectTargetedTaxaSamples:: Was not able to select for the features given. So targeted feature selection was performed. Check to make sure the features are spelled correctly and exist in the abundance file.")
		  return []
	
	  #Select from results
	  return [sSample[0] for sSample in lsTargetedSamples[:iSampleSelectionCount]]
	
	####Group 5## Random
	#Testing: Happy path Tested
	def funcGetRandomSamples(self, lsSamples=None, iNumberOfSamplesToReturn=0):
		"""
		Returns random sample names of the number given. No replacement.
		
		:param	lsSamples:	List of sample names 
		:type:	list	list of strings
		:param	iNumberOfSamplesToReturn:	Number of samples to select
		:type:	integer	integer.
		:return	List:	List of selected samples (strings).
		"""

		#Input matrix sample count
		sampleCount = len(lsSamples)

		#Return the full matrix if they ask for a return matrix where length == original
		if(iNumberOfSamplesToReturn >= sampleCount):
			return lsSamples
	
		#Get the random indices for the sample (without replacement)
		liRandomIndices = random.sample(range(sampleCount), iNumberOfSamplesToReturn)
	
		#Create a boolean array of if indexes are to be included in the reduced array
                return [sSample for iIndex, sSample in enumerate(lsSamples) if iIndex in liRandomIndices]

	#Happy path tested (case 3)
	def funcGetAveragePopulation(self, abndTable, lfCompress):
		"""
		Get the average row per column in the abndtable.

		:param abndTable: AbundanceTable of data to be averaged
		:type: AbudanceTable
		:param lfCompress: List of boolean flags (false means to remove sample before averaging
		:type: List of floats
		:return List of doubles: 
		"""
		if sum(lfCompress) == 0:
			return []

		#Get the average populations
		lAverageRet = []

		for sFeature in abndTable.funcGetAbundanceCopy():
			sFeature = list(sFeature)[1:]
			sFeature=np.compress(lfCompress,sFeature,axis=0)
			lAverageRet.append(sum(sFeature)/float(len(sFeature)))
		return lAverageRet

	#Happy path tested (2 cases)
	def funcGetDistanceFromAverage(self, abndTable,ldAverage,lsSamples,lfSelected):
		"""
		Given an abundance table and an average sample, this returns the distance of each sample
		(measured using brays-curtis dissimilarity) from the average.
		The distances are reduced by needing to be in the lsSamples and being a true in the lfSelected
		(which is associated with the samples in the order of the samples in the abundance table;
		use abundancetable.funcGetSampleNames() to see the order if needed).

		:param abndTable: Abundance table holding the data to be analyzed.
		:type: AbundanceTable
		:param ldAverage: Average population (Average features of the abundance table of samples)
		:type: List of doubles which represent the average population
		:param lsSamples: These are the only samples used in the analysis
		:type: List of strings (sample ids)
		:param lfSelected: Samples to be included in the analysis
		:type: List of boolean (true means include)
		:return: List of distances (doubles)
		"""
		#Get the distance from label 1 of all samples in label0 splitting into selected and not selected lists
		ldSelectedDistances = []

		for sSampleName in [sSample for iindex, sSample in enumerate(lsSamples) if lfSelected[iindex]]:
			#Get the sample measurements
			ldSelectedDistances.append(Metric.funcGetBrayCurtisDissimilarity(np.array([abndTable.funcGetSample(sSampleName),ldAverage]))[0])
		return ldSelectedDistances

	#Happy path tested (1 case)
	def funcMeasureDistanceFromLabelToAverageOtherLabel(self, abndTable, lfGroupOfInterest, lfGroupOther):
		"""
		Get the distance of samples from one label from the average sample of not the label.
		Note: This assumes 2 classes.  

		:param abndTable: Table of data to work out of.
		:type: Abundace Table
		:param lfGroupOfInterest: Boolean indicator of the sample being in the first group.
		:type: List of floats, true indicating an individual in the group of interest.
		:param lfGroupOther:	Boolean indicator of the sample being in the other group.
		:type:	List of floats, true indicating an individual in the 
		:return List of List of doubles: [list of tuples (string sample name,double distance) for the selected population, list of tuples for the not selected population]
		"""
		#Get all sample names
		lsAllSamples = abndTable.funcGetSampleNames()

		#Get average populations
		lAverageOther = self.funcGetAveragePopulation(abndTable=abndTable, lfCompress=lfGroupOther)

		#Get the distance from the average of the other label (label 1)
		ldSelectedDistances = self.funcGetDistanceFromAverage(abndTable=abndTable, ldAverage=lAverageOther,
			lsSamples=lsAllSamples, lfSelected=lfGroupOfInterest)

		return zip([lsAllSamples[iindex] for iindex, fGroup in enumerate(lfGroupOfInterest) if fGroup],ldSelectedDistances)

	#Happy path tested (1 test case)
	def funcPerformDistanceSelection(self, abndTable, iSelectionCount, sLabel, sValueOfInterest):
		"""
		Given metadata, metadata of one value (sValueOfInterest) is measured from the average (centroid) value of another label group.
		An iSelectionCount of samples is selected from the group of interest closest to and furthest from the centroid of the other group.

		:params  abndTable: Abundance of measurements
		:type: AbundanceTable
		:params iSelectionCount: The number of samples selected per sample.
		:type: Integer Integer greater than 0
		:params sLabel: ID of the metadata which is the supervised label
		:type: String
		:params sValueOfInterest: Metadata value in the sLabel metadta row of the abundance table which defines the group of interest.
		:type: String found in the abundance table metadata row indicated by sLabel.
		:return list list of tuples (samplename, distance) [[iSelectionCount of tuples closest to the other centroid], [iSelectionCount of tuples farthest from the other centroid], [all tuples of samples not selected]]
		"""

		lsMetadata = abndTable.funcGetMetadata(sLabel)
		#Other metadata values
		lsUniqueOtherValues = list(set(lsMetadata)-set(sValueOfInterest))

		#Get boolean indicator of values of interest
		lfLabelsInterested = [sValueOfInterest == sValue for sValue in lsMetadata]

                #Get the distances of the items of interest from the other metadata values
		dictDistanceAverages = {}
                for sOtherLabel in lsUniqueOtherValues:
			#Get boolean indicator of labels not of interest 
			lfLabelsOther = [sOtherLabel == sValue for sValue in lsMetadata]

			#Get the distances of data from two different groups to the average of the other
			ldValueDistances = dict(self.funcMeasureDistanceFromLabelToAverageOtherLabel(abndTable, lfLabelsInterested, lfLabelsOther))

			for sKey in ldValueDistances:
				dictDistanceAverages[sKey] = ldValueDistances[sKey] + dictDistanceAverages[sKey] if sKey in dictDistanceAverages else ldValueDistances[sKey]

		#Finish average by dividing by length of lsUniqueOtherValues
		ltpleAverageDistances = [(sKey, dictDistanceAverages[sKey]/float(len(lsUniqueOtherValues))) for sKey in dictDistanceAverages]

                #Sort to extract extremes
                ltpleAverageDistances = sorted(ltpleAverageDistances,key=operator.itemgetter(1))

		#Get the closest and farthest distances
		ltupleDiscriminantSamples = ltpleAverageDistances[:iSelectionCount]
		ltupleDistinctSamples = ltpleAverageDistances[iSelectionCount*-1:]

		#Remove the selected samples from the larger population of distances (better visualization)
		ldSelected = [tpleSelected[0] for tpleSelected in ltupleDiscriminantSamples+ltupleDistinctSamples]

		#Return discriminant tuples, distinct tuples, other tuples
		return [ltupleDiscriminantSamples, ltupleDistinctSamples,
			   [tplData for tplData in ltpleAverageDistances if tplData[0] not in ldSelected]]

	#Run the supervised method surrounding distance from centroids
	#Happy path tested (3 test cases)
	def funcRunSupervisedDistancesFromCentroids(self, abundanceTable, fRunDistinct, fRunDiscriminant,
						xOutputSupFile, xPredictSupFile, strSupervisedMetadata,
						iSampleSupSelectionCount, lsOriginalSampleNames, lsOriginalLabels, fAppendFiles = False):
		"""
		Runs supervised methods based on measuring distances of one label from the centroid of another. NAs are evaluated as theirown group.

		:param	abundanceTable:	AbundanceTable
		:type:	AbudanceTable	Data to analyze
		:param	fRunDistinct:	Run distinct selection method
		:type:	Boolean	boolean (true runs method)
		:param	fRunDiscriminant:	Run discriminant method
		:type:	Boolean	boolean (true runs method)
		:param	xOutputSupFile:	File output from supervised methods detailing data going into the method.
		:type:	String or FileStream
		:param	xPredictSupFile:	File output from supervised methods distance results from supervised methods.
		:type:	String or FileStream
		:param strSupervisedMetadata:	The metadata that will be used to group samples.
		:type:	String
		:param	iSampleSupSelectionCount:	Number of samples to select
		:type:	Integer	int sample selection count
		:param lsOriginalSampleNames:	List of the sample names, order is important and should be preserved from the abundanceTable.
		:type:	List of samples	
		:param	fAppendFiles:	Indicates that output files already exist and appending is occuring.
		:type:	Boolean
		:return	Selected Samples:	A dictionary of selected samples by selection ID
		Dictionary	{"Selection Method":["SampleID","SampleID"...]}
		"""
		#Get labels and run one label against many
		lstrMetadata = abundanceTable.funcGetMetadata(strSupervisedMetadata)
		dictlltpleDistanceMeasurements = {}
		for sMetadataValue in set(lstrMetadata):

			#For now perform the selection here for the label of interest against the other labels
			dictlltpleDistanceMeasurements.setdefault(sMetadataValue,[]).extend(self.funcPerformDistanceSelection(abndTable=abundanceTable,
				iSelectionCount=iSampleSupSelectionCount, sLabel=strSupervisedMetadata, sValueOfInterest=sMetadataValue))

		#Make expected output files for supervised methods
		#1. Output file which is similar to an input file for SVMs
		#2. Output file that is similar to the probabilitic output of a SVM (LibSVM)
		#Manly for making output of supervised methods (Distance from Centroid) similar
		#MicropitaVis needs some of these files
		if xOutputSupFile:
			if fAppendFiles:
				SVM.funcUpdateSVMFileWithAbundanceTable(abndAbundanceTable=abundanceTable, xOutputSVMFile=xOutputSupFile,
					lsOriginalLabels=lsOriginalLabels, lsSampleOrdering=lsOriginalSampleNames)
			else:
				SVM.funcConvertAbundanceTableToSVMFile(abndAbundanceTable=abundanceTable, xOutputSVMFile=xOutputSupFile,
					sMetadataLabel=strSupervisedMetadata, lsOriginalLabels=lsOriginalLabels, lsSampleOrdering=lsOriginalSampleNames)

		#Will contain the samples selected to return
		#One or more of the methods may be active so this is why I am extending instead of just returning the result of each method type
		dictSelectedSamplesRet = dict()
		for sKey, ltplDistances in dictlltpleDistanceMeasurements.items():
			if fRunDistinct:
				dictSelectedSamplesRet.setdefault(ConstantsMicropita.c_strDistinct,[]).extend([ltple[0] for ltple in ltplDistances[1]])
			if fRunDiscriminant:
				dictSelectedSamplesRet.setdefault(ConstantsMicropita.c_strDiscriminant,[]).extend([ltple[0] for ltple in ltplDistances[0]])

		if xPredictSupFile:
			dictFlattenedDistances = dict()
			[dictFlattenedDistances.setdefault(sKey, []).append(tple)
				for sKey, lltple in dictlltpleDistanceMeasurements.items()
				for ltple in lltple for tple in ltple]
			if fAppendFiles:
				self._updatePredictFile(xPredictSupFile=xPredictSupFile, xInputLabelsFile=xOutputSupFile,
					dictltpleDistanceMeasurements=dictFlattenedDistances, abundanceTable=abundanceTable, lsOriginalSampleNames=lsOriginalSampleNames)
			else:
				self._writeToPredictFile(xPredictSupFile=xPredictSupFile, xInputLabelsFile=xOutputSupFile,
					dictltpleDistanceMeasurements=dictFlattenedDistances, abundanceTable=abundanceTable, lsOriginalSampleNames=lsOriginalSampleNames)
		return dictSelectedSamplesRet

	#Two happy path test cases
	def _updatePredictFile(self, xPredictSupFile, xInputLabelsFile, dictltpleDistanceMeasurements, abundanceTable, lsOriginalSampleNames):
		"""
		Manages updating the predict file.

		:param	xPredictSupFile: File that has predictions (distances) from the supervised method.
		:type:	FileStream or String file path
		:param	xInputLabelsFile: File that as input to the supervised methods.
		:type:	FileStream or String file path
		:param	dictltpleDistanceMeasurements: 
		:type:	Dictionary of lists of tuples {"labelgroup":[("SampleName",dDistance)...], "labelgroup":[("SampleName",dDistance)...]}
		"""

		if not isinstance(xPredictSupFile, str):
			xPredictSupFile.close()
			xPredictSupFile = xPredictSupFile.name
		csvr = open(xPredictSupFile,'r')

		f = csv.reader(csvr,delimiter=ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)
		lsHeader = f.next()[1:]
		dictlltpleRead = dict([(sHeader,[]) for sHeader in lsHeader])

		#Read data in 
		iSampleIndex = 0
		for sRow in f:
			sLabel = sRow[0]
			[dictlltpleRead[lsHeader[iDistanceIndex]].append((lsOriginalSampleNames[iSampleIndex],dDistance)) for iDistanceIndex, dDistance in enumerate(sRow[1:])
				if not dDistance == ConstantsMicropita.c_sEmptyPredictFileValue]
			iSampleIndex += 1

		#Combine dictltpleDistanceMeasurements with new data
		#If they share a key then merge keeping parameter data
		#If they do not share the key, keep the full data
		dictNew = {}
		for sKey in dictltpleDistanceMeasurements.keys():
			lsSamples = [tple[0] for tple in dictltpleDistanceMeasurements[sKey]]
			dictNew[sKey] = dictltpleDistanceMeasurements[sKey]+[tple for tple in dictlltpleRead[sKey] if tple[0] not in lsSamples] if sKey in dictlltpleRead.keys() else dictltpleDistanceMeasurements[sKey]
                for sKey in dictlltpleRead:
			if sKey not in dictltpleDistanceMeasurements.keys():
				dictNew[sKey] = dictlltpleRead[sKey]

		#Call writer
		self._writeToPredictFile(xPredictSupFile=xPredictSupFile, xInputLabelsFile=xInputLabelsFile,
			dictltpleDistanceMeasurements=dictNew, abundanceTable=abundanceTable,
			lsOriginalSampleNames=lsOriginalSampleNames, fFromUpdate=True)

	#2 happy path test cases
        def _writeToPredictFile(self, xPredictSupFile, xInputLabelsFile, dictltpleDistanceMeasurements, abundanceTable, lsOriginalSampleNames, fFromUpdate=False):
		"""
		Write to the predict file.

		:param	xPredictSupFile: File that has predictions (distances) from the supervised method.
		:type:	FileStream or String file path
		:param	xInputLabelsFile: File that as input to the supervised methods.
		:type:	FileStream or String file path
		:param	dictltpleDistanceMeasurements: 
		:type:	Dictionary of lists of tuples {"labelgroup":[("SampleName",dDistance)...], "labelgroup":[("SampleName",dDistance)...]}
		:param	abundanceTable: An abundance table of the sample data.
		:type:	AbundanceTable
		:param	lsOriginalSampleNames: Used if the file is being updated as the sample names so that it may be passed in and consistent with other writing.
			Otherwise will use the sample names from the abundance table.
		:type:	List of strings
		:param	fFromUpdate:	Indicates if this is part of an update to the file or not.
		:type:	Boolean
		"""

		xInputLabelsFileName = xInputLabelsFile
		if not isinstance(xInputLabelsFile,str):
			xInputLabelsFileName = xInputLabelsFile.name
		f = csv.writer(open(xPredictSupFile,"w") if isinstance(xPredictSupFile, str) else xPredictSupFile,delimiter=ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)

		lsAllSampleNames = abundanceTable.funcGetSampleNames()
		lsLabels = SVM.funcReadLabelsFromFile(xSVMFile=xInputLabelsFileName, lsAllSampleNames= lsOriginalSampleNames if fFromUpdate else lsAllSampleNames,
						isPredictFile=False)
		dictLabels = dict([(sSample,sLabel) for sLabel in lsLabels.keys() for sSample in lsLabels[sLabel]])

		#Dictionay keys will be used to order the predict file
		lsMeasurementKeys = dictltpleDistanceMeasurements.keys()
		#Make header
		f.writerow(["labels"]+lsMeasurementKeys)

		#Reformat dictionary to make it easier to use
		for sKey in dictltpleDistanceMeasurements:
			dictltpleDistanceMeasurements[sKey] = dict([ltpl for ltpl in dictltpleDistanceMeasurements[sKey]])

		for sSample in lsOriginalSampleNames:
			#Make body of file
			f.writerow([dictLabels.get(sSample,ConstantsMicropita.c_sEmptyPredictFileValue)]+
				[str(dictltpleDistanceMeasurements[sKey].get(sSample,ConstantsMicropita.c_sEmptyPredictFileValue))
				for sKey in lsMeasurementKeys])

	def _funcRunNormalizeSensitiveMethods(self, abndData, iSampleSelectionCount, dictSelectedSamples, lsAlphaMetrics, lsBetaMetrics, lsInverseBetaMetrics,
												fRunDiversity, fRunRepresentative, fRunExtreme, strAlphaMetadata=None,
												istmBetaMatrix=None, istrmTree=None, istrmEnvr=None, fInvertDiversity=False):
		"""
		Manages running methods that are sensitive to normalization. This is called twice, once for the set of methods which should not be normalized and the other
		for the set that should be normalized.
	
		:param	abndData:	Abundance table object holding the samples to be measured.
		:type:	AbundanceTable
		:param	iSampleSelectionCount	The number of samples to select per method.
		:type:	Integer
		:param	dictSelectedSamples	Will be added to as samples are selected {"Method:["strSelectedSampleID","strSelectedSampleID"...]}.
		:type:	Dictionary
		:param	lsAlphaMetrics:	List of alpha metrics to use on alpha metric dependent assays (like highest diversity).
		:type:	List of strings
		:param	lsBetaMetrics:	List of beta metrics to use on beta metric dependent assays (like most representative).
		:type:	List of strings
		:param	lsInverseBetaMetrics:	List of inverse beta metrics to use on inverse beta metric dependent assays (like most dissimilar).
		:type:	List of strings
		:param	fRunDiversity:	Run Diversity based methods (true indicates run).
		:type:	Boolean	
		:param	fRunRepresentative:	Run Representative based methods (true indicates run).
		:type:	Boolean	
		:param	fRunExtreme:	Run Extreme based methods (true indicates run).
		:type:	Boolean	
		:param	istmBetaMatrix:	File that has a precalculated beta matrix
		:type:	File stream or File path string
		:return	Selected Samples:	Samples selected by methods.
				Dictionary	{"Selection Method":["SampleID","SampleID","SampleID",...]}
		"""

		#Sample ids/names
		lsSampleNames = abndData.funcGetSampleNames()
	
		#Generate alpha metrics and get most diverse
		if fRunDiversity:

			#Get Alpha metrics matrix
			internalAlphaMatrix = None
			#Name of technique
			strMethod = [strAlphaMetadata] if strAlphaMetadata else lsAlphaMetrics

			#If given an alpha-diversity metadata
			if strAlphaMetadata:
				internalAlphaMatrix = [[float(strNum) for strNum in abndData.funcGetMetadata(strAlphaMetadata)]]
			else:
				#Expects Observations (Taxa (row) x sample (column))
				#Returns [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
				internalAlphaMatrix = Metric.funcBuildAlphaMetricsMatrix(npaSampleAbundance = abndData.funcGetAbundanceCopy()
							if not abndData.funcIsSummed()
							else abndData.funcGetFeatureAbundanceTable(abndData.funcGetTerminalNodes()).funcGetAbundanceCopy(),
							lsSampleNames = lsSampleNames, lsDiversityMetricAlpha = lsAlphaMetrics)
	
			if internalAlphaMatrix:
				#Invert measurments
				if fInvertDiversity:
					lldNewDiversity = []
					for lsLine in internalAlphaMatrix:
						lldNewDiversity.append([1/max(dValue,ConstantsMicropita.c_smallNumber) for dValue in lsLine])
					internalAlphaMatrix = lldNewDiversity
				#Get top ranked alpha diversity by most diverse
				#Expects [[sample1,sample2,sample3...],[sample1,sample2,sample3..],...]
				#Returns [[sampleName1, sampleName2, sampleNameN],[sampleName1, sampleName2, sampleNameN]]
				mostDiverseAlphaSamplesIndexes = self.funcGetTopRankedSamples(lldMatrix=internalAlphaMatrix, lsSampleNames=lsSampleNames, iTopAmount=iSampleSelectionCount)

				#Add to results
				for index in xrange(0,len(strMethod)):
					strSelectionMethod = self.dictConvertAMetricDiversity.get(strMethod[index],ConstantsMicropita.c_strDiversity+"="+strMethod[index])
					dictSelectedSamples.setdefault(strSelectionMethod,[]).extend(mostDiverseAlphaSamplesIndexes[index])

		logging.info("MicroPITA.funcRunNormalizeSensitiveMethods:: Selected Samples 1b")
		logging.info(dictSelectedSamples)
	
		#Generate beta metrics and 
		if fRunRepresentative or fRunExtreme:

			#Abundance matrix transposed
			npaTransposedAbundance = UtilityMath.funcTransposeDataMatrix(abndData.funcGetAbundanceCopy(), fRemoveAdornments=True)
	
			#Get center selection using clusters/tiling
			#This will be for beta metrics in normalized space
			if fRunRepresentative:

				if istmBetaMatrix:
					#Get representative dissimilarity samples
					medoidSamples=self.funcGetCentralSamplesByKMedoids(npaMatrix=npaTransposedAbundance, sMetric=ConstantsMicropita.c_custom, lsSampleNames=lsSampleNames, iNumberSamplesReturned=iSampleSelectionCount, istmBetaMatrix=istmBetaMatrix, istrmTree=istrmTree, istrmEnvr=istrmEnvr)

					if medoidSamples:
						dictSelectedSamples.setdefault(ConstantsMicropita.c_strRepresentative+"="+ConstantsMicropita.c_custom,[]).extend(medoidSamples)
				else:
					logging.info("MicroPITA.funcRunNormalizeSensitiveMethods:: Performing representative selection on normalized data.")
					for bMetric in lsBetaMetrics:

						#Get representative dissimilarity samples
						medoidSamples=self.funcGetCentralSamplesByKMedoids(npaMatrix=npaTransposedAbundance, sMetric=bMetric, lsSampleNames=lsSampleNames, iNumberSamplesReturned=iSampleSelectionCount, istmBetaMatrix=istmBetaMatrix, istrmTree=istrmTree, istrmEnvr=istrmEnvr)

						if medoidSamples:
							dictSelectedSamples.setdefault(self.dictConvertBMetricToMethod.get(bMetric,ConstantsMicropita.c_strRepresentative+"="+bMetric),[]).extend(medoidSamples)

			#Get extreme selection using clusters, tiling
			if fRunExtreme:
				logging.info("MicroPITA.funcRunNormalizeSensitiveMethods:: Performing extreme selection on normalized data.")
				if istmBetaMatrix:

					#Samples for representative dissimilarity
					#This involves inverting the distance metric,
					#Taking the dendrogram level of where the number cluster == the number of samples to select
					#Returning a repersentative sample from each cluster
					extremeSamples = self.funcSelectExtremeSamplesFromHClust(strBetaMetric=ConstantsMicropita.c_custom, npaAbundanceMatrix=npaTransposedAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount=iSampleSelectionCount, istmBetaMatrix=istmBetaMatrix, istrmTree=istrmTree, istrmEnvr=istrmEnvr)
	
					#Add selected samples
					if extremeSamples:
						dictSelectedSamples.setdefault(ConstantsMicropita.c_strExtreme+"="+ConstantsMicropita.c_custom,[]).extend(extremeSamples)

				else:
					#Run KMedoids with inverse custom distance metric in normalized space
					for bMetric in lsInverseBetaMetrics:

						#Samples for representative dissimilarity
						#This involves inverting the distance metric,
						#Taking the dendrogram level of where the number cluster == the number of samples to select
						#Returning a repersentative sample from each cluster
						extremeSamples = self.funcSelectExtremeSamplesFromHClust(strBetaMetric=bMetric, npaAbundanceMatrix=npaTransposedAbundance, lsSampleNames=lsSampleNames, iSelectSampleCount=iSampleSelectionCount, istmBetaMatrix=istmBetaMatrix, istrmTree=istrmTree, istrmEnvr=istrmEnvr)
	
						#Add selected samples
						if extremeSamples:
							dictSelectedSamples.setdefault(self.dictConvertInvBMetricToMethod.get(bMetric,ConstantsMicropita.c_strExtreme+"="+bMetric),[]).extend(extremeSamples)

		logging.info("MicroPITA.funcRunNormalizeSensitiveMethods:: Selected Samples 2,3b")
		logging.info(dictSelectedSamples)
		return dictSelectedSamples

	def funcRun(self, strIDName, strLastMetadataName, istmInput,
					  ostmInputPredictFile, ostmPredictFile, ostmCheckedFile, ostmOutput,
					  cDelimiter, cFeatureNameDelimiter, strFeatureSelection,
					  istmFeatures, iCount, lstrMethods, strLabel = None, strStratify = None,
					  strCustomAlpha = None, strCustomBeta = None, strAlphaMetadata = None, istmBetaMatrix = None, istrmTree = None, istrmEnvr = None, 
					  iMinSeqs = ConstantsMicropita.c_liOccurenceFilter[0], iMinSamples = ConstantsMicropita.c_liOccurenceFilter[1], fInvertDiversity = False):
		"""
		Manages the selection of samples given different metrics.

		:param	strIDName: Sample Id metadata row
		:type:	String
		:param	strLastMetadataName: The id of the metadata positioned last in the abundance table.
		:type:	String	String metadata id.
		:param	istmInput: File to store input data to supervised methods.
		:type:	FileStream of String file path
		:param	ostmInputPredictFile: File to store distances from supervised methods.
		:type:	FileStream or String file path
		:param	ostmCheckedFile: File to store the AbundanceTable data after it is being checked.
		:type:	FileStream or String file path
		:param	ostmOutPut: File to store sample selection by methods of interest.
		:type:	FileStream or String file path
		:param	cDelimiter: Delimiter of abundance table.
		:type:	Character Char (default TAB).
		:param	cFeatureNameDelimiter: Delimiter of the name of features (for instance if they contain consensus lineages indicating clades).
		:type:	Character (default |).
		:param	stFeatureSelectionMethod: Which method to use to select features in a targeted manner (Using average ranked abundance or average abundance).
		:type:	String (specific values indicated in ConstantsMicropita.lsTargetedFeatureMethodValues).
		:param	istmFeatures: File which holds the features of interest if using targeted feature methodology.
		:type:	FileStream or String file path
		:param	iCount:	Number of samples to select in each methods, supervised methods select this amount per label if possible.
		:type:	Integer	integer.
		:param	lstrMethods: List of strings indicating selection techniques.
		:type:	List of string method names
		:param	strLabel: The metadata used for supervised labels.
		:type:	String
		:param	strStratify: The metadata used to stratify unsupervised data.
		:type:	String
		:param	strCustomAlpha: Custom alpha diversity metric
		:type:	String
		:param	strCustomBeta: Custom beta diversity metric
		:type:	String
		:param	strAlphaMetadata: Metadata id which is a diveristy metric to use in highest diversity sampling
		:type:	String
		:param	istmBetaMatrix: File containing precalculated beta-diversity matrix for representative sampling
		:type:	FileStream or String file path
		:param	istrmTree: File containing tree for phylogentic beta-diversity analysis
		:type:	FileStream or String file path
		:param	istrmEnvr: File containing environment for phylogentic beta-diversity analysis
		:type:	FileStream or String file path
		:param	iMinSeqs: Minimum sequence in the occurence filter which filters all features not with a minimum number of sequences in each of a minimum number of samples.
		:type:	Integer
		:param	iMinSamples: Minimum sample count for the occurence filter.
		:type:	Integer
		:param	fInvertDiversity: When true will invert diversity measurements before using.
		:type:	boolean
		:return	Selected Samples:	Samples selected by methods.
				Dictionary	{"Selection Method":["SampleID","SampleID","SampleID",...]}
		"""

		#Holds the top ranked samples from different metrics
		#dict[metric name] = [samplename,samplename...]
		selectedSamples = dict()
	
		#If a target feature file is given make sure that targeted feature is in the selection methods, if not add
		if ConstantsMicropita.c_strFeature in lstrMethods:
		  if not istmFeatures:
			logging.error("MicroPITA.funcRun:: Did not receive both the Targeted feature file and the feature selection method. MicroPITA did not run.")
			return False

		#Diversity metrics to run
		#Use custom metrics if specified
                #Custom beta metrics set to normalized only, custom alpha metrics set to count only
		diversityMetricsAlpha = [] if strCustomAlpha or strAlphaMetadata else [MicroPITA.c_strInverseSimpsonDiversity]
		diversityMetricsBeta = [] if istmBetaMatrix else [strCustomBeta] if strCustomBeta else [MicroPITA.c_strBrayCurtisDissimilarity]
#		inverseDiversityMetricsBeta = [MicroPITA.c_strInvBrayCurtisDissimilarity]
		diversityMetricsAlphaNoNormalize = [strAlphaMetadata] if strAlphaMetadata else [strCustomAlpha] if strCustomAlpha else []
		diversityMetricsBetaNoNormalize = []
#		inverseDiversityMetricsBetaNoNormalize = []

		#Targeted taxa
		userDefinedTaxa = []
	
		#Perform different flows flags
		c_RUN_MAX_DIVERSITY_1 = ConstantsMicropita.c_strDiversity in lstrMethods
		c_RUN_REPRESENTIVE_DISSIMILARITY_2 = ConstantsMicropita.c_strRepresentative in lstrMethods
		c_RUN_MAX_DISSIMILARITY_3 = ConstantsMicropita.c_strExtreme in lstrMethods
		c_RUN_RANK_AVERAGE_USER_4 = False
		if ConstantsMicropita.c_strFeature in lstrMethods:
			c_RUN_RANK_AVERAGE_USER_4 = True
			if not istmFeatures:
				logging.error("MicroPITA.funcRun:: No taxa file was given for taxa selection.") 
				return False
			#Read in taxa list, break down to lines and filter out empty strings
			userDefinedTaxa = filter(None,(s.strip( ) for s in istmFeatures.readlines()))
		c_RUN_RANDOM_5 = ConstantsMicropita.c_strRandom in lstrMethods
		c_RUN_DISTINCT = ConstantsMicropita.c_strDistinct in lstrMethods
		c_RUN_DISCRIMINANT = ConstantsMicropita.c_strDiscriminant in lstrMethods

		#Read in abundance data
		#Abundance is a structured array. Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
		#Abundance table object to read in and manage data
		totalAbundanceTable = AbundanceTable.funcMakeFromFile(xInputFile=istmInput, lOccurenceFilter = [iMinSeqs, iMinSamples],
								   cDelimiter=cDelimiter, sMetadataID=strIDName, sLastMetadata=strLastMetadataName, cFeatureNameDelimiter=cFeatureNameDelimiter,
								   xOutputFile=ostmCheckedFile)
		if not totalAbundanceTable:
			logging.error("MicroPITA.funcRun:: Could not read in the abundance table. Analysis was not performed."+
				" This often occurs when the Last Metadata is not specified correctly."+
				" Please check to make sure the Last Metadata selection is the row of the last metadata,"+
				" all values after this selection should be microbial measurements and should be numeric.")
			return False

		lsOriginalLabels = SVM.funcMakeLabels(totalAbundanceTable.funcGetMetadata(strLabel)) if strLabel else strLabel

		dictTotalMetadata = totalAbundanceTable.funcGetMetadataCopy()
		logging.debug("MicroPITA.funcRun:: Received metadata=" + str(dictTotalMetadata))
		#If there is only 1 unique value for the labels, do not run the Supervised methods
		if strLabel and ( len(set(dictTotalMetadata.get(strLabel,[]))) < 2 ):
			logging.error("The label " + strLabel + " did not have 2 or more values. Labels found=" + str(dictTotalMetadata.get(strLabel,[])))
			return False

		#Run unsupervised methods###
		#Stratify the data if need be and drop the old data
		lStratifiedAbundanceTables = totalAbundanceTable.funcStratifyByMetadata(strStratify) if strStratify else [totalAbundanceTable]

		#For each stratified abundance block or for the unstratfified abundance
		#Run the unsupervised blocks
		fAppendSupFiles = False
		for stratAbundanceTable in lStratifiedAbundanceTables:
			logging.info("MicroPITA.funcRun:: Running abundance block:"+stratAbundanceTable.funcGetName())

 			###NOT SUMMED, NOT NORMALIZED			
			#Only perform if the data is not yet normalized
			if not stratAbundanceTable.funcIsNormalized( ):
				#Need to first work with unnormalized data
				if c_RUN_MAX_DIVERSITY_1 or c_RUN_REPRESENTIVE_DISSIMILARITY_2 or c_RUN_MAX_DISSIMILARITY_3:

					self._funcRunNormalizeSensitiveMethods(abndData=stratAbundanceTable, iSampleSelectionCount=iCount,
													 dictSelectedSamples=selectedSamples, lsAlphaMetrics=diversityMetricsAlphaNoNormalize,
													 lsBetaMetrics=diversityMetricsBetaNoNormalize,
													 lsInverseBetaMetrics=diversityMetricsBetaNoNormalize,
													 fRunDiversity=c_RUN_MAX_DIVERSITY_1,fRunRepresentative=c_RUN_REPRESENTIVE_DISSIMILARITY_2,
													 fRunExtreme=c_RUN_MAX_DISSIMILARITY_3, strAlphaMetadata=strAlphaMetadata, 
                                                                                                         istrmTree=istrmTree, istrmEnvr=istrmEnvr, fInvertDiversity=fInvertDiversity)


			#Generate selection by the rank average of user defined taxa
			#Expects (Taxa (row) by Samples (column))
			#Expects a column 0 of taxa id that is skipped
			#Returns [(sample name,average,rank)]
			#SUMMED AND NORMALIZED
			stratAbundanceTable.funcSumClades()
			#Normalize data at this point
			stratAbundanceTable.funcNormalize()
			if c_RUN_RANK_AVERAGE_USER_4:
				selectedSamples[ConstantsMicropita.c_strFeature] = self.funcSelectTargetedTaxaSamples(abndMatrix=stratAbundanceTable,
						lsTargetedTaxa=userDefinedTaxa, iSampleSelectionCount=iCount, sMethod=strFeatureSelection)
				logging.info("MicroPITA.funcRun:: Selected Samples Rank")
				logging.info(selectedSamples)

 			###SUMMED AND NORMALIZED analysis block
			#Diversity based metric will move reduce to terminal taxa as needed
			if c_RUN_MAX_DIVERSITY_1 or c_RUN_REPRESENTIVE_DISSIMILARITY_2 or c_RUN_MAX_DISSIMILARITY_3:

				self._funcRunNormalizeSensitiveMethods(abndData=stratAbundanceTable, iSampleSelectionCount=iCount,
												 dictSelectedSamples=selectedSamples, lsAlphaMetrics=diversityMetricsAlpha,
												 lsBetaMetrics=diversityMetricsBeta,
												 lsInverseBetaMetrics=diversityMetricsBeta,
												 fRunDiversity=c_RUN_MAX_DIVERSITY_1,fRunRepresentative=c_RUN_REPRESENTIVE_DISSIMILARITY_2,
												 fRunExtreme=c_RUN_MAX_DISSIMILARITY_3,
                                                                                                 istmBetaMatrix=istmBetaMatrix, istrmTree=istrmTree, istrmEnvr=istrmEnvr, fInvertDiversity=fInvertDiversity)

			#5::Select randomly
			#Expects sampleNames = List of sample names [name, name, name...]
			if(c_RUN_RANDOM_5):
				#Select randomly from sample names
				selectedSamples[ConstantsMicropita.c_strRandom] = self.funcGetRandomSamples(lsSamples=stratAbundanceTable.funcGetSampleNames(), iNumberOfSamplesToReturn=iCount)
				logging.info("MicroPITA.funcRun:: Selected Samples Random")
				logging.info(selectedSamples)

			#Perform supervised selection
			if c_RUN_DISTINCT or c_RUN_DISCRIMINANT:
 				if strLabel:
					dictSelectionRet = self.funcRunSupervisedDistancesFromCentroids(abundanceTable=stratAbundanceTable,
								fRunDistinct=c_RUN_DISTINCT, fRunDiscriminant=c_RUN_DISCRIMINANT,
								xOutputSupFile=ostmInputPredictFile,xPredictSupFile=ostmPredictFile,
								strSupervisedMetadata=strLabel, iSampleSupSelectionCount=iCount,
								lsOriginalSampleNames = totalAbundanceTable.funcGetSampleNames(),
								lsOriginalLabels = lsOriginalLabels,
								fAppendFiles=fAppendSupFiles)

					[selectedSamples.setdefault(sKey,[]).extend(lValue) for sKey,lValue in dictSelectionRet.items()]

					if not fAppendSupFiles:
						fAppendSupFiles = True
					logging.info("MicroPITA.funcRun:: Selected Samples Unsupervised")
					logging.info(selectedSamples)
		return selectedSamples
	
	#Testing: Happy path tested
	@staticmethod
	def funcWriteSelectionToFile(dictSelection,xOutputFilePath):
		"""
		Writes the selection of samples by method to an output file.
		
		:param	dictSelection:	The dictionary of selections by method to be written to a file.
		:type:	Dictionary	The dictionary of selections by method {"method":["sample selected","sample selected"...]}
		:param	xOutputFilePath:	FileStream or String path to file inwhich the dictionary is written.
		:type:	String	FileStream or String path to file
		"""
	
		if not dictSelection:
			return

		#Open file
		f = csv.writer(open(xOutputFilePath,"w") if isinstance(xOutputFilePath, str) else xOutputFilePath, delimiter=ConstantsMicropita.c_outputFileDelim )

		#Create output content from dictionary
		for sKey in dictSelection:
			f.writerow([sKey]+dictSelection[sKey])
			logging.debug("MicroPITA.funcRun:: Selected samples output to file:"+str(dictSelection[sKey]))
	
	#Testing: Happy Path tested
	@staticmethod
	def funcReadSelectionFileToDictionary(xInputFile):
		"""
		Reads in an output selection file from micropita and formats it into a dictionary.
		
		:param	xInputFile:	String path to file or file stream to read and translate into a dictionary.
									{"method":["sample selected","sample selected"...]}
		:type:	FileStream or String Path to file
		:return	Dictionary:	Samples selected by methods.
					Dictionary	{"Selection Method":["SampleID","SampleID","SampleID",...]}
		"""

		#Open file
		istmReader = csv.reader(open(xInputFile,'r') if isinstance(xInputFile, str) else xInputFile, delimiter = ConstantsMicropita.c_outputFileDelim)

		#Dictionary to hold selection data
		return dict([(lsLine[0], lsLine[1:]) for lsLine in istmReader])

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicroPITA.py", 
	description = """Selects samples from abundance tables based on various selection schemes.""" )

args = argp.add_argument_group( "Common", "Commonly modified options" )
args.add_argument(ConstantsMicropita.c_strCountArgument,"--num", dest="iCount", metavar = "samples", default = 10, type = int, help = ConstantsMicropita.c_strCountHelp)
args.add_argument("-m","--method", dest = "lstrMethods", metavar = "method", default = [], help = ConstantsMicropita.c_strSelectionTechniquesHelp,
	choices = ConstantsMicropita.c_lsAllMethods, action = "append")

args = argp.add_argument_group( "Custom", "Selecting and inputing custom metrics" )
args.add_argument("-a","--alpha", dest = "strAlphaDiversity", metavar = "AlphaDiversity", default = None, help = ConstantsMicropita.c_strCustomAlphaDiversityHelp,  choices = Metric.setAlphaDiversities)
args.add_argument("-b","--beta", dest = "strBetaDiversity", metavar = "BetaDiversity", default = None, help = ConstantsMicropita.c_strCustomBetaDiversityHelp,  choices = list(Metric.setBetaDiversities)+[Metric.c_strUnifracUnweighted,Metric.c_strUnifracWeighted])
args.add_argument("-q","--alphameta", dest = "strAlphaMetadata", metavar = "AlphaDiversityMetadata", default = None, help = ConstantsMicropita.c_strCustomAlphaDiversityMetadataHelp)
args.add_argument("-x","--betamatrix", dest = "istmBetaMatrix", metavar = "BetaDiversityMatrix", default = None, help = ConstantsMicropita.c_strCustomBetaDiversityMatrixHelp)
args.add_argument("-o","--tree", dest = "istrmTree", metavar = "PhylogeneticTree", default = None, help = ConstantsMicropita.c_strCustomPhylogeneticTreeHelp)
args.add_argument("-i","--envr", dest = "istrmEnvr", metavar = "EnvironmentFile", default = None, help = ConstantsMicropita.c_strCustomEnvironmentFileHelp)
args.add_argument("-f","--invertDiversity", dest = "fInvertDiversity", action="store_true", default = False, help = ConstantsMicropita.c_strInvertDiversityHelp)

args = argp.add_argument_group( "Miscellaneous", "Row/column identifiers and feature targeting options" )
args.add_argument("-d",ConstantsMicropita.c_strIDNameArgument, dest="strIDName", metavar="sample_id", help= ConstantsMicropita.c_strIDNameHelp)
args.add_argument("-l",ConstantsMicropita.c_strLastMetadataNameArgument, dest="strLastMetadataName", metavar = "metadata_id",
				  help= ConstantsMicropita.c_strLastMetadataNameHelp)
args.add_argument("-r",ConstantsMicropita.c_strTargetedFeatureMethodArgument, dest="strFeatureSelection", metavar="targeting_method", default=ConstantsMicropita.lsTargetedFeatureMethodValues[0],
				  choices=ConstantsMicropita.lsTargetedFeatureMethodValues, help= ConstantsMicropita.c_strTargetedFeatureMethodHelp)
args.add_argument("-t",ConstantsMicropita.c_strTargetedSelectionFileArgument, dest="istmFeatures", metavar = "feature_file", type = argparse.FileType("rU"), help = ConstantsMicropita.c_strTargetedSelectionFileHelp)

args = argp.add_argument_group( "Data labeling", "Metadata IDs for strata and supervised label values" )
args.add_argument("-e",ConstantsMicropita.c_strSupervisedLabelArgument, dest="strLabel", metavar= "supervised_id", help= ConstantsMicropita.c_strSupervisedLabelHelp)
args.add_argument("-s",ConstantsMicropita.c_strUnsupervisedStratifyMetadataArgument, dest="strUnsupervisedStratify", metavar= "stratify_id", 
				  help= ConstantsMicropita.c_strUnsupervisedStratifyMetadataHelp)

args = argp.add_argument_group( "File formatting", "Rarely modified file formatting options" )
args.add_argument("-j",ConstantsMicropita.c_strFileDelimiterArgument, dest="cFileDelimiter", metavar="column_delimiter", default="\t", help=ConstantsMicropita.c_strFileDelimiterHelp) 
args.add_argument("-k",ConstantsMicropita.c_strFeatureNameDelimiterArgument, dest="cFeatureNameDelimiter", metavar="taxonomy_delimiter", default="|", help=ConstantsMicropita.c_strFeatureNameDelimiterHelp) 

args = argp.add_argument_group( "Debugging", "Debugging options - modify at your own risk!" )
args.add_argument("-v",ConstantsMicropita.c_strLoggingArgument, dest="strLogLevel", metavar = "log_level", default="WARNING", 
				  choices=ConstantsMicropita.c_lsLoggingChoices, help= ConstantsMicropita.c_strLoggingHelp)
args.add_argument("-c",ConstantsMicropita.c_strCheckedAbundanceFileArgument, dest="ostmCheckedFile", metavar = "output_qc", type = argparse.FileType("w"), help = ConstantsMicropita.c_strCheckedAbundanceFileHelp)
args.add_argument("-g",ConstantsMicropita.c_strLoggingFileArgument, dest="ostmLoggingFile", metavar = "output_log", type = argparse.FileType("w"), help = ConstantsMicropita.c_strLoggingFileHelp)
args.add_argument("-u",ConstantsMicropita.c_strSupervisedInputFile, dest="ostmInputPredictFile", metavar = "output_scaled", type = argparse.FileType("w"), help = ConstantsMicropita.c_strSupervisedInputFileHelp)
args.add_argument("-p",ConstantsMicropita.c_strSupervisedPredictedFile, dest="ostmPredictFile", metavar = "output_labels", type = argparse.FileType("w"), help = ConstantsMicropita.c_strSupervisedPredictedFileHelp)

argp.add_argument("istmInput", metavar = "input.txt", type = argparse.FileType("rU"), help = ConstantsMicropita.c_strAbundanceFileHelp,
	default = sys.stdin)
argp.add_argument("ostmOutput", metavar = "output.txt", type = argparse.FileType("w"), help = ConstantsMicropita.c_strGenericOutputDataFileHelp,
	default = sys.stdout)

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
	args = argp.parse_args( )

	#Set up logger
	iLogLevel = getattr(logging, args.strLogLevel.upper(), None)
	logging.basicConfig(stream = args.ostmLoggingFile if args.ostmLoggingFile else sys.stderr, filemode = 'w', level=iLogLevel)

	#Run micropita
	logging.info("MicroPITA:: Start microPITA")
	microPITA = MicroPITA()

	#Argparse will append to the default but will not remove the default so I do this here
	if not len(args.lstrMethods):
		args.lstrMethods = [ConstantsMicropita.c_strRepresentative]

	dictSelectedSamples = microPITA.funcRun(
		strIDName		= args.strIDName,
		strLastMetadataName	= args.strLastMetadataName,
		istmInput		= args.istmInput,
		ostmInputPredictFile	= args.ostmInputPredictFile,
		ostmPredictFile		= args.ostmPredictFile,
		ostmCheckedFile		= args.ostmCheckedFile,
		ostmOutput		= args.ostmOutput,
		cDelimiter		= args.cFileDelimiter,
		cFeatureNameDelimiter	= args.cFeatureNameDelimiter,
		istmFeatures		= args.istmFeatures,
		strFeatureSelection	= args.strFeatureSelection,
		iCount			= args.iCount,
		strLabel		= args.strLabel,
		strStratify		= args.strUnsupervisedStratify,
		strCustomAlpha		= args.strAlphaDiversity,
		strCustomBeta		= args.strBetaDiversity,
		strAlphaMetadata	= args.strAlphaMetadata,
		istmBetaMatrix		= args.istmBetaMatrix,
		istrmTree		= args.istrmTree,
		istrmEnvr		= args.istrmEnvr,
		lstrMethods		= args.lstrMethods,
		fInvertDiversity	= args.fInvertDiversity
	)

	if not dictSelectedSamples:
		logging.error("MicroPITA:: Error, did not get a result from analysis.")
		return -1
	logging.info("End microPITA")

	#Log output for debugging
	logging.debug("MicroPITA:: Returned the following samples:"+str(dictSelectedSamples))

	#Write selection to file
	microPITA.funcWriteSelectionToFile(dictSelection=dictSelectedSamples, xOutputFilePath=args.ostmOutput)

if __name__ == "__main__":
	_main( )
