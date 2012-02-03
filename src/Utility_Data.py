#######################################################
#
#	Title:		Utility_Data
#	Author:		Timothy Tickle
#	Date:		12/12/2011
#	Purpose:	Utility class for Data generation.
#
#######################################################

#Import libaries
from Constants import Constants
from CommandLine import CommandLine
from FileIO import FileIO
import math
import numpy as np
import os
import random
import ValidateData

##
#Utility function for data generation
class Utility_Data():

    ##
    #Contructor
    def __init__(): pass

    ##
    #Generate matrix for microPITA
    #@param tempOutPutFile
    @staticmethod
    def generateAbundanceTable(tempFilePath, tempGenerateFigure=True):
        #Matrix descriptors
        #Samples
        #Number of diversity samples
        diversitySampleCount = 8
        #Number of representative dissimilarity samples
        representativeDissimilarityCount = 7
        #Number of extreme dissimilarity samples
        extremeDissimilarityCount = 7
        #Number of taxa-driven samples
        taxaDrivenCount = 2
        #Total sample count
        sampleCount = diversitySampleCount + representativeDissimilarityCount + extremeDissimilarityCount + taxaDrivenCount

        #Taxa
        #Number of taxa
        taxaCount = 56
        #Number of taxa with abundance in diversity samples
        diversityMaximalTaxa = 25
        #Number of taxa with minimal abundance in diversity samples
        diversityMinimalTaxa = 6
        #Number of blocks of dissimlarity
        representiveDiversityBlocks = 7
        #Number of taxa in each block of dissimilarity
        representiveDiversityTaxa = 8
        #Number of blocks of extreme dissimilarity
        extremeDissimilarityBlocks = 7
        #Number of taxa in each block of extreme dissimilarity
        extremeDissimilarityTaxa = 3
        #Number of low abundance taxa in extreme dissimilarity
        extremeDissimilarityLowTaxa = 12

        #Taxa selected for taxa driven
        taxaDriverPositions = set([3,43,19])
        
        ####Create matrix
        #Create the matrix of n taxa (row) by n samples (column)
        dataMatrix = np.zeros((taxaCount,sampleCount))
        #Create taxa names
        taxaNames = []
        for taxa in xrange(0,taxaCount):
            taxaNames.append("Taxa_"+str(taxa))
        #Create sample names
        sampleNames = ["TID"]
        for sample in xrange(0,sampleCount):
            sampleNames.append("Sample_"+str(sample))
        #Create labels
        #TODO


        #Keep tract of which samples are being produced in the matrix
        sampleStart = 0
        sampleStop = diversitySampleCount

        ####Create diversity samples
        for diversitySampleIndex in xrange(sampleStart,sampleStop):
            #Define index populations and set abundance per sample
            population = set(range(taxaCount))
            #Remove taxa drivers
            population = population.difference(taxaDriverPositions)
            #Taxa indices with maximal abundance in diversity samples
            diversityMaximalTaxaPositions = random.sample(population,diversityMaximalTaxa)
            #Set abundance
            for taxonPosition in diversityMaximalTaxaPositions:
                dataMatrix[taxonPosition,diversitySampleIndex] = 20+random.randint(0,5)
            #Remove already selected position from the potential indices population
            population = population.difference(set(diversityMaximalTaxaPositions))
            #Taxa indices with minimal abundance in diversity samples
            diversityMinimalTaxaPositions = random.sample(population,diversityMinimalTaxa)
            #Set abundance
            for taxonPosition in diversityMinimalTaxaPositions:
                dataMatrix[taxonPosition,diversitySampleIndex] = 0+random.randint(1,5)
        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + representativeDissimilarityCount

        #Create blocks of dissimilarity
        #Representative dissimilarity
        #Track taxon blocks
        taxonBlockStart = 0
        for dissimilaritySampleIndex in xrange(sampleStart,sampleStop):
            for blockTaxonPosition in xrange(taxonBlockStart,taxonBlockStart+representiveDiversityTaxa):
                dataMatrix[blockTaxonPosition,dissimilaritySampleIndex] = 50+random.randint(0,5)
            taxonBlockStart = taxonBlockStart+representiveDiversityTaxa

        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + extremeDissimilarityCount

        #Extreme dissimilarity
        taxonBlockStart = 3
        taxonBlockStop = taxonBlockStart + extremeDissimilarityTaxa
        extremeTaxaBlockIncrement = (taxaCount/extremeDissimilarityBlocks)-extremeDissimilarityTaxa
        for extremeSampleIndex in xrange(sampleStart,sampleStop):
            #Set block abundance
            for extremeBlockTaxonPosition in xrange(taxonBlockStart,taxonBlockStop):
                dataMatrix[extremeBlockTaxonPosition,extremeSampleIndex] = 100+random.randint(0,5)
            #Define index populations and set abundance per sample (for noise)
            population = set(range(taxaCount))
            #Remove taxa drivers
            population = population.difference(taxaDriverPositions)
            #Remove current block
            population = population.difference(set(range(taxonBlockStart,taxonBlockStop)))
            #Generate noise population for extreme dissimilarity
            extremeNoisePopulation = random.sample(population,extremeDissimilarityLowTaxa)
            #Set extreme noise
            for extremeNoiseTaxonPosition in extremeNoisePopulation:
                dataMatrix[extremeNoiseTaxonPosition,extremeSampleIndex] = 0+random.randint(1,5)
            #Increment start and stop
            taxonBlockStart = taxonBlockStop+extremeTaxaBlockIncrement
            taxonBlockStop = taxonBlockStart + extremeDissimilarityTaxa

        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + taxaDrivenCount

        #Create the remaining taxa driven samples
        #Update sample bounds
        for drivenSampleIndex in xrange(sampleStart,sampleStop):
            #Define index populations and set abundance per sample (for noise)
            population = set(range(taxaCount))
            #Remove taxa drivers
            population = population.difference(taxaDriverPositions)
            #Set sample abundance
            for drivingTaxonPosition in taxaDriverPositions:
                dataMatrix[drivingTaxonPosition,drivenSampleIndex] = 50+random.randint(0,5)
            #Generate low noise
            noisePositions = random.sample(population,extremeDissimilarityLowTaxa)
            for noiseTaxonPosition in noisePositions:
                dataMatrix[noiseTaxonPosition,drivenSampleIndex] = 0+random.randint(1,5)

        #Write to file
        if(os.path.isfile(tempFilePath)):
            os.remove(tempFilePath)
        #File writer
        outToFile = FileIO(tempFilePath,False,True,True)
        #Data to output
        outputContents = []
        #Output sample id header
        outputContents.append(Constants.TAB.join(sampleNames))
        #Write rows/taxa
        for row in xrange(0,taxaCount):
            taxaData = []
            for sampleAbundance in dataMatrix[row]:
                taxaData.append(str(sampleAbundance))
            taxaData = Constants.TAB.join(taxaData)
            outputContents.append(Constants.TAB.join([taxaNames[row],taxaData]))
        outToFile.writeToFile(Constants.ENDLINE.join(outputContents))
        outToFile.close()

        #
        if(tempGenerateFigure):
            CommandLine().runCommandLine(["./external/hclust/hclust.py", "--in", tempFilePath, "--out", tempFilePath[:-4]+".png", "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1",'--flabel',"1",'--font_size',"5"])

        #Return file name
        return tempFilePath

Utility_Data.generateAbundanceTable(tempFilePath="TestDataGeneration.txt")
