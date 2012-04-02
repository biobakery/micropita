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
from Constants_Arguments import Constants_Arguments
from CommandLine import CommandLine
from Diversity import Diversity
from FileIO import FileIO
import numpy as np
import os
import random
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
    def generateAbundanceTable(tempFilePath):
        #Matrix descriptors
        #Samples
        #Number of diversity samples
        diversitySampleCount = 16#8
        #Number of representative dissimilarity samples
        representativeDissimilarityCount = 14#7
        #Number of extreme dissimilarity samples
        extremeDissimilarityCount = 14#7
        #Number of taxa-driven samples
        taxaDrivenCount = 4#2
        #Total sample count
        sampleCount = diversitySampleCount + representativeDissimilarityCount + extremeDissimilarityCount + taxaDrivenCount

        #TODO remove
        strHCLLoc = "./external/hclust/hclust.py"

        #Taxa
        #Number of blocks of dissimlarity
        representiveDiversityBlocks = representativeDissimilarityCount
        #Number of taxa in each block of dissimilarity
        representiveDiversityTaxa = 8
        #Number of blocks of extreme dissimilarity
        extremeDissimilarityBlocks = representiveDiversityBlocks
        #Number of taxa in each block of extreme dissimilarity
        extremeDissimilarityTaxa = 5#int(representiveDiversityBlocks*.5)
        #Number of low abundance taxa in extreme dissimilarity
        extremeDissimilarityLowTaxa = representiveDiversityTaxa
        #Number of taxa
        taxaCount = representiveDiversityBlocks*representiveDiversityTaxa
        #Number of taxa with abundance in diversity samples
        diversityMaximalTaxa = int(taxaCount*.4)
        #Number of taxa with minimal abundance in diversity samples
        diversityMinimalTaxa = int(taxaCount*.1)
        #Buffer of taxa to skip before startign the first extreme dissimilarity block
        iExtremeDissimilarityStart = representiveDiversityTaxa-extremeDissimilarityTaxa

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
        extremeTaxaBlockIncrement = representiveDiversityTaxa-extremeDissimilarityTaxa
        taxonBlockStart = iExtremeDissimilarityStart
        taxonBlockStop = taxonBlockStart + extremeTaxaBlockIncrement
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
        #Delete current file before writing
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

        #Return file name
        return tempFilePath

    ##
    #Generate matrix for microPITA focused on diversity
    #@param tempOutPutFile
    @staticmethod
    def generateDiversityAbundanceTable(tempFilePath):
        #Matrix descriptors
        #Samples
        #Number of diversity samples
        iSampleCount = 50

        #Taxa
        #Number of taxa
        iTaxaCount = 250

        #Chao1 
        strChao1 = "Chao 1"
        #Inverse Simpson
        strISimpson = "Inverse_Simpson"

        #TODO remove
        strHCLLoc = "./external/hclust/hclust.py"

        #First data row for HCLust
        strFirstDataRow = "3"
        strFirstDataCol = "2"
        
        ####Create matrix
        #Create the matrix of n taxa (row) by n samples (column)
        dataMatrix = np.zeros((iTaxaCount,iSampleCount+1))
        #Create taxa names
        lsTaxaNames = []
        for iTaxaIndex in xrange(0,iTaxaCount):
            lsTaxaNames.append("".join(["Taxa_",str(iTaxaIndex)]))
        #Create sample names
        lsSampleNames = ["ID"]
        for iSampleIndex in xrange(0,iSampleCount+1):
            lsSampleNames.append("".join(["Sample_",str(iSampleIndex)]))
        #Create labels based on Chao 1 diversity
        lsChaoLabels = [strChao1]
        #Create labels based on Inverse Simpson diversity
        lsISLabels = [strISimpson]

        #Keep tract of which samples are being produced in the matrix
        iSampleStart = 0
        iSampleStop = iSampleCount+1

        ####Create diversity samples
        for iDiversitySampleIndex in xrange(iSampleStart,iSampleStop):
            #Define index populations and set abundance per sample
            liPopulation = set(range(iTaxaCount))
            #Taxa indices with maximal abundance in diversity samples
            liSelection = random.sample(liPopulation,int(iDiversitySampleIndex*(iTaxaCount/iSampleCount)))
            #Set abundance
            for iTaxonPosition in liSelection:
                dataMatrix[iTaxonPosition,iDiversitySampleIndex] = 50+random.randint(0,5)
            #Estimate diversity and set as a label
            lsChaoLabels.append(str(Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies=dataMatrix[:,iDiversitySampleIndex])))
            if(sum(dataMatrix[:,iDiversitySampleIndex])==0):
                lsISLabels.append("0")
            else:
                lsISLabels.append(str(Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=dataMatrix[:,iDiversitySampleIndex])))

        #Write to file
        #Delete current file before writing
        if(os.path.isfile(tempFilePath)):
            os.remove(tempFilePath)

        #Data to output
        outputContents = []
        #Output sample id header
        outputContents.append(Constants.TAB.join(lsSampleNames))
        #Output Chao labels
        outputContents.append(Constants.TAB.join(lsChaoLabels))
        #Output InvSimpson labels
        outputContents.append(Constants.TAB.join(lsISLabels))
        #Write rows/taxa
        for row in xrange(0,iTaxaCount):
            taxaData = []
            for sampleAbundance in dataMatrix[row]:
                taxaData.append(str(sampleAbundance))
            taxaData = Constants.TAB.join(taxaData)
            outputContents.append(Constants.TAB.join([lsTaxaNames[row],taxaData]))

        #Write to file
        with open(tempFilePath,'w') as f:
            sFileContents = f.write(Constants.ENDLINE.join(outputContents))
            f.close()

        #Return file name
        return tempFilePath
