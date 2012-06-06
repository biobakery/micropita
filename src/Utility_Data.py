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
import csv
from Diversity import Diversity
from MicroPITA import MicroPITA
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
    def generateAbundanceTable(strOutputFile, strSampleClassification, iScalingFactorForSampleAmount = 1, iRandomMax=5):
        #Matrix descriptors
        #Samples
        #Number of diversity samples
        diversitySampleCount = 16 * iScalingFactorForSampleAmount
        #Number of representative dissimilarity samples
        representativeDissimilarityCount = 14 *iScalingFactorForSampleAmount
        #Number of extreme dissimilarity samples
        extremeDissimilarityCount = 14 * iScalingFactorForSampleAmount
        #Number of taxa-driven samples
        taxaDrivenCount = 4 * iScalingFactorForSampleAmount
        #Total sample count
        sampleCount = diversitySampleCount + representativeDissimilarityCount + extremeDissimilarityCount + taxaDrivenCount

        #Taxa
        #Number of blocks of dissimlarity
        representiveDiversityBlocks = representativeDissimilarityCount
        #Number of taxa in each block of dissimilarity
        representiveDiversityTaxa = 16 * iScalingFactorForSampleAmount
        #Number of blocks of extreme dissimilarity
        extremeDissimilarityBlocks = representiveDiversityBlocks
        #Number of taxa in each block of extreme dissimilarity
        extremeDissimilarityTaxa = 8 * iScalingFactorForSampleAmount
        #Number of low abundance taxa in extreme dissimilarity
        extremeDissimilarityLowTaxa = representiveDiversityTaxa
        #Number of taxa
        taxaCount = representiveDiversityBlocks*representiveDiversityTaxa
        #Number of taxa with abundance in diversity samples
        diversityMaximalTaxa = int(representiveDiversityTaxa*.3)
        #Number of taxa with minimal abundance in diversity samples
#        diversityMinimalTaxa = int(representiveDiversityTaxa*.1)
        #Buffer of taxa to skip before starting the first extreme dissimilarity block
        iExtremeDissimilarityStart = representiveDiversityTaxa-extremeDissimilarityTaxa

        #Abundance levels
        iDiversityAbundanceMax = 50
        iDiversityAbundanceMin = 0
        iRepresentativeAbundanceMax = 50
        iRepresentativeAbundanceMin = 0
        iExtremeAbundanceMax = 100
        iExtremeAbundanceMin = 0
        iTargetedAbundanceMax = 50
        iTargetedAbundanceMin = 0

        #Taxa selected for taxa driven
        taxaDriverPositions = set([3,43,19])
        
        ####Create matrix with noise
        #Create the matrix of n taxa (row) by n samples (column)
#        dataMatrix = np.array([random.randint(0,iRandomMax) for i in xrange(taxaCount*sampleCount)])
        iRandomMax = float(iRandomMax)
        dataMatrix = np.array([random.random()*iRandomMax for i in xrange(taxaCount*sampleCount)])
        dataMatrix.shape = (taxaCount, sampleCount)

        #Create taxa names
        taxaNames = []
        for taxa in xrange(0,taxaCount):
            taxaNames.append("Taxa_"+str(taxa))
        #Create sample names
        sampleNames = ["ID"]

        #Create labels
        lsLabels = ["Label","Class-One","Class-One","Class-One","Class-Two","Class-One","Class-Two","Class-Two","Class-One",
                           "Class-One","Class-Two","Class-One","Class-Two","Class-Two","Class-Two","Class-Two","Class-One",
                           "Class-One","Class-One","Class-One","Class-One","Class-Two","Class-Two","Class-One","Class-Two",
                           "Class-Two","Class-One","Class-Two","Class-Two","Class-Two","Class-Two","Class-One","Class-One",
                           "Class-One","Class-One","Class-Two","Class-Two","Class-One","Class-Two","Class-Two","Class-One",
                           "Class-Two","Class-Two","Class-Two","Class-Two","Class-One","Class-One","Class-One","Class-One"]

        if iScalingFactorForSampleAmount == 2:
          lsLabels = ["Label","Class-Two","Class-Two","Class-Two","Class-One","Class-Two","Class-One","Class-Two","Class-Two",
                    "Class-Two","Class-Two","Class-One","Class-One","Class-One","Class-Two","Class-One","Class-One","Class-One",
                    "Class-Two","Class-One","Class-One","Class-One","Class-One","Class-Two","Class-One","Class-Two","Class-One",
                    "Class-Two","Class-One","Class-Two","Class-Two","Class-One","Class-Two","Class-Two","Class-Two","Class-One",
                    "Class-One","Class-Two","Class-One","Class-Two","Class-One","Class-One","Class-One","Class-Two","Class-One",
                    "Class-Two","Class-Two","Class-Two","Class-One","Class-One","Class-Two","Class-One","Class-Two","Class-Two",
                    "Class-Two","Class-Two","Class-Two","Class-One","Class-Two","Class-One","Class-One","Class-Two","Class-Two",
                    "Class-Two","Class-One","Class-Two","Class-One","Class-Two","Class-One","Class-One","Class-One","Class-Two",
                    "Class-One","Class-Two","Class-Two","Class-Two","Class-Two","Class-One","Class-One","Class-One","Class-Two",
                    "Class-Two","Class-Two","Class-Two","Class-Two","Class-Two","Class-One","Class-Two","Class-One","Class-Two",
                    "Class-Two","Class-Two","Class-Two","Class-Two","Class-Two","Class-Two","Class-Two"]																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																														
        
        #Keep tract of which samples are being produced in the matrix
        sampleStart = 0
        sampleStop = diversitySampleCount

        #Makes the file that holds the classification for each sample and add keys for sample types
        dictActualData = dict()
        dictActualData[MicroPITA.c_DIVERSITY_1] = []
        dictActualData[MicroPITA.c_DIVERSITY_2] = []
        dictActualData[MicroPITA.c_EXTREME_DISSIMILARITY_1] = []
        dictActualData[MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1] = ["Sample_14_D", "Sample_44_T"]
        dictActualData[MicroPITA.c_USER_RANKED] = ["Sample_16_R","Sample_17_R","Sample_18_R","Sample_32_E"]
        dictActualData[MicroPITA.c_SVM_CLOSE] = []
        dictActualData[MicroPITA.c_SVM_FAR] = []

        ####Create diversity samples
        for diversitySampleIndex in xrange(sampleStart,sampleStop):

            #Keep track for taxa as blocks of diversity are created
            iTaxaBlockStart = 0

            #Append Label and name and add to actual classification dictionary
#            lsLabels.append(sLabelSelection[diversitySampleIndex%2])
            strSampleName = "_".join(["Sample",str(diversitySampleIndex),"D"])
            sampleNames.append(strSampleName)
            dictActualData[MicroPITA.c_DIVERSITY_1].append(strSampleName)
            dictActualData[MicroPITA.c_DIVERSITY_2].append(strSampleName)
            dictActualData[MicroPITA.c_SVM_CLOSE].append(strSampleName)

            #Create diversity with in each block
            for diversityBlock in xrange(0,representiveDiversityBlocks):
                #Define index populations and set abundance per sample
                population = set(range(iTaxaBlockStart,iTaxaBlockStart+representiveDiversityTaxa))
                #Remove taxa drivers
                population = population.difference(taxaDriverPositions)
                #Taxa indices with maximal abundance in diversity samples
                diversityMaximalTaxaPositions = random.sample(population,diversityMaximalTaxa)
                #Set abundance
                for taxonPosition in diversityMaximalTaxaPositions:
                    dataMatrix[taxonPosition,diversitySampleIndex] = iDiversityAbundanceMax+random.randint(0,iRandomMax)
                #Remove already selected position from the potential indices population
                population = population.difference(set(diversityMaximalTaxaPositions))
                #Taxa indices with minimal abundance in diversity samples
#                diversityMinimalTaxaPositions = set()
#                if(len(population)>0):
#                    if(len(population)>diversityMinimalTaxa):
#                        diversityMinimalTaxaPositions = random.sample(population,diversityMinimalTaxa)
#                    else:
#                        diversityMinimalTaxaPositions = population

                #Set abundance
#                for taxonPosition in diversityMinimalTaxaPositions:
#                    dataMatrix[taxonPosition,diversitySampleIndex] = iDiversityAbundanceMin+random.randint(1,iRandomMax)

                #Update block
                iTaxaBlockStart = iTaxaBlockStart + representiveDiversityTaxa

        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + representativeDissimilarityCount

        #Create blocks of dissimilarity
        #Representative dissimilarity
        #Track taxon blocks
        taxonBlockStart = 0
        for dissimilaritySampleIndex in xrange(sampleStart,sampleStop):
            #Append Label and name
#            lsLabels.append(sLabelSelection[dissimilaritySampleIndex%2])
            strSampleName = "_".join(["Sample",str(dissimilaritySampleIndex),"R"])
            sampleNames.append(strSampleName)
            dictActualData[MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1].append(strSampleName)
            dictActualData[MicroPITA.c_SVM_FAR].append(strSampleName)

            for blockTaxonPosition in xrange(taxonBlockStart,taxonBlockStart+representiveDiversityTaxa):
                dataMatrix[blockTaxonPosition,dissimilaritySampleIndex] = iRepresentativeAbundanceMax+random.randint(0,iRandomMax)
            taxonBlockStart = taxonBlockStart+representiveDiversityTaxa

        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + extremeDissimilarityCount

        #Extreme dissimilarity
        extremeTaxaBlockIncrement = representiveDiversityTaxa-extremeDissimilarityTaxa
        taxonBlockStart = iExtremeDissimilarityStart
        taxonBlockStop = taxonBlockStart + extremeTaxaBlockIncrement
        for extremeSampleIndex in xrange(sampleStart,sampleStop):
            #Append Label and name
#            lsLabels.append(sLabelSelection[extremeSampleIndex%2])
            strSampleName = "_".join(["Sample",str(extremeSampleIndex),"E"])
            sampleNames.append(strSampleName)
            dictActualData[MicroPITA.c_EXTREME_DISSIMILARITY_1].append(strSampleName)
            dictActualData[MicroPITA.c_SVM_FAR].append(strSampleName)

            #Set block abundance
            for extremeBlockTaxonPosition in xrange(taxonBlockStart,taxonBlockStop):
                dataMatrix[extremeBlockTaxonPosition,extremeSampleIndex] = iExtremeAbundanceMax+random.randint(0,iRandomMax)
            #Define index populations and set abundance per sample (for noise)
            population = set(range(taxaCount))
            #Remove taxa drivers
            population = population.difference(taxaDriverPositions)
            #Remove current block
            population = population.difference(set(range(taxonBlockStart,taxonBlockStop)))
            #Generate noise population for extreme dissimilarity
#            extremeNoisePopulation = random.sample(population,extremeDissimilarityLowTaxa)
            #Set extreme noise
#            for extremeNoiseTaxonPosition in extremeNoisePopulation:
#                dataMatrix[extremeNoiseTaxonPosition,extremeSampleIndex] = iExtremeAbundanceMin+random.randint(1,iRandomMax)
            #Increment start and stop
            taxonBlockStart = taxonBlockStop+extremeTaxaBlockIncrement
            taxonBlockStop = taxonBlockStart + extremeDissimilarityTaxa

        #Update sample bounds
        sampleStart = sampleStop
        sampleStop = sampleStop + taxaDrivenCount

        #Create the remaining taxa driven samples
        #Update sample bounds
        for drivenSampleIndex in xrange(sampleStart,sampleStop):
            #Append Label and name
#            lsLabels.append(sLabelSelection[drivenSampleIndex%2])
            strSampleName = "_".join(["Sample",str(drivenSampleIndex),"T"])
            sampleNames.append(strSampleName)
            dictActualData[MicroPITA.c_USER_RANKED].append(strSampleName)
            dictActualData[MicroPITA.c_SVM_FAR].append(strSampleName)

            #Define index populations and set abundance per sample (for noise)
            population = set(range(taxaCount))
            #Remove taxa drivers
            population = population.difference(taxaDriverPositions)
            #Set sample abundance
            for drivingTaxonPosition in taxaDriverPositions:
                dataMatrix[drivingTaxonPosition,drivenSampleIndex] = iTargetedAbundanceMax+random.randint(0,iRandomMax)
            #Generate low noise
#            noisePositions = random.sample(population,extremeDissimilarityLowTaxa)
#            for noiseTaxonPosition in noisePositions:
#                dataMatrix[noiseTaxonPosition,drivenSampleIndex] = iTargetedAbundanceMin+random.randint(1,iRandomMax)

        #Write to file
        #Delete current file before writing
        if(os.path.isfile(strOutputFile)):
            os.remove(strOutputFile)

        #Data to output
        outputContents = []
        #Output sample id header
        outputContents.append(Constants.TAB.join(sampleNames))
        #Output sample labels
        outputContents.append(Constants.TAB.join(lsLabels))
        #Write rows/taxa
        for row in xrange(0,taxaCount):
            taxaData = []
            for sampleAbundance in dataMatrix[row]:
                taxaData.append(str(sampleAbundance))
            taxaData = Constants.TAB.join(taxaData)
            outputContents.append(Constants.TAB.join([taxaNames[row],taxaData]))
        with open(strOutputFile,'a') as f:
            f.write(Constants.ENDLINE.join(outputContents))
        f.close()

        #Write actual file
        strOutputContent = ""
        for sKey in dictActualData:
            strOutputContent = "".join([strOutputContent,sKey,Constants.COLON,", ".join(dictActualData[sKey]),Constants.ENDLINE])

        #Write to file
        if not strSampleClassification == None:
            if(not strOutputContent == ""):
                fHndlOutput = open(strSampleClassification,'w')
                fHndlOutput.write(str(strOutputContent))
                fHndlOutput.close()

        #Return file name
        return strOutputFile

    ##
    #Generate matrix of random data
    #
    @staticmethod
    def generateRandomMatrix(tempFilePath,iNumberRows,iNumberColumns, iMinValue, iMaxValue, charDelimiter = Constants.TAB):
        dataMatrix = np.zeros([iNumberRows,iNumberColumns])

        #Generate random data
        for row in xrange(0, iNumberRows):
            for col in xrange(0, iNumberColumns):
                dataMatrix[row,col] = random.randint(iMinValue, iMaxValue)
        
        #Write to file
        fileOutput = csv.writer(open(tempFilePath,'w'), delimiter=charDelimiter)
        fileOutput.writerow(["_".join(["Label",str(isample)]) for isample in xrange(0,iNumberColumns)])
        #Write rows/taxa
        for row in dataMatrix:
            fileOutput.writerow(row)
        fileOutput.close()
      
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

        #First data row for HCLust
        strFirstDataRow = "3"
        strFirstDataCol = "2"
        
        ####Create matrix
        #Create the matrix of n taxa /(row) by n samples (column)
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

    @staticmethod
    def funcGenerateCorrelatedFeaturesDataSet(sOutputFile, dFeatureScale, dSampleScale):

        #Make sure scales are ints
        dFeatureScale = int(dFeatureScale)
        dSampleScale = int(dSampleScale)

        lsCorrelation = []

        dSampleCount = 10*dSampleScale
        dFeatureCount = 9*dFeatureScale

        iFeatureIndex = 1

        #Create correlated, anticorrelated, and flat but scaled features
        for feature in xrange(dFeatureCount):
            ldFeature = []
            dCurrentFeatureScale = math.pow(2.0,feature)
            for dSample in xrange(dSampleCount):
                ldFeature.append(dSample*dCurrentFeatureScale)
            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+[str(feature) for feature in ldFeature[:]])
            iFeatureIndex = iFeatureIndex + 1
            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+([str(feature)]*dSampleCount))
            ldFeature.reverse()
            iFeatureIndex = iFeatureIndex + 1
            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+[str(feature) for feature in ldFeature[:]])
            iFeatureIndex = iFeatureIndex + 1
            dCurrentFeatureScale = math.pow(2.0,feature)

        #Update contents to a line
        sContents = "FeatureID"+Constants.TAB+Constants.TAB.join(["".join(["Sample_",str(iSampleIndex)]) for iSampleIndex in xrange(dSampleCount)])
        sContents = Constants.ENDLINE.join([sContents]+[Constants.TAB.join(feature) for feature in lsCorrelation])

        #Write line to a file
        with open(sOutputFile,'w') as f:
            f.write(sContents)
        f.close()
        return

for iRandom in [5]:#(10,15,20,25):
  for i in xrange(1,11):
    Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced96-random"+str(iRandom)+"-"+str(i)+".pcl", strSampleClassification="Unbalanced96-random-"+str(iRandom)+"-"+str(i)+"Actual.txt", iScalingFactorForSampleAmount = 2, iRandomMax = iRandom)
#    Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced-random"+str(iRandom)+"-"+str(i)+".pcl", strSampleClassification="Unbalanced-random-"+str(iRandom)+"-"+str(i)+"Actual.txt", iScalingFactorForSampleAmount = 1, iRandomMax = iRandom)
#Utility_Data.generateRandomMatrix(tempFilePath="Random.txt",iNumberRows=10,iNumberColumns=6, iMinValue=0, iMaxValue=2500)
