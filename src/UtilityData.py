"""
Author: Timothy Tickle
Description: Utility class for data generation.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libaries
from ConstantsMicropita import ConstantsMicropita
from breadcrumbs.CommandLine import CommandLine
import csv
from breadcrumbs.Metric import Metric
import math
#from MicroPITA import MicroPITA
import numpy as np
import os
import random

class UtilityData():
    """
    Utility function for data generation.
    """

    ##
    #Contructor
    def __init__(): pass

    @staticmethod
    def generateAbundanceTable(strOutputFile, strSampleClassification, iScalingFactorForSampleAmount=1, dMaxGeneralNoise = 0.0, dMaxSignalNoise=5.0, dSimpleNoise=0.0):
        """
        dMaxGeneralNoise is an absolute number of reads for noise a number starting with 1
        dMaxSignalNoise is an absolute number of reads for noise a number starting with 1
        dSimpleNoise should be between 0 and 1 (percentage noise)
        """

        #If iSimpleNoise is greater than 0 then the simple noise model will be used over the sparse noise model
        #The simple noise model will uniformly add up to the iSimpleNoise value to all data in the data set
        #The simple noise model will ignore the dMaxGeneralNoise and dMaxSignalNoise variables /  noise structures
        if dSimpleNoise > 0.0:
#            dMaxGeneralNoise = 0.0
            dMaxSignalNoise = 0.0

        #Matrix descriptors
        #Samples
        #Number of diversity samples
        iDiversitySampleCount = 16 * iScalingFactorForSampleAmount
        #Number of representative dissimilarity samples
        iRepresentativeDissimilarityCount = 14 *iScalingFactorForSampleAmount
        #Number of extreme dissimilarity samples
        iExtremeDissimilarityCount = 14 * iScalingFactorForSampleAmount
        #Number of taxa-driven samples
        iTaxaDrivenCount = 4 * iScalingFactorForSampleAmount
        #Total sample count
        iSampleCount = iDiversitySampleCount + iRepresentativeDissimilarityCount + iExtremeDissimilarityCount + iTaxaDrivenCount

        #Taxa
        #Number of blocks of dissimlarity
        iRepresentiveDissimilarityBlocks = iRepresentativeDissimilarityCount
        #Number of taxa in each block of dissimilarity
        iRepresentiveDissimilarityTaxa = 16 * iScalingFactorForSampleAmount
        #Number of blocks of extreme dissimilarity
        iExtremeDissimilarityBlocks = iRepresentiveDissimilarityBlocks
        #Number of taxa in each block of extreme dissimilarity
        iExtremeDissimilarityTaxa = 8 * iScalingFactorForSampleAmount
        #Taxa selected for taxa driven
        liTaxaDriverPositions = set([3,43,19])
        #Number of features in targeted feature
        iTargetedFeature = len(liTaxaDriverPositions)
        #Number of taxa
        iTaxaCount = iRepresentiveDissimilarityBlocks*iRepresentiveDissimilarityTaxa
        #Number of taxa with abundance in diversity samples
        #This is based on teh representative group because the codes is made in a way that
        #The diversity populations of each sample are created per representative block so
        #that per sample the diversity of features associated with certain block are consistent
        iDiversityMaximalTaxa = int(iRepresentiveDissimilarityTaxa*.3)

        #Buffer of taxa to skip before starting the first extreme dissimilarity block
        iExtremeDissimilarityStart = iRepresentiveDissimilarityTaxa-iExtremeDissimilarityTaxa

        #Noise populations
        #Number of low abundance noise in representative dissimilarity samples
        iRepresentativeDissimilarityLowTaxa = int((iTaxaCount-iExtremeDissimilarityTaxa)*.1)
        #Number of taxa with minimal abundance in diversity samples
#        diversityMinimalTaxa = int(iRepresentiveDissimilarityTaxa*.1)
        iDiversityMinimalTaxa = int((iTaxaCount-(iDiversityMaximalTaxa*iRepresentiveDissimilarityBlocks))*.1)
        #Number of low abundance noise in extreme dissimilarity samples
        iExtremeDissimilarityLowTaxa = int((iTaxaCount-iExtremeDissimilarityTaxa)*.1)
        #Number of low abundance noise features in targeted feature samples
        iTargetedFeatureLowTaxa = int((iTaxaCount-iTargetedFeature)*.1)

        #Abundance levels
        iDiversityAbundanceMax = 50
        iDiversityAbundanceMin = 0
        iRepresentativeAbundanceMax = 50
        iRepresentativeAbundanceMin = 0
        iExtremeAbundanceMax = 100
        iExtremeAbundanceMin = 0
        iTargetedAbundanceMax = 50
        iTargetedAbundanceMin = 0

        ####Create matrix, potentially with noise with noise
        #Create the matrix of n taxa (row) by n samples (column)
        npDataMatrix = None
        npDataMatrix = np.array([0.0]*(iTaxaCount*iSampleCount))
        npDataMatrix.shape = (iTaxaCount, iSampleCount)

        #Create taxa names
        lsTaxaNames = []
        for taxa in xrange(0,iTaxaCount):
            lsTaxaNames.append("Taxa_"+str(taxa))
        #Create sample names
        lsSampleNames = ["ID"]

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

        #Makes the file that holds the classification for each sample and add keys for sample types
        strActualData = ""
        if iScalingFactorForSampleAmount == 1:
            strActualData = "".join(["[Class]Diverse:Sample_0_D, Sample_1_D, Sample_2_D, Sample_3_D, Sample_4_D, Sample_5_D, Sample_6_D,",
                                     " Sample_7_D, Sample_8_D, Sample_9_D, Sample_10_D, Sample_11_D, Sample_12_D, Sample_13_D, Sample_14_D,",
                                     " Sample_15_D\n[Class]Dissimilar:Sample_16_R, Sample_17_R, Sample_18_R, Sample_19_R, Sample_20_R,",
                                     " Sample_21_R, Sample_22_R, Sample_23_R, Sample_24_R, Sample_25_R, Sample_26_R, Sample_27_R, Sample_28_R,",
                                     " Sample_29_R\n[Class]Extreme:Sample_30_E, Sample_31_E, Sample_32_E, Sample_33_E, Sample_34_E, Sample_35_E,",
                                     " Sample_36_E, Sample_37_E, Sample_38_E, Sample_39_E, Sample_40_E, Sample_41_E, Sample_42_E, Sample_43_E\n",
                                     "[Class]Taxa:Sample_44_T, Sample_45_T, Sample_46_T, Sample_47_T\n\n",
                                     "Diversity_C:[Class]Diverse\nDiversity_I:[Class]Diverse\nDiscriminant:[Class]Diversen\n",
                                     "Distinct:[Class]Dissimilar,[Class]Extreme,[Class]Taxa\nExtreme_B:[Class]Extreme,[Class]Taxa\n",
                                     "Representative_B:[EVEN]\nTaxa_Defined:[Class]Taxa, Sample_16_R, Sample_17_R, Sample_18_R, Sample_32_E"])
        if iScalingFactorForSampleAmount == 2:
            strActualData = "".join(["[Class]Diverse:Sample_0_D, Sample_1_D, Sample_2_D, Sample_3_D, Sample_4_D, Sample_5_D, Sample_6_D, Sample_7_D,",
                                     " Sample_8_D, Sample_9_D, Sample_10_D, Sample_11_D, Sample_12_D, Sample_13_D, Sample_14_D, Sample_15_D, Sample_16_D,",
                                     " Sample_17_D, Sample_18_D, Sample_19_D, Sample_20_D, Sample_21_D, Sample_22_D, Sample_23_D, Sample_24_D, Sample_25_D,",
                                     " Sample_26_D, Sample_27_D, Sample_28_D, Sample_29_D, Sample_30_D, Sample_31_D\n[Class]Dissimilar:Sample_32_R, Sample_33_R,",
                                     " Sample_34_R, Sample_35_R, Sample_36_R, Sample_37_R, Sample_38_R, Sample_39_R, Sample_40_R, Sample_41_R, Sample_42_R,",
                                     " Sample_43_R, Sample_44_R, Sample_45_R, Sample_46_R, Sample_47_R, Sample_48_R, Sample_49_R, Sample_50_R, Sample_51_R,",
                                     " Sample_52_R, Sample_53_R, Sample_54_R, Sample_55_R, Sample_56_R, Sample_57_R, Sample_58_R, Sample_59_R\n",
                                     "[Class]Extreme:Sample_60_E, Sample_61_E, Sample_62_E, Sample_63_E, Sample_64_E, Sample_65_E, Sample_66_E,",
                                     " Sample_67_E, Sample_68_E, Sample_69_E, Sample_70_E, Sample_71_E, Sample_72_E, Sample_73_E, Sample_74_E, Sample_75_E,",
                                     " Sample_76_E, Sample_77_E, Sample_78_E, Sample_79_E, Sample_80_E, Sample_81_E, Sample_82_E, Sample_83_E, Sample_84_E,",
                                     " Sample_85_E, Sample_86_E, Sample_87_E\n[Class]Taxa:Sample_88_T, Sample_89_T, Sample_90_T, Sample_91_T,",
                                     " Sample_92_T, Sample_93_T, Sample_94_T, Sample_95_T\n\nDiversity_C:[Class]Diverse\nDiversity_I:[Class]Diverse",
                                     "Discriminant:[Class]Diverse\nDistinct:[Class]Dissimilar,[Class]Extreme,[Class]Taxa\nExtreme_B:[Class]Extreme,[Class]Taxa\n",
                                     "Representative_B:[EVEN]\nTaxa_Defined:[Class]Taxa, Sample_32_R, Sample_33_R, Sample_60_E"])

        #Keep tract of which samples are being produced in the matrix
        #Stop after all diversity samples are generated
        iSampleStart = 0
        iSampleStop = iDiversitySampleCount

        ####Create diversity samples
        for iDiversitySampleIndex in xrange(iSampleStart,iSampleStop):

            #Keep track of taxa as blocks of diversity are created
            iTaxaBlockStart = 0

            #Append name and add to actual classification dictionary
            strSampleName = "_".join(["Sample",str(iDiversitySampleIndex),"D"])
            lsSampleNames.append(strSampleName)

            #Features not indicated as max diversity features in the sample
            #Noise samples will be drawn from this population
            setDiversityNoiseFeatures = set()

            #Create max diversity samples
            #Create diversity with in each block
            for diversityBlock in xrange(0,iRepresentiveDissimilarityBlocks):
                #Define index populations and set abundance per sample
                setPopulation = set(range(iTaxaBlockStart,iTaxaBlockStart+iRepresentiveDissimilarityTaxa))
                #Remove targeted features
                setPopulation = setPopulation.difference(liTaxaDriverPositions)
                #Taxa indices with maximal abundance in diversity samples
                liDiversityMaximalTaxaPositions = random.sample(setPopulation,iDiversityMaximalTaxa)
                #Set abundance
                for taxonPosition in liDiversityMaximalTaxaPositions:
                    npDataMatrix[taxonPosition,iDiversitySampleIndex] = iDiversityAbundanceMax+(random.random()*dMaxSignalNoise)
                #Remove all diversityMaximalTaxa already selected position from the potential indices setPopulation
                setDiversityNoiseFeatures = setDiversityNoiseFeatures | setPopulation.difference(set(liDiversityMaximalTaxaPositions))

                #Update block
                iTaxaBlockStart = iTaxaBlockStart + iRepresentiveDissimilarityTaxa

            #Make noisy samples
            setDiversityMinimalTaxaPositions = set()
            if(len(setDiversityNoiseFeatures)>0):
                if(len(setDiversityNoiseFeatures)>iDiversityMinimalTaxa):
                    setDiversityMinimalTaxaPositions = random.sample(setDiversityNoiseFeatures,iDiversityMinimalTaxa)
                else:
                    setDiversityMinimalTaxaPositions = setDiversityNoiseFeatures

            #Set abundance
            for iTaxonPosition in setDiversityMinimalTaxaPositions:
                npDataMatrix[iTaxonPosition,iDiversitySampleIndex] = iDiversityAbundanceMin+(random.random()*dMaxGeneralNoise)

        #Update sample bounds to generate the representative samples
        iSampleStart = iSampleStop
        iSampleStop = iSampleStop + iRepresentativeDissimilarityCount

        #Create blocks of dissimilarity
        #Representative dissimilarity
        #Track taxon blocks
        taxonBlockStart = 0
        for dissimilaritySampleIndex in xrange(iSampleStart,iSampleStop):
            #Append name
            strSampleName = "_".join(["Sample",str(dissimilaritySampleIndex),"R"])
            lsSampleNames.append(strSampleName)

            #Create max representative features in sample
            setBlockFeatures = set()
            for iBlockTaxonPosition in xrange(taxonBlockStart,taxonBlockStart+iRepresentiveDissimilarityTaxa):
                npDataMatrix[iBlockTaxonPosition,dissimilaritySampleIndex] = iRepresentativeAbundanceMax+(random.random()*dMaxSignalNoise)
                setBlockFeatures = setBlockFeatures | set([iBlockTaxonPosition])

            #Create min representative features in sample (noisey samples)
            for iRepresentativeNoiseIndex in random.sample(set(range(iTaxaCount))-setBlockFeatures,iRepresentativeDissimilarityLowTaxa):
                npDataMatrix[iRepresentativeNoiseIndex,dissimilaritySampleIndex] = iRepresentativeAbundanceMin+(random.random()*dMaxGeneralNoise)

            #Update the block range for the next sample
            taxonBlockStart = taxonBlockStart+iRepresentiveDissimilarityTaxa

        #Update sample bounds to make the extreme samples
        iSampleStart = iSampleStop
        iSampleStop = iSampleStop + iExtremeDissimilarityCount

        #Extreme dissimilarity
        extremeTaxaBlockIncrement = iRepresentiveDissimilarityTaxa-iExtremeDissimilarityTaxa
        taxonBlockStart = iExtremeDissimilarityStart
        taxonBlockStop = taxonBlockStart + extremeTaxaBlockIncrement
        for extremeSampleIndex in xrange(iSampleStart,iSampleStop):
            #Append name
            strSampleName = "_".join(["Sample",str(extremeSampleIndex),"E"])
            lsSampleNames.append(strSampleName)

            #Set block abundance
            for extremeBlockTaxonPosition in xrange(taxonBlockStart,taxonBlockStop):
                npDataMatrix[extremeBlockTaxonPosition,extremeSampleIndex] = iExtremeAbundanceMax+(random.random()*dMaxSignalNoise)
            #Define index populations and set abundance per sample (for noise)
            population = set(range(iTaxaCount))
            #Remove taxa drivers
            population = population.difference(liTaxaDriverPositions)
            #Remove current block
            population = population.difference(set(range(taxonBlockStart,taxonBlockStop)))
            #Generate noise population for extreme dissimilarity
            extremeNoisePopulation = random.sample(population,iExtremeDissimilarityLowTaxa)
            #Set extreme noise
            for extremeNoiseTaxonPosition in extremeNoisePopulation:
                npDataMatrix[extremeNoiseTaxonPosition,extremeSampleIndex] = iExtremeAbundanceMin+(random.random()*dMaxGeneralNoise)
            #Increment start and stop
            taxonBlockStart = taxonBlockStop+extremeTaxaBlockIncrement
            taxonBlockStop = taxonBlockStart + iExtremeDissimilarityTaxa

        #Update sample bounds to make the targeted feature samples
        iSampleStart = iSampleStop
        iSampleStop = iSampleStop + iTaxaDrivenCount

        #Create the remaining taxa driven samples
        #Update sample bounds
        for drivenSampleIndex in xrange(iSampleStart,iSampleStop):
            #Append Label and name
            strSampleName = "_".join(["Sample",str(drivenSampleIndex),"T"])
            lsSampleNames.append(strSampleName)

            #Define index populations and set abundance per sample (for noise)
            population = set(range(iTaxaCount))
            #Remove taxa drivers
            population = population.difference(liTaxaDriverPositions)
            #Set sample abundance
            for drivingTaxonPosition in liTaxaDriverPositions:
                npDataMatrix[drivingTaxonPosition,drivenSampleIndex] = iTargetedAbundanceMax+(random.random()*dMaxSignalNoise)
            #Generate low noise
            noisePositions = random.sample(population,iExtremeDissimilarityLowTaxa)
            for noiseTaxonPosition in noisePositions:
                npDataMatrix[noiseTaxonPosition,drivenSampleIndex] = iTargetedAbundanceMin+(random.random()*dMaxGeneralNoise)

        #Apply simple noise structure if indicated to do so
#        if dSimpleNoise > 0.0:
#            for iSample in xrange(iSampleCount):
#                for iFeature in xrange(iTaxaCount):
#                    npDataMatrix[iFeature,iSample] = npDataMatrix[iFeature,iSample] + random.random()*dSimpleNoise

        #Add in noise by shuffling a percentage of the reads to other features
        if dSimpleNoise > 0:
            #For each sample get a percentage of the abundance and add to a random group

            for iSample in xrange(iSampleCount):

                #Holds the noise adjustment and then after all the features in the sample are looked at
                #The noise adjustment is added to the sample.
                #If I added on the noise to features as I went through the sample if the feature with added
                #noise are further down the list the percentage noise of that feature will be changed with the additional
                #Added noise and so skews the noise model
                #Noise is calculated seperately and then added at once.
                npaNoise = np.array([0]*iTaxaCount)

                #Calculate noise and shuffle to other features

                for iFeature in xrange(iTaxaCount):
                    dCurrentAbundance = npDataMatrix[iFeature,iSample]
                    iShuffle = int(dCurrentAbundance*dSimpleNoise)
 
                    #Some of the features will have a measurement of 0 and so this can be skipped
                    if(iShuffle>0):
                        npaNoise[iFeature] = npaNoise[iFeature]-iShuffle
                        
                        for iShuffleInstance in xrange(iShuffle):
                            iLocation = random.randint(0,iTaxaCount-1)
                            npaNoise[iLocation] = npaNoise[iLocation] + 1

                #Add noise to sample
                for iindexNoise, iNoise in enumerate(npaNoise):
                    npDataMatrix[iindexNoise,iSample] = npDataMatrix[iindexNoise,iSample] + iNoise

        #Write to file
        #Delete current file before writing
        if(os.path.isfile(strOutputFile)):
            os.remove(strOutputFile)

        #Data to output
        outputContents = []
        #Output sample id header
        outputContents.append(ConstantsMicropita.TAB.join(lsSampleNames))
        #Output sample labels
        outputContents.append(ConstantsMicropita.TAB.join(lsLabels))
        #Write rows/taxa data to file
        for row in xrange(0,iTaxaCount):
            taxaData = []
            for sampleAbundance in npDataMatrix[row]:
                taxaData.append(str(sampleAbundance))
            taxaData = ConstantsMicropita.TAB.join(taxaData)
            outputContents.append(ConstantsMicropita.TAB.join([lsTaxaNames[row],taxaData]))
        with open(strOutputFile,'a') as f:
            f.write(ConstantsMicropita.ENDLINE.join(outputContents))
        f.close()

        #Write to file the "Actual file" which defines the classes
        if not strSampleClassification == None:
            if(not strActualData == ""):
                fHndlOutput = open(strSampleClassification,'w')
                fHndlOutput.write(str(strActualData))
                fHndlOutput.close()

        #Return file name
        return strOutputFile

    ##
    #Generate matrix of random data
    #
    @staticmethod
    def generateRandomMatrix(tempFilePath,iNumberRows,iNumberColumns, iMinValue, iMaxValue, charDelimiter = ConstantsMicropita.TAB):
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
            lsChaoLabels.append(str(Metric.getChao1DiversityIndex(tempSampleTaxaAbundancies=dataMatrix[:,iDiversitySampleIndex])))
            if(sum(dataMatrix[:,iDiversitySampleIndex])==0):
                lsISLabels.append("0")
            else:
                lsISLabels.append(str(Metric.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=dataMatrix[:,iDiversitySampleIndex])))

        #Write to file
        #Delete current file before writing
        if(os.path.isfile(tempFilePath)):
            os.remove(tempFilePath)

        #Data to output
        outputContents = []
        #Output sample id header
        outputContents.append(ConstantsMicropita.TAB.join(lsSampleNames))
        #Output Chao labels
        outputContents.append(ConstantsMicropita.TAB.join(lsChaoLabels))
        #Output InvSimpson labels
        outputContents.append(ConstantsMicropita.TAB.join(lsISLabels))
        #Write rows/taxa
        for row in xrange(0,iTaxaCount):
            taxaData = []
            for sampleAbundance in dataMatrix[row]:
                taxaData.append(str(sampleAbundance))
            taxaData = ConstantsMicropita.TAB.join(taxaData)
            outputContents.append(ConstantsMicropita.TAB.join([lsTaxaNames[row],taxaData]))

        #Write to file
        with open(tempFilePath,'w') as f:
            sFileContents = f.write(ConstantsMicropita.ENDLINE.join(outputContents))
            f.close()

        #Return file name
        return tempFilePath

#    @staticmethod
#    def funcGenerateCorrelatedFeaturesDataSet(sOutputFile, dFeatureScale, dSampleScale):
#
#        #Make sure scales are ints
#        dFeatureScale = int(dFeatureScale)
#        dSampleScale = int(dSampleScale)
#
#        lsCorrelation = []
#
#        dSampleCount = 10*dSampleScale
#        dFeatureCount = 9*dFeatureScale
#
#        iFeatureIndex = 1
#
#        #Create correlated, anticorrelated, and flat but scaled features
#        for feature in xrange(dFeatureCount):
#            ldFeature = []
#            dCurrentFeatureScale = math.pow(2.0,feature)
#            for dSample in xrange(dSampleCount):
#                ldFeature.append(dSample*dCurrentFeatureScale)
#            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+[str(feature) for feature in ldFeature[:]])
#            iFeatureIndex = iFeatureIndex + 1
#            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+([str(feature)]*dSampleCount))
#            ldFeature.reverse()
#            iFeatureIndex = iFeatureIndex + 1
#            lsCorrelation.append(["Feature_"+str(iFeatureIndex)]+[str(feature) for feature in ldFeature[:]])
#            iFeatureIndex = iFeatureIndex + 1
#            dCurrentFeatureScale = math.pow(2.0,feature)
#
#        #Update contents to a line
#        sContents = "FeatureID"+ConstantsMicropita.TAB+ConstantsMicropita.TAB.join(["".join(["Sample_",str(iSampleIndex)]) for iSampleIndex in xrange(dSampleCount)])
#        sContents = ConstantsMicropita.ENDLINE.join([sContents]+[ConstantsMicropita.TAB.join(feature) for feature in lsCorrelation])
#
#        #Write line to a file
#        with open(sOutputFile,'w') as f:
#            f.write(sContents)
#        f.close()
#        return

    @staticmethod
    def funcGenerateCorrelatedFeaturesDataSet(sOutputFile, dScale=1):

        #Make sure scales are ints
        dScale = int(dScale)

        lsCorrelation = []
        dMaxSampleCount = 10*dScale
        dFeatureCount = dMaxSampleCount
        iIndexFeatureName = 0

        #The mean value the data will center around
        cCenter = math.pow(2.0,dFeatureCount)*dMaxSampleCount/2

        #Create correlated, anticorrelated, and flat but scaled features
        #Need to make sure the sum of the features through the samples are = 1 or you
        #Impose unintented structure
        #Here both the samples and the features sum to same number
        for iFeature in xrange(dFeatureCount):
            ldFeature1 = []
            ldFeature2 = []
            for dSample in xrange(dMaxSampleCount):
                dCurrentFeatureScale = math.pow(2.0,iFeature)*dSample
                ldFeature1.append(cCenter+dCurrentFeatureScale)
                ldFeature2.append(cCenter-dCurrentFeatureScale)
            ldFeature1.reverse()
            ldFeature1 = ldFeature1+ldFeature2
            lsCorrelation.append(["Feature_"+str(iIndexFeatureName)]+[str(dFeature) for dFeature in ldFeature1])
            ldFeature1.reverse()
            lsCorrelation.append(["Feature_"+str(iIndexFeatureName+1)]+[str(dFeature) for dFeature in ldFeature1])
            iIndexFeatureName =iIndexFeatureName + 2

        #Update contents to a line
        sContents = "FeatureID"+ConstantsMicropita.TAB+ConstantsMicropita.TAB.join(["".join(["Sample_",str(iSampleIndex)]) for iSampleIndex in xrange(dMaxSampleCount*2)])
        sContents = ConstantsMicropita.ENDLINE.join([sContents]+[ConstantsMicropita.TAB.join(dFeature) for dFeature in lsCorrelation])

        #Write line to a file
        with open(sOutputFile,'w') as f:
            f.write(sContents)
        f.close()
        return

#for iGeneralRandom in [5.0]:
#  for iSignalRandom in [5.0,10.0,15.0,20.0,25.0]:
#    for i in xrange(1,11):
#      Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced96-GenNoise-"+str(iGeneralRandom)+"-SignalNoise-"+str(iSignalRandom)+"v"+str(i)+".pcl", strSampleClassification="Unbalanced96-GenNoise-"+str(iGeneralRandom)+"-SignalNoise-"+str(iSignalRandom)+"-Actual.txt", iScalingFactorForSampleAmount = 2, dMaxGeneralNoise = iGeneralRandom, dMaxSignalNoise=iSignalRandom)
#      Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced48-GenNoise-"+str(iGeneralRandom)+"-SignalNoise-"+str(iSignalRandom)+"v"+str(i)+".pcl", strSampleClassification="Unbalanced48-GenNoise-"+str(iGeneralRandom)+"-SignalNoise-"+str(iSignalRandom)+"-Actual.txt", iScalingFactorForSampleAmount = 1, dMaxGeneralNoise = iGeneralRandom, dMaxSignalNoise=iSignalRandom)

iGeneralRandom = 0
iSignalRandom = 0
for dSimpleNoise in [.05,.1]:#[.05,.10,.15,.20,.25]:
  for i in xrange(1,11):
    Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced96-SimpleNoise-"+str(int(dSimpleNoise*100))+"v"+str(i)+".pcl", strSampleClassification="Unbalanced96-SimpleNoise-"+str(int(dSimpleNoise*100))+"-Actual.txt", iScalingFactorForSampleAmount = 2, dMaxGeneralNoise = iGeneralRandom, dMaxSignalNoise=iSignalRandom, dSimpleNoise=dSimpleNoise)
    Utility_Data.generateAbundanceTable(strOutputFile="Unbalanced48-SimpleNoise-"+str(int(dSimpleNoise*100))+"v"+str(i)+".pcl", strSampleClassification="Unbalanced48-SimpleNoise-"+str(int(dSimpleNoise*100))+"-Actual.txt", iScalingFactorForSampleAmount = 1, dMaxGeneralNoise = iGeneralRandom, dMaxSignalNoise=iSignalRandom, dSimpleNoise=dSimpleNoise)


#Utility_Data.funcGenerateCorrelatedFeaturesDataSet("TestCor.pcl", 1)
