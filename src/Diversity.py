######################################
# Author: Timothy Tickle
# Description: Calculates diversity metrics
#####################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Update path
import sys
from Constants import Constants
import numpy as np
from ValidateData import ValidateData
if(not Constants.COGENT_SRC in sys.path):
    sys.path.append(Constants.COGENT_SRC)
if(not Constants.QIIME_SRC in sys.path):
    sys.path.append(Constants.QIIME_SRC)

#External libraries
from cogent.maths.stats.alpha_diversity import chao1_uncorrected, chao1_bias_corrected
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode, count_envs
from cogent.parse.tree import DndParser
#from qiime.format import format_unifrac_sample_mapping
from qiime.parse import parse_otu_table
from scipy.spatial.distance import pdist

class Diversity:

    #Diversity metrics Alpha
    c_SHANNON_A_DIVERSITY = "ShannonD"
    c_SIMPSON_A_DIVERSITY = "SimpsonD"
    c_INV_SIMPSON_A_DIVERSITY = "InSimpsonD"
    c_CHAO1_A_DIVERSITY = "Chao1"

    #Diversity metrics Beta
    c_UNIFRAC_B_DIVERSITY = "uUnifrac"
    c_WEIGHTED_UNIFRAC_B_DIVERSITY = "wUnifrac"
    c_BRAY_CURTIS_B_DIVERSITY = "B_Curtis"

    #Addative inverses of beta metrics
    c_INVERSE_BRAY_CURTIS_B_DIVERSITY = "InB_Curtis"
    c_INVERSE_UNIFRAC_B_DIVERSITY = "InuUnifrac"
    c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY = "InwUnifrac"

    #Richness
    c_OBSERVED_COUNT = "Observed_Count"

    #Alpha diversity
    #Testing: Happy Path
    #Calculates the Simpsons diversity index as defined as sum(Pi*Pi)
    #Note***: Assumes that the abundance measurements are already normalized by the total population N
    #@params tempSampleTaxaAbundancies Vector of organisms in a sample
    #@return 1 float is returned
    @staticmethod
    def getSimpsonsDiversityIndex(tempSampleTaxaAbundancies=None):
        #Calculate metric
        return sum((tempSampleTaxaAbundancies)*(tempSampleTaxaAbundancies))

    #Alpha diversity
    #Testing: Happy Path
    #Calculates Inverse Simpsons diversity index 1/sum(Pi*Pi)
    #This is multiplicative inverse which reverses the order of the simpsons diversity index
    #Note***: Assumes that the abundance measurements are already normalized by the total population N
    @staticmethod
    def getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=None):

        simpsons = Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies)
        #Return False if the diversity is 0 before inverting it
        if(simpsons == 0):
            return False
        #If simpsons is false return false, else return inverse
        if(not ValidateData.isFalse(simpsons)):
            simpsons = 1.0/simpsons
        return simpsons

    #Alpha diversity
    #Testing: Happy Path
    #Calculates the Shannon diversity index
    #Note***: Assumes that the abundance measurements are already normalized by the total population N
    #If not normalized, include N in the parameter tempTotalN and it will be
    ## Calculates the Shannon index
    @staticmethod
    def getShannonDiversityIndex(tempSampleTaxaAbundancies=None):

        #Calculate metric
        tempSampleTaxaAbundancies = tempSampleTaxaAbundancies[np.where(tempSampleTaxaAbundancies != 0)]
        tempIntermediateNumber = sum(tempSampleTaxaAbundancies*(np.log(tempSampleTaxaAbundancies)))
        if(tempIntermediateNumber == 0.0):
            return 0.0
        return -1 * tempIntermediateNumber

    #Alpha diversity
    #Testing: Happy Path Tested the no bias option
    #Testing: Need to test the biased option
    #Calculates the Chao1 diversity index
    #Note***: Not normalized by abundance
    #@params tempSampleTaxaAbundance =
    #@params tempCorrectForBias False indicates uncorrected for bias (uncorrected = Chao 1984, corrected = Chao 1987, Eq. 2)
    @staticmethod
    def getChao1DiversityIndex(tempSampleTaxaAbundancies=None, tempCorrectForBias=False):

        #Observed = total number of species observed in all samples pooled
        totalObservedSpecies = len(tempSampleTaxaAbundancies)-len(tempSampleTaxaAbundancies[tempSampleTaxaAbundancies == 0])

        #Singles = number of species that occur in exactly 1 sample
        singlesObserved = len(tempSampleTaxaAbundancies[tempSampleTaxaAbundancies == 1.0])

        #Doubles = number of species that occue in exactly 2 samples
        doublesObserved = len(tempSampleTaxaAbundancies[tempSampleTaxaAbundancies == 2.0])

        #If singles or doubles = 0, return observations so that a divided by zero error does not occur
        if((singlesObserved == 0) or (doublesObserved == 0)):
            return totalObservedSpecies

        #Calculate metric
        if(tempCorrectForBias == True):
            return chao1_bias_corrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)
        else:
            return chao1_uncorrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)

    #Happy Path Tested
    #Count how many bugs / features have a value of greater than 0 or the threshold given
    #Expects a vector of abundances
    #Do not normalize data if using the threshold
    @staticmethod
    def getObservedCount(tempSampleAbundances, dThreshold = 0.0):
        return sum([1 for observation in tempSampleAbundances if observation > dThreshold])

    #Beta diversity
    #Testing: Happy Path
    #Calculates the BrayCurtis Beta diversity index
    #d(u,v)=sum(abs(row1-row2))/sum(row1+row2)
    #This is scale invariant
    #If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order
    #condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)]
    #Note***: Assumes that the abundance measurements are already normalized by the total population N
    #@params tempSampleTaxaAbundancies an np.array of samples (rows) x measurements (columns) in which diversity is measured between rows
    #@return ndarray A condensed distance matrix
    @staticmethod
    def getBrayCurtisDissimilarity(tempSampleTaxaAbundancies=None):

        #Calculate metric
        try:
            return pdist(X=tempSampleTaxaAbundancies, metric='braycurtis')
        except ValueError as error:
            print "".join(["Diversity.getBrayCurtisDissimilarity. Error=",str(error)])
            return False

    #Beta diversity
    #Testing: Happy Path
    #Calculates 1 - the BrayCurtis Beta diversity index
    #d(u,v)=1-(sum(abs(row1-row2))/sum(row1+row2))
    #This is scale invariant
    #If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order
    #condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)]
    #Note***: Assumes that the abundance measurements are already normalized by the total population N
    #@params tempSampleTaxaAbundancies an np.array of samples (rows) x measurements (columns) in which diversity is measured between rows
    #@return ndarray A condensed distance matrix
    @staticmethod
    def getInverseBrayCurtisDissimilarity(tempSampleTaxaAbundancies = None):
        bcValue = Diversity.getBrayCurtisDissimilarity(tempSampleTaxaAbundancies = tempSampleTaxaAbundancies)
        if(not ValidateData.isFalse(bcValue)):
            #TODO Since brays curtis can get larger than 1, need to normalize this with a different value
            #TODO Need all inverses to be inverse in a specific way ... maybe multiplicative inverse is better
            return 1.0-bcValue
        return False

    #Beta diversity
    #Testing: Happy path tested the unweighted option
    #Testing: Need to test the weigthed option
    #Calculates the Unifrac Beta diversity index
    #Note: It seems unifrac takes abundancies not relative abundancies
    #@params tempTaxonomyTree String A rooted (outgroup) Newick format phylogenetics tree
    #@params tempSampleTaxaAbundancies String filename of the Qiime output formatted abundancy matrix
    ##TODO Abundancy matrix Rows = Samples, Columns = Sequences (Taxa,OTU) in a Qiime format
    #@return ndarray A condensed distance matrix
    @staticmethod
    def getUnifracDistance(tempSampleTaxaAbundancies=None, tempTaxonomyTree=None, tempWeighted=True):
        #Validate data
        #Translate abundances into dict for unifrac
        #The following if clause code is from the qiime script convert_otu_table_to_unifrac_sample_mapping.py
        #Used it in this manner to avoid commandline calls and to tie directly into Qiime
        if(ValidateData.isValidFileName(tempSampleTaxaAbundancies)):
            otuReader = open(tempSampleTaxaAbundancies, 'U')
            sample_ids, otu_ids, otu_table_array, lineages = parse_otu_table(otuReader, float)
            envs = format_unifrac_sample_mapping(sample_ids, otu_ids, otu_table_array)
            otuReader.close()

            #Convert to a dictionary for unifrac
            envs_Dict = dict()
            for mapping in envs:
                elements = mapping.split(Constants.TAB)
                if(len(elements) > 1):
                    if(not elements[0] in envs_Dict):
                        envs_Dict[elements[0]] = dict([[elements[1],int(float(elements[2]))]])
                    else:
                        envs_Dict[elements[0]][elements[1]]=int(float(elements[2]))

            #prefunction for the tree
            tr = DndParser(tempTaxonomyTree, UniFracTreeNode)
            #Calculate metric
            #Results of unifrac
            #The distance matrix is res['distance_matrix']
            #The PCoA data is res['pcoa']
            return fast_unifrac(tr,envs_Dict,weighted=tempWeighted)
        else:
            print "".join(["Diversity.getUnifracDistance. Invalid tempSampleTaxaAbundancies filename. Received=",str(tempSampleTaxaAbundancies)])
            return False

    #Testing: Happy Path Tested (4)
    #TODO Need to figure out how to combine the non normalized and normalized metric values going in and going out of metric creation
    #Get alpha abundance of the metric for the vector
    #@params tempAbundancies List of values to compute diversity
    #@params tempMetric Alpha metric to use to define diversity
    #@return float
    @staticmethod
    def getAlphaMetric(tempAbundancies=None, tempMetric=None):
        if(not ValidateData.isValidString(tempMetric)):
            return False
        elif(tempMetric == Diversity.c_SHANNON_A_DIVERSITY):
            return Diversity.getShannonDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == Diversity.c_SIMPSON_A_DIVERSITY):
            return Diversity.getSimpsonsDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == Diversity.c_INV_SIMPSON_A_DIVERSITY):
            return Diversity.getInverseSimpsonsDiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        elif(tempMetric == Diversity.c_OBSERVED_COUNT):
            return Diversity.getObservedCount(tempSampleAbundances=tempAbundancies)
        #Needs NOT Normalized Abundance
        elif(tempMetric == Diversity.c_CHAO1_A_DIVERSITY):
            return Diversity.getChao1DiversityIndex(tempSampleTaxaAbundancies=tempAbundancies)
        else:
            return False

    #Testing: Happy path Testing (3)
    #Build a matrix of alpha diversity metrics for each sample
    #Row = metric, column = sample
    #@params tempSampleAbundance Observations (Taxa (row) x sample (column))
    #@params tempSampleNames List of sample names of samples to measure (do not include the taxa id column name or other column names which should not be read)
    #@params tempDiversityMetricAlpha List of diversity metrics to use in measuring
    #@return A lists of lists. Each internal list is a list of (floats) indicating a specific metric measurement method measuring multiple samples
    #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
    @staticmethod
    def buildAlphaMetricsMatrix(tempSampleAbundance = None, tempSampleNames = None, tempDiversityMetricAlpha = None):

        if not ValidateData.isValidList(tempDiversityMetricAlpha):
            tempDiversityMetricAlpha = [tempDiversityMetricAlpha]

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
                returnMetricsMatrix[metricIndex].append(Diversity.getAlphaMetric(tempAbundancies = sampleAbundance, tempMetric = tempDiversityMetricAlpha[metricIndex]))
        return returnMetricsMatrix
