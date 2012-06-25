######################################
# Author: Timothy Tickle
# Description: Calculates diversity metrics
#####################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Update path
import sys
from Constants import Constants
import numpy as np
from ValidateData import ValidateData
#if(not Constants.COGENT_SRC in sys.path):
#    sys.path.append(Constants.COGENT_SRC)
#if(not Constants.QIIME_SRC in sys.path):
#    sys.path.append(Constants.QIIME_SRC)

#External libraries
from cogent.maths.stats.alpha_diversity import chao1_uncorrected, chao1_bias_corrected
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode, count_envs
from cogent.parse.tree import DndParser
#from qiime.format import format_unifrac_sample_mapping
#from qiime.parse import parse_otu_table
from scipy.spatial.distance import pdist

class Diversity:
    """
    Performs ecological measurements.
    """

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

    #@params ldSampleTaxaAbundancies Vector of organisms in a sample
    #@return 1 float is returned
    @staticmethod
    def funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates the Simpsons diversity index as defined as sum(Pi*Pi).
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	:
        :type
        :return	:
        :type
        """

        #Calculate metric
        return sum((ldSampleTaxaAbundancies)*(ldSampleTaxaAbundancies))

    @staticmethod
    def funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates Inverse Simpsons diversity index 1/sum(Pi*Pi).
        This is multiplicative inverse which reverses the order of the simpsons diversity index.
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	:
        :type
        :return	:
        :type
        """

        simpsons = Diversity.funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies)
        #Return False if the diversity is 0 before inverting it
        if(simpsons == 0):
            return False
        #If simpsons is false return false, else return inverse
        if(not ValidateData.funcIsFalse(simpsons)):
            simpsons = 1.0/simpsons
        return simpsons

    @staticmethod
    def funcGetShannonDiversityIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates the Shannon diversity index.
        Note***: Assumes that the abundance measurements are already normalized by the total population N.
        If not normalized, include N in the parameter tempTotalN and it will be.

        :param	:
        :type
        :return	:
        :type
        """

        #Calculate metric
        ldSampleTaxaAbundancies = ldSampleTaxaAbundancies[np.where(ldSampleTaxaAbundancies != 0)]
        tempIntermediateNumber = sum(ldSampleTaxaAbundancies*(np.log(ldSampleTaxaAbundancies)))
        if(tempIntermediateNumber == 0.0):
            return 0.0
        return -1 * tempIntermediateNumber

    #@params tempSampleTaxaAbundance =
    #@params fCorrectForBias False indicates uncorrected for bias (uncorrected = Chao 1984, corrected = Chao 1987, Eq. 2)
    @staticmethod
    def funcGetChao1DiversityIndex(ldSampleTaxaAbundancies=None, fCorrectForBias=False):
        """
        Calculates the Chao1 diversity index.
        Note***: Not normalized by abundance.

        :param	:
        :type
        :param	:
        :type
        :return	:
        :type
        """

        #Observed = total number of species observed in all samples pooled
        totalObservedSpecies = len(ldSampleTaxaAbundancies)-len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 0])

        #Singles = number of species that occur in exactly 1 sample
        singlesObserved = len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 1.0])

        #Doubles = number of species that occue in exactly 2 samples
        doublesObserved = len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 2.0])

        #If singles or doubles = 0, return observations so that a divided by zero error does not occur
        if((singlesObserved == 0) or (doublesObserved == 0)):
            return totalObservedSpecies

        #Calculate metric
        if(fCorrectForBias == True):
            return chao1_bias_corrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)
        else:
            return chao1_uncorrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)

    @staticmethod
    def funcGetObservedCount(ldSampleAbundances, dThreshold = 0.0):
        """
        Count how many bugs / features have a value of greater than 0 or the threshold given.
        Expects a vector of abundances.
        Do not normalize data if using the threshold.

        :param	:
        :type
        :param	:
        :type
        :return	:
        :type
        """

        return sum([1 for observation in ldSampleAbundances if observation > dThreshold])

    #@params ldSampleTaxaAbundancies an np.array of samples (rows) x measurements (columns) in which diversity is measured between rows
    #@return ndarray A condensed distance matrix
    @staticmethod
    def funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies=None):
        """
        Calculates the BrayCurtis Beta diversity index.
        d(u,v)=sum(abs(row1-row2))/sum(row1+row2).
        This is scale invariant.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	:
        :type
        :return	:
        :type
        """

        #Calculate metric
        try:
            return pdist(X=ldSampleTaxaAbundancies, metric='braycurtis')
        except ValueError as error:
            print "".join(["Diversity.getBrayCurtisDissimilarity. Error=",str(error)])
            return False

    #@params ldSampleTaxaAbundancies an np.array of samples (rows) x measurements (columns) in which diversity is measured between rows
    #@return ndarray A condensed distance matrix
    @staticmethod
    def funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies=None):
        """
        Calculates 1 - the BrayCurtis Beta diversity index.
        d(u,v)=1-(sum(abs(row1-row2))/sum(row1+row2)).
        This is scale invariant.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	:
        :type
        :return	:
        :type
        """

        bcValue = Diversity.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies = ldSampleTaxaAbundancies)
        if(not ValidateData.funcIsFalse(bcValue)):
            return 1.0-bcValue
        return False

    #Beta diversity
    #Testing: Happy path tested the unweighted option
    #Testing: Need to test the weigthed option
    #Calculates the Unifrac Beta diversity index
    #Note: It seems unifrac takes abundancies not relative abundancies
    #@params tempTaxonomyTree String A rooted (outgroup) Newick format phylogenetics tree
    #@params ldSampleTaxaAbundancies String filename of the Qiime output formatted abundancy matrix
    ##TODO Abundancy matrix Rows = Samples, Columns = Sequences (Taxa,OTU) in a Qiime format
    #@return ndarray A condensed distance matrix
#    @staticmethod
#    def getUnifracDistance(ldSampleTaxaAbundancies=None, tempTaxonomyTree=None, tempWeighted=True):
        #Validate data
        #Translate abundances into dict for unifrac
        #The following if clause code is from the qiime script convert_otu_table_to_unifrac_sample_mapping.py
        #Used it in this manner to avoid commandline calls and to tie directly into Qiime
#        if(ValidateData.isValidFileName(ldSampleTaxaAbundancies)):
#            otuReader = open(ldSampleTaxaAbundancies, 'U')
#            sample_ids, otu_ids, otu_table_array, lineages = parse_otu_table(otuReader, float)
#            envs = format_unifrac_sample_mapping(sample_ids, otu_ids, otu_table_array)
#            otuReader.close()
#
#            #Convert to a dictionary for unifrac
#            envs_Dict = dict()
#            for mapping in envs:
#                elements = mapping.split(Constants.TAB)
#                if(len(elements) > 1):
#                    if(not elements[0] in envs_Dict):
#                        envs_Dict[elements[0]] = dict([[elements[1],int(float(elements[2]))]])
#                    else:
#                        envs_Dict[elements[0]][elements[1]]=int(float(elements[2]))

            #prefunction for the tree
#            tr = DndParser(tempTaxonomyTree, UniFracTreeNode)
            #Calculate metric
            #Results of unifrac
            #The distance matrix is res['distance_matrix']
            #The PCoA data is res['pcoa']
#            return fast_unifrac(tr,envs_Dict,weighted=tempWeighted)
#        else:
#            print "".join(["Diversity.getUnifracDistance. Invalid ldSampleTaxaAbundancies filename. Received=",str(ldSampleTaxaAbundancies)])
#            return False

    #TODO Need to figure out how to combine the non normalized and normalized metric values going in and going out of metric creation
    #@params ldAbundancies List of values to compute diversity
    #@params strMetric Alpha metric to use to define diversity
    #@return float
    @staticmethod
    def funcGetAlphaMetric(ldAbundancies=None, strMetric=None):
        """
        Get alpha abundance of the metric for the vector.

        :param	:
        :type
        :param	:
        :type
        :return	:
        :type
        """

        if(not ValidateData.funcIsValidString(strMetric)):
            return False
        elif(strMetric == Diversity.c_SHANNON_A_DIVERSITY):
            return Diversity.funcGetShannonDiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Diversity.c_SIMPSON_A_DIVERSITY):
            return Diversity.funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Diversity.c_INV_SIMPSON_A_DIVERSITY):
            return Diversity.funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Diversity.c_OBSERVED_COUNT):
            return Diversity.funcGetObservedCount(npaSampleAbundances=ldAbundancies)
        #Needs NOT Normalized Abundance
        elif(strMetric == Diversity.c_CHAO1_A_DIVERSITY):
            return Diversity.funcGetChao1DiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        else:
            return False

    #Testing: Happy path Testing (3)
    #@params npaSampleAbundance Observations (Taxa (row) x sample (column))
    #@params lsSampleNames List of sample names of samples to measure (do not include the taxa id column name or other column names which should not be read)
    #@params lsDiversityMetricAlpha List of diversity metrics to use in measuring
    #@return A lists of lists. Each internal list is a list of (floats) indicating a specific metric measurement method measuring multiple samples
    #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
    @staticmethod
    def funcBuildAlphaMetricsMatrix(npaSampleAbundance = None, lsSampleNames = None, lsDiversityMetricAlpha = None):
        """
        Build a matrix of alpha diversity metrics for each sample
        Row = metric, column = sample

        :param	:
        :type
        :param	:
        :type
        :param	:
        :type
        :return	:
        :type
        """

        if not ValidateData.funcIsValidList(lsDiversityMetricAlpha):
            lsDiversityMetricAlpha = [lsDiversityMetricAlpha]

        #Create return
        returnMetricsMatrix = []
        [returnMetricsMatrix.append(list()) for index in lsDiversityMetricAlpha]

        #Get amount of metrics
        metricsCount = len(lsDiversityMetricAlpha)

        #For each sample get all metrics
        #Place in list of lists
        #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
        for sample in lsSampleNames:
            sampleAbundance = npaSampleAbundance[sample]
            for metricIndex in xrange(0,metricsCount):
                returnMetricsMatrix[metricIndex].append(Diversity.funcGetAlphaMetric(ldAbundancies = sampleAbundance, strMetric = lsDiversityMetricAlpha[metricIndex]))
        return returnMetricsMatrix
