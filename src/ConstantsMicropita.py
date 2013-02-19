"""
Author: Timothy Tickle
Description: Constants.
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

from breadcrumbs.Metric import Metric

##
#Used to test the FileIO class
class ConstantsMicropita():
    """
    Class to hold project constants.
    """

    #Character Constants
    COLON = ":"
    COMMA = ","
    FASTA_ID_LINE_START = ">"
    QUOTE = "\""
    TAB = '\t'
    WHITE_SPACE = " "
    PIPE = "|"
    c_outputFileDelim = '\t'

    c_sEmptyPredictFileValue = 'NA'

    #Used to stop divide by zero errors
    c_smallNumber = 0.00000000001

    #SVM related
    c_COST_RANGE_KEY = "range"
    c_lCostRange = [-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]
    c_SCALED_FILE_EXT = ".scaled"
    c_intScaleLowerBound = 0
    #LIBSVM file extensions
    c_SCALING_PARAMETERS = ".range"
    c_CV_FILE_EXT = ".cv.out"
    c_CV_IMAGE_EXT = ".cv.png"
    c_MODEL_FILE_EXT = ".model"
    c_PREDICT_FILE_EXT = ".predict"
    c_fProbabilitistic = True
    c_SCALED_FOR_PREDICTION_FILE_EXT = ".scaledForpredict"

    #SVM output Dictionary keywords for files
    c_strKeywordInputFile = "INPUT"
    c_strKeywordScaledFile = "SCALED"
    c_strKeywordRangeFile = "RANGE"
    c_strKeywordCVOutFile = "CV_OUT"
    c_strKeywordModelFile = "MODEL"
    c_strKeywordScaledPredFile = "SCALED_FOR_PREDICTION"
    c_strKeywordPredFile = "PREDICTION"
    c_strKeywordCostValue = "C"
    c_strKeywordAccuracy = "ACCURACY"

    #Default values for missing data in the Abundance Table
    c_strEmptyAbundanceData = "0"
    c_strEmptyDataMetadata = "NA"
    lNAs = list(set(["NA","na","Na","nA",c_strEmptyDataMetadata]))

    #Occurence filter [min abundance, min samples occuring in]
    #To turn off make == [0,0]
    c_liOccurenceFilter = [0,0]

    #Break ties in targeted feature with diversity
    c_fBreakRankTiesByDiversity = False

    ####Commandline arguments
    #a Custom diversity metrics found in cogent
    c_strCustomAlphaDiversityHelp = "A key word for any PyCogent supplied alpha diveristy metric (Richness, evenness, or diversity). Please supply an unnormalized (counts) abundance table for these metrics. Metrics include "+" ".join(Metric.setAlphaDiversities)+"."

    #b Custom diversity metrics found in cogent
    c_strCustomBetaDiversityHelp = "A key word for any PyCogent supplied beta diversity metric. Metrics include "+" ".join(list(Metric.setBetaDiversities)+[Metric.c_strUnifracUnweighted,Metric.c_strUnifracWeighted])+"."

    #c,checked Checked abundance file
    c_strCheckedAbundanceFileArgument = "--checked"
    c_strCheckedAbundanceFileHelp = "Before analysis abundance files are checked and a new file results which analysis is perfromed on. The name fo the checked file can be specified of the default will will be used (appending a -Checked to the end of the file name)."

    #d,id Name of the sample id row
    c_strIDNameArgument = "--id"
    c_strIDNameHelp = "The row in the abundance file that is the sample name/id row. Should be the sample name/Id in first column of the row."

    #e,label Supervised Label
    c_strSupervisedLabelArgument = "--label"
    c_strSupervisedLabelHelp = "The name of the metadata on which to perform supervised methods"

    #f, invertDiversity
    c_strInvertDiversityHelp = "".join(["When using this flag, the diversity will be inverted (multiplicative inverse) before ranking in the highest diversity method. ",
			       "Recommended to use with dominance, menhinick, reciprocal_simpson, berger_parker_d, mcintosh_e, simpson_e, strong and any metric where 0 indicates most diverse."])

    #g,logging Path of the logging file
    c_strLoggingFileArgument = "--logfile"
    c_strLoggingFileHelp = "File path to save the logging file."

    #h help

    #i,tree
    c_strCustomEnvironmentFileHelp = "File describing the smaple environments; for use with Unifrac distance metrics."

    #j,delim File delimiter
    c_strFileDelimiterArgument = "--delim"
    c_strFileDelimiterHelp = "The delimiter for the abundance table (default = TAB)"

    #k,featdelim Feature delimiter
    c_strFeatureNameDelimiterArgument = "--featdelim"
    c_strFeatureNameDelimiterHelp = "The delimiter for a feature name if it contains a consensus sequence."

    #l,lastmeta The name of the last metadata
    c_strLastMetadataNameArgument = "--lastmeta"
    c_strLastMetadataNameHelp = "The row in the abundance file that is the sample name/id row. Should be the metadata name/Id in first column of the metadta row."

    #m,method
    c_strSelectionTechniquesHelp = "Select techniques listed one after another."

    #n,num The Number of unsupervised sample selection
    c_strCountArgument = "-n"
    c_strCountHelp = "The number of samples to select with unsupervised methodology. (An integer greater than 0.)."

    #o,tree
    c_strCustomPhylogeneticTreeHelp = "Tree for phylogenetic when selecting custom beta-diversities in the representative sampling criteria."

    #p,suppredfile File path fo the predict file for the supervised methods
    c_strSupervisedPredictedFile = "--suppredfile"
    c_strSupervisedPredictedFileHelp = "The file path for the predict file."

    #q,alphameta
    c_strCustomAlphaDiversityMetadataHelp = "Metric in the pcl file which has custom alpha diversity measurements to use with the highest diversity sampling criteria. Should be a number between 0.0 and 1.0 with 1.0 meaning most diverse."

    #r,targetmethod Taxa selection method
    c_strTargetedFeatureMethodArgument = "--feature_method"
    c_strTargetedFeatureMethodHelp = "The ranking method used to select targeted features."

    #s,stratify Unsupervised stratify metadata
    c_strUnsupervisedStratifyMetadataArgument = "--stratify"
    c_strUnsupervisedStratifyMetadataHelp = "The metatdata to stratify unsupervised analysis."

    #t,target Targeted feature file
    c_strTargetedSelectionFileArgument = "--targets"
    c_strTargetedSelectionFileHelp = "A file containing taxa/OTUs/clades to be used in targeted feature sampling criteria."

    #u,supinputfile File path for the input file for the supervised methods
    c_strSupervisedInputFile = "--supinputfile"
    c_strSupervisedInputFileHelp = "The file path for the input file for supervised methods."

    #v,logging String for logging level
    c_strLoggingArgument = "--logging"
    c_strLoggingHelp = "".join(["Logging level which will be logged to a .log file with the",
         " same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL."])
    c_lsLoggingChoices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]

    #x,betamatrix
    c_strCustomBetaDiversityMatrixHelp = "Precalculated beta-diversity matrix to be used in the representative sampling criteria. Should be a number between 0.0 and 1.0 with 1.0 meaning most dissimilar."

    #Order is important, the first is the default
    c_strTargetedRanked = "rank"
    c_strTargetedAbundance = "abundance"
    lsTargetedFeatureMethodValues = [c_strTargetedRanked, c_strTargetedAbundance]

    #Selection methods
    c_strDiversity = "diverse"
    c_strExtreme = "extreme"
    c_strDiscriminant = "discriminant"
    c_strDistinct = "distinct"
    c_strRandom = "random"
    c_strRepresentative = "representative"
    c_strFeature = "features"
    c_custom = "custom"
    c_lsAllUnsupervisedMethods = [c_strRepresentative,c_strDiversity,c_strExtreme,c_strFeature,c_strRandom]
    c_lsAllSupervisedMethods = [c_strDiscriminant,c_strDistinct]
    c_lsAllMethods = c_lsAllUnsupervisedMethods + c_lsAllSupervisedMethods

    #Technique Names
    c_strDiversity2 = c_strDiversity+"_C"

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_strGenericOutputDataFileHelp = "The generated output data file."
