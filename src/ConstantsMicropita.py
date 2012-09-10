"""
Author: Timothy Tickle
Description: Constants.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

##
#Used to test the FileIO class
class ConstantsMicropita():
    """
    Class to hold project constants.
    """

    #References to other projects
#    c_strBreadcrumbsProject = "../breadcrumbs/src/"

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
    c_liOccurenceFilter = [2,2]

    #Break ties in targeted feature with diversity
    c_fBreakRankTiesByDiversity = False

    ####Commandline arguments
    #checked Checked abundance file
    c_strCheckedAbundanceFileArgument = "--checked"
    c_strCheckedAbundanceFileHelp = "Before analysis abundance files are checked and a new file results which analysis is perfromed on. The name fo the checked file can be specified of the default will will be used (appending a -Checked to the end of the file name)."

    #delim File delimiter
    c_strFileDelimiterArgument = "--delim"
    c_strFileDelimiterHelp = "The delimiter for the abundance table (default = TAB)"

    #featdelim Feature delimiter
    c_strFeatureNameDelimiterArgument = "--featdelim"
    c_strFeatureNameDelimiterHelp = "The delimiter for a feature name if it contains a consensus sequence."

    #id Name of the sample id row
    c_strIDNameArgument = "--id"
    c_strIDNameHelp = "The row in the abundance file that is the sample name/id row. Should be the sample name/Id in first column of the row."

    #label Supervised Label
    c_strSupervisedLabelArgument = "--label"
    c_strSupervisedLabelHelp = "The name of the metadata on which to perform supervised methods"

    #lastmeta The name of the last metadata
    c_strLastMetadataNameArgument = "--lastmeta"
    c_strLastMetadataNameHelp = "The row in the abundance file that is the sample name/id row. Should be the metadata name/Id in first column of the metadta row."

    #logging String for logging level
    c_strLoggingArgument = "--logging"
    c_strLoggingHelp = "".join(["Logging level which will be logged to a .log file with the",
         " same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL."])
    c_lsLoggingChoices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]

    #logging Path of the logging file
    c_strLoggingFileArgument = "--logfile"
    c_strLoggingFileHelp = "File path to save the logging file."

    #n The Number of unsupervised sample selection
    c_strCountArgument = "-n"
    c_strCountHelp = "The number of samples to select with unsupervised methodology. (An integer greater than 0.)."

    #stratify Unsupervised stratify metadata
    c_strUnsupervisedStratifyMetadataArgument = "--stratify"
    c_strUnsupervisedStratifyMetadataHelp = "The metatdata to stratify unsupervised analysis."

    #supinputfile File path for the input file for the supervised methods
    c_strSupervisedInputFile = "--supinputfile"
    c_strSupervisedInputFileHelp = "The file path for the input file for supervised methods."

    #suppredfile File path fo the predict file for the supervised methods
    c_strSupervisedPredictedFile = "--suppredfile"
    c_strSupervisedPredictedFileHelp = "The file path for the predict file."

    #targetmethod Taxa selection method
    c_strTargetedFeatureMethodArgument = "--feature_method"
    c_strTargetedFeatureMethodHelp = "The ranking method used to select targeted features."

    #target Targeted feature file
    c_strTargetedSelectionFileArgument = "--targets"
    c_strTargetedSelectionFileHelp = "A file containing taxa/OTUs/clades to be used in targeted feature selection."

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
    c_lsAllUnsupervisedMethods = [c_strRepresentative,c_strDiversity,c_strExtreme,c_strFeature,c_strRandom]
    c_lsAllSupervisedMethods = [c_strDiscriminant,c_strDistinct]
    c_lsAllMethods = c_lsAllUnsupervisedMethods + c_lsAllSupervisedMethods

    #Technique Names
    c_strDiversity2 = c_strDiversity+"_C"

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_strSelectionTechniquesHelp = "Select techniques listed one after another."
    c_strGenericOutputDataFileHelp = "The generated output data file."
