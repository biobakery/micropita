#######################################################
#
#	Title:		Constants
#	Author:		Timothy Tickle
#	Date:		09/16/2011
#	Purpose:	Class to hold Constants
#
#######################################################

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
class Constants():
    """
    Class to hold project constants.
    """

    #Src code locations
    #Source Code for Cogent
#    COGENT_SRC = "./external/PyCogent-1.5.1/"
    #Source Code for QIIME
#    QIIME_SRC = "./external/Qiime-1.3.0/"

    #File Locations
    ##Directory that holds all testing oriented documents but not code
    TEST_DATA_DIRECTORY = "testData/"
    ##Directory that holds ouput files generated during testing
#    TEST_DATA_OUTPUT_DIRECTORY = TEST_DATA_DIRECTORY+"output/"
    ##Directory that is used to create temporary files during testing.
    ##The contents of this directory are subject to deletion at any time.
    TEST_DATA_TEMP_DIRECTORY = TEST_DATA_DIRECTORY+"Temp/"
    ##Abridged versions of documents that represent large data files
    TEST_MICROPITA_DOCUMENTS = TEST_DATA_DIRECTORY+"microPITA/"
    ##Abridged versions of documents that represent large data files
    TEST_ABRIDGED_DOCUMENTS = TEST_DATA_DIRECTORY+"AbridgedDocuments/"
    ##Files that represent the correct output of methods
    ##These are part of the regression suite and should not be moved or changed
    TEST_ANSWER_DOCUMENTS = TEST_DATA_DIRECTORY+"CorrectTestingResults/"
    ##Directory that holds all input data
    INPUT_DATA_DIRECTORY = "data/"
    ##Directory that holds all configure oriented documents but not code
    CONFIG_DIRECTORY = "config/"
    ##Config file for standard mode
    LOG_CONFIGURE_FILE = CONFIG_DIRECTORY + "ConfigureFile.xml"
    ##Config file for testing mode
    LOG_TEST_CONFIGURE_FILE = CONFIG_DIRECTORY + "TestConfigureFile.xml"

    #File Constants
    COLON = ":"
    COMMA = ","
    ENDLINE = "\n"
    FASTA_ID_LINE_START = ">"
    PATH_SEP = "/"
    QUOTE = "\""
    TAB = '\t'
    WHITE_SPACE = " "
    PIPE = "|"

    #Default values for missing data in the Abundance Table
    c_strEmptyAbundanceData = "0"
    c_strEmptyDataMetadata = "NA"

    lNAs = list(set(["NA","na","Na","nA",c_strEmptyDataMetadata]))

    #Testing Variables
    TEMPORARY_TEST_FILE = TEST_DATA_TEMP_DIRECTORY+"TEMPTESTINGFILEEEEEEEANDIMGONE.txt"

    #AbundanceTable
    #Suffix given to a file that is check with the checkRawDataFile method
    OUTPUT_SUFFIX = "-checked.pcl"

    #Actual file details for the confusion matrix
    c_strClassPrefix = "[Class]"
    c_strEvenSelection = "[EVEN]"

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

    #Break microbiota driven selection ties by diversity
    fBreakRankTiesByDiversity = False
