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
class Constants():
    """
    Class to hold project constants.
    """
    #Character Constants
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

    #AbundanceTable
    #Suffix given to a file that is check with the checkRawDataFile method
    OUTPUT_SUFFIX = "-checked.pcl"

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

    #Occurence filter [min abundance, min samples occuring in]
    #To turn off make == None
    c_liOccurenceFilter = [3,3]

    #Break ties in targeted feature with diversity
    c_fBreakRankTiesByDiversity = False
