#######################################################
#
#	Title:		Constants
#	Author:		Timothy Tickle
#	Date:		09/16/2011
#	Purpose:	Class to hold Constants
#
#######################################################

#Import libaries
import Constants

##
#Used to test the FileIO class
class Constants():

    #Src code locations
    #Source Code for Cogent
    COGENT_SRC = "./external/PyCogent-1.5.1/"
    #Source Code for QIIME
    QIIME_SRC = "./external/Qiime-1.3.0/"

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

    #Testing Variables
    TEMPORARY_TEST_FILE = TEST_DATA_TEMP_DIRECTORY+"TEMPTESTINGFILEEEEEEEANDIMGONE.txt"

    #AbundanceTable
    #Suffix given to a file that is check with the checkRawDataFile method
    OUTPUT_SUFFIX = "-checked.txt"

    #Actual file details for the confusion matrix
    c_strClassPrefix = "[Class]"
    c_strEvenSelection = "[EVEN]"
