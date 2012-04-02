#######################################################
#
#	Title:		Constants_Arguments
#	Author:		Timothy Tickle
#	Date:		03/29/2012
#	Purpose:	Class to hold Constants associated with args parsing and arguments
#
#######################################################


class Constants_Arguments():

    #Index of the first data row
    c_strFirstDataRow = "-d"
    c_strFirstDataRowHelp = "".join(["The row in the abundance file that is the first row to contain abundance data. ",
         "This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata."])

    #Flag to record if the colors are inverted
    c_strInvertArgument = "-i"
    c_strInvertHelp = "Invert the image to a black background (default=False)."

    #String for logging level
    c_strLoggingArgument = "-l"
    c_strLoggingHelp = "".join(["Logging level which will be logged to a .log file with the",
         " same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL."])
    #Index of the sample id row
    c_strSampleNameRowArgument = "-n"
    c_strSampleNameRowHelp = "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering."

    #Flag indicator for normalization
    c_strNormalizeArgument = "-r"
    c_strNormalizeHelp = "Normalize the abundance data before working with it (default=False)."

    #Data key to indicate which insilico data set to generate
    c_dataSetKeyHelp = "Key to indicate which data set to generate. Valid values are Diversity, Unbalanced."

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_genericOutputFigureFileHelp = "The generated output figure."
    c_genericOutputDataFileHelp = "The generated output data file."
    c_strMicropitaSelectFileHelp = "A file containing the samples selected which will be visualized."
    c_strSelectionMethodsHelp = "Select techniques listed one after another."

    ####################################
    #Arg choices
    c_lsLoggingChoices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]


