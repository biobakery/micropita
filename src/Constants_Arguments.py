"""
Author: Timothy Tickle
Description: Class to hold commanline argument information.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

class Constants_Arguments():

    #checked Checked abundance file
    c_strCheckedAbundanceFileArgument = "--checked"
    c_strCheckedAbundanceFileHelp = "Before analysis abundance files are checked and a new file results which analysis is perfromed on. The name fo the checked file can be specified of the default will will be used (appending a -Checked to the end of the file name)."

    #delim File delimiter
    c_strFileDelimiterArgument = "--delim"
    c_strFileDelimiterHelp = "The delimiter for the abundance table (default = TAB)"

    #featdelim Feature delimiter
    c_strFeatureNameDelimiterArgument = "--featdelim"
    c_strFeatureNameDelimiterHelp = "The delimiter for a feaure name if it contains a consensus sequence."

    #id Name of the sample id row
    c_strIDNameArgument = "--id"
    c_strIDNameHelp = "The row in the abundance file that is the sample name/id row. Should be the sample name/Id in first column of the row."

    #isnorm flag indicating the original file is normalized as read in
    c_strIsNormalizedArgument = "--isnorm"
    c_strIsNormalizedHelp = "If used this flag indicates the file is already normalized."

    #issum flag indicating the original file is summed as read in
    c_strIsSummedArgument = "--issum"
    c_strIsSummedHelp = "If used this flag indicates the file is already summed."

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

    #nosum Sum data
    c_strSumDataArgument = "--nosum"
    c_strSumDataHelp = "The analysis/plotting should be performed on data with clades which are summed. If data is not provided this way, summation will occur before analysis and plotting as needed. This turns off summing."

    #nsup Supervised label count
    c_strSupervisedLabelCountArgument = "--nsup"
    c_strSupervisedLabelCountHelp = "The number of samples to select per label for supervised methods. (An integer greater than 0)."

    #nun The Number of unsupervised sample selection
    c_strUnsupervisedCountArgument = "--nun"
    c_strUnsupevisedCountHelp = "The number of samples to select with unsupervised methodology. (An integer greater than 0.)."

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
    c_strTargetedFeatureMethodArgument = "--targetmethod"
    c_strTargetedFeatureMethodHelp = "The method used to select targeted features."

    #target Targeted feature file
    c_strTargetedSelectionFileArgument = "--target"
    c_strTargetedSelectionFileHelp = "A file containing taxa/otu/clades to be used in targeted feature selection."

    #Order is important, the first is the default
    c_strTargetedRanked = "Targeted_Rank"
    c_strTargetedAbundance = "Targeted_Abundance"
    lsTargetedFeatureMethodValues = [c_strTargetedRanked, c_strTargetedAbundance]

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_strSelectionTechniquesHelp = "Select techniques listed one after another."
    c_strGenericOutputDataFileHelp = "The generated output data file."
