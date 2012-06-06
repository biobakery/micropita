#######################################################
#
#	Title:		Constants_Arguments
#	Author:		Timothy Tickle
#	Date:		03/29/2012
#	Purpose:	Class to hold Constants associated with args parsing and arguments
#
#######################################################

class Constants_Arguments():
    #afc Abundance filter cut off
    c_strAbundanceFilterCutoffArgument = "-afc"
    c_strAbundanceFilterCutoffHelp = "The percentage of samples inwhich a terminal node must be above the AbundanceFilterPercentile to not be filtered."

    #afp Abundance filter percentile
    c_strAbundanceFilterPercentileArgument = "-afp"
    c_strAbundanceFilterPercentileHelp = "The percentile threshold cuttoff for filtering by abundance."

    #a alpha The threshold used for either pvalues or qvalues in Cladogram feature enrichment
    c_strEnrichmentThresholdArgument = "-a"
    c_strEnrichmentThresholdHelp = "The threshold defining significance for p-value or q-values of feature enrichment."

    #c Circlader highlight clade file
    c_strHighlightCladeFileArgument = "-c"
    c_strHighlightCladeHelp = "The file containing the clades or taxa/OTUs to highlight in the circlader."

    #cfl Circlader clade filter level (when filtering by clades this is the level in the clade which is filtered; the ancestor clade)
    c_strCladeFilterLevelArgument = "-cfl"
    c_strCladeFilterLevelHelp = "The level in the clade lineage to filter (based on your input Taxa/OTU IDs). 1 is the highest clade level (for instance Kingdom)."

    #cfm Circlader clade measure level (when filtering by clades this is the level in the clade which is measured; the more terminal clade)
    c_strCladeMeasureLevelArgument = "-cfm"
    c_strCladeMeasureLevelHelp = "The level in the clade lineage to measure for filtering (based on your input Taxa/OTU IDs)."

    #cfn Circlader clade clade filtering minimum number (when filtering by clades minimum number of measring level clades a filter level clade may have without being removed.)
    c_strCladeFilteringMinLevelArgument = "-cfn"
    c_strCladeFilteringMinLevelHelp = "The minimum clades at CladeFilterLevelMeasure the CladeFilterLevel clade can have without being removed."

    #ct Optional Circlader ticks for the internal dendrogram
    c_strCircladerTicksArgument = "-ct"
    c_strCircladerTicksHelp = "When indicated by the Circlader style file, are the ticks used for the internal dendrogram. First to last = internal to external."

    #d Index of the first data row
    c_strFirstDataRowArgument = "-d"
    c_strFirstDataRowHelp = "".join(["The row in the abundance file that is the first row to contain abundance data. ",
         "This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata."])

    #delim File delimiter
    c_strFileDelimiterArgument = "-delim"
    c_strFileDelimiterHelp = "The delimiter for the abundance table."

    #e Enrichment method for circlader
    c_strEnrichmentMethodArgument = "-e"
    c_strEnrichmentMethodHelp = """The type of result the TAXA/OTU circles are commenting on. 
                        PVALUE shows differentially enriched elements between selected or not selected populations given a Rank sum nonparameteric t-test and a threshold of 0.05,
                        FDR shows the PVALUE option but with a Benjamini and Hochberg FDR correction and a threshold of 0.1"""
    c_strEnrichmentChoices = ["PVALUE","FDR"]

    #featdelim Feature delimiter
    c_strFeatureNameDelimiterArgument = "-featdelim"
    c_strFeatureNameDelimiterHelp = "The delimiter for a feaure name if it contains a consensus sequence."

    #i Flag to record if the colors are inverted
    c_strInvertArgument = "-i"
    c_strInvertHelp = "Invert the image to a black background (default=False)."

    #id Name of the sample id row
    c_strIDNameArgument = "-id"
    c_strIDNameHelp = "The row in the abundance file that is the sample name/id row. Should be the sample name/Id in first column of the row."

    #isnorm flag indicating the original file is normalized as read in
    c_strIsNormalizedArgument = "-isnorm"
    c_strIsNormalizedHelp = "Indicates if the file is normalized when read in (True indicates normalized)."

    #issum flag indicating the original file is summed as read in
    c_strIsSummedArgument = "-issum"
    c_strIsSummedHelp = "Indicates if the file is summed when read in (True indicates summed)."

    #label Supervised Label
    c_strSupervisedLabelArgument = "-label"
    c_strSupervisedLabelHelp = "The name of the phenotype data row on which to perform supervised methods"

    #lastmeta The name of the last metadata
    c_strLastMetadataNameArgument = "-lastmeta"
    c_strLastMetadataNameHelp = "The row in the abundance file that is the sample name/id row. Should be the metadata name/Id in first column of the metadta row."

    #Level Terminal clade level
    c_strTerminalLevelArgument = "-level"
    c_strTerminalLevelHelp = "The terminal clade for analysis defined in the tree made by the consensus lineage of the features."

    #l String for logging level
    c_strLoggingArgument = "-logging"
    c_strLoggingHelp = "".join(["Logging level which will be logged to a .log file with the",
         " same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL."])
    c_lsLoggingChoices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]

    #metric The metric used for evaluation
    c_strMetricArgument = "-metric"
    c_strMetricHelp = "The metric used for measurement."

    #-nosum Sum data
    c_strSumDataArgument = "-nosum"
    c_strSumDataHelp = "The analysis/plotting should be performed on data with clades which are summed. If data is not provided this way, summation will occur before analysis and plotting as needed. This turns off summing."

    #nsup Supervised label count
    c_strSupervisedLabelCountArgument = "-nsup"
    c_strSupervisedLabelCountHelp = "The count of labeled data to select per label. (An integer greater than 0)."

    #nun The Number of unsupervised sample selection
    c_strUnsupervisedCountArgument = "-nun"
    c_strUnsupevisedCountHelp = "The number of samples to select with unsupervised methodology. (An integer greater than 0.)."

    #o Circlader ring order
    c_strRingOrderArgument = "-o"
    c_strRingOrderHelp = "The order to use for the cladogram rings from inner ring to outer ring." 

    #ofe Occurence filtering minimum sequence count
    c_strOccurenceFilterSequenceCountArgument = "-ofe"
    c_strOccurenceFilterSequenceHelp = "The minimum sequence count of a feature to be kept in occurence filtering. (A positive value activates the filter)."

    #ofa Occurence filtering minimum sample count
    c_strOccurenceFilterSampleCountArgument = "-ofa"
    c_strOccurenceFilterSampleHelp = "The minimum sample the c_strOccurenceFilterSequenceCount must occur in for a feature to be kept."

    #p Predict file path
    c_strPredictFilePathArgument = "-p"
    c_strPredictFilePathHelp = "Predict file path used to plot prediction if supervised methods occured."

    #pair Pairing metadata
    c_strPairingMetadataArgument = "-pair"
    c_strPairingMetadataHelp = "Metadata used to pair samples."

    #pltSel PLot the selected methods (used for teh collection curve figure because I need two nargs=* and
    #so need flags for this argument. Is the same as c_strSelectionTechniques
    c_strPlotSelectedArgument = "-pltSel"
    c_strPlotSelectedHelp = "Select techniques listed one after another."

    #r Flag indicator for normalization
    c_strNormalizeArgument = "-r"
    c_strNormalizeHelp = "Normalize the abundance data before working with it (default=False)."

    #rt Flag indicator for the root of the circlader
    c_strRootArgument = "-rt"
    c_strRootHelp = "The Clade at which to root the cladogram. Default None indicates no rooting outside of the structure of the input file."

    #stratify Unsupervised stratify metadata
    c_strUnsupervisedStratifyMetadataArgument = "-stratify"
    c_strUnsupervisedStratifyMetadataHelp = "The metatdata to stratify unsupervised analysis."

    #t Taxa file path
    c_strTaxaFilePathArgument = "-t"
    c_strTaxaFileHelp= "The file containing the clades or taxa/OTUs to highlight in the circlader."

    #targetmethod Taxa selection method
    c_strTargetedFeatureMethodArgument = "-targetmethod"
    c_strTargetedFeatureMethodHelp = "The method used to select targeted features."
    #Order is important, the first is the default
    c_TARGETED_METHOD_RANKED = "Targeted_Rank"
    c_TARGETED_METHOD_ABUNDANCE = "Targeted_Abundance"
    lsTargetedFeatureMethodValues = [c_TARGETED_METHOD_RANKED, c_TARGETED_METHOD_ABUNDANCE]

    #tmp Temporary directory
    c_strTemporaryDirectoryArgument = "-tmp"
    c_genericTMPDirLocationHelp = "Directory to place temporary and intermediate output files."

    #target Targeted feature file
    c_strTargetedSelectionFileArgument = "-target"
    c_strTargetedSelectionFileHelp = "A file containing taxa/otu/clades to be used in targeted feature selection."

    #vin flag indicating the original file is normalized as read in
    c_strValidationIsNormalizedArgument = "-vin"
    c_strValidationIsNormalizedHelp = "Indicates if the validation file is normalized when read in (True indicates normalized)."

    #vis flag indicating the original file is summed as read in
    c_strValidationIsSummedArgument = "-vis"
    c_strValidationIsSummedHelp = "Indicates if the validation file is summed when read in (True indicates summed)."

    #vm The name of the last metadata in the validation data set (when you need both)
    c_strValidationLastMetadataNameArgument = "-vm"
    c_strValidationLastMetadataNameHelp = "The row in the validation abundance file that is the sample name/id row. Should be the metadata name/Id in first column of the metadta row."

    #vn Name of the sample id row in the validation data set (when you need both)
    c_strValidationIDNameArgument = "-vn"
    c_strValidationIDNameHelp = "The row in the validation abundance file that is the sample name/id row. Should be the sample name/Id in first column of the row."


    #Data key to indicate which insilico data set to generate
    c_dataSetKeyHelp = "Key to indicate which data set to generate. Valid values are Diversity, Unbalanced."

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_strCircladerCircleFileHelp = "The name of the input file specifying circle data that is generated and then used by the cladogram program."
    c_strCircladerColorFileHelp = "The name of the input file specifying color that is generated and then used by the cladogram program."
    c_strCircladerHighlightFileHelp = "The name of the input file specifying selection highlighting that is generated and then used by the cladogram program."
    c_strCircladerOutputFigureHelp = "The output cladogram figure."
    c_strCircladerOutputDetailsHelp = "The output file listing the detail of generating the cladogram."
    c_strCircladerSizeFileHelp = "The name of the input file specifying node size that is generated and then used by the cladogram program."
    c_strCircladerStyleFileHelp = "An input file used to specify cladogram syle features."
    c_strCircladerTaxaFileHelp = "The name of the input file specifying taxa that is generated and then used by the cladogram program."
    c_strCircladerTickFileHelp = "The name of the input file specifying levels of the cladogram that is generated and then used by the cladogram program."
    c_strHCLColorFileHelp = "An output file that is used by HClust. This is the color file."
    c_strHCLDataFileHelp = "An output file that is used by HClust. This is the data file."
    c_strHCLLabelFileHelp = "An output file that is used by HClust. This is the label file."
    c_strHCLLocationHelp = "The location and to HClust (for example ./external/hclust/hclust.py)."
    c_strOptionalOutputDataFileHelp = "An optional output file"
    c_strSelectionTechniquesHelp = "Select techniques listed one after another."
    c_genericOutputFigureFileHelp = "The generated output figure."
    c_genericOutputDataFileHelp = "The generated output data file."
    c_strMicropitaActualFileHelp = "A file containing indications of how in silico samples were actually intended on being sampled by sampling methods."
    c_strMicropitaSelectFileHelp = "A file containing the samples selected which will be visualized."
    c_strMicropitaProjectsHelp = "A list of projects whose output will be combined into the meta matrix."
    c_strSelectionMethodsHelp = "Select techniques listed one after another (seperated by whitespace)."
    c_strSelectionMethodsCommaHelp = "Select techniques listed one after another (seperated by commas)."
    c_strOutputFileHelp = "A file to write the output."
    c_strValidationAbundanceFileHelp = "Abundance table wit which to validate the selection"
    c_strSelectionAbundanceFileHelp = "Abundance table where the selection, which is being validated, occured."

