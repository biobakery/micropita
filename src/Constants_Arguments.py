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
    c_strAbundanceFilterCutoff = "-afc"
    c_strAbundanceFilterCutoffHelp = "The percentage of samples inwhich a terminal node must be above the AbundanceFilterPercentile to not be filtered."

    #afp Abundance filter percentile
    c_strAbundanceFilterPercentile = "-afp"
    c_strAbundanceFilterPercentileHelp = "The percentile threshold cuttoff for filtering by abundance."

    #a alpha The threshold used for either pvalues or qvalues in Cladogram feature enrichment
    c_strEnrichmentThreshold = "-a"
    c_strEnrichmentThresholdHelp = "The threshold defining significance for p-value or q-values of feature enrichment."

    #c Circlader highlight clade file
    c_strHighlightCladeFile = "-c"
    c_strHighlightCladeHelp = "The file containing the clades or taxa/OTUs to highlight in the circlader."

    #cfl Circlader clade filter level (when filtering by clades this is the level in the clade which is filtered; the ancestor clade)
    c_strCladeFilterLevel = "-cfl"
    c_strCladeFilterLevelHelp = "The level in the clade lineage to filter (based on your input Taxa/OTU IDs). 1 is the highest clade level (for instance Kingdom)."

    #cfm Circlader clade measure level (when filtering by clades this is the level in the clade which is measured; the more terminal clade)
    c_strCladeMeasureLevel = "-cfm"
    c_strCladeMeasureLevelHelp = "The level in the clade lineage to measure for filtering (based on your input Taxa/OTU IDs)."

    #cfn Circlader clade clade filtering minimum number (when filtering by clades minimum number of measring level clades a filter level clade may have without being removed.)
    c_strCladeFilteringMinLevel = "-cfn"
    c_strCladeFilteringMinLevelHelp = "The minimum clades at CladeFilterLevelMeasure the CladeFilterLevel clade can have without being removed."

    #ct Optional Circlader ticks for the internal dendrogram
    c_strCircladerTicks = "-ct"
    c_strCircladerTicksHelp = "When indicated by the Circlader style file, are the ticks used for the internal dendrogram. First to last = internal to external."

    #d Index of the first data row
    c_strFirstDataRow = "-d"
    c_strFirstDataRowHelp = "".join(["The row in the abundance file that is the first row to contain abundance data. ",
         "This row and after are assumed to be abundance data. The area between the iSampleNameRow and this are assumed to be metadata."])

    #e Enrichment method for circlader
    c_strEnrichmentMethod = "-e"
    c_strEnrichmentMethodHelp = """The type of result the TAXA/OTU circles are commenting on. 
                        PVALUE shows differentially enriched elements between selected or not selected populations given a Rank sum nonparameteric t-test and a threshold of 0.05,
                        FDR shows the PVALUE option but with a Benjamini and Hochberg FDR correction and a threshold of 0.1"""

    #i Flag to record if the colors are inverted
    c_strInvertArgument = "-i"
    c_strInvertHelp = "Invert the image to a black background (default=False)."

    #l String for logging level
    c_strLoggingArgument = "-l"
    c_strLoggingHelp = "".join(["Logging level which will be logged to a .log file with the",
         " same name as the strOutFile (but with a .log extension). Valid values are DEBUG, INFO, WARNING, ERROR, or CRITICAL."])

    #n Index of the sample id row
    c_strSampleNameRowArgument = "-n"
    c_strSampleNameRowHelp = "The row in the abundance file that is the sample name/id row (default 0). 0 Based numbering."

    #o Circlader ring order
    c_strRingOrder = "-o"
    c_StrRingOrderHelp = "The order to use for the cladogram rings from inner ring to outer ring." 

    #p Predict file path
    c_strPredictFilePath = "-p"
    c_strPredictFilePathHelp = "Predict file path used to plot prediction if supervised methods occured."

    #pltSel PLot the selected methods (used for teh collection curve figure because I need two nargs=* and
    #so need flags for this argument. Is the same as c_strSelectionTechniques
    c_strPlotSelectedArgument = "-pltSel"
    c_strPlotSelectedHelp = "Select techniques listed one after another."

    #r Flag indicator for normalization
    c_strNormalizeArgument = "-r"
    c_strNormalizeHelp = "Normalize the abundance data before working with it (default=False)."

    #rt Flag indicator for the root of the circlader
    c_strRoot = "-rt"
    c_strRootHelp = "The Clade at which to root the cladogram. Default None indicates no rooting outside of the structure of the input file."

    #s Supervised Label
    c_strSupervisedLabel = "-s"
    c_strSupervisedLabelCountHelp = "The name of the phenotype data row on which to perform supervised methods"    

    #sc Supervised label count
    c_strSupervisedLabelCount = "-sc"
    c_strSupervisedLabelCountHelp = "The count of labeled data to select per label (default =1)"

    #t Taxa file path
    c_strTaxaFilePath = "-t"
    c_strTaxaFileHelp= "The file containing the clades or taxa/OTUs to highlight in the circlader."

    #u Unsupervised stratify metadata
    c_strUnsupervisedStratifyMetadata = "-u"
    c_strUnsupervisedStratifyMetadataHelp = "The metatdata to stratify unsupervised analysis."

    #Data key to indicate which insilico data set to generate
    c_dataSetKeyHelp = "Key to indicate which data set to generate. Valid values are Diversity, Unbalanced."

    ####################################
    #Arguments without commandline flags
    c_strAbundanceFileHelp = "An abundance table."
    c_strCircladerCircleFile = "The name of the input file specifying circle data that is generated and then used by the cladogram program."
    c_strCircladerColorFile = "The name of the input file specifying color that is generated and then used by the cladogram program."
    c_strCircladerHighlightFile = "The name of the input file specifying selection highlighting that is generated and then used by the cladogram program."
    c_strCircladerOutputFigure = "The output cladogram figure."
    c_strCircladerSizeFile = "The name of the input file specifying node size that is generated and then used by the cladogram program."
    c_strCircladerStyleFile = "An input file used to specify cladogram syle features."
    c_strCircladerTaxaFile = "The name of the input file specifying taxa that is generated and then used by the cladogram program."
    c_strCircladerTickFile = "The name of the input file specifying levels of the cladogram that is generated and then used by the cladogram program."
    c_strCountHelp = "The number of samples to select (An integer greater than 0.)."
    c_strHCLColorFile = "An output file that is used by HClust. This is the color file."
    c_strHCLDataFile = "An output file that is used by HClust. This is the data file."
    c_strHCLLabelFile = "An output file that is used by HClust. This is the label file."
    c_strHCLLocation = "The location and to HClust (for example ./external/hclust/hclust.py)."
    c_strOptionalOutputDataFile = "An optional output file"
    c_strSelectionTechniques = "Select techniques listed one after another."
    c_strTaxaSelectionFile = "A file containing taxa to be used in taxa-directed selection."
    c_genericOutputFigureFileHelp = "The generated output figure."
    c_genericOutputDataFileHelp = "The generated output data file."
    c_genericTMPDirLocationHelp = "Directory to place temporary and intermediate files."
    c_strMicropitaActualFileHelp = "A file containing indications of how in silico samples were actually intended on being sampled by sampling methods."
    c_strMicropitaSelectFileHelp = "A file containing the samples selected which will be visualized."
    c_strSelectionMethodsHelp = "Select techniques listed one after another."

    ####################################
    #Arg choices
    c_lsLoggingChoices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
    c_strEnrichmentChoices = ["PVALUE","FDR"]

