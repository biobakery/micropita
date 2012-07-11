#######################################################
# Author: Timothy Tickle
# Description: Class to Run analysis for the microPITA paper
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
import blist
from CClade import CClade
from Constants import Constants
import csv
import copy
import numpy as np
import os
import re
import scipy.stats
import string
from ValidateData import ValidateData

c_dTarget	= 1.0
c_fRound	= False
c_iSumAllCladeLevels = -1
c_fOutputLeavesOnly = False

class AbundanceTable:
    """
    Represents an abundance table and contains common function to perform on the object.

    This class is made from an abundance data file. What is expected is a text file delimited by
    a character (which is given to the object). The first column is expected to be the id column
    for each of the rows. Metadata is expected before measurement data. Columns are samples and
    rows are features (bugs). 
    """

    def __init__(self, npaAbundance, dictMetadata, strName, fIsNormalized, fIsSummed, cFileDelimiter = Constants.TAB, cFeatureNameDelimiter="|"):
      """
      Averages feature abundance.

      :param	npaAbundance:	Structured Array of abundance data (Row=Features, Columns=Samples)
      :type	Numpy Structured Array:	Structured Array of abundance data (Row=Features, Columns=Samples)
      :param	dictMetadata:	Structured Array of abundance data (Row=Features, Columns=Samples)
      :type	Dictionary:	Dictionary of metadata {"String ID":["strValue","strValue","strValue","strValue","strValue"]}
      :param	strName:	The name of the metadata that serves as the ID for the columns (For example a sample ID)
      :type	String:	Structured Array of abundance data (Row=Features, Columns=Samples)
      :param	fIsNormalized:	Indicates if the data is already normalized upon reading
      :type	Boolean:	Boolean indicator of normalization (True=Already Normalized)
      :param	fIsSummed:	Indicates if the data is already summed upon reading
      :type	Boolean:	Boolean indicator of already being summed (True=Summed)
      :param	cFileDelimiter:	Character used as the delimiter of the file that is read in to create the abundance table.
                                Will also be used to write the abudance table file to a file to keep file consistency.
      :type	Character:	Character delimiter for reading the data in (default = TAB)
      :param	cFeatureNameDelimiter:	Character used as the delimiter of the feature names (column 1). This is useful if the name are complex, for instance consensus lineages in metagenomics.
      :type	Character:	Character delimiter for feature names (default = |)
      """

      #The abundance data
      self._npaFeatureAbundance = None

      #The metdata
      self._dictTableMetadata = None

      #The name of the object relating to the file it was read from or would have been read from if it exists
      #Keeps tract of changes to the file through the name
      #Will be used to wrtie out the object to a file as needed
      self._strOriginalName = strName

      #The original number of features in the table
      self._iOriginalFeatureCount = -1

      #The original number of samples in the table
      self._iOriginalSampleCount = -1

      #Current normalization state
      self._fIsNormalized = fIsNormalized

      #Indicates if the internal clades are summed
      self._fIsSummed = fIsSummed

      #Indicates if the table has been filtered and how
      self._iCurrentFilterState = ""

      #The feature name delimiter
      self._cFeatureDelimiter = cFeatureNameDelimiter

      #The delimiter from the source file
      self._cDelimiter = cFileDelimiter

      #If contents is not a false then set contents to appropriate objects
      if (not npaAbundance == None) and dictMetadata:
        self._npaFeatureAbundance = npaAbundance
        self._dictTableMetadata = dictMetadata
        self._iOriginalFeatureCount = self._npaFeatureAbundance.shape[0]
        self._iOriginalSampleCount = len(self.funcGetSampleNames())
      else:
        print "Abundance or metadata was None, should be atleast an empty object"

    @staticmethod
    def funcMakeFromFile(strInputFile, fIsNormalized, fIsSummed, cDelimiter = Constants.TAB, sMetadataID = None, sLastMetadata = None, cFeatureNameDelimiter="|"):
        """
        Creates an abundance table from a table file.

        :param	strInputFile:	Path to input file.
        :type	String		String path.
        :param	fIsNormalized:	Indicates if the data is already normalized on read.
        :type	Boolean		True indicates the data IS normalized.
        :param	fIsSummed:	Indicates if the data is already summed.
        :type	Boolean		True indicates the data IS summed at the clade levels.
        :param	cDelimiter:	Delimiter for parsing the input file.
        :type	Character	Character.
        :param	sMetadataID:	String ID that is a metadata row ID (found on the first column) and used as an ID for samples
        :type	String		String ID
        :param	sLastMetadata:	The ID of the metadata that is the last metadata before measurement or feature rows.
        :type	String		String ID
        :param	cFeatureNameDelimiter:	Used to parse Feature (bug) names if they are complex.
                                        For example if they are consensus lineages and contain parent clade information.
        :type	Character	Delimiting letter
        :return	AbundanceTable:	Will return an AbundanceTable object on no error. Returns False on error.
        :type	AbundanceTable or False
        """


        #Read in from text file to create the abundance and metadata structures
        lContents = AbundanceTable._funcTextToStructuredArray(strInputFile=strInputFile, cDelimiter=cDelimiter,
                                                        sMetadataID = sMetadataID, sLastMetadata = sLastMetadata)

        #If contents is not a false then set contents to appropriate objects
        if lContents:
            return AbundanceTable(npaAbundance=lContents[0], dictMetadata=lContents[1], strName=strInputFile, 
                              fIsNormalized=fIsNormalized, fIsSummed=fIsSummed, cFileDelimiter=cDelimiter, cFeatureNameDelimiter=cFeatureNameDelimiter)
        return False

    def __repr__(self):
        """
        Represent or print object
        """
        return "AbundanceTable"

    def __str__(self):
      """
      Create a string representation of the Abundance Table.
      """

      return "".join(["Sample count:", str(len(self._npaFeatureAbundance.dtype.names[1:])),
      "\nFeature count:", str(len(self._npaFeatureAbundance[self._npaFeatureAbundance.dtype.names[0]])),
      "\nId Metadata:", self._npaFeatureAbundance.dtype.names[0],
      "\nMetadata ids:", str(self._dictTableMetadata.keys()),
      "\nMetadata count:", str(len(self._dictTableMetadata.keys())),
      "\nOriginating source:",self._strOriginalName,
      "\nOriginal feature count:", str(self._iOriginalFeatureCount),
      "\nOriginal sample count:", str(self._iOriginalSampleCount),
      "\nIs normalized:", str(self._fIsNormalized),
      "\nIs summed:", str(self._fIsSummed),
      "\nCurrent filtering state:", str(self._iCurrentFilterState),
      "\nFeature delimiter:", self._cFeatureDelimiter,
      "\nFile delimiter:",self._cDelimiter])

    #Private Methods

    #Testing Status: Light happy path testing
    @staticmethod
    def _funcTextToStructuredArray(strInputFile = None, cDelimiter = Constants.TAB, sMetadataID = None, sLastMetadata = None):
        """
        Private method
        Used to read in a file that is samples (column) and taxa (rows) into a structured array.

        :param	strInputFile:	Path to input file.
        :type	String		String path.
        :param	cDelimiter:	Delimiter for parsing the input file.
        :type	Character	Character.
        :param	sMetadataID:	String ID that is a metadata row ID (found on the first column) and used as an ID for samples
        :type	String		String ID
        :param	sLastMetadata:	The ID of the metadata that is the last metadata before measurement or feature rows.
        :type	String		String ID
        :return	[taxData,metadata]:	Numpy Structured Array of abundance data and dictionary of metadata.
                                        Metadata is a dictionary as such {"ID", [value,value,values...]}
                                        Values are in the order thety are read in (and the order of the sample names).
                                        ID is the first column in each metadata row.
        :type [Numpy structured Array, Dictionary]
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strInputFile)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, input file not valid. File:",str(strInputFile)])
            return False
        if(not ValidateData.funcIsValidStringType(cDelimiter)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, cDelimiter was invalid. Value =",str(cDelimiter)])
            return False
        if(not ValidateData.funcIsValidString(sMetadataID)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, sMetadataID was invalid. Value =",str(sMetadataID)])
            return False
        if(not ValidateData.funcIsValidString(sLastMetadata)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, sLastMetadata was invalid. Value =",str(sLastMetadata)])
            return False

        #Read in file
        with open(strInputFile,'r') as f:
            contents = f.read()

        #Turn to lines of the file
        contents = contents.replace("\"","")
        contents = filter(None,contents.split(Constants.ENDLINE))

        #Get metadata and sample ids
        iFirstDataRow = -1
        namesRow = None
        metadata = dict()
        for iIndex, line in enumerate(contents):
            lsLineElements = line.split(cDelimiter)
            sIdElement = lsLineElements[0].strip()
            if sIdElement == sMetadataID:
                namesRow = lsLineElements
            metadata[sIdElement]=lsLineElements[1:]
            if sIdElement == sLastMetadata:
                iFirstDataRow = iIndex + 1
                break

        #Check to make sure the aLastMetadata was found
        if iFirstDataRow == -1:
            print "".join(["AbundanceTable:textToStructuredArray::Error, sLastMetadata was not found. Value =",str(sLastMetadata)])
            return False

        #Check to make sure there is abundance
        if len(contents) <= iFirstDataRow:
            return [np.array([]),metadata]

        #Make sure the names are found
        if namesRow == None:
            print "".join(["AbundanceTable:textToStructuredArray::Error, did not find the row for the unique sample/column. File:",str(strInputFile)," Identifier:",str(sMetadataID)])
            return False

        #Build data type object (data name,data type)
        incompleteDataTypeVector = []
        for sampleColumn in xrange(1,len(namesRow)):
            incompleteDataTypeVector.append((namesRow[sampleColumn],'f8'))

        #Extract data
        #Create the tuple with the first data row
        rowData = contents[iFirstDataRow]
        rowData = rowData.split(cDelimiter)
        taxId = rowData[0]
        longestTaxId = len(taxId)
        sampleReads = rowData[1:]
        tempSampleReads = [taxId]
        for reads in sampleReads:
            tempSampleReads.append(float(reads))
        sampleReads = tempSampleReads
        dataMatrix = [tuple(sampleReads)]

        #Add the rest of the rows
        for rowData in contents[iFirstDataRow+1:]:
            rowData = rowData.split(cDelimiter)
            taxId = rowData[0]
            tempTaxIdLen = len(taxId)
            if(tempTaxIdLen > longestTaxId):
                longestTaxId = len(taxId)
            sampleReads = rowData[1:]
            tempSampleReads = [taxId]
            for reads in sampleReads:
                tempSampleReads.append(float(reads))
            sampleReads = tempSampleReads
            #Validate data before adding to matrix
            if(ValidateData.funcIsValidString(taxId)):
                if(ValidateData.funcIsValidList(sampleReads)):
                    if(len(sampleReads)>0):
                        dataMatrix.append(tuple(sampleReads))

        #Now we know the longest taxId we can define the first column holding the tax id
        dataTypeVector = [(namesRow[0],'a'+str(longestTaxId*2))]
        dataTypeVector.extend(incompleteDataTypeVector)
        dataTypeVector = np.dtype(dataTypeVector)

        #Create structured array
        taxData = np.array(dataMatrix,dtype=dataTypeVector)

        return [taxData,metadata]

    #Happy path tested
    def funcGetSampleNames(self):
        """
        Returns the sample names (IDs) contained in the abundance table.

        :return	Sample Name:	A List of sample names indicated by the metadata associated with the sMetadataId given in table creation.
        :type	List of strings	A list tof string names or empty list on error as well as no underlying table.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.dtype.names[1:]
        else:
            return []

    #Happy Path Tested
    def funcGetIDMetadataName(self):
        """
        Returns the metadata id.

        :return	ID:	The metadata id (the sample Id).
        :type	String	Returns none on error.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.dtype.names[0]
        else:
            return None

    #Happy path tested
    def funcGetAbundanceCopy(self):
        """
        Returns a deep copy of the abundance table.

        :return	Numpy Structured Array:	The measurement data in the Abundance table. Can use sample names to access each column of measurements.
        :type	Numpy Structured Array	Returns none on error.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.copy()
        return None

    #Happy path tested
    def funcGetAverageAbundancePerSample(self, lsTargetedFeatures):
	"""
	Averages feature abundance.

	:param	lsTargetedFeatures:	String names of features to average
	:type	list	list of string names of features which are measured
        :return	List: of lists or boolean:	List of lists or False on error. One internal list per sample indicating the sample and the feature's average abudance
        :type	list	[[sample,average abundance of selected taxa]] or False on error
	"""

        #Sample rank averages [[sample,average abundance of selected taxa]]
        #Returned
        sampleAbundanceAverages = []
        
        #Get sample names
        sampleNames = self.funcGetSampleNames()
        #Get taxa names
        allTaxaNames = self.funcGetFeatureNames()

        #Get an abundance table compressed to features of interest
        abndReducedTable = self.funcGetFeatureAbundanceTable(lsTargetedFeatures)
        if abndReducedTable == None:
            return False  

        #If the taxa to be selected are not in the list
        #Return nothing and log
        lsMissing = []
        for sFeature in lsTargetedFeatures:
            if not sFeature in allTaxaNames:
                lsMissing.append(sFeature)
            else:
                #Check to make sure the taxa of interest is not average abundance of 0
                iFeatureSum = abndReducedTable.funcGetFeatureSumAcrossSamples(sFeature)
                if (iFeatureSum == None) or (iFeatureSum == 0):
                    lsMissing.append(sFeature)

        if len(lsMissing) > 0:
            return False        

        #For each sample name get average abundance
        for sName in sampleNames:
            npaFeaturesSample = abndReducedTable.funcGetSample(sName)
            sampleAbundanceAverages.append([sName,sum(npaFeaturesSample)/float(len(npaFeaturesSample))])

        #Sort based on average
        sampleAbundanceAverages = sorted(sampleAbundanceAverages, key = lambda sampleData: sampleData[1], reverse = True)

        #return
        return sampleAbundanceAverages

    #Happy Path Tested
    def funcGetFeatureAbundanceTable(self, lsFeatures):
        """
        Returns a copy of the current abundance table with the abundance of just the given features.

        :param	lsFeatures:	String Feature IDs that are kept in the compressed abundance table.
        :type	List of strings	Feature IDs (found as the first entry of a filter in the input file.
        :return	AbundanceTable:	A compressed version of the abundance table.
        :type	AbundanceTable	On an error None is returned.
        """

        if (not self._npaFeatureAbundance == None) and lsFeatures:
            #Get a list of boolean indicators that the row is from the features list
            lfFeatureData = [sRowID.strip() in lsFeatures for sRowID in self.funcGetFeatureNames()]

            #compressed version as an Abundance table
            strFeatureName = "".join([os.path.splitext(self._strOriginalName)[0],"-",str(len(lsFeatures)),"-Features.txt"])
            abndTableRet = AbundanceTable(npaAbundance=np.compress(lfFeatureData, self._npaFeatureAbundance, axis = 0),
                           dictMetadata = self.funcGetMetadataCopy(), strName = strFeatureName,
                           fIsNormalized = self.funcIsNormalized(), fIsSummed = self.funcIsSummed(),
                           cFileDelimiter = self.funcGetFileDelimiter(), cFeatureNameDelimiter= self.funcGetFeatureDelimiter())

            #Make sure the features are all found
            if ",".join(sorted(lsFeatures)) == ",".join(sorted(abndTableRet.funcGetFeatureNames())):
                return abndTableRet
        return None

    #Happy path tested
    def funcGetFeatureDelimiter(self):
        """
        The delimiter of the feature names (For example to use on concensus lineages).

        :return	Delimiter:	Delimiter for the feature name pieces if it is complex.
        :type	Character	Delimiter.
        """

        return self._cFeatureDelimiter

    #Happy path tested
    def funcGetFeatureCount(self):
        """
        Returns the current feature count.

        :return	Count:	Returns the int count of features in the abundance table.
        :type	Int	Returns None on error.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.shape[0]
        else:
            return None

    #Happy path tested
    def funcGetFeatureSumAcrossSamples(self,sFeatureName):
        """
        Returns float sum of feature values across the samples.

        :param	sFeatureName:
        :type	String.
        :return	Sum:	Sum of one feature across samples.
        :type	Double
        """

        if (not self._npaFeatureAbundance == None):
            for sFeature in self._npaFeatureAbundance:
                if sFeature[0] == sFeatureName:
                    return sum(list(sFeature)[1:])
        else:
            return None

    #Happy path tested
    def funcGetFeatureNames(self):
        """
        Return the feature names as a list.

        :return	Feature Names:	List of feature names (or IDs) as strings.
        :type	List of strings.	As an error returns empty list.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[self.funcGetIDMetadataName()]
        return []

    #Happy path tested
    def funcGetFileDelimiter(self):
        """
        The delimiter of the file the data was read from and which is also the delimiter which would be used to write the data to a file.

        :return	Delimiter:	Delimiter for the parsing and writing the file.
        :type	Character	Delimiter.
        """

        return self._cDelimiter

    #Happy path tested
    def funcGetSample(self,sSampleName):
        """
        Return a copy of the feature measurements of a sample.

        :param	sSampleName:	Name of sample to return.	
        :type	String	
        :return	Sample: Measurements	Feature measurements of a sample.
        :type	Numpy Array	Empty numpy array returned on error.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[sSampleName].copy()
        return np.array([])

    #Happy path tested
    def funcGetMetadata(self, strMetadataName):
        """
        Returns a list of metadata that is associated with the given metadata name (id).

        :param	strMetadataName:	String metadata ID to be returned
        :type	String	ID
        :return	Metadata:	List of metadata
        :type	List	
        """

        if self._dictTableMetadata:
            retValue = self._dictTableMetadata.get(strMetadataName,None)
            if retValue:
                retValue = copy.deepcopy(retValue)
            return retValue
        return None

    #Happy path tested
    def funcGetMetadataCopy(self):
        """
        Returns a deep copy of the metadata.

	:return	Metadata copy:	{"ID":[value,value...]}
        :type	Dictionary
        """

        if self._dictTableMetadata:
            return copy.deepcopy(self._dictTableMetadata)
        return None

    #Happy path tested
    def funcGetName(self):
        """
        Returns the name of the object which is the file name that generated it.
        If the object was generated from an Abundance Table (for instance through stratification)
        the name is still in the form of a file that could be written to which is informative
        of the changes that have occured on the dataset.
        """

        return self._strOriginalName

    #Happy path tested. could do more
    def funcGetTerminalNodes(self):
        """
        Returns the terminal nodes given the current feature names in the abundance table. The 
        features must contain a consensus lineage or all will be returned.
        """
        return AbundanceTable.funcGetTerminalNodesFromList(lsNames=self.funcGetFeatureNames(),cNameDelimiter=self.funcGetFeatureDelimiter())

    @staticmethod
    def funcGetTerminalNodesFromList(lsNames,cNameDelimiter):
        """
        Returns the terminal nodes given the current feature names in the abundance table. The 
        features must contain a consensus lineage or all will be returned.

        :param	lsNames:	The list of string names to parse and filter.
        :type	List of strings
        :param	cNameDelimiter:	The delimiter for the name of the features.
        :type	Character	Delimiter
        """

        #Return list
        lsRetList = list()

        #Build hash
        dictCounts = dict()
        for strTaxaName in lsNames:
            #Split into the elements of the clades
            lsClades = filter(None,strTaxaName.split(cNameDelimiter))
            #Count clade levels
            iCladeLength = len(lsClades)
    
            #Make sure there is data to work with
            if iCladeLength < 0:
                pass

            #Evaluate first element
            sClade = lsClades[0]
            if sClade in dictCounts:
                dictCounts[sClade] = False
            else:
                dictCounts[sClade] = True

            #Evaluate the rest of the elements
            if iCladeLength > 1:
                for iIndex in xrange(1,iCladeLength):
                    prevClade = sClade
                    sClade = cNameDelimiter.join([sClade,lsClades[iIndex]])
                    if sClade in dictCounts:
                        dictCounts[sClade] = False
                        dictCounts[prevClade] = False
                    else:
                        dictCounts[sClade] = True
                        dictCounts[prevClade] = False

        #Return only the elements that were of count 1
        for sName in dictCounts:
            if dictCounts[sName]==True:
                lsRetList.append(sName)
        return lsRetList

    #Happy path tested
    def funcIsNormalized(self):
        """
        Returns if the data has been normalized.

        :return	Boolean:	Indicates if the data is normalized.
        :type	Boolean	True indicates it the data is normalized.
        """

        return self._fIsNormalized

    #Happy path tested
    def funcIsPrimaryIdMetadata(self,sMetadataName):
        """
        Checks the metadata data associatd with the sMetadatName and returns if the metadata is unique.
        This is important to some of the functions in the Abundance Table specifically when translating from one metadata to another.
        
        :param	sMetadataName:	ID of metadata to check for uniqueness.
        :type	String	Metadata ID.
        :return	Boolean:	Returns indicator of uniqueness.
        :type	Boolean	True indicates unique.
        """

        lMetadata = self.funcGetMetadata(sMetadataName)
        if not lMetadata:
            return False
        return (len(lMetadata) == len(set(lMetadata)))

    #Happy path tested
    def funcIsSummed(self):
        """
        Return is the data is summed.

        :return	Boolean:	Indicator of being summed. 
        :type	Boolean	True indicates summed.
        """

        return self._fIsSummed

    #Happy path tested
    def funcFilterAbundanceByPercentile(self, dPercentileCutOff = 95.0, dPercentageAbovePercentile=1.0):
        """
        Filter on features.
        A feature is removed if it's abundance is not found in the top X percentile a certain percentage of the samples.

        :param	dPercentileCutOff:	The percentile used for filtering.
        :type	double	A double between 0.0 and 100.0
        :param	dPercentageAbovePercentile:	The percentage above the given percentile (dPercentileCutOff) that must exist to keep the feature.
        :type	double	Between 0.0 and 100.0
        :return	Boolean:	Indicator of filtering occuring without error.	
        :type	Boolean	True indicates filtering occuring.
        """

        #No need to do anything
        if(dPercentileCutOff==0.0) or (dPercentageAbovePercentile==0.0):
            return True

        #Sample names
        lsSampleNames = self.funcGetSampleNames()

        #Scale percentage out of 100
        dPercentageAbovePercentile = dPercentageAbovePercentile/100.0

        #Sample count
        iSampleCount = len(lsSampleNames)

        #Get a threshold score of the value at the specified percentile for each sample
        #In the order of the sample names
        ldScoreAtPercentile = [scipy.stats.scoreatpercentile(self._npaFeatureAbundance[lsSampleNames[iIndex]],dPercentileCutOff) for iIndex in xrange(iSampleCount)]

        #Record how many entries for each feature have a value equal to or greater than the dPercentileCutOff
        #If the percentile of entries passing the criteria are above the dPercentageAbovePercentile put index in list to keep
        liKeepIndices = []
        iSampleCount = float(iSampleCount)
        for iRowIndex, npaRow in enumerate(self._npaFeatureAbundance):
            iCountPass = sum([1 if dValue >= ldScoreAtPercentile[iValueIndex] else 0 for iValueIndex, dValue in enumerate(list(npaRow)[1:])])
            if (iCountPass / iSampleCount) >= dPercentageAbovePercentile:
                liKeepIndices.append(iRowIndex)

        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepIndices,:]

        #Update filter state
        self._iCurrentFilterState = "".join([self._iCurrentFilterState, ":", "dPercentileCutOff", "=", str(dPercentileCutOff),
                                             ",", "dPercentageAbovePercentile", "=", str(dPercentageAbovePercentile)])

        return True

    #Happy path tested
    def funcFilterAbundanceBySequenceOccurence(self, iMinSequence = 2, iMinSamples = 2):
        """
        Filter abundance by requiring features to have a minimum sequence occurence in a minimum number of samples.
        Will evaluate greater than or equal to the iMinSequence and iMinSamples.

        :param	iMinSequence:	Minimum sequence to occur.
        :type	Integer	Number Greater than 1.
        :param	iMinSamples:	Minimum samples to occur in.
        :type	Integer	Number greater than 1.
        :return	Boolean:	Indicator of the filter running without error.
        :type	Boolean	False indicates error.
        """

        #No need to do anything
        if(iMinSequence==0) or (iMinSamples==0):
            return True

        #This normalization requires the data to be reads
        if self._fIsNormalized:
            print "Could not filter by sequence occurence because the data is already normalized."
            return False

        #Holds which indexes are kept
        liKeepFeatures = []
        for iRowIndex, dataRow in enumerate(self._npaFeatureAbundance):
            #See which rows meet the criteria and keep the index if needed.
            if sum([1 if fRead >= iMinSequence else 0 for fRead in list(dataRow)[1:]]) >= iMinSamples:
                liKeepFeatures.append(iRowIndex)

        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]

        #Update filter state
        self._iCurrentFilterState = "".join([self._iCurrentFilterState, ":", "iMinSequence", "=", str(iMinSequence),
                                             ",", "iMinSamples", "=", str(iMinSamples)])

        return True
   
    #1 Happy path test
    def funcFilterFeatureBySD(self, dMinSDCuttOff = 0.0):
        """
        A feature is removed if it's abundance is not found to have standard deviation more than the given dMinSDCutoff.

        :param	dMinSDCuttOff:	Standard deviation threshold.
        :type	Double	A double greater than 0.0.
        :return	Boolean:	Indicator of success.
        :type	Boolean	False indicates error.
        """

        #No need to do anything
        if(dMinSDCuttOff==0.0):
            return True

        #Holds which indexes are kept
        liKeepFeatures = []

        #Evaluate each sample
        for iRowIndex, dataRow in enumerate(self._npaFeatureAbundance):
            if(np.std(list(dataRow)[1:])>=dMinSDCuttOff):
                liKeepFeatures.append(iRowIndex)
        
        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]

        #Update filter state
        self._iCurrentFilterState = "".join([self._iCurrentFilterState, ":", "dMinSDCuttOff", "=", str(dMinSDCuttOff)])

        return True

    #Happy path tested
    def funcNormalize(self):
        """
        Convenience method which will call which ever normalization is approriate on the data.
        """

        if self._fIsSummed:
            return self.funcNormalizeColumnsWithSummedClades()
        else:
            return self.funcNormalizeColumnsBySum()

    #Testing Status: Light happy path testing
    def funcNormalizeColumnsBySum(self):
        """
        Normalize the data in a manner that is approrpiate for NOT summed data.
        Normalize the columns (samples) of the abundance table.
        Normalizes as a fraction of the total (number/(sum of all numbers in the column)).
        Will not act on summed tables.

        :return	Boolean:	Indicator of success.
        :type	Boolean	False indicates error.
        """

        if self._fIsNormalized:
            print "This table is already normalized, did not perform new normalization request."
            return False

        if self._fIsSummed:
            print "This table has clades summed, this normalization is not appropriate. Did not perform."
            return False

        #Normalize
        for columnName in self.funcGetSampleNames():
            column = self._npaFeatureAbundance[columnName]
            columnTotal = sum(column)
            if(columnTotal > 0.0):
                column = column/columnTotal
            self._npaFeatureAbundance[columnName] = column

        #Indicate normalization has occured
        self._fIsNormalized = True

        return True

    #Happy path tested
    def funcNormalizeColumnsWithSummedClades(self):
        """
        Normalizes a summed Abundance Table.
        If this is called on a dataset which is not summed and not normalized.
        The data will be summed first and then normalized.
        If already normalized, the current normalization is kept.

        :return	Boolean:	Indicator of success.
        :type	Boolean	False indicates error.
        """

        if self._fIsNormalized:
            print "This table is already normalized, did not perform new normalization request."
            return False

        if not self._fIsSummed:
            print "This table does not have clades summed, this normalization is not appropriate until the clades are summed. The clades are being summed now before normalization."
            self.funcSumClades()

        #Load a hash table with root data {sKey: npaAbundances}
        hashRoots = {}
        for npaRow in self._npaFeatureAbundance:

            curldAbundance = np.array(list(npaRow)[1:])
            curFeatureNameLength = len(npaRow[0].split(self._cFeatureDelimiter))
            curlRootData = hashRoots.get(npaRow[0].split(self._cFeatureDelimiter)[0])

            if not curlRootData:
                hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]
            elif curlRootData[0] > curFeatureNameLength:
                hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]

        #Normalize each feature by thier root feature
        dataMatrix = list()
        for npaRow in self._npaFeatureAbundance:

            curHashRoot = list(hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]][1])
            dataMatrix.append(tuple([npaRow[0]]+[npaRow[i+1]/curHashRoot[i] if curHashRoot[i] > 0 else 0 for i in xrange(len(curHashRoot))]))

        self._npaFeatureAbundance = np.array(dataMatrix,self._npaFeatureAbundance.dtype)

        #Indicate normalization has occured
        self._fIsNormalized = True

        return True

    #1 Happy path test
    def funcRankAbundance(self):
        """
        Rank abundances of features with in a sample.

        :return	AbundanceTable:	Abundance table data ranked (Features with in samples).
        :type	AbundanceTable	None is returned on error.
        """

        if not self._npaFeatureAbundance == None:

            #Sample rank averages [[sample,average rank of selected taxa]]
            #Returned
            sampleRankAverages = []
        
            #Get sample names
            lsSampleNames = self.funcGetSampleNames()
            #Get taxa names
            allTaxaNames = self.funcGetFeatureNames()
            #Get abundance copy
            npRankAbundance = self.funcGetAbundanceCopy()
            #If there are no features then return self.
            if len(allTaxaNames)==0:
               return self 

            liRanks = []
            #For each sample name get the ranks
            for sName in lsSampleNames:
            
                #Enumerate for order and sort abundances
                lfSample = list(enumerate(npRankAbundance[sName]))
                lfSample = sorted(lfSample, key = lambda sampleData: sampleData[1], reverse = True)

                #Replace abundance with rank
                fPreviousAbundance = None
                iRank = 0
                iRankedList = []
                for iOrder, fAbundance in lfSample:
                    if not fPreviousAbundance == fAbundance:
                        iRank = iRank + 1
                    iRankedList.append((iOrder,iRank))
                    fPreviousAbundance = fAbundance

                #Sort back to original order
                iRankedList = sorted(iRankedList, key = lambda sampleData: sampleData[0], reverse = False)
                npRankAbundance[sName] = np.array([tuplData[1] for tuplData in iRankedList])

            lsPathPieces = os.path.splitext(self.funcGetName())
            return AbundanceTable(npaAbundance=npRankAbundance, dictMetadata=self.funcGetMetadataCopy(),
                  strName="".join([lsPathPieces[0],"-Ranked",lsPathPieces[1]]), fIsNormalized=self.funcIsNormalized(),
                  fIsSummed=self.funcIsSummed(), cFileDelimiter=self.funcGetFileDelimiter(),
                  cFeatureNameDelimiter=self.funcGetFeatureDelimiter())

        return None

    #Happy Path Tested
    def funcReduceFeaturesToCladeLevel(self, iCladeLevel):
        """
        Reduce the current table to a certain clade level.

        :param	iCladeLevel:	The level of the clade to trim the features to.
        :type	Integer	The higher the number the more clades are presevered in the consensus lineage contained in the feauture name.
        :return	Boolean:	Indicator of success. False indicates error.
        """

        if iCladeLevel < 1: return False
        if not self._npaFeatureAbundance == None:
            liFeatureKeep = []
            [liFeatureKeep.append(tplFeature[0]) if (len(tplFeature[1][0].split(self.funcGetFeatureDelimiter())) <= iCladeLevel) else 0
             for tplFeature in enumerate(self._npaFeatureAbundance)]
            #Compress array
            self._npaFeatureAbundance = self._npaFeatureAbundance[liFeatureKeep,:]

            #Update filter state
            self._iCurrentFilterState = "".join([self._iCurrentFilterState, ":", "iCladeLevel", "=", str(iCladeLevel)])
            return True
        else:
            return False

    #Happy path tested
    def funcRemoveSamples(self,lsSampleNames):
        """
        Removes the samples given in the list.

        :param	lsSampleNames:	A list of string names of samples to remove.
        :type	List of strings	Unique values
        :return Boolean:	Indicator of success (True = success, no error)
        """

        #Samples to remove
        setSamples = set(lsSampleNames)

        #Get orignal sample count
        iOriginalCount  = self._iOriginalSampleCount

        #The samples to keep
        lsKeepSamples = [sSample for sSample in self.funcGetSampleNames() if not sSample in setSamples]
        #The sample to keep as boolean flags for compressing the metadata
        lfKeepSamples = [not sSample in setSamples for sSample in self.funcGetSampleNames()]
        
        #Reduce the abundance data and update
        self._npaFeatureAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+lsKeepSamples]

        #Reduce the metadata and update
        for sKey in self._dictTableMetadata:
            self._dictTableMetadata[sKey] = [value for iindex, value in enumerate(self._dictTableMetadata[sKey]) if lfKeepSamples[iindex]]

        #Update sample number count
        self._iOriginalSampleCount = len(self.funcGetSampleNames())

        return self._iOriginalSampleCount == (iOriginalCount-len(setSamples))

    #Happy path tested
    def funcRemoveSamplesByMetadata(self, sMetadata, lValuesToRemove):
        """
        Removes samples from the abundance table based on values of a metadata.
        If a metadata has any value given the associated sample is removed.

        :param	sMetadata:	ID of the metdata to check the given values.
        :type	String	Metadata ID
        :param	lValuesToRemove:	A list of values which if equal to a metadata entry indicate to remove the associated sample.
        :type	List of values:	List
        :return	Boolean:	Indicator of success (True = success, no error)
        """

        lsSampleNames = self.funcGetSampleNames()
        return self.funcRemoveSamples([lsSampleNames[iindex] for iindex, sValue in enumerate(self.funcGetMetadata(sMetadata)) if sValue in lValuesToRemove])

    #Happy path testing
    def funcSumClades(self):
        """
        Sums abundance data by clades indicated in the feature name (as consensus lineages).

        :return	Boolean:	Indicator of success.
        :type	Boolean	False indicates an error.
        """

        if not self.funcIsSummed():

            #Read in the data
            #Find the header column (iCol) assumed to be 1 or 2 depending on the location of "NAME"
            #Create a list (adSeq) that will eventually hold the sum of the columns of data
            astrHeaders = iCol = None
            adSeqs = np.array([0] * len(self.funcGetSampleNames()))
            pTree = CClade( )
            aastrRaw = []

            #For each row in the npaAbundance
            #Get the feature name, feature abundances, and sum up the abudance columns
            #Keep the sum for later normalization
            #Give a tree the feature name and abundance
            for dataRow in self._npaFeatureAbundance:
                
                sFeatureName = dataRow[0]
                ldAbundances = list(dataRow)[1:]

                #Add to the sum of the columns (samples)
                adSeqs = adSeqs + np.array(list(dataRow)[1:])

                #Build tree
                pTree.get( sFeatureName.split(self._cFeatureDelimiter) ).set( ldAbundances )

            #Create tree of data
            #Input missing data
            #Fill hashFeatures with the clade name (key) and a blist of values (value) of the specified level interested.
            pTree.impute( )
            hashFeatures = {}
            pTree.freeze( hashFeatures, c_iSumAllCladeLevels, c_fOutputLeavesOnly )
            setstrFeatures = hashFeatures.keys( )

            #Remove parent clades that are identical to child clades
            for strFeature, adCounts in hashFeatures.items( ):
                    astrFeature = strFeature.strip( ).split( "|" )
                    while len( astrFeature ) > 1:
                        astrFeature = astrFeature[:-1]
                        strParent = "|".join( astrFeature )
                        adParent = hashFeatures.get( strParent )
                        if adParent == adCounts:
                            del hashFeatures[strParent]
                            setstrFeatures.remove( strParent )

            #Sort features to be nice
            astrFeatures = sorted( setstrFeatures )

            #Change the hash table to an array
            dataMatrix = list()
            for sFeature in astrFeatures:
                dataMatrix.append(tuple([sFeature]+list(hashFeatures[sFeature])))
            self._npaFeatureAbundance=np.array(dataMatrix,self._npaFeatureAbundance.dtype)

            #Indicate summation has occured
            self._fIsSummed = True

        return True

    #Happy path tested
    def funcStratifyByMetadata(self, strMetadata, fWriteToFile=False):
        """
        Stratifies the AbundanceTable by the given metadata.
        Will write each stratified abundance table to file
        if fWriteToFile is True the object will used it's internally stored name as a file to write to
        if fWriteToFile is a string then it should be a directory and end with "." This will rebase the file
        and store it in a different directory but with an otherwise unchanged name.
        Note: If the metadata used for stratification has NAs, they will be segregated to thier own table and returned.

        :param	strMetadata:	Metadata ID to stratify data with.
        :type	String	ID for a metadata.
        :param	fWriteToFile:	Indicator to write to file.
        :type	Boolean	True indicates to write to file.
        :return	List:	List of AbundanceTables which are deep copies of the original.
        :type	List	List of AbundanceTables. Empty list on error.
        """

        retlAbundanceTables = []

        if (not self._npaFeatureAbundance == None) and self._dictTableMetadata :
            dictAbundanceBlocks = dict()
            #Get unique metadata values to stratify by
            lsMetadata = self._dictTableMetadata.get(strMetadata,[])
            setValues = set(lsMetadata)

            #If there is only one metadata value then no need to stratify so return the original in the list (and write if needed)
            if len(setValues) == 0:
              return retlAbundanceTables
            elif len(setValues) == 1:
              retlAbundanceTables.append(self)
              #Write to file if need be before returning data
              if fWriteToFile:
                  if isinstance(fWriteToFile, basestring):
                      self.funcWriteToFile("".join([fWriteToFile,os.path.split(self.getName())[1]]))
                  else:
                      self.funcWriteToFile(self.getName())
              return retlAbundanceTables

            #Given here there are multiple metadata values, continue to stratify
            lsNames = self.funcGetSampleNames()
            #Get index of values to break up
            for value in setValues:
                fDataIndex = [sData==value for sData in lsMetadata]
                #Get abundance data for the metadata value
                #The true is added to keep the first column which should be the feature id
                npaStratfiedAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+list(np.compress(fDataIndex,lsNames))]

                #Get metadata for the metadata value
                dictStratifiedMetadata = dict()
                for metadataType in self._dictTableMetadata:
                    dictValues = self.funcGetMetadata(metadataType)
                    dictStratifiedMetadata[metadataType] = np.compress(fDataIndex,dictValues).tolist()

                #Make abundance table
                #Add abundance table to the list
                strStratifiedName = "".join([os.path.splitext(self._strOriginalName)[0],"-StratBy-",value,".txt"])
                objStratifiedAbundanceTable = AbundanceTable(npaAbundance=npaStratfiedAbundance, dictMetadata=dictStratifiedMetadata, strName=strStratifiedName,
                                              fIsNormalized=self._fIsNormalized, fIsSummed=self._fIsSummed, cFeatureNameDelimiter=self._cFeatureDelimiter)

                #Write to file if needed
                if fWriteToFile:
                    if isinstance(fWriteToFile, basestring):
                        objStratifiedAbundanceTable.funcWriteToFile("".join([fWriteToFile,os.path.split(objStratifiedAbundanceTable.funcGetName())[1]]))
                    else:
                        objStratifiedAbundanceTable.funcWriteToFile(objStratifiedAbundanceTable.funcGetName())

                #Append abundance table to returning list
                retlAbundanceTables.append(objStratifiedAbundanceTable)

        return retlAbundanceTables

    #Happy Path Tested
    def funcTranslateIntoMetadata(self, lsValues, sMetadataFrom, sMetadataTo, fFromPrimaryIds=True):
        """
        Takes the given data values in one metadata and translates it to values in another
        metadata of the sample samples holding the values of the first metadata
        FPrimaryIds, if true the sMetadataFrom are checked for unique values,
        If FPrimaryIds is not true, duplicate values can stop the preservation of order
        Or may cause duplication in the "to" group. This is not advised.
        if the sMetadataFrom has any duplicates the function fails and return false.

        :param	lsValues:	Values to translate.
        :type	List	List of values.
        :param	sMetadataFrom:	The metadata the lsValues come from.
        :type	String	ID for the metadata.
        :param	sMetadataTo:	The metadata the lsValues will be translated into keeping the samples the same.
        :type	String	ID for the metadata.
        :param	fFromPrimaryIds:	The metadata that are in the from metadata list must be unique in each sample.
        :type	Boolean	True indicates the metadata list should be unique in each sample. Otherwise a false will return.
        """

        #Get metadata
        lFromMetadata = self.funcGetMetadata(sMetadataFrom)
        if not lFromMetadata:
                print "Abundancetable::funcTranlateIntoMetadata. Did not receive lFromMetadata."
                return False

        lToMetadata = self.funcGetMetadata(sMetadataTo)
        if not lToMetadata:
                print "Abundancetable::funcTranlateIntoMetadata. Did not receive lToMetadata."
                return False

        #Check to see if the values are unique if indicated to do so
        if fFromPrimaryIds:
            if not len(lFromMetadata) == len(set(lFromMetadata)):
                print "Abundancetable::funcTranlateIntoMetadata. sMetadataFrom did not have unique values."
                return False

        #Translate over
        if lFromMetadata and lToMetadata:
            return [lToMetadata[iIndex] for iIndex in [lFromMetadata.index(value) for value in lsValues]]

        return False

    #Happy path tested
    def funcToArray(self):
        """
        Returns a numpy array of the current Abundance Table.
        Removes the first ID head column and the numpy array is
        Made of lists, not tuples.

        :return Numpy Array:	np.array([[float,float,...],[float,float,...],[float,float,...]])
        :type	Numpy Array	None is returned on error.
        """

        if not self._npaFeatureAbundance == None:
            return np.array([list(tplRow)[1:] for tplRow in self._npaFeatureAbundance],'float')
        return None

    #Happy Path tested
    def funcWriteToFile(self, strOutputFile, cDelimiter=None):
        """
        Writes the AbundanceTable to a file strOutputFile.
        Will rewrite over a file as needed.
        Will use the cDelimiter to delimit columns if provided.

        :param	strOutputFile:	File path to write the file to.
        :type	String	File Path
        :param	cDelimiter:	Delimiter for the output file.
        :type	Character	If cDlimiter is not specified, the internally stored file delimiter is used.
        """

        with open(strOutputFile, 'w') as f:
            #Check delimiter argument
            if not cDelimiter:
                cDelimiter = self._cDelimiter 
            #Write Ids
            f.write(cDelimiter.join([self.funcGetIDMetadataName()]+list(self.funcGetSampleNames()))+Constants.ENDLINE)
            #Write metadata
            lsKeys = list(set(self._dictTableMetadata.keys())-set([self.funcGetIDMetadataName()]))
            f.write(Constants.ENDLINE.join([cDelimiter.join([sMetaKey]+self.funcGetMetadata(sMetaKey)) for sMetaKey in lsKeys])+Constants.ENDLINE)
            #Write abundance
            lsOutput = list()
            curAbundance = self._npaFeatureAbundance.tolist()
            for curAbundanceRow in curAbundance:
                lsOutput.append(cDelimiter.join([str(curAbundanceElement) for curAbundanceElement in curAbundanceRow]))
            f.write(Constants.ENDLINE.join(lsOutput))

    #Testing Status: 1 Happy path test
    @staticmethod
    def funcPairTables(strFileOne, strFileTwo, strIdentifier, cDelimiter, strOutFileOne, strOutFileTwo, lsIgnoreValues=None):
        """
        This method will read in two files and abridge both files (saved as new files)
        to just the samples in common between the two files given a common identifier.
        ***If the identifier is not unique in each data set, the first sample with the pairing id is taken so make sure the ID is unique.
        Expects the files to have the sample delimiters.

        :param	strFileOne:	Path to file one to be paired.
        :type	String	File path.
        :param	strFileTwo:	Path to file two to be paired.
        :type	String	File path.
        :param	strIdentifier:	Metadata ID that is used for pairing.
        :type	String	Metadata ID.
        :param	cDelimiter:	Character delimiter to read the files.
        :type	Character	Delimiter.
        :param	strOutFileOne:	The output file for the paired version of the first file.
        :type	String	File path.
        :param	strOutFileTwo:	The output file for the paired version of the second file.
        :type	String	File path.
        :param	lsIgnoreValues:	These values are ignored even if common IDs between the two files.
        :type	List	List of strings.
        :return	Boolean:	Indicator of no errors.
        :type	Boolean	False indicates errors.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strFileOne)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strFileOne)
            return False
        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strFileTwo)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strFileTwo)
            return False

        #Make file one
        #Read in file
        with open(strFileOne,'r') as f:
            sContentsOne = f.read()

        #Get the file identifier for file one
        fileOneIdentifier = None
        for sLine in filter(None, sContentsOne.split(Constants.ENDLINE)):
            lsLineContents = sLine.split(cDelimiter)
            if lsLineContents[0] == strIdentifier:
                fileOneIdentifier = lsLineContents
                break

        #Make file two
        #Read in file
        with open(strFileTwo,'r') as f:
            sContentsTwo = f.read()

        #Get the file identifier for file two
        fileTwoIdentifier = None
        for sLine in filter(None, sContentsTwo.split(Constants.ENDLINE)):
            lsLineContents = sLine.split(cDelimiter)
            if lsLineContents[0] == strIdentifier:
                fileTwoIdentifier = lsLineContents
                break

        #Get what is in common between the identifiers
        #And find which columns to keep in the tables based on the common elements
        setsCommonIdentifiers = set(fileOneIdentifier) & set(fileTwoIdentifier)
        if lsIgnoreValues:
            setsCommonIdentifiers = setsCommonIdentifiers - set(lsIgnoreValues)

        #Get positions of common identifiers in each data set, if the identifier is not unique in a date set just take the first index
        lfFileOneIDIndexes = [fileOneIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]
        lfFileTwoIDIndexes = [fileTwoIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]

        #Convert index list to list of boolean
        lfFileOneElements = [iIndex in lfFileOneIDIndexes for iIndex, sIdentifier in enumerate(fileOneIdentifier)]
        lfFileTwoElements = [iIndex in lfFileTwoIDIndexes for iIndex, sIdentifier in enumerate(fileTwoIdentifier)]

        #Write out file one
        with open(strOutFileOne, 'w') as f:
            f.write(Constants.ENDLINE.join([cDelimiter.join(np.compress(lfFileOneElements,sLine.split(cDelimiter)))
                                           for sLine in filter(None, sContentsOne.split(Constants.ENDLINE))]))

        #Write out file two
        with open(strOutFileTwo, 'w') as f:
            f.write(Constants.ENDLINE.join([cDelimiter.join(np.compress(lfFileTwoElements,sLine.split(cDelimiter)))
                                           for sLine in filter(None, sContentsTwo.split(Constants.ENDLINE))]))

        return True

    #Testing Status: Light happy path testing
    @staticmethod
    def funcCheckRawDataFile(strReadDataFileName, iFirstDataIndex = -1, sLastMetadataName = None, lOccurenceFilter = None, strOutputFileName = "", cDelimiter = Constants.TAB):
        """
        Check the input otu or phlotype abundance table.
        Currently reduces the features that have no occurence.
        Also inserts a NA for blank metadata and a 0 for blank abundance data.
        Gives the option to filter features through an occurence filter (a feature must have a level of abundance in a minimal number of samples to be included).
        Either iFristDataIndex or sLastMetadataName must be given

        :param	strReadDataFileName:	File path of file to read and check.
        :type	String	File path.
        :param	iFirstDataIndex:	First (row) index of data not metadata in the abundance file.
        :type	Integer	Index starting at 0.
        :param	sLastMetadataName:	The ID of the last metadata in the file. Rows of measurements should follow this metadata.
        :param	lOccurenceFilter:	The lowest number of occurences in the lowest number of samples needed for a feature to be kept
        :type	List[2]	List length 2 [lowest abundance (not normalized), lowest number of samples to occur in] (eg. [2.0,2.0])
        :type	String	Matadata ID.
        :param	strOutputFileName:	File path of out put file.
        :type	String	File path.
        :param	cDelimiter:	Character delimiter for reading and writing files.
        :type	Character	Delimiter.
        :return	Output Path:	Output path for written checked file.
        :type	String	File Path.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strReadDataFileName)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strReadDataFileName)
            return False
        if(not ValidateData.funcIsValidStringType(cDelimiter)):
            print "AbundanceTable:checkRawDataFile::Error, Delimiter is not a valid string/char type. Delimiter ="+str(cDelimiter)+"."
            return False
        if (iFirstDataIndex == -1) and (sLastMetadataName == None):
            print "AbundanceTable:checkRawDataFile::Error, either iFirstDataIndex or sLastMetadataNamemust be given."
            return False

        #Get output file and remove if existing
        outputFile = strOutputFileName
        if not strOutputFileName:
            outputFile = os.path.splitext(strReadDataFileName)[0]+Constants.OUTPUT_SUFFIX
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Read input file lines
        #Drop blank lines
        readData = ""
        with open(strReadDataFileName,'r') as f:
            readData = f.read()
        readData = filter(None,readData.split(Constants.ENDLINE))

        #Read the length of each line and make sure there is no jagged data
        #Also hold row count for the metadata
        iLongestLength = len(readData[0].split(cDelimiter))
        iMetadataRow = -1
        if not sLastMetadataName:
            sLastMetadataName = "None"
        for iIndex, strLine in enumerate(readData):
            sLineElements = strLine.split(cDelimiter)
            if sLineElements[0] == sLastMetadataName:
                iMetadataRow = iIndex
            iLongestLength = max(iLongestLength, len(sLineElements))

        #If not already set, set iFirstDataIndex
        if iFirstDataIndex < 0:
            iFirstDataIndex = iMetadataRow + 1

        #Used to substitute . to -
        reSubPeriod = re.compile('\.')

        #File writer
        with open(outputFile,'a') as f:

            #Write metadata
            #Empty data is changed to a default
            #Jagged ends are filled with a default
            for strDataLine in readData[0:iFirstDataIndex]:
                lsLineElements = strDataLine.split(cDelimiter)
                for iindex, sElement in enumerate(lsLineElements):
                    if not sElement.strip():
                        lsLineElements[iindex] = Constants.c_strEmptyDataMetadata
                if len(lsLineElements) < iLongestLength:
                    lsLineElements = lsLineElements + ([Constants.c_strEmptyDataMetadata]*(iLongestLength-len(lsLineElements)))
                f.write(cDelimiter.join(lsLineElements)+Constants.ENDLINE)

            #For each data line in the table
            for line in readData[iFirstDataIndex:]:
                writeToFile = False
                cleanLine = list()
                #Break line into delimited elements
                lineElements = line.split(cDelimiter)

                #Clean feature name
                sCleanFeatureName = reSubPeriod.sub("-",lineElements[0])

                #For each element but the first (taxa name)
                #Element check to see if not == zero
                #If so add to output
                for element in lineElements[1:]:
                    if(element.strip() in string.whitespace):
                        cleanLine.append(Constants.c_strEmptyAbundanceData)
                    #Set abundance of 0 but do not indicate the line should be saved
                    elif(element == "0"):
                        cleanLine.append(element)
                    #If an abundance is found set the line to be saved.
                    else:
                        cleanLine.append(element)
                        writeToFile = True

                #Occurence filtering
                #Removes features that do not have a given level iLowestAbundance in a given amount of samples iLowestSampleOccurence
                if lOccurenceFilter:
                    iLowestAbundance, iLowestSampleOccurence = lOccurenceFilter
                    if iLowestSampleOccurence > sum([1 if float(sEntry) >= iLowestAbundance else 0 for sEntry in cleanLine]):
                        writeToFile = False

                #Write to file
                if writeToFile:    
                    f.write(sCleanFeatureName+cDelimiter+cDelimiter.join(cleanLine)+Constants.ENDLINE)
        return outputFile

    #Testing Status: Light happy path testing
    @staticmethod
    def funcStratifyAbundanceTableByMetadata(strInputFile = None, strDirectory = "", cDelimiter = Constants.TAB, iStratifyByRow = 1, llsGroupings = []):
        """
        Splits an abundance table into multiple abundance tables stratified by the metadata

        :param	strInputFile:	String file path to read in and stratify.
        :type	String	File path.
        :param	strDirectory:	Output directory to write stratified files.
        :type	String	Output directory path.
        :param	cDelimiter:	The delimiter used in the adundance file.
        :type	Character	Delimiter.
        :param	iStratifyByRow:	The row which contains the metadata to use in stratification.
        :type	Integer	Positive integer index.
        :param	llsGroupings:	A list of string lists where each string list holds values that are equal and should be grouped together.
                                So for example, if you wanted to group metadata "1", "2", and "3" seperately but "4" and "5" together you would
                                Give the following [["4","5"]].
                                If you know what "1" and "3" also together you would give [["1","3"],["4","5"]]
        :type	List	List of list of strings
        :return	Boolean:	Indicator of NO error.
        :type	Booelan	False indicates an error.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strInputFile)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, file not valid. File:"+str(strInputFile)
            return False
        if(not ValidateData.funcIsValidStringType(cDelimiter)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Delimiter is not a valid string/char type. Delimiter ="+str(cDelimiter)+"."
            return False
        if(not ValidateData.funcIsValidPositiveInteger(iStratifyByRow, tempZero = True) and (not ValidateData.funcIsValidString(iStratifyByRow))):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Stratify by row is not a positive integer or string keyword. Row ="+str(iStratifyByRow)+"."
            return False

        #Get the base of the file path
        #This is dependent on the given output directory and the prefix of the file name of the input file
        #If no output file is given then the input file directory is used.
        baseFilePath = strDirectory
        lsFilePiecesExt = os.path.splitext(strInputFile)
        if baseFilePath:
            baseFilePath = baseFilePath + os.path.splitext(os.path.split(strInputFile)[1])[0]
        else:
            baseFilePath = lsFilePiecesExt[0]

        #Read in file
        sFileContents = None
        with open(strInputFile,'r') as f:
            sFileContents = f.read()
        sFileContents = filter(None,re.split(Constants.ENDLINE,sFileContents))

        #Collect metadata
        metadataInformation = dict()

        #If the tempStratifyRow is by key word than find the index
        if ValidateData.funcIsValidString(iStratifyByRow):
            for iLineIndex, strLine in enumerate(sFileContents):
                if strLine.split(cDelimiter)[0].strip("\"") == iStratifyByRow:
                    iStratifyByRow = iLineIndex
                    break

        #Stratify by metadata row
        #Split metadata row into metadata entries
        #And put in a dictionary containing {"variable":[1,2,3,4 column index]}
        stratifyByRow = sFileContents[iStratifyByRow].split(cDelimiter)
        for metaDataIndex in xrange(1,len(stratifyByRow)):
            metadata = stratifyByRow[metaDataIndex]
            #Put all wierd categories, none, whitespace, blank space metadata cases into one bin
            if not metadata or metadata in string.whitespace:
                metadata = "Blank"
            #Remove any extraneous formatting
            metadata = metadata.strip(string.whitespace)
            #Store processed metadata with column occurence in dictionary
            if(not metadata in metadataInformation):
                metadataInformation[metadata] = []
            metadataInformation[metadata].append(metaDataIndex)

        #For each of the groupings
        #Use the first value as the primary value which the rest of the values in the list are placed into
        #Go through the dict holding the indices and extend the list for the primary value with the secondary values
        #Then set the secondary value list to empty so that it will be ignored.
        if llsGroupings:
            for lSKeyGroups in llsGroupings:
                if len(lSKeyGroups) > 1:
                    for sGroup in lSKeyGroups[1:]:
                        if sGroup in metadataInformation:
                            metadataInformation[lSKeyGroups[0]].extend(metadataInformation[sGroup])
                            metadataInformation[sGroup] = []

        #Stratify data
        stratifiedAbundanceTables = dict()
        for tableRow in sFileContents:
            row = tableRow.split(cDelimiter) 
            if(len(row)> 1):
                for metadata in metadataInformation:
                    #[0] includes the taxa line
                    columns = metadataInformation[metadata]
                    if columns:
                        columns = [0] + columns
                        lineList = list()
                        for column in columns:
                            lineList.append(row[column])
                        if(not metadata in stratifiedAbundanceTables):
                            stratifiedAbundanceTables[metadata] = list()
                        stratifiedAbundanceTables[metadata].append(cDelimiter.join(lineList))

        #Write to file
        lsFilesWritten = []
        for metadata in stratifiedAbundanceTables:
            sOutputFile = baseFilePath+"-by-"+metadata.strip("\"")+lsFilePiecesExt[1]
            with open(sOutputFile,'w') as f:
                sFileContents = f.write(Constants.ENDLINE.join(stratifiedAbundanceTables[metadata]))
                lsFilesWritten.append(sOutputFile)

        return lsFilesWritten
