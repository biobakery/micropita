#######################################################
# Author: Timothy Tickle
# Description: Class to Run analysis for the microPITA paper
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from Constants import Constants
import csv
import copy
import numpy as np
import os
import re
import scipy.stats
import string
from ValidateData import ValidateData

class AbundanceTable:

    #Constructors..ish

    def __init__(self, npaAbundance, dictMetadata, strName, fIsNormalized, fIsSummed, cFileDelimiter = Constants.TAB, cFeatureNameDelimiter="|"):

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
      self._iCurrentFilterState = 0

      #The feature name delimiter
      self._cFeatureDelimiter = cFeatureNameDelimiter

      #The delimiter from the source file
      self._cDelimiter = cFileDelimiter

      #If contents is not a false then set contents to appropriate objects
      if (not npaAbundance == None) and (not dictMetadata == None):
        self._npaFeatureAbundance = npaAbundance
        self._dictTableMetadata = dictMetadata
        self._iOriginalFeatureCount = self._npaFeatureAbundance.shape[0]
        self._iOriginalSampleCount = len(self.funcGetSampleNames())
      else:
        print "Abundance or metadata was None, should be atleast an empty object"

    #Allows an Abundance Table to be directly made from a file
    @staticmethod
    def makeFromFile(strInputFile, fIsNormalized, fIsSummed, cDelimiter = Constants.TAB, iNameRow = 0, iFirstDataRow = 0, cFeatureNameDelimiter="|"):

      #Read in from text file to create the abundance and metadata structures
      lContents = AbundanceTable._textToStructuredArray(strInputFile=strInputFile, cDelimiter=cDelimiter,
                                              iNameRow=iNameRow, iFirstDataRow=iFirstDataRow)

      #If contents is not a false then set contents to appropriate objects
      if lContents:
        return AbundanceTable(npaAbundance=lContents[0], dictMetadata=lContents[1], strName=strInputFile, 
                              fIsNormalized=fIsNormalized, fIsSummed=fIsSummed, cFileDelimiter=cDelimiter, cFeatureNameDelimiter=cFeatureNameDelimiter)

    #Private Methods

    #Testing Status: Light happy path testing
    #Reads in a file that is Samples (column) by Taxa (rows) into a structured array.
    #@params strInputFile File to read in abundacy data. Tested against Qiime output.
    #@params cDelimiter The delimiter used in the file being read.
    #@params iNameRow The index that is the id row.
    #@params iFirstDataRow The index that is the first data to be read.
    #@returns [A structured array of all data read, A structured array of the metadata].
    #Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
    #Metadata include the taxa ID column name if one is given
    @staticmethod
    def _textToStructuredArray(strInputFile = None, cDelimiter = Constants.TAB, iNameRow = 0, iFirstDataRow = 0):
    #Validate parameters
        if(not ValidateData.isValidFileName(strInputFile)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, input file not valid. File:",str(strInputFile)])
            return False
        if(not ValidateData.isValidStringType(cDelimiter)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, cDelimiter was invalid. Value =",str(cDelimiter)])
            return False
        if(not ValidateData.isValidPositiveInteger(iNameRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, iNameRow was invalid. Value =",str(iNameRow)])
            return False
        if(not ValidateData.isValidPositiveInteger(iFirstDataRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, iFirstDataRow was invalid. Value =",str(iFirstDataRow)])
            return False

        #Read in file
        with open(strInputFile,'r') as f:
            contents = f.read()
        f.close()

        #Turn to lines of the file
        contents = contents.replace("\"","")
        contents = filter(None,contents.split(Constants.ENDLINE))
        namesRow = contents[iNameRow]
        namesRow = namesRow.split(cDelimiter)

        #Get metadata lines
        metadata = dict()
        if((iFirstDataRow-iNameRow)>1):
            for line in contents[iNameRow+1:iFirstDataRow]:
                asmetadataElements = line.split(cDelimiter)
                metadata[asmetadataElements[0]]=asmetadataElements[1:]

        #Check to make sure there is abundance
        if len(contents) <= iFirstDataRow:
            return [np.array([]),metadata]

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
            if(ValidateData.isValidString(taxId)):
                if(ValidateData.isValidList(sampleReads)):
                    if(len(sampleReads)>0):
                        dataMatrix.append(tuple(sampleReads))

        #Now we know the longest taxId we can define the first column holding the tax id
        dataTypeVector = [(namesRow[0],'a'+str(longestTaxId*2))]
        dataTypeVector.extend(incompleteDataTypeVector)
        dataTypeVector = np.dtype(dataTypeVector)

        #Create structured array
        taxData = np.array(dataMatrix,dtype=dataTypeVector)

        return [taxData,metadata]

    #Public methods

    #Returns the samples of the Abundance table
    #Returns an empty list on error or no underlying table
    def funcGetSampleNames(self):
        if not self._npaFeatureAbundance == None:
            return self._npaFeatureAbundance.dtype.names[1:]
        else:
            return []

    #Returns the metadata name which links to the metadata used
    #as IDs (for instance Sample Ids)
    def funcGetIDMetadataName(self):
        if not self._npaFeatureAbundance == None:
            return self._npaFeatureAbundance.dtype.names[0]
        else:
            return []

    #Returns a deep copy of the abundance data
    def funcGetAbundanceCopy(self):
        if not self._npaFeatureAbundance == None:
            return self._npaFeatureAbundance.copy()
        return None

    #The delimiter of the feature names (concensus lineages)
    def funcGetFeatureDelimiter(self):
        return self._cFeatureDelimiter

    #Returns the current feature count
    def funcGetFeatureCount(self):
        if not self._npaFeatureAbundance == None:
            return 0
        else:
            self._npaFeatureAbundance.shape[0]

    #The delimiter of the file the data was read from and which is also the
    #delimiter which would be used to write the data to a file
    def funcGetFileDelimiter(self):
        return self._cDelimiter

    #Returns a specific list of metadata associated with the
    #metadata name (key) that is given
    def funcGetMetadata(self, strMetadataName):
        if not self._dictTableMetadata == None:
            retValue = self._dictTableMetadata.get(strMetadataName,None)
            if retValue:
                retValue = copy.deepcopy(retValue)
            return retValue
        return None

    #Returns a deep copy of the metadata
    def funcGetMetadataCopy(self):
        if not self._dictTableMetadata == None:
            return copy.deepcopy(self._dictTableMetadata)
        return None

    #Returns the name of the object which is the file name that generated it.
    #If the object was generated from an Abundance Table (for instance through stratification)
    #The name is still in the form of a file that could be written to which is informative
    #Of the changes that have occured on the dataset.
    def funcGetName(self):
        return self._strOriginalName

    #True indicates it is normalized
    def funcIsNormalized(self):
        return self._fIsNormalized

    #True indicates the clades are summed
    def funcIsSummed(self):
        return self._fIsSummed

    #Filter abundance by percentile.
    #This method of filtering is inspired from Nicola
    #A feature is removed if it's abundance is not found in the top X percetile a certain percentage of the samples
    def funcFilterAbundanceByPercentile(self, dPercentileCutOff = 95.0, dPercentageAbovePercentile=1.0):
        """
        A feature is removed if it's abundance is not found in the top X percentile a certain percentage of the samples.
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
        self._npaFeatureAbundance =self._npaFeatureAbundance[liKeepIndices,:]
        return True

    #Filter abundance by requiring features to have a minimum sequence occurence in a minimum number of samples
    #Will evaluate greater than or equal to the iMinSequence and iMinSamples
    def funcFilterAbundanceBySequenceOccurence(self, iMinSequence = 2, iMinSamples = 2):
        """
        Filter abundance by requiring features to have a minimum sequence occurence in a minimum number of samples
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
        return True

    #Convenience method which will call which ever normalization is approriate on the data.
    def funcNormalize(self):
        if self._fIsSummed:
            return self.funcNormalizeColumnsByHeirarchy()
        else:
            return self.funcNormalizeColumnsBySum()

    #Testing Status: Light happy path testing
    #Normalize the columns (samples) of the abundance table
    #Normalizes as a fraction of the total (number/(sum of all numbers in the column))
    #@returns Normalized structured array or False on error. 
    #All columns are returned, this could be no normalization, all normalization or mixed normalized columns.
    def funcNormalizeColumnsBySum(self):

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

    #TODO
    #Normalizes a summed Abundance Table
    #If this is called on a dataset which is not summed and not normalized
    #The data will be summed first and then normalized
    #If already normalized, the current normalization is kept
    def funcNormalizeColumnsByHeirarchy(self):
        if self._fIsNormalized:
            print "This table is already normalized, did not perform new normalization request."
            return False

        if not self._fIsSummed:
            print "This table does not have clades summed, this normalization is not appropriate until the clades are summed. The clades are being summed now before normalization."
            self.funcSumClades()

        print "funcNormalizeColumnsByHeirarchy IS NOT IMPLEMENTED"
        #Indicate normalization has occured
        self._fIsNormalized = True

        return False

    #TODO
    #Sums abundance data by clades indicated in the feature name (as consensus lineages)
    def funcSumClades(self):
        if not self._fIsSummed:
          print "funcSumClades IS NOT IMPLEMENTED"
          return False

        #Indicate summation has occured
        self._fIsSummed = True

        return True

    #Happy path tested
    #Expectes a structured array for npData (rows = Taxa/OTU)
    #Expectes the first entry of every row to be an id that is ignored but carried forward
    #Metadata is used to collapse by
#    def funcStratifyDataByMetadata(self,lsMetadata, npData):
#        dictAbundanceBlocks = dict()
#        setValues = set(lsMetadata)
#        lsNames = npData.dtype.names
#        #Get index of values to break up
#        for value in setValues:
#            fDataIndex = [sData==value for sData in lsMetadata]
#            #The true is added to keep the first column which should be the feature id
#            dictAbundanceBlocks[value] = npData[np.compress([True]+fDataIndex,lsNames)]
#        return dictAbundanceBlocks
#
#    def funcStratifyMetadataByMetadata(self, lsMetadata, dictMetadataToStratify):
#        dictMetadataBlocks = dict()
#        setValues = set(lsMetadata)
#        #Get index of values to break up
#        for value in setValues:
#            fDataIndex = [sData==value for sData in lsMetadata]
#            dictBrokenMetadata = dict()
#            for metadataType in dictMetadataToStratify:
#                dictValues = dictMetadataToStratify[metadataType]
#                dictBrokenMetadata[metadataType] = np.compress(fDataIndex,dictValues).tolist()
#            #The true is added to keep the first column which should be the feature id
#            dictMetadataBlocks[value] = dictBrokenMetadata
#        return dictMetadataBlocks

    #Stratifies the AbundanceTable by the given metadata.
    #Will write each stratified abundance table to file
    #if fWriteToFile is True the object will used it's internally stored name as a file to write to
    #if fWriteToFile is a string then it should be a dirctory and end with /. This will rebase the file
    #and store it in a different directory but with an otherwise unchanged name
    #Returns a list of AbundanceTables which are deep copies of the original
    def funcStratifyByMetadata(self,strMetadata,xWriteToFile=False):
        retlAbundanceTables = []

        if (not self._npaFeatureAbundance == None) and (not self._dictTableMetadata == None):
            dictAbundanceBlocks = dict()
            #Get unique metadata values to stratify by
            setValues = set(self._dictTableMetadata.get(strMetadata,[]))

            #If there is only one metadata value then no need to stratify so return the original in the list (and write if needed)
            if len(setValues) == 0:
              return retlAbundanceTables
            elif len(setValues) == 1:
              retlAbundanceTables.append(self)
              #Write to file if need be before returning data
              if xWriteToFile:
                  if isinstance(xWriteToFile, basestring):
                      self.funcWriteToFile("".join([xWriteToFile,os.path.split(self.getName())[1]]))
                  else:
                      self.funcWriteToFile(self.getName())
              return retlAbundanceTables

            #Given here there are multiple metadata values, continue to stratify
            lsNames = self.getSampleNames()
            #Get index of values to break up
            for value in setValues:
                fDataIndex = [sData==value for sData in lsMetadata]
                #Get abundance data for the metadata value
                #The true is added to keep the first column which should be the feature id
                npaStratfiedAbundance = self._npaFeatureAbundance[np.compress([True]+fDataIndex,lsNames)]

                #Get metadata for the metadata value
                dictStratifiedMetadata = dict()
                for metadataType in self._dictTableMetadata:
                    dictValues = dictMetadataToStratify[metadataType]
                    dictStratifiedMetadata[metadataType] = np.compress(fDataIndex,dictValues).tolist()

                #Make abundance table
                #Add abundance table to the list
                strStratifiedName = "".join([os.path.splitext(self._strOriginalName)[0],"-StratBy-",value,".txt"])
                objStratifiedAbundanceTable = AbundanceTable(npaAbundance=npaStratfiedAbundance, dictMetadata=dictStratifiedMetadata, strName=strStratifiedName,
                                              fIsNormalized=self._fIsNormalized, fIsSummed=self._fIsSummed, cFeatureNameDelimiter=self._cFeatureDelimiter)

                #Write to file if needed
                if xWriteToFile:
                    if isinstance(xWriteToFile, basestring):
                        objStratifiedAbundanceTable.funcWriteToFile("".join([xWriteToFile,os.path.split(objStratifiedAbundanceTable.getName())[1]]))
                    else:
                        objStratifiedAbundanceTable.funcWriteToFile(objStratifiedAbundanceTable.getName())

                #Append abundance table to returning list
                retlAbundanceTables.append(objStratifiedAbundanceTable)

        return retlAbundanceTables


    #Returns a numpy array of the current Abundance Table.
    #Removes the first ID head column and the numpy array is
    #Made of lists, not tuples.
    def funcToArray(self):
        if not self._npaFeatureAbundance == None:
            return np.array([list(tplRow)[1:] for tplRow in self._npaFeatureAbundance],'float')
        return None

    #Writes the AbundanceTable to a file strOutputFile.
    #Will rewrite over a file as needed.
    #Will use the cDelimiter to delimit columns if provided.
    def funcWriteToFile(self, strOutputFile, cDelimiter=self._cDelimiter):
        with open(strOutputFile, 'w') as f:
            #Write Ids
            f.write(cDelimiter.join([self.funcGetIDMetadataName()]+self.getSampleNames())+Constants.ENDLINE)
            #Write metadata
            f.write(Constants.ENDLINE.join([cDelimiter.join([sMetaKey]+dictCurMetadata[sMetaKey]) for sMetaKey in self._dictMetadata]))
            #Write abundance
            lsOutput = list()
            curAbundance = self._npaFeatureAbundance.tolist()
            for curAbundanceRow in curAbundance:
                lsOutput.append(cDelimiter.join([str(curAbundanceElement) for curAbundanceElement in curAbundanceRow]))
                f.write(Constants.ENDLINE.join(lsOutput))
            f.close()

    #Static methods

    #Testing Status: Light happy path testing
    #Splits an abundance table into multiple abundance tables stratified by the metadata
    #@params strInputFile String file path to read in and stratify
    #These assumptions are based on the current default formating of an abundance file
    #@cDelimiter The delimiter used in the adundance file
    #@iStratifyByRow The row which contains the metadata to use in stratification.
    #Can be the index of the row starting with 0 or a keyword found in the first column of the row.
    #llsGroupings is a list of string lists where each string list holds values that are equal and should be grouped together.
    #So for example, if you wanted to group metadata "1", "2", and "3" seperately but "4" and "5" together you would
    #Give the following [["4","5"]].
    #If you know what "1" and "3" also together you would give [["1","3"],["4","5"]]
    #@return Returns boolean true or false. True indicating completion without detected errors
    @staticmethod
    def funcStratifyAbundanceTableByMetadata(strInputFile = None, strDirectory = "", cDelimiter = Constants.TAB, iStratifyByRow = 1, llsGroupings = []):

        #Validate parameters
        if(not ValidateData.isValidFileName(strInputFile)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, file not valid. File:"+str(strInputFile)
            return False
        if(not ValidateData.isValidStringType(cDelimiter)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Delimiter is not a valid string/char type. Delimiter ="+str(cDelimiter)+"."
            return False
        if(not ValidateData.isValidPositiveInteger(iStratifyByRow, tempZero = True) and (not ValidateData.isValidString(iStratifyByRow))):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Stratify by row is not a positive integer or string keyword. Row ="+str(iStratifyByRow)+"."
            return False

        #Get the base of the file path
        #This is dependent on the given output directory and the prefix of the file name of the input file
        #If no output file is given then the input fiel directory is used.
        baseFilePath = strDirectory
        if baseFilePath:
            baseFilePath = baseFilePath + os.path.splitext(os.path.split(strInputFile)[1])[0]
        else:
            baseFilePath = os.path.splitext(strInputFile)[0]

        #Read in file
        sFileContents = None
        with open(strInputFile,'r') as f:
            sFileContents = f.read()
        f.close()
        sFileContents = filter(None,re.split(Constants.ENDLINE,sFileContents))

        #Collect metadata
        metadataInformation = dict()

        #If the tempStratifyRow is by key word than find the index
        if ValidateData.isValidString(iStratifyByRow):
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
            sOutputFile = baseFilePath+"-by-"+metadata.strip("\"")+".txt"
            with open(sOutputFile,'w') as f:
                sFileContents = f.write(Constants.ENDLINE.join(stratifiedAbundanceTables[metadata]))
                lsFilesWritten.append(sOutputFile)
            f.close()

        return lsFilesWritten

    #Testing Status: Light happy path testing
    #Check the input otu or phlotype abundance table
    #Currently reduces the taxa that have no occurence
    #Also inserts a NA for blank metadta and a 0 for blank abundance data
    #@params strReadDataFileName String. Data file name
    #@params strOutputFileName String. The file path to save the checked file as.
    #If left empty (defualt) a path will be created from the input file name
    #Placing the file in the same directory as the current file.
    #@params cDelimiter Character. Delimiter for the data
    #@return Return string file path
    @staticmethod
    def funcCheckRawDataFile(strReadDataFileName, iFirstDataIndex, strOutputFileName = "", cDelimiter = Constants.TAB):
        #Validate parameters
        if(not ValidateData.isValidFileName(strReadDataFileName)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strReadDataFileName)
            return False
        if(not ValidateData.isValidStringType(cDelimiter)):
            print "AbundanceTable:checkRawDataFile::Error, Delimiter is not a valid string/char type. Delimiter ="+str(cDelimiter)+"."
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
        f.close()
        readData = filter(None,readData.split(Constants.ENDLINE))

        #Read the length of each line and make sure there is no jagged data
        iLongestLength = len(readData[0].split(cDelimiter))
        for strLine in readData:
            iLongestLength = max(iLongestLength, len(strLine.split(cDelimiter)))

        #File writer
        with open(strOutputFileName,'a') as f:

            #Write metadata
            #Empty data is changed to a default
            #Jagged ends are filled with a default
            for strDataLine in readData[0:iFirstDataIndex]:
                lsLineElements = strDataLine.split(cDelimiter)
                for iindex, sElement in enumerate(lsLineElements):
                    if not sElement.strip():
                        lsLineElements[iindex] = Constants.c_strEmptyDataMetadata
                if len(lsLineElements) < iLongestLength:
                    lsLineElements = lsLineElements + (["NA"]*(iLongestLength-len(lsLineElements)))
                f.write(cDelimiter.join(lsLineElements)+Constants.ENDLINE)

            #For each data line in the table
            for line in readData[iFirstDataIndex:]:
                writeToFile = False
                cleanLine = list()
                #Break line into delimited elements
                lineElements = line.split(cDelimiter)

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
                #Write to file
                if(writeToFile):    
                    f.write(lineElements[0]+cDelimiter+cDelimiter.join(cleanLine)+Constants.ENDLINE)
        f.close()
        return outputFile


