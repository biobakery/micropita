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
      if (not npaAbundance == None) and dictMetadata:
        self._npaFeatureAbundance = npaAbundance
        self._dictTableMetadata = dictMetadata
        self._iOriginalFeatureCount = self._npaFeatureAbundance.shape[0]
        self._iOriginalSampleCount = len(self.funcGetSampleNames())
      else:
        print "Abundance or metadata was None, should be atleast an empty object"

    #Allows an Abundance Table to be directly made from a file
    @staticmethod
    def makeFromFile(strInputFile, fIsNormalized, fIsSummed, cDelimiter = Constants.TAB, sMetadataID = None, sLastMetadata = None, cFeatureNameDelimiter="|"):

      #Read in from text file to create the abundance and metadata structures
      lContents = AbundanceTable._textToStructuredArray(strInputFile=strInputFile, cDelimiter=cDelimiter,
                                                        sMetadataID = sMetadataID, sLastMetadata = sLastMetadata)

      #If contents is not a false then set contents to appropriate objects
      if lContents:
        return AbundanceTable(npaAbundance=lContents[0], dictMetadata=lContents[1], strName=strInputFile, 
                              fIsNormalized=fIsNormalized, fIsSummed=fIsSummed, cFileDelimiter=cDelimiter, cFeatureNameDelimiter=cFeatureNameDelimiter)
      return False

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
    def _textToStructuredArray(strInputFile = None, cDelimiter = Constants.TAB, sMetadataID = None, sLastMetadata = None):
    #Validate parameters
        if(not ValidateData.isValidFileName(strInputFile)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, input file not valid. File:",str(strInputFile)])
            return False
        if(not ValidateData.isValidStringType(cDelimiter)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, cDelimiter was invalid. Value =",str(cDelimiter)])
            return False
        if(not ValidateData.isValidString(sMetadataID)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, sMetadataID was invalid. Value =",str(sMetadataID)])
            return False
        if(not ValidateData.isValidString(sLastMetadata)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, sLastMetadata was invalid. Value =",str(sLastMetadata)])
            return False

        #Read in file
        with open(strInputFile,'r') as f:
            contents = f.read()
        f.close()

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
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.dtype.names[1:]
        else:
            return []

    #Returns the metadata name which links to the metadata used
    #as IDs (for instance Sample Ids)
    def funcGetIDMetadataName(self):
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.dtype.names[0]
        else:
            return []

    #Returns a deep copy of the abundance data
    def funcGetAbundanceCopy(self):
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.copy()
        return None

    #Returns a copy of the current abundance table with the abundance of just the given features
    def funcGetFeatureAbundanceTable(self, lsFeatures):
        if (not self._npaFeatureAbundance == None) and lsFeatures:
            #Get a list of boolean indicators that the row is from the features list
            lfFeatureData = [sRowID.strip() in lsFeatures for sRowID in self.funcGetFeatureNames()]

            #compressed version as an Abundance table
            abndTableRet = AbundanceTable(npaAbundance=np.compress(lfFeatureData, self._npaFeatureAbundance, axis = 0),
                           dictMetadata = self.funcGetMetadataCopy(), strName = self.funcGetIDMetadataName(),
                           fIsNormalized = self.funcIsNormalized(), fIsSummed = self.funcIsSummed(),
                           cFileDelimiter = self.funcGetFileDelimiter(), cFeatureNameDelimiter= self.funcGetFeatureDelimiter())

            #Make sure the features are all found
            if ",".join(sorted(lsFeatures)) == ",".join(sorted(abndTableRet.funcGetFeatureNames())):
                return abndTableRet
        return None

    #The delimiter of the feature names (concensus lineages)
    def funcGetFeatureDelimiter(self):
        return self._cFeatureDelimiter

    #Returns the current feature count
    def funcGetFeatureCount(self):
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance.shape[0]
        else:
            return 0

    def funcGetFeatureNames(self):
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[self.funcGetIDMetadataName()]
        return []

    #The delimiter of the file the data was read from and which is also the
    #delimiter which would be used to write the data to a file
    def funcGetFileDelimiter(self):
        return self._cDelimiter

    def funcGetSample(self,sSampleName):
        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[sSampleName].copy()
        return []

    #Returns a specific list of metadata associated with the
    #metadata name (key) that is given
    def funcGetMetadata(self, strMetadataName):
        if self._dictTableMetadata:
            retValue = self._dictTableMetadata.get(strMetadataName,None)
            if retValue:
                retValue = copy.deepcopy(retValue)
            return retValue
        return None

    #Returns a deep copy of the metadata
    def funcGetMetadataCopy(self):
        if self._dictTableMetadata:
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

    #True indicates the metadata exists and are of unique values
    def funcIsPrimaryIdMetadata(self,sMetadataName):
        lMetadata = self.funcGetMetadata(sMetadataName)
        if not lMetadata:
            return False
        return (len(lMetadata) == len(set(lMetadata)))
        

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
        print "AbundanceTable:funcNormalize called"
        if self._fIsSummed:
            return self.funcNormalizeColumnsWithSummedClades()
        else:
            return self.funcNormalizeColumnsBySum()

    #Testing Status: Light happy path testing
    #Normalize the columns (samples) of the abundance table
    #Normalizes as a fraction of the total (number/(sum of all numbers in the column))
    #@returns Normalized structured array or False on error. 
    #All columns are returned, this could be no normalization, all normalization or mixed normalized columns.
    def funcNormalizeColumnsBySum(self):
        print "AbundanceTable:funcNormalizeColumnsBySum called"

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

    #Normalizes a summed Abundance Table
    #If this is called on a dataset which is not summed and not normalized
    #The data will be summed first and then normalized
    #If already normalized, the current normalization is kept
    def funcNormalizeColumnsWithSummedClades(self):
        print "AbundanceTable:funcNormalizeColumnsWithSummedClades called"

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

    #Sums abundance data by clades indicated in the feature name (as consensus lineages)
    def funcSumClades(self):
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

    #Stratifies the AbundanceTable by the given metadata.
    #Will write each stratified abundance table to file
    #if fWriteToFile is True the object will used it's internally stored name as a file to write to
    #if fWriteToFile is a string then it should be a dirctory and end with /. This will rebase the file
    #and store it in a different directory but with an otherwise unchanged name
    #Returns a list of AbundanceTables which are deep copies of the original
    def funcStratifyByMetadata(self,strMetadata,xWriteToFile=False):
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
              if xWriteToFile:
                  if isinstance(xWriteToFile, basestring):
                      self.funcWriteToFile("".join([xWriteToFile,os.path.split(self.getName())[1]]))
                  else:
                      self.funcWriteToFile(self.getName())
              return retlAbundanceTables

            #Given here there are multiple metadata values, continue to stratify
            lsNames = self.funcGetSampleNames()
            #Get index of values to break up
            for value in setValues:
                ffuncGetIDMetadataNameDataIndex = [sData==value for sData in lsMetadata]
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
                if xWriteToFile:
                    if isinstance(xWriteToFile, basestring):
                        objStratifiedAbundanceTable.funcWriteToFile("".join([xWriteToFile,os.path.split(objStratifiedAbundanceTable.funcGetName())[1]]))
                    else:
                        objStratifiedAbundanceTable.funcWriteToFile(objStratifiedAbundanceTable.funcGetName())

                #Append abundance table to returning list
                retlAbundanceTables.append(objStratifiedAbundanceTable)

        return retlAbundanceTables

    #Takes the given data values in one metadata and translates it to values in another
    #metadata of the sample samples holding the values of the first metadata
    #FPrimaryIds, if true the sMetadataFrom are checked for unique values,
    #if the sMetadataFrom has any duplicates the function fails and return false
    def funcTranslateIntoMetadata(self, lsValues, sMetadataFrom, sMetadataTo, fFromPrimaryIds=True):

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
    def funcWriteToFile(self, strOutputFile, cDelimiter=None):
        with open(strOutputFile, 'w') as f:
            #Check delimiter argument
            if not cDelimiter:
                cDelimiter = self._cDelimiter 
            #Write Ids
            f.write(cDelimiter.join([self.funcGetIDMetadataName()]+list(self.funcGetSampleNames()))+Constants.ENDLINE)
            #Write metadata
            f.write(Constants.ENDLINE.join([cDelimiter.join([sMetaKey]+self.funcGetMetadata(sMetaKey)) for sMetaKey in self._dictTableMetadata])+Constants.ENDLINE)
            #Write abundance
            lsOutput = list()
            curAbundance = self._npaFeatureAbundance.tolist()
            for curAbundanceRow in curAbundance:
                lsOutput.append(cDelimiter.join([str(curAbundanceElement) for curAbundanceElement in curAbundanceRow]))
            f.write(Constants.ENDLINE.join(lsOutput))
            f.close()

    #Static methods

    #Testing Status: 1 Happy path test
    #This method will read in two files and abridge both files (saved as new files)
    #to just the samples in common between the two files given a common identifier.
    #Expects the files to have the sample delimiters
    @staticmethod
    def funcPairTables(strFileOne, strFileTwo, strIdentifier, cDelimiter, strOutFileOne, strOutFileTwo, lsIgnoreValues=None):
        #Validate parameters
        if(not ValidateData.isValidFileName(strFileOne)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strFileOne)
            return False
        #Validate parameters
        if(not ValidateData.isValidFileName(strFileTwo)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strFileTwo)
            return False

        #Make file one
        #Read in file
        with open(strFileOne,'r') as f:
            sContentsOne = f.read()
        f.close()

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
        f.close()

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
        lfFileOneElements = [sIdentifier in setsCommonIdentifiers for sIdentifier in fileOneIdentifier]
        lfFileTwoElements = [sIdentifier in setsCommonIdentifiers for sIdentifier in fileTwoIdentifier]

        #Write out file one
        with open(strOutFileOne, 'w') as f:
            f.write(Constants.ENDLINE.join([cDelimiter.join(np.compress(lfFileOneElements,sLine.split(cDelimiter)))
                                           for sLine in filter(None, sContentsOne.split(Constants.ENDLINE))]))
        f.close()

        #Write out file two
        with open(strOutFileTwo, 'w') as f:
            f.write(Constants.ENDLINE.join([cDelimiter.join(np.compress(lfFileTwoElements,sLine.split(cDelimiter)))
                                           for sLine in filter(None, sContentsTwo.split(Constants.ENDLINE))]))
        f.close()

        return True

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
    def funcCheckRawDataFile(strReadDataFileName, iFirstDataIndex= -1, sLastMetadataName = None, strOutputFileName = "", cDelimiter = Constants.TAB):
        #Validate parameters
        if(not ValidateData.isValidFileName(strReadDataFileName)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(strReadDataFileName)
            return False
        if(not ValidateData.isValidStringType(cDelimiter)):
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
        f.close()
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
            sOutputFile = baseFilePath+"-by-"+metadata.strip("\"")+lsFilePiecesExt[1]
            with open(sOutputFile,'w') as f:
                sFileContents = f.write(Constants.ENDLINE.join(stratifiedAbundanceTables[metadata]))
                lsFilesWritten.append(sOutputFile)
            f.close()

        return lsFilesWritten
