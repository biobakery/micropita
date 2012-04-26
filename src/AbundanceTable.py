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
from FileIO import FileIO
import numpy as np
import os
import scipy.stats
from ValidateData import ValidateData

class AbundanceTable:

    #Testing Status: Light happy path testing
    #Check the input otu or phlotype abundance table
    #Currently reduces the taxa that have no occurence
    #@params tempReadDataFileName String. Data file name
    #@params tempDelimiter Character. Delimiter for the data
    #@return Return string file path
    @staticmethod
    def checkRawDataFile(tempReadDataFileName, tempDelimiter = Constants.TAB, tempFirstDataIndex = 2):
        #Validate parameters
        if(not ValidateData.isValidFileName(tempReadDataFileName)):
            print "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+str(tempReadDataFileName)
            return False
        if(not ValidateData.isValidStringType(tempDelimiter)):
            print "AbundanceTable:checkRawDataFile::Error, Delimiter is not a valid string/char type. Delimiter ="+str(tempDelimiter)+"."
            return False

        #Get output file and remove if existing
        outputFile = os.path.splitext(tempReadDataFileName)[0]+Constants.OUTPUT_SUFFIX
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Read input file lines
        readData = ""
        with open(tempReadDataFileName,'r') as f:
            readData = f.read()
            f.close()
        readData = filter(None,readData.split(Constants.ENDLINE))

        #File writer
        writeData = FileIO(outputFile,False,True,True)

        #Write metadata without change
        for strDataLine in readData[0:tempFirstDataIndex]:
            writeData.writeToFile("".join([strDataLine,Constants.ENDLINE]))

        #For each data line in the table
        for line in readData[tempFirstDataIndex:]:
            writeToFile = False
            cleanLine = list()
            #Break line into delimited elements
            lineElements = filter(None,line.split(tempDelimiter))

            #For each element but the first (taxa name)
            #Element check to see if not == zero
            #If so add to output
            for element in lineElements[1:]:
                if(element.rstrip() == ''):
                    cleanLine.append("0")
                elif(element == "0"):
                    cleanLine.append("0")
                #If an abundance is found set the line to be saved.
                else:
                    cleanLine.append(element)
                    writeToFile = True
            #Write to file
            if(writeToFile):    
                writeData.writeToFile(lineElements[0]+tempDelimiter+tempDelimiter.join(cleanLine)+Constants.ENDLINE)
        writeData.close()
        return outputFile

    #Filter an abundance file by abundance
    def filterByAbundance(self, npaAbundance, dPercentileCutOff = 0.95, dPercentageAbovePercentile=0.1):
      """
      Filter by abundance. Specifically this version requires elements of
      the tree to have a certain percentage of a certain percentile in samples.

      """

      print(npaAbundance)
      if(npaAbundance is not None):
          #Sample names
          lsSampleNames = npaAbundance.dtype.names[1:]
          print("lsSampleNames")
          print(lsSampleNames)
          #List of indexes of taxa to keep
          liTaxa = []
          #Hold the cuttoff score (threshold) for the percentile of interest {SampleName(string):score(double)}
          dictPercentiles = dict()
          for index in xrange(0,len(lsSampleNames)):
              print("npaAbundance[lsSampleNames[index]]")
              print(npaAbundance[lsSampleNames[index]])
              print("dPercentileCutOff")
              print(dPercentileCutOff)
              print("scipy.stats.scoreatpercentile(npaAbundance[lsSampleNames[index]],dPercentileCutOff)")
              print(scipy.stats.scoreatpercentile(npaAbundance[lsSampleNames[index]],dPercentileCutOff))
              dScore = scipy.stats.scoreatpercentile(npaAbundance[lsSampleNames[index]],dPercentileCutOff)
              dictPercentiles[lsSampleNames[index]] = dScore

          print("dictPercentiles")
          print(dictPercentiles)
          #Sample count (Ignore sample id [position 0] which is not a name)
          dSampleCount = float(len(lsSampleNames[1:]))

          #Check each taxa
          iTaxaIndex = 0
          for rowTaxaData in npaAbundance:
              sCurTaxaName = rowTaxaData[0]
              print("sCurTaxaName")
              print(sCurTaxaName)
              dCountAbovePercentile = 0.0
              ldAbundanceMeasures = list(rowTaxaData)[1:]
              print("ldAbundanceMeasures")
              print(ldAbundanceMeasures)
              #Check to see if the abundance score meets the threshold and count if it does
              for iScoreIndex in xrange(0,len(ldAbundanceMeasures)):
                  if(ldAbundanceMeasures[iScoreIndex] >= dictPercentiles[lsSampleNames[iScoreIndex]]):
                      dCountAbovePercentile = dCountAbovePercentile + 1.0
              dPercentOverPercentile = dCountAbovePercentile / dSampleCount
              if(dPercentOverPercentile >= (dPercentageAbovePercentile/100.0)):
                  liTaxa.append(iTaxaIndex)
              #Increment taxa
              iTaxaIndex = iTaxaIndex + 1
          print("liTaxa")
          print(liTaxa)

    #Testing Status: Light happy path testing
    #Normalize the given columns of a structured array given the structured array and a list of the column names
    #Normalizes as a fraction of the total (number/(sum of all numbers in the column))
    #@params tempStructuredArray Structured array where columns will be normalized
    #@params tempColumnNames List of column names to be normalized
    #@returns Normalized structured array or False on error. 
    #All columns are returned, this could be no normalization, all normalization or mixed normalized columns.
    def normalizeColumns(self, tempStructuredArray = None, tempColumns = None):
        #Validate parameters
        if(not ValidateData.isValidStructuredArray(tempStructuredArray)):
            print "AbundanceTable:normalizeColumns::Error, structured array is not valid. Structured array:"+str(tempStructuredArray)
            return False
        if(not ValidateData.isValidList(tempColumns)):
            print "AbundanceTable:normalizeColumns::Error, Column name not a valid list. List ="+str(tempColumns)+"."
            return False

        for columnName in tempColumns:
            column = tempStructuredArray[columnName]
            columnTotal = sum(column)
            if(columnTotal > 0.0):
                column = column/columnTotal
            tempStructuredArray[columnName] = column
        return tempStructuredArray

    #Testing Status: Light happy path testing
    #Splits an abundance table into multiple abundance tables stratified by the metadata
    #@params tempInputFile String file path to read in and stratify
    #These assumptions are based on the current default formating of an abundance file
    #@tempDelimiter The delimiter used in the adundance file
    #@tempStratifyByRow The row which contains the metadata to use in stratification.
    #Can be the index of the row starting with 0 or a keyword found in the first column of the row.
    #@return Returns boolean true or false. True indicating completion without detected errors
    def stratifyAbundanceTableByMetadata(self, tempInputFile = None, tempDelimiter = Constants.TAB, tempStratifyByRow = 1):
        #Validate parameters
        if(not ValidateData.isValidFileName(tempInputFile)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, file not valid. File:"+str(tempInputFile)
            return False
        if(not ValidateData.isValidStringType(tempDelimiter)):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Delimiter is not a valid string/char type. Delimiter ="+str(tempDelimiter)+"."
            return False
        if(not ValidateData.isValidPositiveInteger(tempStratifyByRow, tempZero = True) and (not ValidateData.isValidString(tempStratifyByRow))):
            print "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Stratify by row is not a positive integer or string keyword. Row ="+str(tempStratifyByRow)+"."
            return False

        #Get the base of the file path
        baseFilePath = os.path.splitext(tempInputFile)[0]

        #Read in file
        inputRead = FileIO(tempInputFile,True, False, False)
        contents = inputRead.readFullFile()
        inputRead.close()
        contents = filter(None,contents.split(Constants.ENDLINE))

        #Collect metadata
        metadataInformation = dict()

        #If the tempStratifyRow is by key word than find the index
        if ValidateData.isValidString(tempStratifyByRow):
            for iLineIndex, strLine in enumerate(contents):
                if strLine.split(tempDelimiter)[0].strip("\"") == tempStratifyByRow:
                    tempStratifyByRow = iLineIndex
                    break

        #Stratify by metadata row
        stratifyByRow = filter(None,contents[tempStratifyByRow].split(tempDelimiter))
        for metaDataIndex in xrange(1,len(stratifyByRow)):
            metadata = stratifyByRow[metaDataIndex]
            if(not metadata in metadataInformation):
                #Include the taxa line
                metadataInformation[metadata] = [0]
            metadataInformation[metadata].append(metaDataIndex)

        #Stratify data
        stratifiedAbundanceTables = dict()
        for tableRow in contents:
            row = filter(None,tableRow.split(tempDelimiter))
            if(len(row)> 1):
                for metadata in metadataInformation:
                    columns = metadataInformation[metadata]
                    lineList = list()
                    for column in columns:
                        lineList.append(row[column])
                    if(not metadata in stratifiedAbundanceTables):
                        stratifiedAbundanceTables[metadata] = list()
                    stratifiedAbundanceTables[metadata].append(tempDelimiter.join(lineList))

        #Write to file
        for metadata in stratifiedAbundanceTables:
            write = FileIO(baseFilePath+"-by-"+metadata.strip("\"")+".txt",False,True,False)
            write.writeToFile(Constants.ENDLINE.join(stratifiedAbundanceTables[metadata]))
            write.close()

        return True

    #Testing Status: Light happy path testing
    #Reads in a file that is Samples (column) by Taxa (rows) into a structured array.
    #@params tempInputFile File to read in abundacy data. Tested against Qiime output.
    #@params tempDelimiter The delimiter used in the file being read.
    #@params tempNameRow The index that is the id row.
    #@params tempFirstDataRow The index that is the first data to be read.
    #@params tempNormalize If the abundancy data needs to be normalized. If so, normalized by column (sample) columnElement/sum(column).
    #@returns [A structured array of all data read, A structured array of the metadata].
    #Samples (column) by Taxa (rows) with the taxa id row included as the column index=0
    #Metadata include the taxa ID column name if one is given
    def textToStructuredArray(self, tempInputFile = None, tempDelimiter = Constants.TAB, tempNameRow = 0, tempFirstDataRow = 0, tempNormalize = True):
    #Validate parameters
        if(not ValidateData.isValidFileName(tempInputFile)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, input file not valid. File:",str(tempInputFile)])
            return False
        if(not ValidateData.isValidStringType(tempDelimiter)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempDelimiter was invalid. Value =",str(tempDelimiter)])
            return False
        if(not ValidateData.isValidPositiveInteger(tempNameRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempNameRow was invalid. Value =",str(tempNameRow)])
            return False
        if(not ValidateData.isValidPositiveInteger(tempFirstDataRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempFirstDataRow was invalid. Value =",str(tempFirstDataRow)])
            return False
        if(not ValidateData.isValidBoolean(tempNormalize)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempNormalize was invalid. Value =",str(tempNormalize)])
            return False

        #Read in file
        inputRead = FileIO(tempInputFile,True, False, False)
        contents = inputRead.readFullFile()
        inputRead.close()

        #Turn to lines of the file
        contents = contents.replace("\"","")
        contents = filter(None,contents.split(Constants.ENDLINE))
        namesRow = contents[tempNameRow]
        namesRow = filter(None,namesRow.split(tempDelimiter))

        #Get metadata lines
        metadata = dict()
        if((tempFirstDataRow-tempNameRow)>1):
            for line in contents[tempNameRow+1:tempFirstDataRow]:
                asmetadataElements = filter(None,line.split(tempDelimiter))
                metadata[asmetadataElements[0]]=asmetadataElements[1:]

        #Build data type object (data name,data type)
        incompleteDataTypeVector = []
        for sampleColumn in xrange(1,len(namesRow)):
            incompleteDataTypeVector.append((namesRow[sampleColumn],'f8'))

        #Extract data
        #Create the tuple with the first data row
        rowData = contents[tempFirstDataRow]
        rowData = filter(None,rowData.split(tempDelimiter))
        taxId = rowData[0]
        longestTaxId = len(taxId)
        sampleReads = rowData[1:]
        tempSampleReads = [taxId]
        for reads in sampleReads:
            tempSampleReads.append(float(reads))
        sampleReads = tempSampleReads
        dataMatrix = [tuple(sampleReads)]

        #Add the rest of the rows
        for rowData in contents[tempFirstDataRow+1:]:
            rowData = filter(None,rowData.split(tempDelimiter))
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

        #Normalize the sample reads by the total sample reads
        #Normalize to 1
        if(tempNormalize):
            for columnName in taxData.dtype.names[1:]:
                column = taxData[columnName]
                columnTotal = sum(column)
                if(columnTotal > 0.0):
                    column = column/columnTotal
                taxData[columnName] = column

        #Successful return
        return [taxData,metadata]

    #Not tested
    #Reads in a file that is Samples (column) by Taxa (rows) into an numpy array of floats.
    #@params tempInputFile File to read in abundacy data. Tested against Qiime output.
    #@params tempDelimiter The delimiter used in the file being read.
    #@params tempIgnoreIdColumn boolean indicator that the first column is an id column which will be ignored.
    #@params tempNameRow The index that is the id row.
    #@params tempFirstDataRow The index that is the first data to be read.
    #@params tempNormalize If the abundancy data needs to be normalized. If so, normalized by column (sample) columnElement/sum(column).
    #@returns [A structured array of all data read, A structured array of the metadata].
    #Samples (column) by Taxa (rows) without the column at index=0 if tempIgnoreIdColumn=True
    #Metadata include the taxa ID column name if one is given
    def textToArray(self, tempInputFile = None, tempDelimiter = Constants.TAB, tempIgnoreIdColumn=True, tempNameRow = 0, tempFirstDataRow = 0, tempNormalize = True):
    #Validate parameters
        if(not ValidateData.isValidFileName(tempInputFile)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, input file not valid. File:",str(outputFile)])
            return False
        if(not ValidateData.isValidStringType(tempDelimiter)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempDelimiter was invalid. Value =",str(tempDelimiter)])
            return False
        if(not ValidateData.isValidPositiveInteger(tempNameRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempNameRow was invalid. Value =",str(tempNameRow)])
            return False
        if(not ValidateData.isValidPositiveInteger(tempFirstDataRow, tempZero = True)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempFirstDataRow was invalid. Value =",str(tempFirstDataRow)])
            return False
        if(not ValidateData.isValidBoolean(tempNormalize)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempNormalize was invalid. Value =",str(tempNormalize)])
            return False
        if(not ValidateData.isValidBoolean(tempIgnoreIdColumn)):
            print "".join(["AbundanceTable:textToStructuredArray::Error, tempIgnoreIdColumn was invalid. Value =",str(tempIgnoreIdColumn)])
            return False

        #Read in file
        inputRead = FileIO(tempInputFile,True, False, False)
        contents = inputRead.readFullFile()
        inputRead.close()

        #Turn to lines of the file
        contents = contents.replace("\"","")
        contents = filter(None,contents.split(Constants.ENDLINE))
        namesRow = contents[tempNameRow]
        namesRow = filter(None,namesRow.split(tempDelimiter))

        #Get metadata lines
        metadata = dict()
        if((tempFirstDataRow-tempNameRow)>1):
            for line in contents[tempNameRow+1:tempFirstDataRow]:
                asMetadataElement=filter(None,line.split(tempDelimiter))
                metadata[asMetadataElement[0]]=asMetadataElement[1:]

        #Build data type object (data name,data type)
        dataTypeVector = []
        #If ignoring the first column, do not include it
        if(tempIgnoreIdColumn==True):
            for sampleColumn in xrange(1,len(namesRow)):
                dataTypeVector.append((namesRow[sampleColumn],'f8'))
        else:
            for sampleColumn in xrange(0,len(namesRow)):
                dataTypeVector.append((namesRow[sampleColumn],'f8'))

        #Extract data
        #Create the list with the first data row
        rowData = contents[tempFirstDataRow]
        rowData = filter(None,rowData.split(tempDelimiter))
        #Ignore first column if needed
        if(tempIgnoreIdColumn==True):
            rowData = rowData[1:]
        #Change to floats
        sampleReads = list()
        for reads in rowData:
            sampleReads.append(float(reads))
        dataMatrix = [sampleReads]

        #Add the rest of the rows
        for rowData in contents[tempFirstDataRow+1:]:
            rowData = filter(None,rowData.split(tempDelimiter))
            #Ignore first column if needed
            if(tempIgnoreIdColumn==True):
                rowData = rowData[1:]
            #Change to floats
            sampleReads = list()
            for reads in rowData:
                sampleReads.append(float(reads))
            #Validate data before adding to matrix
            if(len(sampleReads)>0):
                dataMatrix.append(sampleReads)

        #Set the data types
        dataTypeVector = np.dtype(dataTypeVector)

        #Create structured array
        taxaData = np.array(dataMatrix,'float')

        #Normalize the sample reads by the total sample reads
        #Normalize to 1
        if(tempNormalize):
            for rowIndex in xrange(0,len(taxaData[0])):
                column = taxaData[:,rowIndex]
                columnTotal = sum(column)
                if(columnTotal > 0.0):
                    column = column/columnTotal
                taxaData[:,rowIndex] = column

        #Successful return
        return [taxaData,metadata]

    #Testing Status: Light happy path testing
    #Transposes a matrix.
    #Removes the first column before transposing if tempRemoveAdornments = True
    #@params tempMatrix Strucutred array to transpose
    #@params tempRemoveAdornments Remove first column before transposing
    #@return Returns Transposed structured array (list of lists, list=previous columns) with potentially the first column removed.
    def transposeDataMatrix(self,tempMatrix, tempRemoveAdornments=False):
        #Validate parameters
        if(not ValidateData.isValidStructuredArray(tempMatrix)):
            print "".join(["AbundanceTable:transposeDataMatrix::Error, transposeDataMatrix was an invalid structured array. Value =",str(tempMatrix)])
            return False
        if(not ValidateData.isValidBoolean(tempRemoveAdornments)):
            print "".join(["AbundanceTable:transposeDataMatrix::Error, tempRemoveAdornments was an invalid boolean. Value =",str(tempRemoveAdornments)])
            return False

        #Change to samples x taxa as is needed for the compute method below
        #Also remove the first row which is taxa identification
        startColumnElement = 0
        if(tempRemoveAdornments == True):
            startColumnElement = 1
        conversionMatrix = list()
        for row in tempMatrix:
            conversionMatrix.append(list(row)[startColumnElement:])
        tempMatrix = None
        conversionMatrix = np.array(conversionMatrix)
        conversionMatrix = conversionMatrix.transpose()
        return conversionMatrix
