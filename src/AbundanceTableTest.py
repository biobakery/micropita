#######################################################
# Author: Timothy Tickle
# Description: Class to test the RunMicroPITA class
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Import libraries
import AbundanceTable
import Constants
import FileIO
import os
import unittest
import Utility_File

##
#Tests the Blog object
class AbundanceTableTest(unittest.TestCase):

    def testCheckRawDataFileForGoodCase(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking.txt"
        outputFile = Utility_File.Utility_File.getFileNamePrefix(inputFile)+Constants.Constants.OUTPUT_SUFFIX
        delimiter = Constants.Constants.TAB

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\"\t700098986\t700098984\t700098982\t700098980\t700098988\t700037470\t700037472\t700037474\t700037476\t700037478\n\"STSite\"\t\"L_Antecubital_fossa\"\t\"R_Retroauricular_crease\"\t\"L_Retroauricular_crease\"\t\"Subgingival_plaque\"\t\"R_Antecubital_fossa\"\t\"L_Retroauricular_crease\"\t\"R_Retroauricular_crease\"\t\"L_Antecubital_fossa\"\t\"R_Antecubital_fossa\"\t\"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\"\t1\t0\t0\t12\t0\t6\t0\t2\t1\t0\n\"Bacteria|unclassified|4904\"\t0\t10\t0\t43\t6\t0\t23\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\"\t3\t0\t0\t29\t0\t45\t0\t1\t1\t0\n\"Bacteria|3417\"\t0\t45\t0\t34\t3\t0\t0\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\"\t5\t0\t0\t2\t0\t6\t0\t1\t1\t0\n"

        #Call method
        AbundanceTable.AbundanceTable.checkRawDataFile(tempReadDataFileName=inputFile, tempDelimiter=delimiter)

        #Get answer
        readFile = FileIO.FileIO(outputFile,True,False,False)
        result = readFile.readFullFile()
        readFile.close()

        #Remove output file after the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testCheckRawDataFileForGoodCaseDelimterSpace(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking_Space.txt"
        outputFile = Utility_File.Utility_File.getFileNamePrefix(inputFile)+Constants.Constants.OUTPUT_SUFFIX
        delimiter = Constants.Constants.WHITE_SPACE

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\" 700098986 700098984 700098982 700098980 700098988 700037470 700037472 700037474 700037476 700037478\n\"STSite\" \"L_Antecubital_fossa\" \"R_Retroauricular_crease\" \"L_Retroauricular_crease\" \"Subgingival_plaque\" \"R_Antecubital_fossa\" \"L_Retroauricular_crease\" \"R_Retroauricular_crease\" \"L_Antecubital_fossa\" \"R_Antecubital_fossa\" \"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\" 1 0 0 12 0 6 0 2 1 0\n\"Bacteria|unclassified|4904\" 0 10 0 43 6 0 23 0 1 0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\" 3 0 0 29 0 45 0 1 1 0\n\"Bacteria|3417\" 0 45 0 34 3 0 0 0 1 0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\" 5 0 0 2 0 6 0 1 1 0\n"

        #Call method
        AbundanceTable.AbundanceTable.checkRawDataFile(tempReadDataFileName=inputFile, tempDelimiter=delimiter)

        #Get answer
        readFile = FileIO.FileIO(outputFile,True,False,False)
        result = readFile.readFullFile()
        readFile.close()

        #Remove output file after the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testNormalizeColumnsForGoodCaseNoNormalize(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        columnNames = []

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0)\n ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.normalizeColumns(tempStructuredArray=result[0], tempColumns=columnNames)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testNormalizeColumnsForGoodCaseNormalize1(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        columnNames = ["700098986"]

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0)\n ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.normalizeColumns(tempStructuredArray=result[0], tempColumns=columnNames)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testNormalizeColumnsForGoodCaseNormalize14(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        columnNames = ["700098986","700098980"]

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 0.35833333333333334, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0)\n ('Bacteria|3417', 0.0, 45.0, 0.0, 0.2833333333333333, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.normalizeColumns(tempStructuredArray=result[0], tempColumns=columnNames)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testNormalizeColumnsForGoodCaseNormalizeAll(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        columnNames = ["700098986","700098984","700098982","700098980","700098988","700037470","700037472","700037474","700037476","700037478"]

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)\n ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.normalizeColumns(tempStructuredArray=result[0], tempColumns=columnNames)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTextToStructuredArrayForGoodCaseNoNormalize(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        microPITA = AbundanceTable.AbundanceTable()

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), [['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']]]"

        #Call method
        result = microPITA.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTextToStructuredArrayForGoodCaseNormalize(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = True
        microPITA = AbundanceTable.AbundanceTable()

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0),\n       ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), [['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']]]"

        #Call method
        result = microPITA.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=",str(answer),".\nReceived=",str(result),"."]))

    def testTextToStructuredArrayForGoodCaseSpaceDelimiter(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged_Space.txt"
        delimiter = Constants.Constants.WHITE_SPACE
        nameRow = 0
        firstDataRow = 2
        normalize = False
        microPITA = AbundanceTable.AbundanceTable()

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), [['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']]]"

        #Call method
        result = microPITA.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testStratifyAbundanceTableByMetadataForGoodCase(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        stratifyRow = 1
        table = AbundanceTable.AbundanceTable()

        #Correct Answer is blank indicating no error
        answer = ""

        #Should generate the following files
        anteriorNaresFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"
        anteriorNaresFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"
        lAntecubitalFossaFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"
        lAntecubitalFossaFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"
        lRetroauricularCreaseFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"
        lRetroauricularCreaseFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"
        rAntecubitalFossaFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"
        rAntecubitalFossaFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"
        rRetroauricularCreaseFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"
        rRetroauricularCreaseFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"
        subgingivalPlaqueFileName = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"
        subgingivalPlaqueFileNameAnswer = "./testData/CorrectTestingResults/hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"

        #Check file creation
        error = ""
        if(os.path.exists(anteriorNaresFileName)):
            read = FileIO.FileIO(anteriorNaresFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(anteriorNaresFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+anteriorNaresFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+anteriorNaresFileName

        if(os.path.exists(lAntecubitalFossaFileName)):
            read = FileIO.FileIO(lAntecubitalFossaFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(lAntecubitalFossaFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+lAntecubitalFossaFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+lAntecubitalFossaFileName

        if(os.path.exists(lRetroauricularCreaseFileName)):
            read = FileIO.FileIO(lRetroauricularCreaseFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(lRetroauricularCreaseFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+lRetroauricularCreaseFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+lRetroauricularCreaseFileName

        if(os.path.exists(rAntecubitalFossaFileName)):
            read = FileIO.FileIO(rAntecubitalFossaFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(rAntecubitalFossaFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+rAntecubitalFossaFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+rAntecubitalFossaFileName

        if(os.path.exists(rRetroauricularCreaseFileName)):
            read = FileIO.FileIO(rRetroauricularCreaseFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(rRetroauricularCreaseFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+rRetroauricularCreaseFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+rRetroauricularCreaseFileName

        if(os.path.exists(subgingivalPlaqueFileName)):
            read = FileIO.FileIO(subgingivalPlaqueFileName,True,False,False)
            contents = read.readFullFile()
            read.close()
            read = FileIO.FileIO(subgingivalPlaqueFileNameAnswer,True,False,False)
            contentsAnswer = read.readFullFile()
            read.close()
            if(not contents == contentsAnswer):
                error = error + "\nFile: "+subgingivalPlaqueFileName+"\nExpected:"+contentsAnswer+".\nReceived:"+contents+"."
        else:
            error = error + "\nFile count not be found. Path:"+subgingivalPlaqueFileName

        #Call method
        table.stratifyAbundanceTableByMetadata(tempInputFile = inputFile, tempDelimiter = delimiter, tempStratifyByRow = stratifyRow)
        result = error

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))


    def testTransposeDataMatrixForGoodCase(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        removeAdornment = False

        #Correct Answer
        answer = "[[ 'Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72'\n  'Bacteria|unclassified|4904'\n  'Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361'\n  'Bacteria|3417'\n  'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368']\n ['1.0' '0.0' '3.0' '0.0' '5.0']\n ['0.0' '10.0' '0.0' '45.0' '0.0']\n ['0.0' '0.0' '0.0' '0.0' '0.0']\n ['12.0' '43.0' '29.0' '34.0' '2.0']\n ['0.0' '6.0' '0.0' '3.0' '0.0']\n ['6.0' '0.0' '45.0' '0.0' '6.0']\n ['0.0' '23.0' '0.0' '0.0' '0.0']\n ['2.0' '0.0' '1.0' '0.0' '1.0']\n ['1.0' '1.0' '1.0' '1.0' '1.0']\n ['0.0' '0.0' '0.0' '0.0' '0.0']]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.transposeDataMatrix(tempMatrix=result[0], tempRemoveAdornments=removeAdornment)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTransposeDataMatrixForGoodCaseRemoveAdornments(self):
        
        #Inputs
        inputFile = "./testData/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"
        delimiter = Constants.Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.AbundanceTable()
        removeAdornment = True

        #Correct Answer
        answer = "[[  1.   0.   3.   0.   5.]\n [  0.  10.   0.  45.   0.]\n [  0.   0.   0.   0.   0.]\n [ 12.  43.  29.  34.   2.]\n [  0.   6.   0.   3.   0.]\n [  6.   0.  45.   0.   6.]\n [  0.  23.   0.   0.   0.]\n [  2.   0.   1.   0.   1.]\n [  1.   1.   1.   1.   1.]\n [  0.   0.   0.   0.   0.]]"

        #Call method
        result = data.textToStructuredArray(tempInputFile=inputFile, tempDelimiter=delimiter, tempNameRow=nameRow, tempFirstDataRow=firstDataRow, tempNormalize=normalize)
        result = data.transposeDataMatrix(tempMatrix=result[0], tempRemoveAdornments=removeAdornment)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(AbundanceTableTest)
    return suite
