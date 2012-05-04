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
from AbundanceTable import AbundanceTable
from Constants import Constants
from Constants_Testing import Constants_Testing
import os
import re
import unittest

##
#Tests the Blog object
class AbundanceTableTest(unittest.TestCase):

    def testCheckRawDataFileForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,os.path.splitext(os.path.split(inputFile)[1])[0],Constants.OUTPUT_SUFFIX])
        delimiter = Constants.TAB

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\"\t700098986\t700098984\t700098982\t700098980\t700098988\t700037470\t700037472\t700037474\t700037476\t700037478\n\"STSite\"\t\"L_Antecubital_fossa\"\t\"R_Retroauricular_crease\"\t\"L_Retroauricular_crease\"\t\"Subgingival_plaque\"\t\"R_Antecubital_fossa\"\t\"L_Retroauricular_crease\"\t\"R_Retroauricular_crease\"\t\"L_Antecubital_fossa\"\t\"R_Antecubital_fossa\"\t\"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\"\t1\t0\t0\t12\t0\t6\t0\t2\t1\t0\n\"Bacteria|unclassified|4904\"\t0\t10\t0\t43\t6\t0\t23\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\"\t3\t0\t0\t29\t0\t45\t0\t1\t1\t0\n\"Bacteria|3417\"\t0\t45\t0\t34\t3\t0\t0\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\"\t5\t0\t0\t2\t0\t6\t0\t1\t1\t0\n"

        #Call method
        AbundanceTable.funcCheckRawDataFile(strReadDataFileName=inputFile, iFirstDataIndex=2, strOutputFileName = outputFile, cDelimiter=delimiter)

        #Get answer
        result = ""
        with open(outputFile,'r') as f:
            result = f.read()
            f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testCheckRawDataFileForReduceByZero(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking_ReduceZeros.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,os.path.splitext(os.path.split(inputFile)[1])[0],Constants.OUTPUT_SUFFIX])
        delimiter = Constants.TAB

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\"\t700098986\t700098984\t700098982\t700098980\t700098988\t700037470\t700037472\t700037474\t700037476\t700037478\n\"STSite\"\t\"L_Antecubital_fossa\"\t\"R_Retroauricular_crease\"\t\"L_Retroauricular_crease\"\t\"Subgingival_plaque\"\t\"R_Antecubital_fossa\"\t\"L_Retroauricular_crease\"\t\"R_Retroauricular_crease\"\t\"L_Antecubital_fossa\"\t\"R_Antecubital_fossa\"\t\"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\"\t1\t0\t0\t12\t0\t6\t0\t2\t1\t0\n\"Bacteria|unclassified|4904\"\t0\t10\t0\t43\t6\t0\t23\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\"\t3\t0\t0\t29\t0\t45\t0\t1\t1\t0\n\"Bacteria|3417\"\t0\t45\t0\t34\t3\t0\t0\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\"\t5\t0\t0\t2\t0\t6\t0\t1\t1\t0\n"

        #Call method
        AbundanceTable.funcCheckRawDataFile(strReadDataFileName=inputFile, iFirstDataIndex=2, strOutputFileName = outputFile, cDelimiter=delimiter)

        #Get answer
        result = ""
        with open(outputFile,'r') as f:
            result = f.read()
            f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testCheckRawDataFileForReduceByZeroAndMissingData(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking_MissingData.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,os.path.splitext(os.path.split(inputFile)[1])[0],Constants.OUTPUT_SUFFIX])
        delimiter = Constants.TAB

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\"\t700098986\t700098984\tNA\t700098980\t700098988\t700037470\t700037472\t700037474\tNA\t700037478\n\"STSite\"\t\"L_Antecubital_fossa\"\t\"R_Retroauricular_crease\"\t\"L_Retroauricular_crease\"\t\"Subgingival_plaque\"\t\"R_Antecubital_fossa\"\t\"L_Retroauricular_crease\"\t\"R_Retroauricular_crease\"\t\"L_Antecubital_fossa\"\t\"R_Antecubital_fossa\"\t\"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\"\t1\t0\t0\t12\t0\t6\t0\t2\t1\t0\n\"Bacteria|unclassified|4904\"\t0\t10\t0\t43\t6\t0\t23\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\"\t3\t0\t0\t29\t0\t45\t0\t1\t1\t0\n\"Bacteria|3417\"\t0\t45\t0\t34\t3\t0\t0\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\"\t5\t0\t0\t2\t0\t6\t0\t1\t1\t0\n"

        #Call method
        AbundanceTable.funcCheckRawDataFile(strReadDataFileName=inputFile, iFirstDataIndex=2, strOutputFileName = outputFile, cDelimiter=delimiter)

        #Get answer
        result = ""
        with open(outputFile,'r') as f:
            result = f.read()
            f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testCheckRawDataFileForReduceByZeroAndMissingDataOfDifferentLength(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking_MissingDataVary.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,os.path.splitext(os.path.split(inputFile)[1])[0],Constants.OUTPUT_SUFFIX])
        delimiter = Constants.TAB

        #Remove output file before the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Correct Answer
        answer = "\"TID\"\t700098986\t700098984\tNA\t700098980\t700098988\t700037470\t700037472\t700037474\tNA\t700037478\n\"STSite\"\t\"L_Antecubital_fossa\"\t\"R_Retroauricular_crease\"\t\"L_Retroauricular_crease\"\t\"Subgingival_plaque\"\t\"R_Antecubital_fossa\"\t\"L_Retroauricular_crease\"\t\"R_Retroauricular_crease\"\t\"L_Antecubital_fossa\"\t\"R_Antecubital_fossa\"\t\"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\"\t1\t0\t0\t12\t0\t6\t0\t2\t1\t0\n\"Bacteria|unclassified|4904\"\t0\t10\t0\t43\t6\t0\t23\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\"\t3\t0\t0\t29\t0\t45\t0\t1\t1\t0\n\"Bacteria|3417\"\t0\t45\t0\t34\t3\t0\t0\t0\t1\t0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\"\t5\t0\t0\t2\t0\t6\t0\t1\t1\t0\n"

        #Call method
        AbundanceTable.funcCheckRawDataFile(strReadDataFileName=inputFile, iFirstDataIndex=2, strOutputFileName = outputFile, cDelimiter=delimiter)

        #Get answer
        result = ""
        with open(outputFile,'r') as f:
            result = f.read()
            f.close()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testCheckRawDataFileForGoodCaseDelimiterSpace(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.ForChecking_Space.txt"])
        outputFile = "".join([Constants_Testing.c_strTestingTMP,os.path.splitext(os.path.split(inputFile)[1])[0],Constants.OUTPUT_SUFFIX])
        delimiter = Constants.WHITE_SPACE

        #Correct Answer
        answer = "\"TID\" 700098986 700098984 700098982 700098980 700098988 700037470 700037472 700037474 700037476 700037478\n\"STSite\" \"L_Antecubital_fossa\" \"R_Retroauricular_crease\" \"L_Retroauricular_crease\" \"Subgingival_plaque\" \"R_Antecubital_fossa\" \"L_Retroauricular_crease\" \"R_Retroauricular_crease\" \"L_Antecubital_fossa\" \"R_Antecubital_fossa\" \"Anterior_nares\"\n\"Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\" 1 0 0 12 0 6 0 2 1 0\n\"Bacteria|unclassified|4904\" 0 10 0 43 6 0 23 0 1 0\n\"Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\" 3 0 0 29 0 45 0 1 1 0\n\"Bacteria|3417\" 0 45 0 34 3 0 0 0 1 0\n\"Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\" 5 0 0 2 0 6 0 1 1 0\n"

        #Call method
        AbundanceTable.funcCheckRawDataFile(strReadDataFileName=inputFile, iFirstDataIndex=2, strOutputFileName = outputFile, cDelimiter=delimiter)

        #Get answer
        result = ""
        with open(outputFile,'r') as f:
            result = f.read()
            f.close()

        #Remove output file after the test
        if(os.path.exists(outputFile)):
            os.remove(outputFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testStratifyAbundanceTableByMetadataForGoodCaseByIndex(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        stratifyRow = 1
        strOutputDirectory = Constants_Testing.c_strTestingTMP

        #Correct Answer is blank indicating no error
        answer = ""

        #Should generate the following files
        anteriorNaresFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"])
        anteriorNaresFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"])
        lAntecubitalFossaFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"])
        lAntecubitalFossaFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"])
        lRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"])
        lRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"])
        rAntecubitalFossaFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"])
        rAntecubitalFossaFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"])
        rRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"])
        rRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"])
        subgingivalPlaqueFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"])
        subgingivalPlaqueFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"])

        dictCreatedFiles = {anteriorNaresFileName:anteriorNaresFileNameAnswer,
                            lAntecubitalFossaFileName:lAntecubitalFossaFileNameAnswer,
                            lRetroauricularCreaseFileName:lRetroauricularCreaseFileNameAnswer,
                            rAntecubitalFossaFileName:rAntecubitalFossaFileNameAnswer,
                            rRetroauricularCreaseFileName:rRetroauricularCreaseFileNameAnswer,
                            subgingivalPlaqueFileName:subgingivalPlaqueFileNameAnswer}

        #Delete files in abridged documents if they exist
        for strFile in dictCreatedFiles:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Call method
        AbundanceTable.funcStratifyAbundanceTableByMetadata(strInputFile = inputFile, strDirectory = strOutputDirectory, cDelimiter = delimiter, iStratifyByRow = stratifyRow)

        #Check file creation
        error = ""
        for strFile in dictCreatedFiles:
            if(os.path.exists(strFile)):

                contents = list()
                contentsAnswer = list()
                with open(strFile) as f:
                    contents = f.read()
                    contents = filter(None,re.split("\n",contents))
                    f.close()
                with open(dictCreatedFiles[strFile]) as f:
                    contentsAnswer = f.read()
                    contentsAnswer = filter(None,re.split("\n",contentsAnswer))
                    f.close()
                if(not contents == contentsAnswer):
                    error = error + "\nFile: "+strFile+"\nExpected:"+",".join(contentsAnswer)+".\nReceived:"+",".join(contents)+"."
            else:
                error = error + "\nFile count not be found. Path:"+strFile
        result = error

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testStratifyAbundanceTableByMetadataForGoodCaseByKeyWord(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        stratifyRow = "STSite"
        strOutputDirectory = Constants_Testing.c_strTestingTMP

        #Correct Answer is blank indicating no error
        answer = ""

        #Should generate the following files
        anteriorNaresFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"])
        anteriorNaresFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Anterior_nares.txt"])
        lAntecubitalFossaFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"])
        lAntecubitalFossaFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"])
        lRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"])
        lRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"])
        rAntecubitalFossaFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"])
        rAntecubitalFossaFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Antecubital_fossa.txt"])
        rRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"])
        rRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"])
        subgingivalPlaqueFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"])
        subgingivalPlaqueFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-Subgingival_plaque.txt"])

        dictCreatedFiles = {anteriorNaresFileName:anteriorNaresFileNameAnswer,
                            lAntecubitalFossaFileName:lAntecubitalFossaFileNameAnswer,
                            lRetroauricularCreaseFileName:lRetroauricularCreaseFileNameAnswer,
                            rAntecubitalFossaFileName:rAntecubitalFossaFileNameAnswer,
                            rRetroauricularCreaseFileName:rRetroauricularCreaseFileNameAnswer,
                            subgingivalPlaqueFileName:subgingivalPlaqueFileNameAnswer}

        #Delete files if they exist
        for strFile in dictCreatedFiles:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Call method
        AbundanceTable.funcStratifyAbundanceTableByMetadata(strInputFile = inputFile, strDirectory = strOutputDirectory, cDelimiter = delimiter, iStratifyByRow = stratifyRow)

        #Check file creation
        error = ""
        for strFile in dictCreatedFiles:
            if(os.path.exists(strFile)):

                contents = list()
                contentsAnswer = list()
                with open(strFile) as f:
                    contents = f.read()
                    contents = filter(None,re.split("\n",contents))
                    f.close()
                with open(dictCreatedFiles[strFile]) as f:
                    contentsAnswer = f.read()
                    contentsAnswer = filter(None,re.split("\n",contentsAnswer))
                    f.close()
                if(not contents == contentsAnswer):
                    error = error + "\nFile: "+strFile+"\nExpected:"+",".join(contentsAnswer)+".\nReceived:"+",".join(contents)+"."
            else:
                error = error + "\nFile count not be found. Path:"+strFile
        result = error

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testStratifyAbundanceTableByMetadataForGoodCaseByKeyWordAndGroup(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        stratifyRow = "STSite"
        strOutputDirectory = Constants_Testing.c_strTestingTMP

        llsGroupings = [["\"L_Antecubital_fossa\"","\"Anterior_nares\"","\"Subgingival_plaque\""],["\"R_Retroauricular_crease\"","\"R_Antecubital_fossa\""]]

        #Correct Answer is blank indicating no error
        answer = ""

        #Should generate the following files
        lAntecubitalFossaFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa.txt"])
        lAntecubitalFossaFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Antecubital_fossa-GROUPED.txt"])
        lRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease.txt"])
        lRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-L_Retroauricular_crease-GROUPED.txt"])
        rRetroauricularCreaseFileName = "".join([Constants_Testing.c_strTestingTMP,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease.txt"])
        rRetroauricularCreaseFileNameAnswer = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-by-R_Retroauricular_crease-GROUPED.txt"])

        dictCreatedFiles = {lAntecubitalFossaFileName:lAntecubitalFossaFileNameAnswer,
                            lRetroauricularCreaseFileName:lRetroauricularCreaseFileNameAnswer,
                            rRetroauricularCreaseFileName:rRetroauricularCreaseFileNameAnswer}

        #Delete files if they exist
        for strFile in dictCreatedFiles:
            if os.path.exists(strFile):
                os.remove(strFile)

        #Call method
        AbundanceTable.funcStratifyAbundanceTableByMetadata(strInputFile = inputFile, strDirectory = strOutputDirectory, cDelimiter = delimiter, iStratifyByRow = stratifyRow, llsGroupings = llsGroupings)

        #Check file creation
        error = ""
        for strFile in dictCreatedFiles:
            if(os.path.exists(strFile)):

                contents = list()
                contentsAnswer = list()
                with open(strFile) as f:
                    contents = f.read()
                    contents = filter(None,re.split("\n",contents))
                    f.close()
                with open(dictCreatedFiles[strFile]) as f:
                    contentsAnswer = f.read()
                    contentsAnswer = filter(None,re.split("\n",contentsAnswer))
                    f.close()
                if(not contents == contentsAnswer):
                    error = error + "\nFile: "+strFile+"\nExpected:"+",".join(contentsAnswer)+".\nReceived:"+",".join(contents)+"."
            else:
                error = error + "\nFile count not be found. Path:"+strFile
        result = error

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncNormalizeColumnsBySumForGoodCaseNormalizeAll(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)\n ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)]"

        #Call method
        data.funcNormalizeColumnsBySum()
        result = data.funcGetAbundanceCopy()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTextToStructuredArrayForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), {'STSite': ['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']}]"

        #Call method
        result = data._textToStructuredArray(strInputFile=inputFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testTextToStructuredArrayForGoodCaseSpaceDelimiter(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged_Space.txt"])
        delimiter = Constants.WHITE_SPACE
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), {'STSite': ['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']}]"

        #Call method
        result = data._textToStructuredArray(strInputFile=inputFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFilterAbundanceBySequenceOccurenceForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering-Answer.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        iMinSampleThreshold = 2
        iMinSequenceThreshold = 2

        #Correct Answer
        abndDataAnswer, metadataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceBySequenceOccurence(iMinSequence = iMinSequenceThreshold, iMinSamples = iMinSampleThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer),"".join([str(self),"::Expected=",str(abndDataAnswer),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer1.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 50
        dPercentageAboveThreshold = 60

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCase2(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer2.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 70
        dPercentageAboveThreshold = 70

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCase3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer3.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 20

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCase3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer3.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 20

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCaseAllGone(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-AnswerAllGone.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 100

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFilterAbundanceByPercentileForGoodCaseAll(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 50
        dPercentageAboveThreshold = 30

        #Correct Answer
        abndDataAnswer = data._textToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, iNameRow=nameRow, iFirstDataRow=firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    def testFuncToArrayForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = 0
        firstDataRow = 2
        normalize = False

        #Make array
        data = AbundanceTable.makeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             iNameRow = nameRow, iFirstDataRow = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcToArray()

        #Correct Answer
        answer = "[[  1.   0.   0.  12.   0.   6.   0.   2.   1.   0.]\n [  0.  10.   0.  43.   6.   0.  23.   0.   1.   0.]\n [  3.   0.   0.  29.   0.  45.   0.   1.   1.   0.]\n [  0.  45.   0.  34.   3.   0.   0.   0.   1.   0.]\n [  5.   0.   0.   2.   0.   6.   0.   1.   1.   0.]]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(AbundanceTableTest)
    return suite
