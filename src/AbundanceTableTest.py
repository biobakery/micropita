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
import numpy as np
import os
import re
import unittest

##
#Tests the Blog object
class AbundanceTableTest(unittest.TestCase):

    #Test constructors
    def testInitForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        result = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        #Correct Answer
        answer = "".join(["Sample count:10",
                 "\nFeature count:5",
                 "\nId Metadata:TID",
                 "\nMetadata ids:['TID', 'STSite']",
                 "\nMetadata count:2",
                 "\nOriginating source:",
                 Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt",
                 "\nOriginal feature count:5",
                 "\nOriginal sample count:10",
                 "\nIs normalized:",str(fIsNormalized),
                 "\nIs summed:",str(fIsSummed),
                 "\nCurrent filtering state:",
                 "\nFeature delimiter:|",
                 "\nFile delimiter:",Constants.TAB])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testInitForGoodCase2(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = True
        fIsNormalized = True

        result = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        #Correct Answer
        answer = "".join(["Sample count:10",
                 "\nFeature count:5",
                 "\nId Metadata:TID",
                 "\nMetadata ids:['TID', 'STSite']",
                 "\nMetadata count:2",
                 "\nOriginating source:",
                 Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt",
                 "\nOriginal feature count:5",
                 "\nOriginal sample count:10",
                 "\nIs normalized:",str(fIsNormalized),
                 "\nIs summed:",str(fIsSummed),
                 "\nCurrent filtering state:",
                 "\nFeature delimiter:|",
                 "\nFile delimiter:",Constants.TAB])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testInitForGoodCase3(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False

        result = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        #Correct Answer
        answer = "".join(["Sample count:10",
                 "\nFeature count:5",
                 "\nId Metadata:TID",
                 "\nMetadata ids:['TID', 'STSite']",
                 "\nMetadata count:2",
                 "\nOriginating source:",
                 Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt",
                 "\nOriginal feature count:5",
                 "\nOriginal sample count:10",
                 "\nIs normalized:",str(fIsNormalized),
                 "\nIs summed:",str(fIsSummed),
                 "\nCurrent filtering state:",
                 "\nFeature delimiter:|",
                 "\nFile delimiter:",Constants.TAB])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test private methods
    #TestfuncTextToStructuredArray
    def testfuncTextToStructuredArrayForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), {'TID': ['700098986', '700098984', '700098982', '700098980', '700098988', '700037470', '700037472', '700037474', '700037476', '700037478'], 'STSite': ['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']}]"

        #Call method
        result = data._funcTextToStructuredArray(strInputFile=inputFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testfuncTextToStructuredArrayForGoodCaseSpaceDelimiter(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged_Space.txt"])
        delimiter = Constants.WHITE_SPACE
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[array([ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 1.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 2.0, 1.0, 0.0),\n       ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0),\n       ('Bacteria|3417', 0.0, 45.0, 0.0, 34.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0),\n       ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)], \n      dtype=[('TID', '|S158'), ('700098986', '<f8'), ('700098984', '<f8'), ('700098982', '<f8'), ('700098980', '<f8'), ('700098988', '<f8'), ('700037470', '<f8'), ('700037472', '<f8'), ('700037474', '<f8'), ('700037476', '<f8'), ('700037478', '<f8')]), {'TID': ['700098986', '700098984', '700098982', '700098980', '700098988', '700037470', '700037472', '700037474', '700037476', '700037478'], 'STSite': ['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']}]"

        #Call method
        result = data._funcTextToStructuredArray(strInputFile=inputFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #Test public methods
    #Test funcGetSampleNames
    def testFuncGetSampleNames(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        result = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = result.funcGetSampleNames()

        #Correct Answer
        answer = "('700098986', '700098984', '700098982', '700098980', '700098988', '700037470', '700037472', '700037474', '700037476', '700037478')"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetIDMetadataName
    def testFuncGetIDMetadataName(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        result = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = result.funcGetIDMetadataName()

        #Correct Answer
        answer = sNameRow

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetAbundanceCopy()
    def testFuncGetAbundanceCopyEquals(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetAbundanceCopy()

        #Correct Answer
        answer = abndData._npaFeatureAbundance

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetAbundanceCopyIsACopy(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetAbundanceCopy()
        result[0][2] = 99999999.9

        #Correct Answer
        answer = abndData._npaFeatureAbundance

        #Check result against answer
        self.assertNotEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetAverageAbundancePerSampleForGoodCase1Feature(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"]
        answer= [["700098980",12.0],["700037470",6.0],["700098986",1.0],["700098988",1.0],["700098982",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = abndData.funcGetAverageAbundancePerSample(liFeatures)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetAverageAbundancePerSampleForGoodCase2Feature(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361"]
        answer= [["700037470",25.5],["700098980",20.5],["700098986",2.0],["700098988",2.0],["700098982",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = abndData.funcGetAverageAbundancePerSample(liFeatures)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetAverageAbundancePerSampleForGoodCaseAllFeature(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                      "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                      "Bacteria|unclassified|4904","Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368",
                      "Bacteria|3417"]
        answer= [["700098980",24.0],["700037470",11.4],["700098988",3.0],["700098986",1.8],["700098982",0.0]]

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = abndData.funcGetAverageAbundancePerSample(liFeatures)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetAverageAbundancePerSampleForBadCaseBadFeature(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/TestFeatureAverages.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        liFeatures = [""]
        answer= False

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)

        result = abndData.funcGetAverageAbundancePerSample(liFeatures)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureAbundanceTable
    def testFuncGetFeatureAbundanceTableForGoodCase(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        #Correct Answer
        answer1 = "".join(["[ ('Bacteria|unclassified|4904', 0.0, 10.0, 0.0, 43.0, 6.0, 0.0, 23.0, 0.0, 1.0, 0.0)\n",
                  " ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 3.0, 0.0, 0.0, 29.0, 0.0, 45.0, 0.0, 1.0, 1.0, 0.0)\n",
                  " ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 5.0, 0.0, 0.0, 2.0, 0.0, 6.0, 0.0, 1.0, 1.0, 0.0)]"])

        answer2 = "".join(["Sample count:10",
                 "\nFeature count:3",
                 "\nId Metadata:TID",
                 "\nMetadata ids:['TID', 'STSite']",
                 "\nMetadata count:2",
                 "\nOriginating source:",
                 Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-3-Feature.txt",
                 "\nOriginal feature count:5",
                 "\nOriginal sample count:10",
                 "\nIs normalized:",str(fIsNormalized),
                 "\nIs summed:",str(fIsSummed),
                 "\nCurrent filtering state:",
                 "\nFeature delimiter:|",
                 "\nFile delimiter:",Constants.TAB])

        #Check result against answer
        self.assertNotEqual(str(result.funcGetAbundanceCopy())+str(result),str(answer1+answer2),"".join([str(self),"::\nExpected=\n",str(answer1+answer2),". \nReceived=\n",str(result.funcGetAbundanceCopy())+str(result),"."]))

    def testFuncGetFeatureAbundanceTableForGoodCaseCheckMetadata(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        result = result.funcGetMetadataCopy()
        result = [result[key] for key in sorted(result.keys())]

        #Correct Answer
        answer = ""

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::\nExpected=\n",answer,". \nReceived=\n",result,"."]))

    def testFuncGetFeatureAbundanceTableForGoodCaseCheckOtheData(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        result = str(result)

        #Correct Answer
        answer = ""

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::\nExpected=\n",answer,". \nReceived=\n",result,"."]))

    def testFuncGetFeatureAbundanceTableForGoodCaseCheckAbudanceData(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        result = str(result.funcGetAbundanceCopy())

        #Correct Answer
        answer = ""

        #Check result against answer
        self.assertEqual(result,answer,"".join([str(self),"::\nExpected=\n",answer,". \nReceived=\n",result,"."]))

    #Test funcGetFeatureDelimiter
    def testFuncGetFeatureDelimiter(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureDelimiter()

        #Correct Answer
        answer = cFeatureDelimiter

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetFeatureDelimiterForPow(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "!"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureDelimiter()

        #Correct Answer
        answer = cFeatureDelimiter

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureCount
    def testFuncGetFeatureCount(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureCount()

        #Correct Answer
        answer = 5

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetFeatureCountFor3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        result = abndData.funcGetFeatureCount()

        #Correct Answer
        answer = 3

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureSumAcrossSamples
    def testFuncGetFeatureSumAcrossSamplesGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-FeatureSum.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sFeature = "Bacteria|3417"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureSumAcrossSamples(sFeature)

        #Correct Answer
        answer = 83.0

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureSumAcrossSamples
    def testFuncGetFeatureSumAcrossSamplesGoodCase2(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-FeatureSum.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sFeature = "Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureSumAcrossSamples(sFeature)

        #Correct Answer
        answer = 22.0

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureSumAcrossSamples
    def testFuncGetFeatureSumAcrossSamplesGoodCase3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-FeatureSum.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sFeature = "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureSumAcrossSamples(sFeature)

        #Correct Answer
        answer = 15.0

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureSumAcrossSamples
    def testFuncGetFeatureSumAcrossSamplesGoodCase0(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-FeatureSum.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sFeature = "Bacteria|Lost"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureSumAcrossSamples(sFeature)

        #Correct Answer
        answer = 0.0

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFeatureNames
    def testFuncGetFeatureNames(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFeatureNames()

        #Correct Answer
        answer = str("".join(["[ \'Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72\'\n ",
                  "\'Bacteria|unclassified|4904\'\n ",
                  "\'Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\'\n ",
                  "\'Bacteria|3417\'\n ",
                  "\'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\']"]))

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetFeatureNamesFor3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData = abndData.funcGetFeatureAbundanceTable(["Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                                                        "Bacteria|unclassified|4904",
                                                        "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"])

        result = abndData.funcGetFeatureNames()

        #Correct Answer
        answer = str("".join(["[\'Bacteria|unclassified|4904\'\n ",
                  "\'Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361\'\n ",
                  "\'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368\']"]))

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetFileDelimiter
    def testFuncGetFileDelimiter(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetFileDelimiter()

        #Correct Answer
        answer = delimiter

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetSample
    def testFuncGetSampleForGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700098986"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetSample(sSample)

        #Correct Answer
        answer = str(np.array([1.0,0.0,3.0,0.0,5.0]))

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetSampleForGoodCase2(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037470"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetSample(sSample)

        #Correct Answer
        answer = str(np.array([6.0,0.0,45.0,0.0,6.0]))

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetSampleForGoodCase3(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetSample(sSample)

        #Correct Answer
        answer = str(np.array([0.0,0.0,0.0,0.0,0.0]))

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetMetadata
    def testFuncGetMetadataForGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetMetadata(sLastMetadata)

        #Correct Answer
        answer = str(['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease', 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease', 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares'])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetMetadataForBadCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetMetadata("error")

        #Correct Answer
        answer = None

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetMetadataForBadCase2(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetMetadata(None)

        #Correct Answer
        answer = None

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncGetMetadataCopyForGoodCase(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetMetadataCopy()

        #Correct Answer
        answer = "".join(["{'TID': ['700098986', '700098984', '700098982', '700098980', '700098988',",
                          " '700037470', '700037472', '700037474', '700037476', '700037478'],",
                          " 'STSite': ['L_Antecubital_fossa', 'R_Retroauricular_crease', 'L_Retroauricular_crease',",
                          " 'Subgingival_plaque', 'R_Antecubital_fossa', 'L_Retroauricular_crease',",
                          " 'R_Retroauricular_crease', 'L_Antecubital_fossa', 'R_Antecubital_fossa', 'Anterior_nares']}"])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetMetadataCopy
    def testFuncGetMetadataCopyForGoodCaseISACopy(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetMetadataCopy()
        result["TID"]= 1

        #Correct Answer
        answer = abndData._npaFeatureAbundance

        #Check result against answer
        self.assertNotEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetName
    def testFuncGetNameForGoodCase(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetName()

        #Correct Answer
        answer = inputFile

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetTerminalNode
    def testFuncGetTerminalNodesForGoodCase(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcGetTerminalNodes()

        #Correct Answer
        answer = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                  "Bacteria|unclassified|4904",
                  "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                  "Bacteria|3417",
                  "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"]
        result.sort()
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcGetTerminalNode
    def testFuncGetTerminalNodesForGoodCase2(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        abndData.funcSumClades()
        result = abndData.funcGetTerminalNodes()

        #Correct Answer
        answer = ["Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72",
                  "Bacteria|unclassified|4904",
                  "Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361",
                  "Bacteria|3417",
                  "Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368"]
        result.sort()
        answer.sort()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))


    #Test funcIsNormalized
    def testFuncIsNormalizedForGoodCaseTrue(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = True
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsNormalized()

        #Correct Answer
        answer = fIsNormalized

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncIsNormalizedForGoodCaseFalse(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsNormalized()

        #Correct Answer
        answer = fIsNormalized

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcIsPrimaryMetadata
    def testFuncIsPrimaryIdMetadataForGoodCaseTrue(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsPrimaryIdMetadata(sNameRow)

        #Correct Answer
        answer = True

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncIsPrimaryIdMetadataForGoodCaseFalse(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsPrimaryIdMetadata(sLastMetadata)

        #Correct Answer
        answer = False

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test funcIsSummed
    def testFuncIsSummedForGoodCaseFalse(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = False
        fIsNormalized = False
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsSummed()

        #Correct Answer
        answer = False

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    def testFuncIsSummedForGoodCaseTrue(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        sNameRow = "TID"
        sLastMetadata = "STSite"
        cFeatureDelimiter = "|"
        fIsSummed = True
        fIsNormalized = False
        sSample = "700037478"

        abndData = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=fIsNormalized, fIsSummed=fIsSummed,
                                             cDelimiter = delimiter, sMetadataID = sNameRow,
                                             sLastMetadata = sLastMetadata, cFeatureNameDelimiter=cFeatureDelimiter)
        result = abndData.funcIsSummed()

        #Correct Answer
        answer = True

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::\nExpected=\n",str(answer),". \nReceived=\n",str(result),"."]))

    #Test filterAbundanceByPercentile
    def testFilterAbundanceByPercentileForGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestPercentileFiltering-Answer1.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 50
        dPercentageAboveThreshold = 60

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

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
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 70
        dPercentageAboveThreshold = 70

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

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
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 20

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

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
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 20

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

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
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 90
        dPercentageAboveThreshold = 100

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

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
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dPercentileCutOffThreshold = 50
        dPercentageAboveThreshold = 30

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

        #Call method
        data.funcFilterAbundanceByPercentile(dPercentileCutOff=dPercentileCutOffThreshold, dPercentageAbovePercentile=dPercentageAboveThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    #Test filterAbundanceBySequenceOccurence
    def testFilterAbundanceBySequenceOccurenceForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering-Answer.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        iMinSampleThreshold = 2
        iMinSequenceThreshold = 2

        #Correct Answer
        abndDataAnswer, metadataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

        #Call method
        data.funcFilterAbundanceBySequenceOccurence(iMinSequence = iMinSequenceThreshold, iMinSamples = iMinSampleThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer),"".join([str(self),"::Expected=",str(abndDataAnswer),". Received=",str(result),"."]))

    def testFilterAbundanceBySequenceOccurenceForGoodCaseNoFiltering(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestOccurenceFiltering.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        iMinSampleThreshold = 0
        iMinSequenceThreshold = 0

        #Correct Answer
        abndDataAnswer, metadataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

        #Call method
        data.funcFilterAbundanceBySequenceOccurence(iMinSequence = iMinSequenceThreshold, iMinSamples = iMinSampleThreshold)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer),"".join([str(self),"::Expected=",str(abndDataAnswer),". Received=",str(result),"."]))

    #Test filterFeatureBySD
    def testFilterFeatureBYSDForGoodCase1(self):
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        strAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-TestSDFiltering-Answer.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Filter options
        dMinSD = 5

        #Correct Answer
        abndDataAnswer = data._funcTextToStructuredArray(strInputFile=strAnswerFile, cDelimiter=delimiter, sMetadataID = nameRow, sLastMetadata = firstDataRow)  

        #Call method
        data.funcFilterFeatureBySD(dMinSDCuttOff=dMinSD)
        result = data.funcGetAbundanceCopy()
        
        #Check result against answer
        self.assertEqual(str(result),str(abndDataAnswer[0]),"".join([str(self),"::Expected=",str(abndDataAnswer[0]),". Received=",str(result),"."]))

    #Test funcNormalize
    def testFuncNormalizeForGoodCaseSummed(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = Constants.ENDLINE.join(["[('Bacteria', 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)",
                   " ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes', 1.0, 0.0, 0.0, 0.35833333333333334, 0.0, 1.0, 0.0, 1.0, 0.6, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli', 0.8888888888888888, 0.0, 0.0, 0.25833333333333336, 0.0, 0.8947368421052632, 0.0, 0.5, 0.4, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)",
                   " ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)]"])


        #Call method
        data.funcSumClades()
        data.funcNormalize()
        result = data.funcGetAbundanceCopy()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncNormalizeForGoodCaseNotSummed(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = Constants.ENDLINE.join(["[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)",
                 " ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)",
                 " ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)",
                 " ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)",
                 " ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)]"])

        #Call method
        data.funcNormalize()
        result = data.funcGetAbundanceCopy()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #Test funcNormalizeColumnsBySum
    def testFuncNormalizeColumnsBySumForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = "[ ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)\n ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)\n ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)\n ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)]"

        #Call method
        data.funcNormalizeColumnsBySum()
        result = data.funcGetAbundanceCopy()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #Test funcNormalizeColumnsWithSummed
    def testFuncNormalizeColumnsWithSummedForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Correct Answer
        answer = Constants.ENDLINE.join(["[('Bacteria', 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)",
                   " ('Bacteria|3417', 0.0, 0.8181818181818182, 0.0, 0.2833333333333333, 0.3333333333333333, 0.0, 0.0, 0.0, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes', 1.0, 0.0, 0.0, 0.35833333333333334, 0.0, 1.0, 0.0, 1.0, 0.6, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli', 0.8888888888888888, 0.0, 0.0, 0.25833333333333336, 0.0, 0.8947368421052632, 0.0, 0.5, 0.4, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|1368', 0.5555555555555556, 0.0, 0.0, 0.016666666666666666, 0.0, 0.10526315789473684, 0.0, 0.25, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes|Bacilli|Lactobacillales|Lactobacillaceae|Lactobacillus|1361', 0.3333333333333333, 0.0, 0.0, 0.24166666666666667, 0.0, 0.7894736842105263, 0.0, 0.25, 0.2, 0.0)",
                   " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Clostridiaceae|Clostridium|72', 0.1111111111111111, 0.0, 0.0, 0.1, 0.0, 0.10526315789473684, 0.0, 0.5, 0.2, 0.0)",
                   " ('Bacteria|unclassified|4904', 0.0, 0.18181818181818182, 0.0, 0.35833333333333334, 0.6666666666666666, 0.0, 1.0, 0.0, 0.2, 0.0)]"])


        #Call method
        data.funcSumClades()
        data.funcNormalizeColumnsWithSummedClades()
        result = data.funcGetAbundanceCopy()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    def testFuncNormalizeColumnsWithSummedCladesForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/test1.otu"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "#OTU ID"
        firstDataRow = "#OTU ID"

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        data.funcSumClades()
        data.funcNormalizeColumnsWithSummedClades()
        result = data._npaFeatureAbundance

        #Correct Answer
        answer = Constants.ENDLINE.join(["[('Bacteria', 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)",
                  " ('Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales', 0.8958333333333334, 1.0, 1.0, 0.0, 1.0, 0.3333333333333333, 0.0, 0.9210526315789473, 0.6, 0.0)",
                  " ('Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|1008', 0.8958333333333334, 0.9850746268656716, 0.9375, 0.0, 0.0, 0.0, 0.0, 0.9210526315789473, 0.6, 0.0)",
                  " ('Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Prevotellaceae|Prevotella', 0.0, 0.014925373134328358, 0.0625, 0.0, 1.0, 0.3333333333333333, 0.0, 0.0, 0.0, 0.0)",
                  " ('Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Prevotellaceae|Prevotella|1010', 0.0, 0.014925373134328358, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)",
                  " ('Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Prevotellaceae|Prevotella|1013', 0.0, 0.0, 0.0625, 0.0, 0.0, 0.3333333333333333, 0.0, 0.0, 0.0, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia', 0.10416666666666667, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.3333333333333333, 0.013157894736842105, 0.3, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae', 0.09375, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.3333333333333333, 0.013157894736842105, 0.2, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|1000', 0.052083333333333336, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|1029', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Ruminococcus|1023', 0.041666666666666664, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.013157894736842105, 0.0, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|unclassified|101', 0.0, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.0, 0.0, 0.0, 0.0)",
                  " ('Bacteria|Firmicutes|Clostridia|unclassified|Lachnospiraceae|1034', 0.010416666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0)",
                  " ('Bacteria|Tenericutes', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.06578947368421052, 0.1, 1.0)",
                  " ('Bacteria|Tenericutes|Erysipelotrichi|Erysipelotrichales|Erysipelotrichaceae|Clostridium|1026', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6666666666666666, 0.06578947368421052, 0.0, 1.0)",
                  " ('Bacteria|Tenericutes|Mollicutes|RF39|unclassified|1', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0)",
                  " ('unclassified|1035', 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0)]"])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #TestFuncRankAbundance
    def testFuncRankAbundanceForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-ForRanked.txt"])

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        data = data.funcRankAbundance()
        data._strOriginalName = ""
        resultStr = str(data)+str(data.funcGetAbundanceCopy())
        

        #Correct Answer
        answer = AbundanceTable.funcMakeFromFile(strInputFile=sAnswerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        answer._strOriginalName = ""
        answerStr = str(answer)+str(answer.funcGetAbundanceCopy())

        #Check result against answer
        self.assertEqual(resultStr,answerStr,"".join([str(self),"::Expected=",answerStr,". Received=",resultStr,"."]))

    #Test funcReduceFeaturesToCladeLevelForLevel1
    def testFuncReduceFeaturesToCladeLevelForLevel1(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeClades.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        iCladeLevel = 1
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeCladesAnswerClade1.txt"])

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Reduce Clades
        fError = data.funcReduceFeaturesToCladeLevel(iCladeLevel)
        resultStr = str(data.funcGetAbundanceCopy())

        #Correct Answer
        answer = AbundanceTable.funcMakeFromFile(strInputFile=sAnswerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        answerStr = str(answer.funcGetAbundanceCopy())

        #Check result against answer
        self.assertEqual(resultStr,answerStr,"".join([str(self),"::Expected=",answerStr,". Received=",resultStr,"."]))

    #Test funcReduceFeaturesToCladeLevelForLevel3
    def testFuncReduceFeaturesToCladeLevelForLevel3(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeClades.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        iCladeLevel = 3
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeCladesAnswerClade3.txt"])

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Reduce Clades
        fError = data.funcReduceFeaturesToCladeLevel(iCladeLevel)
        resultStr = str(data.funcGetAbundanceCopy())

        #Correct Answer
        answer = AbundanceTable.funcMakeFromFile(strInputFile=sAnswerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        answerStr = str(answer.funcGetAbundanceCopy())

        #Check result against answer
        self.assertEqual(resultStr,answerStr,"".join([str(self),"::Expected=",answerStr,". Received=",resultStr,"."]))

    #Test funcReduceFeaturesToCladeLevelForLevel5
    def testFuncReduceFeaturesToCladeLevelForLevel5(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeClades.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        iCladeLevel = 5
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeCladesAnswerClade5.txt"])

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Reduce Clades
        fError = data.funcReduceFeaturesToCladeLevel(iCladeLevel)
        resultStr = str(data.funcGetAbundanceCopy())

        #Correct Answer
        answer = AbundanceTable.funcMakeFromFile(strInputFile=sAnswerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        answerStr = str(answer.funcGetAbundanceCopy())

        #Check result against answer
        self.assertEqual(resultStr,answerStr,"".join([str(self),"::Expected=",answerStr,". Received=",resultStr,"."]))

    #Test funcReduceFeaturesToCladeLevelForLevel10
    def testFuncReduceFeaturesToCladeLevelForLevel10(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeClades.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        iCladeLevel = 10
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-FakeCladesAnswerClade10.txt"])

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")

        #Reduce Clades
        fError = data.funcReduceFeaturesToCladeLevel(iCladeLevel)
        resultStr = str(data.funcGetAbundanceCopy())

        #Correct Answer
        answer = AbundanceTable.funcMakeFromFile(strInputFile=sAnswerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        answerStr = str(answer.funcGetAbundanceCopy())

        #Check result against answer
        self.assertEqual(resultStr,answerStr,"".join([str(self),"::Expected=",answerStr,". Received=",resultStr,"."]))

    #Test funcSumClades
    def testFuncSumCladesForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/test1.otu"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "#OTU ID"
        firstDataRow = "#OTU ID"

        #Make array
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        data.funcSumClades()
        result = data._npaFeatureAbundance

        #Correct Answer
        answerFile = "".join([Constants_Testing.c_strTestingTruth,"test1otu-summed.txt"])
        answer = AbundanceTable.funcMakeFromFile(strInputFile=answerFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")._npaFeatureAbundance

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcStratifyByMetadata(self,strMetadata,xWriteToFile=False)
    def testFuncStratifyByMetadataForGoodCaseNoFile(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        fWriteToFile = False

        #Remove any test files and check to see if this generates files
        testFile1= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Retroauricular_crease.txt"])
        testFile2= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Retroauricular_crease.txt"])
        testFile3= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Antecubital_fossa.txt"])
        testFile4= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Antecubital_fossa.txt"])
        testFile5= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Anterior_nares.txt"])
        testFile6= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Subgingival_plaque.txt"])

        for sFile in [testFile1,testFile2,testFile3,testFile4,testFile5,testFile6]:
            if os.path.exists(sFile):
                os.remove(sFile)

        #Get result
        table = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = table.funcStratifyByMetadata(strMetadata=firstDataRow,fWriteToFile=fWriteToFile)
        result = Constants.ENDLINE.join([str(abndTable) for abndTable in result])

        answer = ""
        #Check to see if files were generated
        for sFile in [testFile1,testFile2,testFile3,testFile4,testFile5]:
            if os.path.exists(sFile):
                answer = "".join(["The following file was created and should not have been :",sFile,"."])

        #Answer
        answer = "".join([answer,"Sample count:2\nFeature count:5\nId Metadata:TID\nMetadata ids:['TID', 'STSite']\n"
                          "Metadata count:2\nOriginating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Retroauricular_crease.txt\n",
                          "Original feature count:5\nOriginal sample count:2\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n"
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB,"\nSample count:2\nFeature count:5\nId Metadata:TID\n",
                          "Metadata ids:['TID', 'STSite']\nMetadata count:2\n",
                          "Originating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Retroauricular_crease.txt\n",
                          "Original feature count:5\nOriginal sample count:2\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n"
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB,"\nSample count:2\nFeature count:5\nId Metadata:TID\n"
                          "Metadata ids:['TID', 'STSite']\nMetadata count:2\n",
                          "Originating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Antecubital_fossa.txt\n",
                          "Original feature count:5\nOriginal sample count:2\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n",
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB,"\nSample count:2\nFeature count:5\nId Metadata:TID\n",
                          "Metadata ids:['TID', 'STSite']\nMetadata count:2\n",
                          "Originating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Antecubital_fossa.txt\n",
                          "Original feature count:5\nOriginal sample count:2\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n",
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB,"\nSample count:1\nFeature count:5\nId Metadata:TID\n",
                          "Metadata ids:['TID', 'STSite']\nMetadata count:2\n",
                          "Originating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Anterior_nares.txt\n",
                          "Original feature count:5\nOriginal sample count:1\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n",
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB,"\nSample count:1\nFeature count:5\nId Metadata:TID\nMetadata ids:['TID', 'STSite']\n",
                          "Metadata count:2\nOriginating source:input/micropita/src/Testing/Data/AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Subgingival_plaque.txt\n",
                          "Original feature count:5\nOriginal sample count:1\nIs normalized:False\nIs summed:False\nCurrent filtering state:\n",
                          "Feature delimiter:|\nFile delimiter:",Constants.TAB])

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcStratifyByMetadata(self,strMetadata,xWriteToFile=False)
    def testFuncStratifyByMetadataForGoodCaseCreateFile(self):

        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        fWriteToFile = True

        #Remove any test files and check to see if this generates files
        testFile1= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Retroauricular_crease.txt"])
        testFile2= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Retroauricular_crease.txt"])
        testFile3= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Antecubital_fossa.txt"])
        testFile4= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Antecubital_fossa.txt"])
        testFile5= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Anterior_nares.txt"])
        testFile6= "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Subgingival_plaque.txt"])

        answerFile1= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Retroauricular_crease.txt"])
        answerFile2= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Retroauricular_crease.txt"])
        answerFile3= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-R_Antecubital_fossa.txt"])
        answerFile4= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-L_Antecubital_fossa.txt"])
        answerFile5= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Anterior_nares.txt"])
        answerFile6= "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-StratBy-Subgingival_plaque.txt"])

        #Delete any previous files if exist
        for sFile in [testFile1,testFile2,testFile3,testFile4,testFile5,testFile6]:
            if os.path.exists(sFile):
                os.remove(sFile)

        #Get result
        table = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = table.funcStratifyByMetadata(strMetadata=firstDataRow,fWriteToFile=fWriteToFile)
        result = Constants.ENDLINE.join([str(abndTable) for abndTable in result])

        result = ""
        answer = ""
        #Check to see if files were generated
        for sFile, answerFile in [(testFile1,answerFile1),(testFile2,answerFile2),(testFile3,answerFile3),(testFile4,answerFile4),(testFile5,answerFile5),(testFile6,answerFile6)]:
            if not os.path.exists(sFile):
                result = result+"".join(["Did not generate the following file:",sFile,". "])

            sResultString = ""
            sAnswerString = ""

            with open(sFile, 'r') as f:
                sResultString = f.read()+Constants.ENDLINE
            with open(answerFile, 'r') as f:
                sAnswerString = f.read()

            if not sResultString == sAnswerString:
                result = "".join([result,"Expected=",sAnswerString," Received=",sResultString," "])

        #Clean up Delete any previous files if exist
        for sFile in [testFile1,testFile2,testFile3,testFile4,testFile5]:
            if os.path.exists(sFile):
                os.remove(sFile)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcTranslateIntoMetadata
    def testFuncTranslateIntoMetadataForGoodCase1(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"

        sFrom = nameRow
        sTo = firstDataRow
        fPrimary = False
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcTranslateIntoMetadata(lsValues=data.funcGetMetadata(sFrom),
                                                sMetadataFrom=sFrom, sMetadataTo=sTo, fFromPrimaryIds=fPrimary)

        #Correct Answer
        answer = data.funcGetMetadata(sTo)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcTranslateIntoMetadata
    def testFuncTranslateIntoMetadataForGoodCase2(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"

        sFrom = firstDataRow
        sTo = nameRow
        fPrimary = False
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcTranslateIntoMetadata(lsValues=data.funcGetMetadata(sFrom),
                                                sMetadataFrom=sFrom, sMetadataTo=sTo, fFromPrimaryIds=fPrimary)
        #Correct Answer
        answer = "['700098986', '700098984', '700098982', '700098980', '700098988', '700098982', '700098984', '700098986', '700098988', '700037478']"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcTranslateIntoMetadata
    def testFuncTranslateIntoMetadataForGoodCase3(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"

        sFrom = nameRow
        sTo = firstDataRow
        fPrimary = True
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcTranslateIntoMetadata(lsValues=data.funcGetMetadata(sFrom),
                                                sMetadataFrom=sFrom, sMetadataTo=sTo, fFromPrimaryIds=fPrimary)

        #Correct Answer
        answer = data.funcGetMetadata(sTo)

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcTranslateIntoMetadata
    def testFuncTranslateIntoMetadataForGoodCase4(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"

        sFrom = firstDataRow
        sTo = nameRow
        fPrimary = True
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcTranslateIntoMetadata(lsValues=data.funcGetMetadata(sFrom),
                                                sMetadataFrom=sFrom, sMetadataTo=sTo, fFromPrimaryIds=fPrimary)

        #Correct Answer
        answer = False

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcTranslateIntoMetadata
    def testFuncTranslateIntoMetadataForGoodCase5Reverse(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"

        sFrom = nameRow
        sTo = firstDataRow
        fPrimary = False

        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        lsRevValues = data.funcGetMetadata(sFrom)
        lsRevValues.reverse()
        result = data.funcTranslateIntoMetadata(lsValues=lsRevValues, sMetadataFrom=sFrom,
                                                sMetadataTo=sTo, fFromPrimaryIds=fPrimary)

        #Correct Answer
        answer = data.funcGetMetadata(sTo)
        answer.reverse()

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #Test funcToArray
    def testFuncToArrayForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        result = data.funcToArray()

        #Correct Answer
        answer = "[[  1.   0.   0.  12.   0.   6.   0.   2.   1.   0.]\n [  0.  10.   0.  43.   6.   0.  23.   0.   1.   0.]\n [  3.   0.   0.  29.   0.  45.   0.   1.   1.   0.]\n [  0.  45.   0.  34.   3.   0.   0.   0.   1.   0.]\n [  5.   0.   0.   2.   0.   6.   0.   1.   1.   0.]]"

        #Check result against answer
        self.assertEqual(str(result),str(answer),"".join([str(self),"::Expected=",str(answer),". Received=",str(result),"."]))

    #funcWriteToFile
    def testFuncWriteToForGoodCase(self):
        
        #Inputs
        inputFile = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/hq.otu_04-nul-nul-mtd-trn-flt-abridged.txt"])
        delimiter = Constants.TAB
        normalize = False
        nameRow = "TID"
        firstDataRow = "STSite"
        sWriteToTempFile = "".join([Constants_Testing.c_strTestingTMP,"testFuncWriteToForGoodCase.txt"])
        sAnswerFile = "".join([Constants_Testing.c_strTestingTruth,"hq.otu_04-nul-nul-mtd-trn-flt-abridged-ForWrite.txt"])
        data = AbundanceTable.funcMakeFromFile(strInputFile=inputFile, fIsNormalized=False, fIsSummed=False, cDelimiter = delimiter,
                                             sMetadataID = nameRow, sLastMetadata = firstDataRow, cFeatureNameDelimiter="|")
        data.funcWriteToFile(strOutputFile=sWriteToTempFile, cDelimiter=data.funcGetFileDelimiter())

        #Correct Answer
        sContents = None
        with open(sWriteToTempFile, 'r') as f:
            sContents = f.read()+Constants.ENDLINE

        sAnswer = None
        with open(sAnswerFile, 'r') as f:
            sAnswer = f.read()
        
        #Check result against answer
        self.assertEqual(str(sContents),str(sAnswer),"".join([str(self),"::Expected=",str(sAnswer),". Received=",str(sContents),"."]))

    #Test funcPairTables
    def testFuncPairTablesForGoodCase(self):
        
        #Inputs
        sInputFileOne = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/PairTables1.txt"])
        sInputFileTwo = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/PairTables2.txt"])
        cDelimiter = Constants.TAB
        sIdentifier = "ID"
        sOutputFileOne = "".join([Constants_Testing.c_strTestingTMP,"PairTables1with2.txt"])
        sOutputFileTwo = "".join([Constants_Testing.c_strTestingTMP,"PairTables2with1.txt"])

        #Answer files
        sAnswerFile1 = "".join([Constants_Testing.c_strTestingTruth,"PairedTables1with2-Correct.txt"])
        sAnswerFile2 = "".join([Constants_Testing.c_strTestingTruth,"PairedTables2with1-Correct.txt"])

        #Collects all the errors to display
        sError = ""

        #Pair files
        result = AbundanceTable.funcPairTables(strFileOne=sInputFileOne, strFileTwo=sInputFileTwo,
                                      strIdentifier=sIdentifier, cDelimiter=cDelimiter,
                                      strOutFileOne=sOutputFileOne, strOutFileTwo=sOutputFileTwo)
        if not result:
            sError = "AbundanceTableTest.testFuncPairTablesForGoodCase::Received false when calling the AbundanceTable.funcPairTables()"

        #Check answers
        if result:
            sCorrect = None
            sResult = None
            with open(sOutputFileOne, 'r') as f:
                sResult = f.read()
            f.close()
            with open(sAnswerFile1, 'r') as f:
                sCorrect = f.read()
            f.close()
            if not sCorrect.strip() == sResult:
                sError = " ".join([sError, "Did not receive the correct output for file 1. Expected:",sCorrect,". Received:",sResult,"."])

            sCorrect = None
            sResult = None
            with open(sOutputFileTwo, 'r') as f:
                sResult = f.read()
            f.close()
            with open(sAnswerFile2, 'r') as f:
                sCorrect = f.read()
            f.close()
            if not sCorrect.strip() == sResult:
                sError = " ".join([sError, "Did not receive the correct output for file 2. Expected:",sCorrect,". Received:",sResult,"."])

        #Check result against answer
        self.assertEqual("",sError,sError)

    #Test funcPairTables
    def testFuncPairTablesForBadCaseDuplicates(self):
        
        #Inputs
        sInputFileOne = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/PairTablesDup1.txt"])
        sInputFileTwo = "".join([Constants_Testing.c_strTestingData,"AbridgedDocuments/PairTablesDup2.txt"])
        cDelimiter = Constants.TAB
        sIdentifier = "ID"
        sOutputFileOne = "".join([Constants_Testing.c_strTestingTMP,"PairTablesDup1with2.txt"])
        sOutputFileTwo = "".join([Constants_Testing.c_strTestingTMP,"PairTablesDup2with1.txt"])

        #Answer files
        sAnswerFile1 = "".join([Constants_Testing.c_strTestingTruth,"PairedTablesDup1with2-Correct.txt"])
        sAnswerFile2 = "".join([Constants_Testing.c_strTestingTruth,"PairedTablesDup2with1-Correct.txt"])

        #Collects all the errors to display
        sError = ""

        #Pair files
        result = AbundanceTable.funcPairTables(strFileOne=sInputFileOne, strFileTwo=sInputFileTwo,
                                      strIdentifier=sIdentifier, cDelimiter=cDelimiter,
                                      strOutFileOne=sOutputFileOne, strOutFileTwo=sOutputFileTwo)
        if not result:
            sError = "AbundanceTableTest.testFuncPairTablesForGoodCase::Received false when calling the AbundanceTable.funcPairTables()"

        #Check answers
        if result:
            sCorrect = None
            sResult = None
            with open(sOutputFileOne, 'r') as f:
                sResult = f.read()
            f.close()
            with open(sAnswerFile1, 'r') as f:
                sCorrect = f.read()
            f.close()
            if not sCorrect.strip() == sResult:
                sError = " ".join([sError, "Did not receive the correct output for file 1. Expected:",sCorrect,". Received:",sResult,"."])

            sCorrect = None
            sResult = None
            with open(sOutputFileTwo, 'r') as f:
                sResult = f.read()
            f.close()
            with open(sAnswerFile2, 'r') as f:
                sCorrect = f.read()
            f.close()
            if not sCorrect.strip() == sResult:
                sError = " ".join([sError, "Did not receive the correct output for file 2. Expected:",sCorrect,". Received:",sResult,"."])

        #Check result against answer
        self.assertEqual("",sError,sError)

    #Test funcCheckRawDataFile
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

    #Test funcStratifyAbundanceTableByMetadata
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

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(AbundanceTableTest)
    return suite
