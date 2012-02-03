#######################################################
#
#	Title:		FileIOTest
#	Author:		Timothy Tickle
#	Date:		01/29/2010
#	Purpose:	Class to test the FileIO Class
#
#######################################################

#Import libaries
import Constants
import FileIO
import os
import unittest

##
#Used to test the FileIO class
class FileIOTest(unittest.TestCase):

    ##FileIO to test
    testedFileIO = None

    ##Test file name to write to
    TestFileName = Constants.Constants.TEMPORARY_TEST_FILE

    ##
    #Test the constructor for initialization FileName = valid Read = True, Write = True Append = True
    def testConstructorForInitialConstructionVTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertTrue(True,self.methodName+" No error occured passed.")

    ##
    #Test the constructor for initialization FileName = valid Read = True, Write = True Append = False
    def testConstructorForInitialConstructionVTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,False)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertTrue(True,self.methodName+" No error occured passed.")

    ##
    #Test the constructor for initialization FileName = valid Read = True, Write = False Append = True
    def testConstructorForInitialConstructionVTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,True)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertTrue(True,self.methodName+" No error occured passed.")

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = True Append = True
    def testConstructorForInitialConstructionVFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,True)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertTrue(True,self.methodName+" No error occured passed.")

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = False Append = True
    def testConstructorForInitialConstructionVFFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,False,False,True)
        except IOError, e:
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(True,self.methodName+" Error was throw and this is correct.")
        else:
            #Close FileIO (incase the no error is thrown)
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Fail if no error is thrown
            self.fail(self.methodName+" No read or write function should have created an error.")

    ##
    #Test the constructor for initialization FileName = valid Read = True, Write = False Append = False
    def testConstructorForInitialConstructionVTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = True Append = False
    def testConstructorForInitialConstructionVFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,False)

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertTrue(True,self.methodName+" No error occured passed.")

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = False Append = False
    def testConstructorForInitialConstructionVFFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,False,False,False)
        except IOError, e:
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(True,self.methodName+" Error was throw and this is correct.")
        else:
            #Close FileIO (incase the no error is thrown)
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Fail if no error is thrown
            self.fail(self.methodName+" No read or write function should have created an error.")

    ##
    #Test the constructor for initialization FileName = None Read = True, Write = True Append = True
    def testConstructorForInitialConstructionNTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,True,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = True, Write = True Append = True
    def testConstructorForInitialConstructionETTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionETTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",True,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = True, Write = True Append = True
    def testConstructorForInitialConstructionBTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",True,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = list() Read = True, Write = True Append = True
    def testConstructorForInitialConstructionITTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionITTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(list(),True,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = True, Write = True Append = False
    def testConstructorForInitialConstructionNTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,True,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = True, Write = True Append = False
    def testConstructorForInitialConstructionETTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionETTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",True,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = True, Write = True Append = False
    def testConstructorForInitialConstructionBTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",True,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = True, Write = False Append = True
    def testConstructorForInitialConstructionNTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,True,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = True, Write = False Append = True
    def testConstructorForInitialConstructionETFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionETFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",True,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = True, Write = False Append = True
    def testConstructorForInitialConstructionBTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",True,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = False, Write = True Append = True
    def testConstructorForInitialConstructionNFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,False,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = False, Write = True Append = True
    def testConstructorForInitialConstructionEFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionEFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",False,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = False, Write = True Append = True
    def testConstructorForInitialConstructionBFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",False,True,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = False, Write = False Append = True
    def testConstructorForInitialConstructionNFFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNFFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,False,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = False, Write = False Append = True
    def testConstructorForInitialConstructionEFFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionEFFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",False,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = False, Write = False Append = True
    def testConstructorForInitialConstructionBFFT(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBFFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",False,False,True)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = True, Write = False Append = False
    def testConstructorForInitialConstructionNTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,True,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = True, Write = False Append = False
    def testConstructorForInitialConstructionETFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionETFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",True,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = True, Write = False Append = False
    def testConstructorForInitialConstructionBTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",True,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = False, Write = True Append = False
    def testConstructorForInitialConstructionNFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,False,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = False, Write = True Append = False
    def testConstructorForInitialConstructionEFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionEFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",False,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = False, Write = True Append = False
    def testConstructorForInitialConstructionBFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",False,True,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = None Read = False, Write = False Append = False
    def testConstructorForInitialConstructionNFFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionNFFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(None,False,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "" Read = False, Write = False Append = False
    def testConstructorForInitialConstructionEFFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionEFFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("",False,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = "    " Read = False, Write = False Append = False
    def testConstructorForInitialConstructionBFFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionBFFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO("    ",False,False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = valid Read = list(), Write = False Append = False
    def testConstructorForInitialConstructionVIFF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVIFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,list(),False,False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = list() Append = False
    def testConstructorForInitialConstructionVFIF(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFIF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,False,list(),False)
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")

    ##
    #Test the constructor for initialization FileName = valid Read = False, Write = False Append = list()
    def testConstructorForInitialConstructionVFFI(self):
        #Method Name
        self.methodName = "FileIOTest.testConstructorForInitialConstructionVFFI"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,False,False,list())
        except ValueError, ve:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.assertTrue(True,self.methodName+" Error occured passed.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            #Pass if no error is thrown
            self.fail(self.methodName+" Failed, an error should have been thrown.")


    ##
    #Test the get file name FileName = valid Read = True, Write = True Append = True
    def testGetFileNameForValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testGetFileNameForValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Actual Answer
        self.actualAnswer = self.testedFileIO.getFileName()

        #Correct answer
        self.correctAnswer = self.TestFileName

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the string function FileName = valid Read = True, Write = True Append = True
    def testStringForValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testStringForValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Actual Answer
        self.actualAnswer = str(self.testedFileIO)
        #Correct answer
        self.correctAnswer = "FileIO:: File name ="+self.TestFileName+" ReadMode:True WriteMode:True File Object:<open file 'FileIOTestWriteFile.txt', mode 'a+' at "

        #Close FileIO
        if not self.testedFileIO == None:
            self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer[0:len(self.correctAnswer)],self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the close function FileName = valid Read = True, Write = True Append = True
    def testCloseForValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testCloseForValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.correctAnswer = True

        #Actual Answer
        #Close FileIO
        if not self.testedFileIO == None:
            self.actualAnswer = self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the multiple close function FileName = valid Read = True, Write = True Append = True
    def testMultCloseForValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testCloseForValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.correctAnswer = True

        #Actual Answer
        #Close FileIO
        if not self.testedFileIO == None:
            self.actualAnswer = self.testedFileIO.close()
            self.actualAnswer = self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the none close function FileName = valid Read = True, Write = True Append = True
    def testNoneCloseForValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testCloseForValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.correctAnswer = None

        #Actual Answer
        #Close FileIO
        if not self.testedFileIO == None:
            self.actualAnswer = self.testedFileIO.close()
            self.testedFileIO._FileIO__file = None
            self.actualAnswer = self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")


    ##
    #Test the close function FileName = valid Read = False, Write = True Append = False
    def testCloseForValidFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testCloseForValidFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,False)

        #Correct answer
        self.correctAnswer = True

        #Actual Answer
        #Close FileIO
        if not self.testedFileIO == None:
            self.actualAnswer = self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the close function FileName = valid Read = True, Write = False Append = False
    def testCloseForValidTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testCloseForValidTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)

        #Correct answer
        self.correctAnswer = True

        #Actual Answer
        #Close FileIO
        if not self.testedFileIO == None:
            self.actualAnswer = self.testedFileIO.close()
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)

        #Pass if no error is thrown
        self.assertEqual(self.actualAnswer,self.correctAnswer,self.methodName+" Should have received:"+str(self.correctAnswer)+". but received:"+str(self.actualAnswer)+".")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = True Append = True
    def testWriteToFileForOneLineValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = False Append = True
    def testWriteToFileForOneLineValidTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.fail(self.methodName+" Did create error.")

        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(not returnValue,self.methodName+" Did not create error but returned that it did write when it should not have.")
    ##
    #Test the writeToFile function FileName = valid Read = True, Write = True Append = False
    def testWriteToFileForOneLineValidTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = False, Write = True Append = False
    def testWriteToFileForOneLineValidFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Did not create error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = False Append = False
    def testWriteToFileForOneLineValidTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.fail(self.methodName+" Did create error.")

        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(not returnValue,self.methodName+" Did not create error but returned that it did write when it should not have.")

    ##
    #Test the writeToFile function FileName = valid Read = False, Write = True Append = True
    def testWriteToFileForOneLineValidFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneLineValidFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False, True, True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.fail(self.methodName+" Did create error.")

        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = True Append = True
    def testWriteToFileForThreeLineValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = False Append = True
    def testWriteToFileForThreeLineValidTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.fail(self.methodName+" Did create error.")

        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(not returnValue,self.methodName+" Did not create error but returned that it did write when it should not have.")
    ##
    #Test the writeToFile function FileName = valid Read = True, Write = True Append = False
    def testWriteToFileForThreeLineValidTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = False, Write = True Append = False
    def testWriteToFileForThreeLineValidFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Did create error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the writeToFile function FileName = valid Read = True, Write = False Append = False
    def testWriteToFileForThreeLineValidTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w').close()
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.fail(self.methodName+" Did not create error.")

        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)
            self.assertTrue(not returnValue,self.methodName+" Did not create error but returned that it did write when it should not have.")

    ##
    #Test the writeToFile function FileName = valid Read = False, Write = True Append = True
    def testWriteToFileForOneNoneLineValidFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForOneNoneLineValidFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,True)

        #Correct answer
        self.lineContents = None
        returnValue = None
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)

        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Did create error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(not returnValue,self.methodName+" Did not create error but returned that it did write when it should not have.")

    ##
    #Test the writeToFile function FileName = valid Read = False, Write = True Append = True
    def testWriteToFileForThreeLineValidFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileForThreeLineValidFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = False
        #Attempt method
        try:
            returnValue = self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
            returnValue = returnValue and self.testedFileIO.writeToFile(self.lineContents)
        except (IOError, ValueError):
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Did create error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(returnValue,self.methodName+" Did not create error but returned that it did not write.")

    ##
    #Test the readFullFile function FileName = valid Read = True, Write = True Append = True
    def testReadFullFileForThreeLineValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test."
        returnValue = None
        self.testedFileIO.writeToFile(self.lineContents)
        self.testedFileIO.writeToFile(self.lineContents)
        self.testedFileIO.writeToFile(self.lineContents)

        #Attempt method
        try:
            returnValue = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(returnValue,self.lineContents+self.lineContents+self.lineContents,self.methodName+" Did not create error.")


    ##
    #Test the readFullFile function FileName = valid Read = True, Write = True Append = True
    def testReadFullFileForNoFileValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForNoFileValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Attempt method
        try:
            returnValue = self.testedFileIO.readFullFile()

        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertTrue(True,self.methodName+" Created an error.")

    ##
    #Test the readFullFile function FileName = valid Read = True, Write = True Append = True
    def testReadFullFileForNoneValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForNoneValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Attempt method
        try:
            self.testedFileIO._FileIO__file = None
            returnValue = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False, returnValue, self.methodName+" Created an error.")

    ##
    #Test the readFullFile function FileName = valid Read = True, Write = True Append = False
    def testReadFullFileForThreeLineValidTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidTTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        self.testedFileIO.writeToFile(self.lineContents)
        self.testedFileIO.writeToFile(self.lineContents)
        self.testedFileIO.writeToFile(self.lineContents)

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(self.actualAnswer,self.lineContents+self.lineContents+self.lineContents,self.methodName+" Did not create error.")

    ##
    #Test the readFullFile function FileName = valid Read = True, Write = False Append = True
    def testReadFullFileForThreeLineValidTFT(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidTFT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,True)

        #Correct answer
        self.lineContents = "This is a test."
        writer = open(self.TestFileName, 'w')
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.close()

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(self.actualAnswer,self.lineContents+self.lineContents+self.lineContents,self.methodName+" Did not create error.")

    ##
    #Test the readFullFile function FileName = valid Read = True, Write = False Append = False
    def testReadFullFileForThreeLineValidTFF(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidTFF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)

        #Correct answer
        self.lineContents = "This is a test."
        writer = open(self.TestFileName, 'w')
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.close()

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(self.actualAnswer,self.lineContents+self.lineContents+self.lineContents,self.methodName+" Did not create error.")

    #Test the readFullFile function FileName = valid Read = False, Write = True Append = False
    def testReadFullFileForThreeLineValidFTF(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidFTF"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,False)

        #Correct answer
        self.lineContents = "This is a test."
        writer = open(self.TestFileName, 'w')
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.close()

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not create error.")

    #Test the readFullFile function FileName = valid Read = False, Write = True Append = True
    def testReadFullFileForThreeLineValidFTT(self):
        #Method Name
        self.methodName = "FileIOTest.testReadFullFileForThreeLineValidFTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')
        self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,True)

        #Correct answer
        self.lineContents = "This is a test."
        writer = open(self.TestFileName, 'w')
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.write(self.lineContents)
        writer.close()

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not create error.")

    #Test the readline function FileName = valid Read = True, Write = False Append = False
    def testReadline1linesTFF1(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline1linesTFF1"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 1\n"
        writer = open(self.TestFileName, 'w')
        writer.write(answer)
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.actualAnswer = str(self.testedFileIO.readline())
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(answer,self.actualAnswer,self.methodName+" Did not create error. received:"+self.actualAnswer+": instead of:"+answer+".")

    #Test the readline function FileName = valid Read = True, Write = False Append = False
    def testReadline3linesTFF1(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesTFF1"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 1\n"
        writer = open(self.TestFileName, 'w')
        writer.write(answer)
        writer.write("line 2 	a b 	c\n")
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.actualAnswer = str(self.testedFileIO.readline())
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(answer,self.actualAnswer,self.methodName+" Did not create error. received:"+self.actualAnswer+": instead of:"+answer+".")

    #Test the readline function FileName = valid Read = True, Write = False Append = False
    def testReadline3linesTFF2(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesTFF2"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 2 	a b 	c\n"
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write(answer)
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.testedFileIO.readline()
            self.actualAnswer = str(self.testedFileIO.readline())
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(answer,self.actualAnswer,self.methodName+" Did not create error. received:"+self.actualAnswer+": instead of:"+answer+".")

    #Test the readline function after closed file FileName = valid Read = True, Write = False Append = False
    def testReadline3linesClosedTFF2(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesClosedTFF2"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 2 	a b 	c\n"
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write(answer)
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.testedFileIO.readline()
            self.testedFileIO.close()
            self.actualAnswer = self.testedFileIO.readline()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not create error. received:"+str(self.actualAnswer)+": instead of:False.")

    #Test the readline function from none file FileName = valid Read = True, Write = False Append = False
    def testReadline3linesNoneTFF2(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesNoneTFF2"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 2 	a b 	c\n"
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write(answer)
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.testedFileIO.readline()
            self.testedFileIO._FileIO__file = None
            self.actualAnswer = self.testedFileIO.readline()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not create error. received:"+str(self.actualAnswer)+": instead of:False.")


    #Test the readline function FileName = valid Read = True, Write = False Append = False
    def testReadline3linesTFF3(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesTFF3"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\""
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write("line 2 	a b 	c\n")
        writer.write(answer)
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,False,False)
            self.testedFileIO.readline()
            self.testedFileIO.readline()
            self.actualAnswer = str(self.testedFileIO.readline())
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(answer,self.actualAnswer,self.methodName+" Did not create error. received:"+self.actualAnswer+": instead of:"+answer+".")

    #Test the readline function FileName = valid Read = True, Write = True Append = True
    def testReadline3linesTTT2(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesTTT2"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = "line 2 	a b 	c\n"
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write(answer)
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)
            self.testedFileIO.readline()
            self.actualAnswer = str(self.testedFileIO.readline())
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(answer,self.actualAnswer,self.methodName+" Did not create error. received:"+self.actualAnswer+": instead of:"+answer+".")

    #Test the readline function FileName = valid Read = False, Write = True Append = True
    def testReadline3linesFTT2(self):
        #Method Name
        self.methodName = "FileIOTest.testReadline3linesFTT2"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        open(self.TestFileName, 'w')

        #Correct answer
        answer = False
        writer = open(self.TestFileName, 'w')
        writer.write("line 1\n")
        writer.write("line 2 	a b 	c\n")
        writer.write("line 3 hhgfhyrjumfdgtyrtuoupl,mn;.,m;oij\"")
        writer.close()

        #Attempt method
        try:
            self.testedFileIO = FileIO.FileIO(self.TestFileName,False,True,True)
            self.testedFileIO.readline()
            self.actualAnswer = self.testedFileIO.readline()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(str(answer),str(self.actualAnswer),self.methodName+" Did not create error. received:"+str(self.actualAnswer)+": instead of:"+str(answer)+".")

    ##
    #Test the WriteToFile and readFullFile function with one line FileName = valid Read = True, Write = True Append = True
    def testWriteToFileAndReadFullFileForOneLineValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileAndReadFullFileForOneLineValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test.This is a test.This is a test."
        self.correctAnswer = self.lineContents

        #Attempt method
        try:
            self.testedFileIO.writeToFile(self.lineContents)
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(self.correctAnswer,self.actualAnswer,self.methodName+" Did not create error but received:"+str(self.actualAnswer)+" when it should have been "+str(self.correctAnswer))

    ##
    #Test the WriteToFile and readFullFile function from no file FileName = valid Read = True, Write = True Append = True
    def testWriteToFileAndReadFullFileForNoFileValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileAndReadFullFileForNoFileValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test.This is a test.This is a test."
        self.correctAnswer = self.lineContents

        #Attempt method
        try:
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual("",self.actualAnswer,self.methodName+" Did not indicate a blank read")

    ##
    #Test the WriteToFile and readFullFile function with none FileName = valid Read = True, Write = True Append = True
    def testWriteToFileAndReadFullFileForNoneValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileAndReadFullFileForNoneValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test.This is a test.This is a test."
        self.correctAnswer = self.lineContents

        #Attempt method
        try:
            self.testedFileIO._FileIO__file = None
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not indicate a false read")

    ##
    #Test the WriteToFile and readFullFile function on a closed file FileName = valid Read = True, Write = True Append = True
    def testWriteToFileAndReadFullFileForClosedValidTTT(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileAndReadFullFileForClosedValidTTT"

        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,True)

        #Correct answer
        self.lineContents = "This is a test.This is a test.This is a test."
        self.correctAnswer = self.lineContents

        #Attempt method
        try:
            self.testedFileIO.close()
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(False,self.actualAnswer,self.methodName+" Did not indicate a false read")

    ##
    #Test the WriteToFile and readFullFile function with three line FileName = valid Read = True, Write = True Append = False
    def testWriteToFileAndReadFullFileForThreeLineValidTTF(self):
        #Method Name
        self.methodName = "FileIOTest.testWriteToFileAndReadFullFileForThreeLineValidTTF"
        self.TestFileName = "trythis.txt"
        #CreateObject
        if os.path.isfile(self.TestFileName):
            os.remove(self.TestFileName)
        self.testedFileIO = FileIO.FileIO(self.TestFileName,True,True,False)

        #Correct answer
        self.lineContents = "This is a test.This is a test.This is a test."
        self.correctAnswer = self.lineContents

        #Attempt method
        try:
            self.testedFileIO.writeToFile(self.lineContents)
            self.actualAnswer = self.testedFileIO.readFullFile()
        except IOError:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.fail(self.methodName+" Created error.")
        else:
            #Close FileIO
            if not self.testedFileIO == None:
                self.testedFileIO.close()
            if os.path.isfile(self.TestFileName):
                os.remove(self.TestFileName)

            self.assertEqual(self.correctAnswer,self.actualAnswer,self.methodName+" Did not create error but received:"+str(self.actualAnswer)+" when it should have been "+str(self.correctAnswer))

##
#Creates a suite of tests
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(FileIOTest)
    return suite
