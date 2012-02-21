#######################################################
#
#	Title:		FileIO
#	Author:		Timothy Tickle
#	Date:		01/28/2010
#	Purpose:	Class to read and write text data to file
#
#######################################################

#Import libaries
from ValidateData import ValidateData
import os

##
#Used to read and write text file data.
class FileIO():

    ##Current file name
    __fileName = None

    ##Current file
    __file = None

    ##Indicator if the file is for writing
    __write = None

    ##Indicator if the file is for reading
    __read = None

    ##
    #Contructor
    #Created 1/29/2010
    #Tested 1/30/2010
    #Reviewed 5/20/10
    #@param fileName File to read or write to
    #@param tempRead Boolean indicator of if data will be read to the file
    #@param tempWrite Boolean indicator of if data will be write to the file
    #@param tempFileAppend Boolean indicator of if data will be appended when writing
    def __init__(self, fileName, tempRead, tempWrite, tempFileAppend):
        #Check parameters
        if not ValidateData.isValidString(fileName):
            raise ValueError, "Error while opening file. File name is invalid :"+str(fileName)+"."
        if not ValidateData.isValidBoolean(tempRead):
            raise ValueError, "Error while opening file. Read file is invalid :"+str(tempRead)+"."
        if not ValidateData.isValidBoolean(tempWrite):
            raise ValueError, "Error while opening file. Write file is invalid :"+str(tempWrite)+"."
        if not ValidateData.isValidBoolean(tempFileAppend):
            raise ValueError, "Error while opening file. File Append is invalid :"+str(tempFileAppend)+"."

        #Set variables
        self.__fileName = fileName
        self.__read = tempRead
        self.__write = tempWrite

        #Attempt open
        try:
            if tempRead:
                if tempWrite:
                    if tempFileAppend:
                        self.__file = open(fileName, 'a+')
                        return
                    else:
                        self.__file = open(fileName, 'w+')
                        return
                else:
                    self.__file = open(fileName, 'r')
                    return
            #Attempt write
            if tempWrite:
                    if tempFileAppend:
                        self.__file = open(fileName, 'a')
                        return
                    else:
                        self.__file = open(fileName, 'w')
                        return
            if (not tempRead) and (not tempWrite):
                raise IOError, "File IO was not read or write, no reason to create it to interact with "+fileName
        except IOError, e:
            raise IOError, "Error while opening file :"+fileName+"."+str(e)

    ##
    #Writes data to a file given the FileIO has write capacity
    #Does not append a "\n" to end of line, must add it
    #Created 1/29/2010
    #Tested 1/30/2010
    #Reviewed 5/20/10
    #@param tempData Data to write
    #@return boolean indicator of write, None is given on bad parameter
    def writeToFile(self, tempData):
        if not self.__write:
            return False
        if (not ValidateData.isValidString(tempData)) or (self.__file == None):
            return False

        #Write to file
        self.__file.write(tempData)

        return True

    ##
    #Reads the contents of a file into one string
    #Created 1/29/2010
    #Tested 1/31/2010
    #Reviewed 5/20/10
    #@params returns the file contents
    def readFullFile(self):
        if not self.__read:
            return False
        if self.__file == None:
            return False
        if self.__file.closed:
            return False
        self.__file.seek(0,0)
        return self.__file.read()

    ##
    #Reads the contents of one line of a file into one string
    #Created 2/07/2010
    #Tested 2/08/2010
    #Reviewed 5/20/10
    #@params returns the file line contents
    def readline(self):
        if not self.__read:
            return False
        if self.__file == None:
            return False
        if self.__file.closed:
            return False
        return self.__file.readline()

    ##
    #Returns the current file name
    #Created 1/29/2010
    #Tested 1/30/2010
    #Reviewed 5/20/10
    #@return string fileName
    def getFileName(self):
        return self.__fileName

    ##
    #Closes the file
    #Created 1/29/2010
    #Tested 1/30/2010
    #Reviewed 5/20/10
    #@return Indicator of file closing
    def close(self):
        if not self.__file == None:
          #Close file
          self.__file.close()
          return self.__file.closed
        return None

    ##
    #Returns a string representation of the object
    #Created 1/29/2010
    #Tested 1/30/2010
    #Reviewed 5/20/10
    def __str__(self):
        return "FileIO:: File name ="+str(self.__fileName)+" ReadMode:"+str(self.__read)+" WriteMode:"+str(self.__write)+" File Object:"+str(self.__file)
