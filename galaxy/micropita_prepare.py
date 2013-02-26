#!/usr/bin/env python

"""
Author: George Weingart
Description: Prepare parameters to call micropita
"""

#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################

__author__ = "George Weingart"
__copyright__ = "Copyright 2012"
__credits__ = ["George Weingart"]
__license__ = "MIT"
__maintainer__ = "George Weingart"
__email__ = "george.weingart@gmail.com"
__status__ = "Development"

import argparse
from cStringIO import StringIO
import sys,string,time
import os
from time import gmtime, strftime
from pprint import pprint

##################################################################################
#  Decode Parms                                                                  #
##################################################################################
def read_params(x):
	parser = argparse.ArgumentParser(description='Micropita Annotate Argparser')
	parser.add_argument('--input', action="store",dest='inputname')
	parser.add_argument('--output', action="store",dest='outputname')
	parser.add_argument('-m', action="store",dest='MParameter')
	parser.add_argument('-n', action="store",dest='NSamples')
	parser.add_argument('--lastmeta', action="store",dest='lastmeta')
	parser.add_argument('--stratify_value', action="store",dest='stratify_value')

	try:
		parser.add_argument('--feature_method', action="store",dest='feature_method')
	except:
		pass
	try:
		parser.add_argument('--targets', action="store",dest='targets')
	except:
		pass
	try:
		parser.add_argument('--label_value', action="store",dest='label_value')
	except:
		pass

	

	return  parser
	

##################################################################################
#  Main Program                                                                  #
##################################################################################
parser = read_params( sys.argv )
results = parser.parse_args()

fname =  results.inputname
input_file = open(fname,'rU')
input_lines = input_file.readlines()
input_file.close()
table_lines = []
for x in input_lines:
	first_column = x.split('\t')[0]
	table_lines.append(first_column)



FileTimeStamp =  strftime("%Y%m%d%H%M%S", gmtime())
LastMetaInt = 0
if results.lastmeta and not results.lastmeta == "None":
	LastMetaInt = int(results.lastmeta) - 1

StratifyValueInt = 0
if  results.stratify_value and not   results.stratify_value == "None":
	StratifyValueInt = int(results.stratify_value) - 2 

LabelValueInt = 0
if results.label_value and not results.label_value == "None":
	LabelValueInt = int(results.label_value) - 1

stratify_string = ""
q = '"'
if  not results.stratify_value == '1':
	stratify_string = " --stratify " + q + table_lines[StratifyValueInt] + q + " "

if results.MParameter == "features":
	TBTargets = list()
	TableTargets = results.targets.split(',')
	for t in TableTargets:
		tb_entry = int(t) + LastMetaInt 
		TBTargets.append(int(tb_entry))


	TempTargetsFileName = "/tmp/micropita_targets" + FileTimeStamp
	OutTargetsFile = open(TempTargetsFileName,"w")
	indx = -1
	for  c in table_lines:
		indx+=1
		if  indx in TBTargets:
			OutputString = table_lines[indx] + "\n"
			OutTargetsFile.write(OutputString)
	OutTargetsFile.close()
	os_command = "python " + \
        "/usr/local/galaxy-dist/tools/micropita/" + \
	"MicroPITA.py "+\
	"--lastmeta " + table_lines[LastMetaInt]+ " " +\
        "--feature_method " + results.feature_method + " " + \
	"--target " + TempTargetsFileName + " " +\
	"-m " + results.MParameter + " " + \
	"-n " + results.NSamples + " " +\
	stratify_string + " " +\
	results.inputname + " " +\
	results.outputname
	#print os_command
	os.system(os_command)

if results.MParameter == "representative"\
or results.MParameter == "diverse"\
or results.MParameter == "extreme": 
	os_command = "python " + \
        "/usr/local/galaxy-dist/tools/micropita/" + \
	"MicroPITA.py "+\
	"--lastmeta " + table_lines[LastMetaInt]+ " " +\
	"-m " + results.MParameter + " " + \
	"-n " + results.NSamples + " " +\
        stratify_string + " " + \
	results.inputname + " " +\
	results.outputname
	#print os_command
	os.system(os_command)

if results.MParameter == "distinct"\
or results.MParameter == "discriminant": 
	os_command = "python " + \
        "/usr/local/galaxy-dist/tools/micropita/" + \
	"MicroPITA.py "+\
	"--lastmeta " + table_lines[LastMetaInt]+ " " +\
	"--label " + table_lines[LastMetaInt]+ " " +\
	"-m " + results.MParameter + " " + \
	"-n " + results.NSamples + " " +\
        stratify_string + " " + \
	results.inputname + " " +\
	results.outputname
	#print os_command
	os.system(os_command)
