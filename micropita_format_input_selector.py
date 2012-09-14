#!/usr/bin/env python
import sys,string,time
from pprint import pprint





def red(st,l):
	if len(st) <= l: return st 
	l1,l2 = l/2,l/2
	return st[:l1]+".."+st[len(st)-l2:]


def get_cols(data,full_names):
	if data == "": return []
	max_len =32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)

	opt = []
	rc = ''
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines



def get_cols_add_line(data,full_names,lastmeta):
	if data == "": return []
	display_to = 1
	try:
		display_to = int(lastmeta)
	except:		
		pass

	max_len = 32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)
	table_lines.insert(0,'-')
	if  not display_to == 1:
		del  table_lines[display_to + 1:]


	opt = []
	rc = ''
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines


def get_cols_features(data,full_names,lastmeta):
	if data == "": return []
	display_from = 1
	try:
		display_from = int(lastmeta)
	except:		
		pass
	max_len = 32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)
	
	opt = []
	rc = ''
	del table_lines[:display_from]
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines





