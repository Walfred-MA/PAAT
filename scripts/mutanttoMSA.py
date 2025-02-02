#!/usr/bin/env python3

import re
import collections as cl
import pandas as pd
import numpy as np
import os 
import argparse

def sizefiltration(cigar):
	
	cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
	
	newcigars = []
	for cigarstr in cigars:
		
		thename,thecigar = cigarstr.split(":")
		
		thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
		
		if thesize >0 and thesize < 200:
			
			fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
			
			thecigar = "{}H".format(fullsize )
			
		newcigars.append(thename+":"+thecigar)
		
	newcigars= "".join(newcigars)
	
	return newcigars

def cigartoseq(index, cigarstr, refseq, globalposi, allinserts):
	
	segments = re.findall(r'\d+[HMX=ID][ACTGNactgn]*', cigarstr)

	queryseq = ""
	posi = 0
	for segment in segments:
		
		size = re.findall(r'\d+', segment)[0]
		rest = segment[len(size):]
		size = int(size)
		thetype = rest[0]
		theseq = rest[1:]
		
		if thetype == 'M':
			
			queryseq += refseq[posi:(posi+size)]
			posi += size
		
		elif thetype == 'H' or thetype == 'D':
			
			queryseq += '-' * size
			posi += size
		
		elif thetype == 'X':
			
			queryseq += theseq
			posi += size
		
		elif thetype == 'I':
			
			allinserts[globalposi + posi][index]= theseq
	
	return queryseq, globalposi + posi, allinserts


def toMSA(inputfile, outputfile):

	highlight = {}	
	with open(inputfile, mode = 'r') as f:
		
		lines = [ line.split('\t') for line in f.read().splitlines() if len(line) and line[0] == 'L'] 
		
	table = pd.DataFrame.from_records(lines)
	
	allseginfo = dict()
	pathsizes = {}
	allpathes = {}
	allgenes = []
	with open(inputfile, mode = 'r') as f:
		
			for line in f:
				
				if len(line.strip()) == 0:
					continue
				
				line = line.strip().split('\t')
				if line[0][0] == "S":
					name = line[1][2:]
					allseginfo[name] = line[3:]
					allpathes[name] = line[-4]
					pathsizes[name] = int([x for x in line if x.startswith("LN:i:")][0].split(':')[-1])
				elif line[0][0] == "L":
					allgenes.append(line)
	

	names = list(table[1])
	cigars = list(table[5])
	sizes = [int(x.split(':')[-1]) for x in list(table[6])]
	
	sizes_used = []
	elements = []
	counter = 0 
	alltypes = []
	
	allnames = []
	allinserts = cl.defaultdict(lambda :cl.defaultdict(str))
	allfasta = []
	for i,(genename,cigar) in enumerate(zip(names,cigars)):
		
		cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
		
		#cigars = [cigar for cigar in cigars if "_" not in cigar]
				
		globalposi = 0
		fasta = ""
		for cigar in cigars:
			
			name,cigar = cigar.split(":")
			
			name = name[1:]
			
			
			segmentseq,globalposi,allinserts  = cigartoseq(i,cigar, allpathes[name.split('_')[0]], globalposi, allinserts)
			
			fasta += segmentseq
		
		allnames.append(genename)
		allfasta.append(list(fasta))
		
	allinserts = sorted([(posi, insertrecord) for posi, insertrecord in allinserts.items()])
	
	for index, (posi, insertrecord) in enumerate(allinserts):
		
		insertforsample = []
		largestsize = max([len(x) for x in insertrecord.values()])
		
		for i, fasta in enumerate(allfasta):
			
			insert = insertrecord[i] + '-'*(largestsize - len(insertrecord[i]))
			allfasta[i].insert(posi+index, insert)
		
	for i, fasta in enumerate(allfasta):
		
		allfasta[i] = "".join(allfasta[i])
		
	
	with open(outputfile, mode = 'w') as f:
		
		for index, (genename, fasta) in enumerate(zip(allnames, allfasta)):
			
			f.write(">{}\n{}\n".format(genename, fasta))
	
	
	
			
	
def main(args):
	
	toMSA(args.input, args.output)
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
