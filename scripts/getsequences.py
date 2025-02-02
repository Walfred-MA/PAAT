#!/usr/bin/env python3


import os
import argparse
import collections as cl


def filesplit(inputpath, outpath, samplelist):
	
	outputset = set()
	if len(samplelist):
		
		with open(samplelist, mode = 'r') as f:
			for line in f:
				outputset.add(line.split()[0].strip())


	ifwrite = 0 
	
	wfile = open(outpath, mode = 'w')
	
	with open(inputpath, mode = "r") as f:
		
		for line in f:
			
			if len(line) == 0:
				continue
			
			if line[0] == ">":
				
				header = line[1:].strip()
				name = header.split()[0]
				if name in outputset:
						
					ifwrite = 1
				else:
					ifwrite = 0
					
			if ifwrite:
				wfile.write(line)
				
	wfile.close()
			
			
			
def main(args):
	
	inputpath = args.input
	outpath = args.output
	samplelist = args.samplelist
	
	filesplit(inputpath, outpath, samplelist)
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	parser.add_argument("-l", "--samplelist", help="path to output file", dest="samplelist", type=str,default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
