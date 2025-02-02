#!/usr/bin/env python3


import os
import argparse


def blastnNamefix(inputpath, querypath):
	
	querylist = []
	with open(querypath, mode = 'r') as f:
		for line in f:
			if len(line) and line[0] == ">":
				querylist.append(line.split()[0].strip()[1:])
	
	ifwrite = 0 
	
	wfile = open(inputpath+"_fixtemp", mode = 'w')
	
	with open(inputpath, mode = "r") as f:
		
		for line in f:
			
			if len(line) == 0:
				continue
			
			line = line.split()
			
			queryname = line[2]
			queryname = querylist[int(queryname.split('_')[-1])-1]
			
			line[2] = queryname 
			
			wfile.write("\t".join(line)+"\n")
				
	wfile.close()
	
	os.system("mv {} {}".format(inputpath+"_fixtemp", inputpath))
			
			
			
def main(args):
	
	inputpath = args.input
	querypath = args.query
	
	blastnNamefix(inputpath, querypath)
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
