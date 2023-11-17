#!/usr/bin/env python3

import os
import argparse
import collections as cl


def filesplit(inputpath, outpath, ifsample=1):
	
	wfiles= {}
	wfile =None	
	with open(inputpath, mode = "r") as f:
		
		for line in f:
			
			if len(line) == 0:
				continue
			
			if line[0] == ">":
				
				header = line[1:].strip()
				name = header.split()[0].replace("#","_h").replace(":","_").replace("-","_").replace(".","v")
		
				if ifsample == 0:
					
					if wfile is not None:
						wfile.close()
					
					wfile = open(outpath+name+".fa",mode = 'w') 
		
				elif ifsample == 1:
					
					haplo = "_h".join(header.split()[1].split("#")[:2])
					if ":" in haplo:
						if  "NC_0609" in haplo:
							haplo = "CHM13_h1"
						else:
							haplo = "HG38_h1"
					
					wfilename = outpath+haplo+".fa"
					
					if wfilename not in wfiles:
						
						wfiles[wfilename] = open(outpath+haplo+".fa",mode = 'w')
					
					wfile = wfiles[wfilename]
					
			if wfile is not None:
				wfile.write(line)
	
	if ifsample == 0:
		
		if wfile is not None:
			wfile.close()
					
	elif ifsample == 1:
	
		for wfilename, wfile in wfiles.items():
		
			wfile.close()
		
					
def main(args):
	
	inputpath = args.input
	outpath = args.output
	ifsample = args.ifsample
	
	filesplit(inputpath, outpath, ifsample)
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str,required=True)
	parser.add_argument("-s", "--ifsample", help="path to output file", dest="ifsample",
						type=int,default = 1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
