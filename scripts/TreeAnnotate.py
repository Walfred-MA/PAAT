#!/usr/bin/env python3

import os 
import argparse
import sys
import numpy as np
import math
import collections as cl

	
def normtotdm(matrix):
	
	
	matrix = np.matrix(matrix)
	
	dm = np.ones((matrix.shape))
	vars = [matrix[i,i] for i in range(len(matrix))]
	for i in range(len(matrix)):
		
		var1 = max(1.0,vars[i])
		for j in range(i+1):
			
			var2 = max(1.0,vars[j])
			
			dm [i,j] = 1 - matrix[i,j]/math.sqrt((var1*var2))
			dm [j,i] =  dm [i,j]
			
	return dm


def main(args):
	
	
	inputfile = args.input
	outputfile = args.output
	
	normfile = args.input + "_norm.txt"
	treefile = args.input + "_tree.ph"
			
	normmatrix = np.genfromtxt(normfile, delimiter=',')
	
	dm = normtotdm(normmatrix)
	
	header_info = dict()	
	headers = dict()
	with open(args.input, mode = 'r') as f:
		for line in f:
			if len(line) ==0 or line[0] != ">":
				continue
			lines = line.strip().split()
			headers[lines[0][1:]] = lines[1]
			header_info[lines[0][1:]] = lines[-1]
	nonrefindexes = []
	refindexes = []
	for i, header in enumerate(headers.keys()):
		
		if "chr" in header and "v" not in header:
			refindexes.append(i)
		else:
			nonrefindexes.append(i)

	headers_list = list(headers.keys())
	
	allpairs = []
	for i in nonrefindexes:
		
		for j in refindexes:
			
			allpairs.append( ( dm[i,j], i,j) )
		
	allpairs = sorted(allpairs)
	
	finish_haplos = set()
	best_matches = {}
	
	for (distance, i, j) in allpairs:
	
		haplo1 = "_".join(headers_list[i].split('_')[:-1])
		haplo2 = headers_list[j]
		
		if distance > 0.2:
			break
		
		if headers_list[i] in best_matches or (haplo1,haplo2) in finish_haplos:
			continue
		
		best_matches[headers_list[i]] = (distance, haplo2 , 0)
		
		finish_haplos.add((haplo1,haplo2))
		
	for (distance, i, j) in allpairs:
		
		if distance > 0.2:
			break
		
		if headers_list[i] in best_matches:
			continue
		
		best_matches[headers_list[i]] = (distance, headers_list[j], 1)
		
	tags = ["Pri","Dup"]
	
	with open(outputfile, mode = 'w') as f:
		
		for header in headers:
		
			if header in best_matches:
		
				distance, ref, iffirst = best_matches[header]
			
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, ref, headers[header], headers[ref], int("ENSE" in header_info[header]) , distance,tags[iffirst]))
			
			else:
				if "chr" not in header:	
					f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, "NA", headers[header], "NA", int("ENSE" in header_info[header]),0.0,"Novel"))
				else:
					f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, "NA", headers[header], "NA", int("ENSE" in header_info[header]),0.0,"Ref"))
	
	
	
	return 0

def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine genes")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
