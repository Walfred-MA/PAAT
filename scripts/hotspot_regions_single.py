#!/usr/bin/env python3

import pandas as pd
import os 
import re
import argparse
import sys
import multiprocessing as mul
import collections as cl

def querycontig_overlap(allaligns, allow_gap = 0):
	
	all_coordinates = [int(x) for y in allaligns for x in y[:2]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group = set([])
	current_group_size = 0 
	current_group_start = 0
	allgroups_sizes = []
	
	number_curr_seq = 0
	coordinate = 0
	
	last_coordinate = -allow_gap-1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if index %2 == 0 :
			
			number_curr_seq += 1
			
			if number_curr_seq == 1:
				
				current_group_start = int(coordinate)
				
				if coordinate - last_coordinate>allow_gap :
					
					allgroups_sizes.append(current_group_size)
					allgroups.append(current_group)
					
					current_group = set([])
					current_group_size = 0
					
					
			current_group.add(index//2)
			
		else:
			
			number_curr_seq -= 1
			
			if number_curr_seq == 0:
				
				current_group_size += coordinate - current_group_start			
				
		last_coordinate = coordinate
		
		
	allgroups_sizes.append(current_group_size)	
	allgroups.append(current_group)
	
	allgroups_sizes = allgroups_sizes[1:]	
	allgroups = allgroups[1:]
	
	results = []
	for i,group in enumerate(allgroups):
		
		aligns = [allaligns[x] for x in group]
		
		coordinates = [x for y in aligns for x in y[:2]]
		start = min(coordinates)
		end = max(coordinates)
	
		start = max(0,start-allow_gap//2)
		end = end+allow_gap//2

		results.append([start, end])
		
		
	return results


def main(args):
	
	inputpath, querypath ,outputpath= args.input, args.query , args.output
	

	if querypath[-3:]==".fa"  or querypath[-6:]==".fasta":
		querypaths = querypath

	else:
		with open(querypath, mode = 'r') as f:
		
			querypaths = {x.split()[0]:x.split()[1] for x in f.read().splitlines()}

	try:
		table = pd.read_csv(inputpath, header = None , sep = "\t")
	except:
		print("warnning, no kmer found in ", inputpath)
		exit(0)
	ncol = len(table.columns.values.tolist())
	
	eachcontig = cl.defaultdict(list)
	for index, row in table.iterrows():
		
		contig = row.iat[0].strip().split()[0]
		
		eachcontig[contig].append(index)

	
	allqueryloci = []
	
	for contig, indexes in eachcontig.items():
		
		subdata = table.iloc[indexes]
		
		usecols = subdata[[ncol-2,ncol-1]].values.tolist()
		
		loci = querycontig_overlap(usecols, 20000)
	
		loci = [[contig]+x for x in loci]
		
		allqueryloci.extend(loci)

	try:
		os.remove(outputpath)	
	except:
		pass

	for i,raw in enumerate(allqueryloci):
		
		contig, start, end = raw
		
		samplename = "_h".join(contig.split("#")[:2])	
		
		if type(querypaths) == type({}):
			queryfile = querypaths[samplename]
		else:
			queryfile = querypaths
		
		title = ">{}_{}_{}\t{}:{}-{}".format(contig, start, end, contig, start, end)
		
		cm = "echo \"{}\"  >> {}".format(title, outputpath)
		os.system(cm)
		
		cm = "samtools faidx {} {}:{}-{}  | grep -v '>' >> {}".format(queryfile,contig, start+1, end+2, outputpath)

		os.system(cm)
		
		cm = "echo $\'\\n\' >> {}".format(outputpath)
		os.system(cm)
			
	
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-q", "--query", help="path to output file", dest="query",
						type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
