#!/usr/bin/env python3

import pandas as pd
import os 
import re
import argparse
import sys
import multiprocessing as mul
import collections as cl

def readfile(querypath):

	read = ""
	header = ""
	query = {}
	with open(querypath, mode = 'r') as f:
		for line in f:
			if len(line) == 0:
				continue
			if line[0]==">":
				query[header] = read
				read = ""
				header = line[1:].split()[0]
			else:
				read+=line.strip()

	query[header] = read
	read = ""
	

	return query


def Ngap_trim(start, end, seq):

	beginN = len(re.findall(r'^N+',seq))
	endN = len(re.findall(r'N+$',seq))

	seq = seq[beginN:(len(seq)-endN)]

	return beginN, endN, seq

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

def process(inputpath, outputpath, query):
	
	try:
		table = pd.read_csv(inputpath, header = None , sep = "\t")
	except:
		print("warnning, no kmer found in ", inputpath)
		return
	
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
	
	write = open(outputpath,mode = "w")
	for i,raw in enumerate(allqueryloci):
		
		contig, start, end = raw
		
		samplename = "_h".join(contig.split("#")[:2])   
	
		seq = query[contig.split()[0]][start:end]

	
		title = ">{}_{}_{}\t{}:{}-{}".format(contig, start , end , contig, start , end )
		
		write.write(title+"\n")
		write.write(seq+"\n")
		
	write.close()

def main(args):
	
	inputpath, querypath ,outputpath= args.input, args.query , args.output
	megainput, megaoutput, prefix = args.megainput, args.megaoutput,args.prefix
	substart,num = args.sub,args.num

	query = readfile(querypath)	
	
	with open(megainput,mode='r') as f:
		inputfolders = ["/".join(x.split("/")[:-1]) for x in f.read().splitlines()[substart:(substart+num)]]
		
	
	if len(megainput):
		#allinputs = [inputfolder+"/"+x for inputfolder in inputfolders for x in os.listdir(inputfolder) if prefix in x and x[-12:]== "_hotspot.txt"]
		allinputs = [inputfolder+"/"+prefix+"_hotspot.txt" for inputfolder in inputfolders]
		alloutputs = [x+".fa" for x in allinputs]

	else:
		allinputs = [inputpath]
		alloutputs = [outputpath]
		

	for inputpath, outputpath in zip(allinputs, alloutputs):

		process(inputpath, outputpath, query)

		
		
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-I", "--megainput", help="path to input data file",
						dest="megainput", type=str, default = '')
	parser.add_argument("-q", "--query", help="path to output file", dest="query",
						type=str, required=True)
	parser.add_argument("-p", "--prefix", help="prefix of query", dest="prefix",
						type=str, default = '')
	parser.add_argument("-O", "--megaoutput", help="path to output file", dest="megaoutput",
						type=str, default = '')
					
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, default = '')
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str, default = '')
	parser.add_argument("-s", "--subsetstart", help="path to output file", dest="sub",
                        type=int, default = 0)
	parser.add_argument("-n", "--num", help="path to output file", dest="num",
                        type=int, default = 100)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

	
