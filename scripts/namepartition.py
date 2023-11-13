#!/usr/bin/env python3

import os
import collections as cl
import pandas as pd
import argparse
import re

class GenodeDB:
	
	def __init__(self, path):
		
		self.pattern = re.compile(r'(?:gene_name=)([^;]+)')
		self.path = path
		self.gene = cl.defaultdict(list)
		
	
		with open(self.path, mode = 'r') as r:

			for line in r:
				if len(line) == 0 or line[0]=="#":
					continue
			
				row = line.split("\t")

			
				if row[2] != "gene":
					continue
			
				gene_name = re.findall( self.pattern ,row[8])
		
				if len(gene_name) == 1:
			
					geneinfo = [int(row[3]), int(row[4])] + gene_name
			
					self.gene[row[0]].append(geneinfo)
	
	def OverlapChr(self, chr , queries):
		
		results = [[] for x in queries]
		
		refs = self.gene[chr]

		length = len(refs) + len(queries)
		lgenes = len(refs)
		
		coordinates = [x for y in refs + queries for x in y[:2]] 
		names = [x[2] for x in refs]
		
		coordinates_sortindex = sorted( range(length * 2), key = lambda x: coordinates[x])
		
		current_refs = []
		current_queries = []
		for index in coordinates_sortindex:

			index2 = index // 2
			
			if index % 2 == 0:
				
				if index2 < lgenes:
					
					current_refs.append(index2)
					
					for query in current_queries:
						
						results[query].append(index2)
			
				else:
					
					current_queries.append(index2 - lgenes)
					
					results[index2 - lgenes].extend(current_refs)
			
			else:
				
				if index2 < lgenes:
					
					current_refs.remove(index2)
				
				else:
					
					current_queries.remove(index2 - lgenes)
	
			
		outputs = cl.defaultdict(list)
		for i,result in enumerate(results):
			
			queryname = queries[i][2]
			
			outputs[queryname].extend([refs[refindex][2] for refindex in result])
		
		return outputs
		
	def Overlap(self, locations):
		
		results = cl.defaultdict(list)
		for chr, locations_bychr in locations.items():
			
			result_chr = self.OverlapChr(chr, locations_bychr)
			
			for queryname,genes in result_chr.items():
				
				results[queryname].extend(genes)
		
		return results

def HG38Locations(file, allLocations):
	
	with open(file, mode ='r') as f:
		for line in f:
			if len(line) == 0 or line[0] != ">":
				continue
			region = line.split()[1]
			chr,generange = region.split(":")[0], region.split(":")[1][:-1]
			start = max(0,int("-".join(generange.split("-")[:-1])))
			end = int(generange.split("-")[-1])
			allLocations[chr].append([start,end,file])
	
def main(args):
	
	with open( args.input , mode = 'r' ) as f:
		
		infiles = f.read().splitlines()
	
	genecode = GenodeDB(args.gene)
	
	allLocations = cl.defaultdict(list)
	for file in infiles:
		
		HG38Locations(file, allLocations)
	
	results = genecode.Overlap(allLocations)

	for file in infiles:

		result = list(set(results[file]))

		if len(result):

			prefix = file.split("/")[-1].split("_")[0]

			allgenes = ",".join(result)+","
			
			oldfile = file
			if file.split("/")[-1].count("_")>1:
				file = "/".join(file.split("/")[:-1])+'/' +  "_".join(file.split("/")[-1].split("_")[:2]) + ".fa"
				os.system("mv {} {}".format(oldfile,file))
				os.system("mv {}_kmer.list {}_kmer.list".format(oldfile,file))
			
			tag = "OOO".join((sorted([ x for x in  result if prefix in x])+ sorted([ x for x in  result if prefix not in x ]))[:3]).replace(".","v")+".fasta"

			newfilename =  ".".join(file.split(".")[:-1]) +"_"+tag
		
			print(file + "\t" + newfilename + "\t" +allgenes )


def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-g", "--gene", help="refgene file path", dest="gene",
                        type=str,required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
