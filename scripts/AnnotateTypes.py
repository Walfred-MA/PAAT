#!/usr/bin/env python3

import os 
import argparse
import sys
import numpy as np
import math
import collections as cl
import pandas as pd
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
	
	
	typefile = args.type
	annofile = args.anno
	
	
	refs = dict()
	annotations = dict()
	with open(annofile, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip())==0:
				continue
			line = line.strip().split('\t')
			annotations[line[0]] = [line[1], line[6], line[4]]
			
			
			if line[6] == "Ref":
				refs[line[0]] = line[2]
				annotations[line[0]] = [line[0], line[6], line[4]]
				
	genecode = GenodeDB(args.gene)
	
	alltypes = []
	with open(typefile, mode = 'r') as f:
		
		for iline,line in enumerate(f):
			
			line = line.strip()	
			if len(line) == 0:
				continue
			
			
			names = [x for x in line.split(",") if len(x)]
			annos = [annotations[name] for name in names]
			
			type_annos = cl.Counter([x[0] for x in annos]).most_common(1)[0][0]
			
			
			type_tags =  [x[1] for x in annos]
			if "Ref" in type_tags:
				type_tags = "Ref" 
			else:
				type_tags = cl.Counter(type_tags).most_common(1)[0][0]
				
			type_exon = [ "Intron" , "Exon" ] [ int("1" in [x[2] for x in annos]) ]
			
			type_refs = [name for name in names if annotations[name][1] == "Ref"]
			
			if len(type_refs) == 0 and type_annos in refs:
				type_refs  = [type_annos] 


			if len(type_refs) > 0:


				locations = [refs[name][:-1].split(":") for name in type_refs]

				locations_dic = cl.defaultdict(list)

				for i,loc in enumerate(locations):
					locations_dic[loc[0]].append(  list(map(int,loc[1].split("-"))) + [str(i)] )

				if len(locations):
					type_locs = list(set(sum ( [loc for name, loc in genecode.Overlap(locations_dic).items()],[])))

					if len(args.pref):
						type_locs = [x for x in list(type_locs) if x.startswith(args.pref)]

				
				type_refs = ";".join([name+":"+refs[name] for name in type_refs])

				shortnames = set(type_locs)
				newtype_locs = []
				for gene in type_locs:
					if '-' in gene and gene.split('-')[0] in shortnames:
						continue
					newtype_locs.append(gene)
					
				type_locs=newtype_locs
				type_locs = ";".join(type_locs)

				if len(type_locs) == 0:
					type_locs = "NA"


			else:
				type_locs = "NA"
				type_refs = "NA"	
							
			alltypes.append([str(iline),type_tags, type_locs, type_refs, type_exon,line])
		
			
		
		
	with open(args.output, mode = 'w') as f:
		
		for thetype in alltypes:
			print(thetype)	
			f.write("\t".join(thetype).replace('\n','')+'\n')
			
			
			
	return 0

def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine genes")
	parser.add_argument("-t", "--type", help="path to input data file", dest="type", type=str, required=True)
	parser.add_argument("-a", "--anno", help="path to input data file", dest="anno", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
	parser.add_argument("-g", "--genecode", help="path to output file", dest="gene", type=str, required=True)
	parser.add_argument("-p", "--pref", help="path to output file", dest="pref", type=str, default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()	
