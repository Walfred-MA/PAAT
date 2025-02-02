#!/usr/bin/env python3

import os 
import argparse
import sys
import numpy as np
import math
import collections as cl
import pandas as pd
import re


def overlap_eachcontig(allaligns, allow_gap = 0):
	
	all_coordinates = [x for y in allaligns for x in y[1:3]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group_index = []
	current_group_coordi = []
	
	number_curr_seq = 0
	last_coordinate = -allow_gap - 10
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		if index %2 == 0 :
			number_curr_seq += 1
			
			if number_curr_seq == 1 and coordinate - last_coordinate>allow_gap :
				
				allgroups.append([current_group_index, current_group_coordi])
				
				current_group_coordi = [coordinate, coordinate]
				current_group_index = [index//2]
				
			else:
				current_group_coordi[-1] = last_coordinate
				current_group_index.append(index//2)
				
		else:
			number_curr_seq -= 1
			if number_curr_seq == 0:
				current_group_coordi[-1] = last_coordinate
				
		last_coordinate = coordinate
		
	current_group_coordi[-1] = last_coordinate
	allgroups.append([current_group_index, current_group_coordi])
	
	return allgroups[1:]

def overlap(allaligns, allow_gap = 0):
	
	contigs = cl.defaultdict(list)
	for region in allaligns:
		contigs[region[0]].append(region)
		
	names= []
	allregioncombines = []
	for contig, regions in contigs.items():
		
		regioncombines = overlap_eachcontig(regions)
		
		for regionindex, region in regioncombines:
			
			#name = ";".join(list(set(sum([regions[i][4].split(";") for i in regionindex],[]))))
			#names.append(name)
			
			strds = [regions[i][3].split(";") for i in regionindex]
			strd = '+' if strds.count('+') >= strds.count('-') else '-'
			
			allregioncombines.append("{}:{}-{}{}".format(contig, region[0], region[1],strd))
			
	return ",".join(allregioncombines)

def filtersegment_bysize(cigar):
	
	segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
	
	curr_path = ""
	for segment_cigar in segment_cigars:
		
		thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', segment_cigar)]+[0])
		if thesize > 300:
			curr_path += segment_cigar.split(":")[0]
			
	return curr_path


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
	
	
	annofile = args.input
	
	liftinfo = cl.defaultdict(lambda: "NA")
	reflocs = cl.defaultdict(lambda: "NA")
	annotations = cl.defaultdict(list)
	alltypes = cl.defaultdict(list)
	all_pathes = cl.defaultdict(list)
	
	chunkfix = cl.defaultdict(list)
	with open(annofile, mode = 'r') as f:
		
		for line in f:
			name, typename, gene_overlap, locusindex, gene_count, category, ifexon, qlocation, closeref, similar, liftlocus, liftref, kmernum,note ,sv,exons, path, pathcigar,liftcigar = line.split('\t')
			
			if ":" in liftlocus:
				contig = liftlocus.split(":")[-2]
				therange = liftlocus.split(":")[-1][:-1]
				start = max(0, int("-".join(therange.split("-")[:-1])))
				end = int(therange.split("-")[-1])
				liftinfo[name] = [contig, start,end , liftlocus[-1]]
				
				
			all_pathes[name] = filtersegment_bysize(pathcigar)
			
			gene_counts = [  x.split(":")[1] + (re.search(r'\(.*?\)', x).group() if re.search(r'\(.*?\)', x) else '')  if float(x.split(":")[-1]) > 90.0 else x.split(":")[1]+(re.search(r'\(.*?\)', x).group() if re.search(r'\(.*?\)', x) else '')+"(Partial)"  for x in gene_count.split(";")[:-1]] 
			
			gene_counts = dict(cl.Counter(gene_counts))
			
			genenames = sorted(list(gene_counts.keys()), key = lambda x: gene_counts[x], reverse = 1)
			
			counts = ";".join([genename+":"+str(gene_counts[genename]) for genename in genenames])
			
			annotations[name] = [closeref,counts,category,ifexon, gene_overlap,sv,kmernum,note]
			if closeref == name:
				
				reflocs[name] = qlocation
				annotations[name] = [name,counts,category,ifexon, gene_overlap,sv,kmernum,note]
				
			alltypes[typename.split(":")[0]].append(name)
			
			if ":" in typename:
				chunkfix[name] = [typename.split(":")[1].split("&")[0]]
				chunkfix[typename.split(":")[0]].append(typename.split(":")[1])
				
	if len(args.type):
		
		alltypes = cl.defaultdict(list)
		with open(args.type, mode = 'r') as f:
			for line in f:
				if len(line.strip()) == 0:
					continue
				
				line = line.strip().split()
				typename, members = line[0],line[-1].split(",")
				alltypes[typename] = [x for x in members if len(x)]
				
				
				
	alltype_pathes = {typename:list(set([all_pathes[name] for name in members]))   for typename, members in alltypes.items()}
	
	allmajortypes = dict()
	majorfix = dict()

	typenames_sorted = sorted(list(alltypes.keys()), key = lambda x: int(x.split("_")[-1]))
	for typename in typenames_sorted:
	
		names = alltypes[typename]	
		alltypes_path = alltype_pathes[typename]
		alltypes_path = cl.Counter(alltypes_path).most_common(1)[0][0]
		if alltypes_path not in allmajortypes:
			allmajortypes[ alltypes_path] = len(allmajortypes)
		majorfix[typename] = typename + "_"+str(allmajortypes[ alltypes_path])
		
	
	chunkfix = {name:[";".join([majorfix[y] for y in x.split(";")]) for x in values] for name,values in chunkfix.items()}	
	

	external_include = cl.defaultdict(list)
	for typename, names in alltypes.items():

		for name in names:
			if name in chunkfix:
				blocks = chunkfix[name][0]
				for block in blocks.split(";"):
					external_include[block].append(name+"@"+blocks)

	
		
	outputs= []
	for typename, names in alltypes.items():
		
		
		majortype = majorfix[typename]
		
		type_chunkfix = chunkfix.get(typename,[])		
		type_chunkfix = cl.Counter(type_chunkfix).most_common(1)[0][0] if 2*len(type_chunkfix) > len(names)  else ""
		
		annos = [annotations[name] for name in names]
		
		type_genecount = cl.Counter([x[1] for x in annos]).most_common(1)[0][0]
		
		if len(type_genecount) == 0:
			type_genecount = "NA"
			
			
		type_gene = [x[4].replace(";",',') for x in annos ]
		type_gene_counts = dict(cl.Counter(type_gene))
		type_gene = ";".join(sorted(list(type_gene_counts.keys()), key = lambda x: type_gene_counts[x], reverse = 1)[:1] )
		
		
		
		type_tags =  [x[2] for x in annos]
		counts = cl.Counter(type_tags)
		type_tags = sorted(list(set(type_tags)), key=lambda x: counts[x], reverse=True)
		type_tags =  ";".join(type_tags)
		
		type_exon = [ "Intron" , "Exon" ] [ int("1" in [x[3] for x in annos]) ]
		
		type_exon = [x[3] for x in annos]
		type_exon = "Exon" if "Exon" in type_exon else "Decoy" if "Decoy" in type_exon else "Intron"
		
		type_refs = [x[0] for x in annos]
		type_ref_counts = dict(cl.Counter(type_refs))
		type_ref = sorted(list(type_ref_counts.keys()), key = lambda x: type_ref_counts[x], reverse = 1)[0]
		if type_ref != "NA":
			type_refloc = type_ref+":"+reflocs[type_ref.split("(")[0]]
		else:
			type_refloc = "NA"
			
		type_sv = [x[5].split(";") for x in annos if x[0] == type_ref ]
		
		svs_dict = cl.defaultdict(list)
		for svs in type_sv:
			
			svs_sort = sorted([(sv.split('_')[0], int(sv.split('_')[-1])) if "_" in sv else ("NA",0) for sv in svs ])
			svs = tuple(sorted([sv[0] for sv in svs_sort ]))
			svs_dict[svs].append([sv[1] for sv in svs_sort])
			
		svs_mostcommon = sorted(list(svs_dict.keys()), key = lambda x: svs_dict[x], reverse = 1)[0]
		if list(svs_mostcommon) != ["NA"]:
			
			svs_modesizes = [max(set(x), key=x.count) for x in zip(*svs_dict[svs_mostcommon])]
			svs_mostcommon = ";".join(["{}({})".format(x, y) for x,y in zip(svs_mostcommon,svs_modesizes)])
		else:
			svs_mostcommon = "NA"
			
		type_note = cl.Counter([x[-1] for x in annos]).most_common(1)[0][0]
		
		
		usekmer_min = min([int(annotations[name][6]) for name in names])
		
		prefix = "_".join(line.split("_")[:2])
		outputs.append([majortype,type_tags, type_gene,type_genecount, type_refloc, type_exon,str(usekmer_min),svs_mostcommon,type_chunkfix+"|"+type_note,",".join([name +":" + chunkfix[name][0] if name in chunkfix else name for name in names] + external_include[majortype])])
		
		
		
		
	with open(args.output, mode = 'w') as f:
		
		for thetype in outputs:
			
			f.write("\t".join(thetype).replace('\n','')+'\n')
			
			
			
	return 0

def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine genes")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-t", "--type", help="path to output file", dest="type", type=str, default = "")
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
