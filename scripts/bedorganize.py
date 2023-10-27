#!/usr/bin/env python3

import os 
import argparse
import sys
import collections as cl

def mergeregions(regions):

	allcoordi = regions
	allcoordi_index = sorted(range(len(allcoordi)), key = lambda x: allcoordi[x])

	allregions = []
	start = 0
	coverage = 0
	for index in allcoordi_index:
		coordi = allcoordi[index]
		if index % 2 ==0:
			coverage += 1
			if coverage == 1:
				start = coordi
		else:
			coverage -= 1
			if coverage == 0:
				allregions.append([start, coordi])
	
	return allregions	


def main(args):
	
	inputpath = args.input
	outputpath = args.output
	
	gene_regions = cl.defaultdict(list)
	with open(inputpath, mode = 'r') as f:
		
		for line in f:
			
			chr, start, end, strand, name, gene = line.split()
		
			gene_regions[(name,chr,strand,gene)].append([ int(start), int(end)])
	
	names = cl.defaultdict(int)
	with open(outputpath, mode = 'w') as f:
		
		for (name, chr, strand, gene),regions in gene_regions.items():
			allcoordi = [x for y in regions for x in y]
		
			if name in names:
				names[name] += 1
				name += "00"+str(names[name])
			
			names[name] = 1
	
			f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,min(allcoordi),max(allcoordi),strand ,name,gene))
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str,required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
