#!/usr/bin/env python3

import pandas as pd
import collections as cl
import argparse
import pandas as pd


def segment_merge(segments):

	coordis = [x for y in segments for x in y[:2]]
	coordis_sortindex = sorted(range(len(coordis)), key = lambda x: coordis[x])

	merges = []
	start = 0
	depth = 0
	for x in coordis_sortindex:
		coordi = coordis[x]
		if x % 2 :
			depth -= 1
			if depth == 0:
				merges.append([start, coordi])
		else:
			depth += 1
			if depth == 1:
				start = coordi

	return merges
		


def main(args):
	
	inputfile = args.input
	region = args.region.replace('/rc','')
		
	outpath = args.output
	
	contig = region.split(':')[0]
	contigstart = int(region.split(':')[1].split('-')[0])
	contigend = int(region.split(':')[1].split('-')[1])
	contigsize = abs(contigend - contigstart)
			
	table = pd.read_csv(inputfile, header = None, sep = '\t')
	table = table[(table[0].str.contains(contig)) & (table[5].str.contains('ENSE') )]
	table = table.sort_values(by=[9], ascending = 0)
	
	aligns = table[[0, 5, 6,2,3,4, 7,8,9,12]].values.tolist()

	aligned_onquery = []
	for row in aligns:
	
	
		querystart = min(int(row[0].split('_')[-2]), int(row[0].split('_')[-1]))
		
		contigstart_onquery = contigstart - querystart
		
		contigend_onquery = contigend - querystart
		
		if row[4] - row[3]  + contigsize - max(row[4],contigend_onquery ) + min(row[3], contigstart_onquery) <= 0:
			continue
		if row[-2] < 0.9 * row[2]:
			continue
	
	
		if '.' in row [1]:
			row[1] = ".".join(row[1].split(".")[:-1])
			
		aligned_onquery.append(  [ max(0,row[3]-contigstart_onquery),  min(row[4]-contigstart_onquery, contigsize) ])
	
	merged = segment_merge(aligned_onquery)

	outtable = ";".join(["{}_{}".format(x[0],x[1]) for x in merged])

	with open(outpath, mode = 'w') as f:
		f.write(outtable+"\n")
	
		
	return 

def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-r", "--region", help="path to output file", dest="region",type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
