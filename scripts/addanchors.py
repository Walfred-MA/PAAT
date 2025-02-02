#!/usr/bin/env python3

import argparse
import os

def main(args):
	
	coordi_file = args.input
	path_file = args.query
	output = args.output

	if len(output) == 0:
		output = args.input + "_largeanchor.fa"

	with open(path_file,mode='r') as f:
		
		pathes = {line.split()[0]: line.split()[1] for line in f.read().splitlines()}
		
	f.close()

	pathes["hg38"], pathes["chm13"] = args.refs.split(",")

	with open(coordi_file,mode='r') as f:
		
		coordinates = [line for line in f.read().splitlines() if len(line) and line[0]==">"]
		
	f.close()

	try:
		os.remove(output)
	except:
		pass

	stdflag = {"+":"",'-':""}
	
	for i,line in enumerate(coordinates):
		
		name, info = line[1:].split()[:2]

		contig,coordi = info.split(":")
			
		strand = coordi[-1]
		coordi = coordi[:-1]

		start = coordi[0]+coordi[1:].split("-")[0]
		end = coordi[1:].split("-")[1]
		
		
		start = int(start) 
		end =int(end) 
		haplo = "_h".join(contig.split("#")[:2])

		
		path = pathes.get(haplo,pathes["hg38"])
		if "NC_0609" in contig:
			path = pathes.get(haplo,pathes["chm13"])
		
		header = name + " {}:{}-{}{}".format(contig, int(start) - args.Anchor, int(end) + args.Anchor  , strand)
		
		cm = "echo \>{} >> {}".format(header,output)
		
		os.system(cm)
		
		cm = "samtools faidx {} {}:{}-{} {} | grep -v \">\" >> {}".format(path, contig, max(0,int(start)-args.Anchor), int(end)+args.Anchor, stdflag[strand] ,output)
		
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
						type=str, default="")
	parser.add_argument("-A", "--Anchor", help="new anchor", dest="Anchor",
						type=int, default=100000)
	parser.add_argument("-r", "--refs", help="reference", dest="refs",type=str)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
