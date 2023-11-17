#!/usr/bin/env python3

import os
import sys
import argparse
import collections as cl

class error_region:
	def __init__(self, loadfile):
		
		self.error_regions = cl.defaultdict(list)

		if len(loadfile) == 0:
			return

		with open(loadfile,mode = 'r') as f:
			
			for line in f:
				
				if len(line.strip())==0 or line[0] in ["@","#","=","\t"," "]:
					continue
				contig, start, end = line.split()[:3]
				self.error_regions[contig].append(sorted([int(start),int(end)]))
				
	def findoverlap(self, contig, region):
		
		error_oncontig = self.error_regions[contig]
		
		region_size = region[1]-region[0]
		for error_region in error_oncontig:
			
			allcordi = error_region + region
			
			overlap =   max(allcordi) - min(allcordi)  -   (error_region[1]-error_region[0]) - region_size
			if overlap <= 0:
				
				return True
			
		return False

def runfilter(input, output, scafsizes, min_edge,error_regions, cutoff):
	
	passnum = 0
	totalnum = 0
	iffilter = 0		
	w = open(output,mode ='w') 
	
	with open(input, mode = 'r') as f:
		
		for line in f:
			if len(line) ==0:
				continue
			if line[0] == ">":
				
				iffilter = 1
				header = line[1:]
				name = header.split()[0]
				region = header.split()[1][:-1]
				contig = region.split(":")[0]
				scafsize = scafsizes.get(contig, 1000000000)
				
				coordi = region.split(":")[1].split("-")
				coordi = ["-".join(coordi[:-1]),coordi[-1]]
				coordi = sorted([int(x) for x in coordi])
				
				edge = min(coordi[0], scafsize - coordi[1]) 
				size = coordi[1]-coordi[0]
				
				scafcutoff = cutoff
				if "chr" in name or "NC_0609" in name or "CN000" in name or "HQ100" in name:
					scafcutoff = 0
				elif "cluster" in contig:
					scafcutoff = cutoff
				elif "#1" in name or "#2" in name:
					scafcutoff = cutoff
					
				scafcutoff = min(scafcutoff, min_edge/2)				
				
				totalnum += 1
				if error_regions.findoverlap(contig,coordi) == False and edge >= scafcutoff:
					iffilter = 0
					passnum += 1
					
			if not iffilter:
				w.write(line)
				
	w.close()

	return passnum , totalnum 
	
	
def main(args):
	
	#args.error = "/project/mchaisso_100/cmb-16/walfred/projects/distract_genes/hprc_error.bed"
	
	scafsizes = {}
	if len(args.ref):
		refs = args.ref.split(",")
		for ref in refs:
			with open(ref+".fai", mode = 'r') as f:
				for info in f:
					if len(info):
						info = info.strip().split()
						scafsizes[info[0]] = int(info[1])
						
	if len(args.query):
		with open(args.query, mode = 'r') as f:
			for line in f:
				line = line.split()
				
				with open(line[1]+".fai", mode = 'r') as f:
					for info in f:
						if len(info):
							info = info.strip().split()
							scafsizes[info[0]] = int(info[1])
							
							
	error_regions = error_region(args.error)
	
	min_edge = 1000000000
	with open(args.input, mode = 'r') as f:
		for line in f:
			if len(line) ==0 or line[0] != ">" or ("chr" not in line and "NC_0609" not in line):
				continue
			
			header = line[1:]
			region = header.split()[1][:-1]
			contig = region.split(":")[0]
			
			if "NC_0609" in contig or ( "chr" in contig and "_" not in contig and "v"  not in contig ) or "CN000" in contig :			
				
				coordi = region.split(":")[1].split("-")
				coordi = ["-".join(coordi[:-1]),coordi[-1]]
				coordi = sorted([int(x) for x in coordi])
				
				scafsize = scafsizes.get(contig, 1000000000)
				
				edge = min(coordi[0], scafsize - coordi[1])
				
				min_edge = min(min_edge, edge)
				
	passnum , totalnum = runfilter(args.input, args.output, scafsizes, min_edge, error_regions, 20000)
	
	if passnum < 0.5 * totalnum:
		pass
		#runfilter(args.input, args.output, scafsizes, min_edge, error_regions, 100)
	
	
	
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to input data file", dest="output", type=str, required=True)
	parser.add_argument("-e", "--error", help="path to input data file",dest="error", type=str,default = "")
	parser.add_argument("-q", "--querypath", help="path to input data file",dest="query", type=str, default = "")
	parser.add_argument("-s", "--scaffolds", help="path to scaffold file",dest="scaf", type=str, default = "")
	parser.add_argument("-r", "--ref", help="path to scaffold file",dest="ref", type=str, default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
