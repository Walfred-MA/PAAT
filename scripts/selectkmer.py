#!/usr/bin/env python3

import os
import argparse
import collections as cl
import multiprocessing as mul


def kmerencode(kmer):
	
	kmerint = 0
	for base in kmer:
		
		kmerint *= 4
		if base=='a' or base =='A':
			
			kmerint += 0
			
		elif base=='t' or base =='T':
			
			kmerint += 3
			
		elif base=='c' or base =='C':
			
			kmerint += 1
		elif base=='g' or base =='G':
			
			kmerint += 2
			
	return kmerint

def kmerdecode(kmerint, size = 31):
	
	if type(kmerint) == type(""):
		return kmerint
	
	index = 0
	text = ['A' for a in range(size)]
	
	while kmerint:
		
		text[index] = "ACGT"[kmerint % 4]
		kmerint //= 4
		index += 1
		
	return "".join(text[::-1])

def select_kmers(kmerlist1, kmerlist2):
	include = {}
	exclude = {}
	
	allkmers = set(kmerlist1.keys()).union(kmerlist2.keys())
	
	for x in allkmers:
		
		count1 = kmerlist1.get(x, 0)
		count2 = kmerlist2.get(x, 0)
		if count1 == count2 and count1 < 255:
			
			include[x] = count1
		else:
			exclude[x] = abs(count2 - count1)
			
	return include, exclude

def compare_kmers(outputpath, pair):
	
	if len(pair) != 2:
		return 
	
	savememo = 1
	if max(os.path.getsize(pair[0]), os.path.getsize(pair[1])) < 100000000:
		savememo = 0
	
	with open(pair[0], mode = 'r') as f:
		
		if savememo:
			genekmers = {kmerencode(x.split()[0]):int(x.split()[1]) for x in f.readlines() if len(x) and x[0] != "@"}
		else:
			genekmers = {x.split()[0]:int(x.split()[1]) for x in f.readlines() if len(x) and x[0] != "@"}
		
	with open(pair[1],mode = 'r') as f:
		
		if savememo:
			assemkmers = {kmerencode(x.split()[0]):int(x.split()[1]) for x in f.readlines() if len(x) and x[0] != "@"}
		else:
			assemkmers = {x.split()[0]:int(x.split()[1]) for x in f.readlines() if len(x) and x[0] != "@"}
		
		
	includelist, excludelist = select_kmers(genekmers, assemkmers)
	
	with open(outputpath+"_kmer.txt", mode = 'w') as f:
		f.write("\n".join(["{}\t{}".format(kmerdecode(key), count) for key, count in includelist.items()]))
		
	with open(outputpath+"_exclude.txt", mode = 'w') as f:
		f.write("\n".join(["{}\t{}".format(kmerdecode(key), count) for key, count in excludelist.items()]))
		
	return  


def main(args):
	
	
	if os.path.isdir(args.gene):
		genekmersfiles = [args.gene+"/"+thefile for thefile in os.listdir(args.gene) if ".txt" in thefile]
	else:
		genekmersfiles = [args.gene]
		
	if os.path.isdir(args.assem):
		assemkmersfiles = [args.assem+"/"+thefile for thefile in os.listdir(args.assem) if ".txt" in thefile]
	else:
		assemkmersfiles = [args.assem]
		
		
	haplo_pair = cl.defaultdict(list)
	for thefile in genekmersfiles:
		
		haplo = "_".join(thefile.split("/")[-1].split(".")[0].split("_")[:2])
		haplo_pair[haplo].append(thefile)
		
	for thefile in assemkmersfiles:
		
		haplo = "_".join(thefile.split("/")[-1].split(".")[0].split("_")[:2])
		haplo_pair[haplo].append(thefile)
		
		
	p=mul.Pool(processes=args.threads)
	
	for haplo, pair in haplo_pair.items():
		
		#compare_kmers(args.output+"/"+haplo, pair)
		p.apply_async(compare_kmers,(args.output+"/"+haplo, pair))
		
	p.close()
	p.join()
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-g", "--gene", help="path to input data file", dest="gene", type=str, required=True)
	parser.add_argument("-a", "--assem", help="path to input data file", dest="assem", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	parser.add_argument("-n", "--nthreads", help="number of threads", dest="threads", type=int,default = 1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
