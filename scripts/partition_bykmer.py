#!/usr/bin/python

import collections as cl
import math
import numpy as np
import scipy 
import pandas as pd
import argparse
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
import kmer_cluster32 as kmer_cluster



def group_matrix(matrix, size, similar, common):

	groups = [[] for x in range(size)]
	group_edges = [x for x in range(size)]
	for i in range(size):
	
		indexstart = i*size - (i*(i+1)//2)
	
		#diag_cutoff = matrix[indexstart + i ]* similar
		
		for j in range(i+1, size):
			
			edge = matrix[indexstart + j ]
			if edge > common :
				
				if group_edges[i] < group_edges[j]:
					group_edges[j] = group_edges[i]
				else:
					group_edges[i] = group_edges[j]
				
	for i, group_index in enumerate(group_edges):
		
		last_group_index = i
		
		while group_index < last_group_index:
			
			last_group_index = group_index
			group_index = group_edges[group_index]
			
		group_edges[i] = group_index
		groups[group_index].append(i)
		
	return groups, group_edges

def readfasta(reffile):
	
	with open(reffile, mode = 'r') as f:
		
		readtext = f.read()
		
		readtext = readtext.split(">")[1:]
		
		#readtext = [read for read in readtext if samplename not in read.splitlines()[0]]
		
		headers = [read.splitlines()[0] for i,read in enumerate(readtext)]
		text = ["".join(read.splitlines()[1:]) for i,read in enumerate(readtext)]
		names = [header.split()[1].split("#")[0] for header in headers]
		
	return text, headers, names 


def makereverse(seq):
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)



def similar(table1, table2):
	
	
	if len(table1)>len(table2):
		larger, smaller = table1, table2
	else:
		larger, smaller = table2, table1
		
	total_counts = sum(list(smaller.values()))
	
	
	concordance = 0
	for kmer, count in smaller.items():
		
		larger_count = larger[kmer]
		
		concordance += min(count, larger_count)
		
	return concordance/total_counts


class KmerData:
	
	def __init__(self,samplelist, selectlist, matrix):
		
		self.selectkmers = selectlist
		self.samplelist = samplelist
		
		self.kmerlist = dict()
		self.tran=str.maketrans('ATCGatcg', 'TAGCtagc') 
		self.matrix = matrix
		self.matchscores = 0.0  
		self.weightscores = 0.0
		
	def getkmerindex(self, kmer):
		
		kmer = max(kmer, kmer[::-1].translate(self.tran))
		
		return self.kmerlist.get(kmer, -1)       
	
	def getsampleindex(self,kmer):
		
		kmerindex = self.getkmerindex(kmer)
		
		if kmerindex != -1:
		
			samples = set(self.matrix.getkmersamples(kmerindex).keys())
		
		else:
			samples = set()
			
		return samples 
	
	def countkmer(self, seq, index ,kmertype = 31):
		
		kmer_counter = 0
		for posi in range(len(seq)-kmertype+1):
			
			kmer = seq[posi:posi+kmertype]  
			
			kmer = max(kmer.upper(), kmer[::-1].translate(self.tran).upper())
			
			if len(self.selectkmers) and kmer not in self.selectkmers:
				continue
			
			kmerindex = self.kmerlist.get(kmer, -1) 
			
			if kmerindex < 0:
				
				kmerindex = len(self.kmerlist)
				self.kmerlist[kmer] = kmerindex
				self.matrix.addkmer()
				
			self.matrix.add(kmerindex, index)
			kmer_counter += 1
			
		return kmer_counter
	
	
def main(args):
	
	
	infile = args.input
	
	outfile = args.output
	
	kmerfile = args.kmer
	
	selectlist = set()
	
	if len(kmerfile):
		with open(kmerfile, mode = 'r') as f:
			
			for line in f:
				
				if len(line) and line[0] != ">":
					selectlist.add(line.strip().split()[0])
	
	samplename = infile.split("/")[-1].split(".")[0]
	
	seqs, headers, names = readfasta(infile)
	numseq = len(seqs)
	
	KmerMatrix = kmer_cluster.SparseKmerMartrix(numseq, 0)
	
	kmer_data = KmerData(names, selectlist, KmerMatrix)
	
	sizes = []
	
	for index,seq in enumerate(seqs):
		
		kmer_counter = kmer_data.countkmer(seq, index)
		#print(kmer_counter, name)
		
		sizes.append(kmer_counter)
		
	#matrix = np.reshape(KmerMatrix.MatchMatrix(), (numseq, numseq))

	#matrix = np.where(matrix<=0.1, 0, matrix)	
	#matrix_frac = matrix/(np.diag(matrix)+1)

	#matrix_samegroup = np.where((matrix<500) & (matrix_frac<0.1), 0, 1)

	groups,sample_groupindex = group_matrix(KmerMatrix.MatchMatrix_LM(), numseq, args.similar, args.common)
	
	group_kmers = [set() for group in groups]
	if len(selectlist):
		
		for kmer in selectlist:
			
			samplesfound = kmer_data.getsampleindex(kmer)
		
			groupsfound = set([sample_groupindex[sampleindex] for sampleindex in samplesfound ])
			
			if len(groupsfound) == 1:
				
				group_kmers[list(groupsfound)[0]].add(kmer)
	
	
	for group_index, group in enumerate(groups):
		
		if len(group) == 0:
			continue
		
		group_kmer = group_kmers[group_index]
		
		with open(outfile+"{}.fa".format(group_index+1), mode = "w") as f:
			
			for index in group:
				
				f.write(">{}\n{}\n".format(headers[index],seqs[index]))
				
		with open(outfile+"{}.fa_kmer.list".format(group_index+1), mode = "w") as f:
			
			for kmer in group_kmer:
				
				f.write(">\n{}\n".format(kmer))
				
				
				
				
				
				
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str, required=True)
	parser.add_argument("-k", "--kmer", help="path to kmer file", dest="kmer",
						type=str, default="")
	parser.add_argument("-s", "--similar", help="path to kmer file", dest="similar",
                        type=float, default=0.1)
	parser.add_argument("-c", "--common", help="path to kmer file", dest="common",
                        type=int, default=500)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
