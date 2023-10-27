#!/usr/bin/env python3

import collections as cl
import math
import numpy as np
import scipy 
import pandas as pd
import argparse

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
	
	index = 0
	text = ['A' for a in range(size)]
	
	while kmerint:
		
		text[index] = "ACGT"[kmerint % 4]
		kmerint //= 4
		index += 1
		
	return "".join(text[::-1])

class KmerData:
	
	def __init__(self, ref,  kmersize = 31):
		
		self.kmersize = kmersize
		self.ref= ref
		
		self.current_k = 0
		self.reverse_k = 0
		self.kmerslist = dict()
		self.intoperator = 4**(kmersize-1)
		
		self.translator = str.maketrans('ATCGatcg', 'TAGCtagc')
		
		samples = []
		with open(self.ref,mode = 'r') as r:
			for line in r:
				
				line = line.strip()
				
				if line[0]==">":
					samples.append(line.split()[0][1:])
					
					
		self.sampleslist = list(samples)
		
		self.matrix = kmer_cluster.SparseKmerMartrix(len(self.sampleslist), 0)
		
		
	def addkmer(self,kmer, index):
		
		kmerindex = self.kmerslist.get(kmer, -1) 
		
		if kmerindex < 0:
			
			kmerindex = len(self.kmerslist)
			self.kmerslist[kmer] = kmerindex
			self.matrix.addkmer()
			
		self.matrix.add(kmerindex, index)
		
	def getKmers(self):
		
		currentline = ""
		
		with open(self.ref,mode = 'r') as r:
			
			index = -1
			for line in r:
				
				line = line.strip()
				
				if line[0]==">":
					index += 1
					
				if len(line) ==0 or line[0]==">":
					
					self.current_k=0
					self.reverse_k=0
					self.current_size = 0
					
					continue
				for char in line:
					
					if char =="\n":
						continue
					
					if char not in ['A','T','C','G']:
						
						self.current_k = 0
						self.reverse_k = 0
						self.current_size = 0
						continue
					
					if self.current_size >= self.kmersize:
						
						self.current_k %= self.intoperator
						self.current_k <<= 2
						self.current_k += ['A','C','G','T'].index(char)
						
						self.reverse_k >>= 2
						self.reverse_k += (3-['A','C','G','T'].index(char)) * self.intoperator
						
					else:
						
						
						self.current_k <<= 2
						self.current_k += ['A','C','G','T'].index(char)
						self.reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*self.current_size))
						
						self.current_size += 1
						
					if self.current_size >= self.kmersize:  
						
						kmer = max(self.current_k, self.reverse_k)
						
						self.addkmer(kmer,index)


class Kmer_annotation:
	
	def __init__(self, kmerlist, kmersize =31):
		
		self.tran=str.maketrans('ATCGatcg', 'TAGCtagc') 
		self.kmerlist = kmerlist
		self.kmersize = kmersize
		self.intoperator = 4**(kmersize-1)
		
	def annotate(self, seq, kmersize = 31, hardmask = 0):
		
		outseq = seq
		
		foundposi = []
		
		self.current_k = 0
		self.reverse_k = 0
		self.current_size = 0
		
		lastposi = -10000000
	
		
		for posi,char in enumerate(seq.upper()):
			
			if char not in ['A','T','C','G']:
				
				self.current_k = 0
				self.reverse_k = 0
				self.current_size = 0
				continue
			
			if self.current_size >= kmersize:
				
				self.current_k %= self.intoperator
				self.current_k <<= 2
				self.current_k += ['A','C','G','T'].index(char)
				
				self.reverse_k >>= 2
				self.reverse_k += (3-['A','C','G','T'].index(char)) * self.intoperator
				
			else:
				
				
				self.current_k <<= 2
				self.current_k += ['A','C','G','T'].index(char)
				self.reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*self.current_size))
				
				self.current_size += 1
				
			if self.current_size >= kmersize:  
				
				kmer = max(self.current_k, self.reverse_k)
				
				if kmer in self.kmerlist:
					
					if posi - lastposi <= kmersize:
						
						foundposi[-1] = (foundposi[-1][0], posi)

						
					else:
						
						foundposi.append((posi, posi))
						
					lastposi = posi
	
		length = len(seq)	

		if hardmask == 0:

			seq = bytearray(seq.encode())

			for region in foundposi:
	
				for posi in range(region[0],region[1]):
				#posi = region[1]	
					if posi < length and  seq[posi] >= ord('a'):
					
						seq[posi] -= 32
	
			seq = seq.decode()

		else:

			seq_o = seq
			seq_o = bytearray(seq_o.upper().encode())

			seq = bytearray(seq.encode())

			for region in foundposi:
				for posi in range(region[0],region[1]):
					if posi < length:
						seq_o[posi] = seq[posi]
			
			seq = seq_o.decode()
				

		return seq
	
	
	
def main(args):
	
	
	infile = args.input
	
	outfile = args.output
	
	
	
	kmerfile = args.kmer
	
	kmerlist = set([])
	with open(kmerfile, mode = 'r') as f:
		
		for line in f:
			
			if len(line) and line[0] != ">":
				
				kmerlist.add(kmerencode(line.strip().split()[0]))
				
				
	annotator = Kmer_annotation(kmerlist)

	if args.hardmask:	
		tran=str.maketrans('atcg', 'nnnn')
	else:
		tran=str.maketrans('', '') 
	header = ""
	seq = ""
	with open(infile, mode = 'r') as read: 
		
		with open(outfile, mode ='w') as write: 
			
			for line in read:
			
				if len(line) and line[0] not in [">", "@"]:
					
					if args.unmask and args.hardmask == 0: 
						seq += line.strip().lower()
					else:
						seq += line.strip()
				else:
				
					seq = annotator.annotate(seq,hardmask=args.hardmask).translate(tran)
					
					write.write(header + "\n"+ seq + "\n")
					
					header = line.strip()
					seq = ""
			
			seq = annotator.annotate(seq,hardmask=args.hardmask).translate(tran)


			write.write(header + "\n"+ seq)
			
			
			
			
			
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
						type=str, required=True)
	parser.add_argument("-u", "--unmask", help="if unmask first", dest="unmask",
			type=int, default = 1)
	parser.add_argument("-H", "--hardmask", help="if unmask first", dest="hardmask",
			type=int, default = 0)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
