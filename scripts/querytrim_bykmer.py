#!/usr/bin/env python3

import os
import argparse
import collections as cl
import pandas as pd
import numpy as np


def trim(header, seq, start, end):

	
	headers = header.split()
	
	header_loc = headers[1]
	header_strd = header_loc[-1]
	header_contig,header_range = header_loc.split(":")
	header_range = header_range[:-1]
	header_coordi = header_range.split('-')
	header_start = max(0, int('-'.join(header_coordi[:-1])))
	header_end = int(header_coordi[-1])
	
	if header_strd == '+':
		newstart = header_start + start
		newend = header_start + len(seq) - end -1
	else:
		newstart = header_start + end
		newend = header_start + len(seq) - start - 1
		
	seq = seq[start:len(seq)-end]
	headers[1] = "{}:{}-{}{}".format(header_contig,newstart,newend, header_strd)
	
	header = "\t".join(headers)
	
	return header, seq

class kmer_database:
	
	def __init__(self):
		
		self.tran = str.maketrans('ATCGatcg', 'TAGCtagc')
		self.kmer_db = set()
		self.exclude = set()
		
	def load(self, kmerfile):
		
		with open(kmerfile, mode = 'r') as f:
			
			line = f.readline()
			if len(line) and line[0] != ">" and line[0]!= "@":
				self.kmer_db.add(line.strip().split()[0])
				
			for line in f:
				if len(line) and line[0] != ">":
					self.kmer_db.add(line.strip().split()[0])
		
	
	def countkmer(self, name, seq , kmertype = 31):
		
		self.foundkmers = [""] * (len(seq)-kmertype+1)

		if len(self.kmer_db) == 0:
			return 0
			
		for posi in range(len(seq)-kmertype+1):
			
			kmer = seq[posi:posi+kmertype].upper()
			kmer = max(kmer, kmer[::-1].translate(self.tran))
			
			if kmer in self.kmer_db:
								
				self.foundkmers[posi] = kmer
		
		count = len(self.foundkmers) - self.foundkmers.count("")
		
		return count
	
	def trimkmer(self, seq, cutoff, rcutoff):

		if len(self.kmer_db) == 0:
			return -len(seq),len(seq),len(seq)

	

		nonkmerscore = -1 * rcutoff/(1 - rcutoff)
		kmerscore = 1
		
		start = 0
		score = 0
		for pos, ifkmer in enumerate(self.foundkmers):
			
			if ifkmer:
				score += kmerscore
			else:
				score += nonkmerscore
			
			if score < 0:
				start = pos
				score = 0
			elif score >= cutoff:
				break
		
		end = 0
		score = 0
		for pos, ifkmer in enumerate(self.foundkmers[::-1]):
			
			if ifkmer:
				score += kmerscore
			else:
				score += nonkmerscore
				
			if score < 0:
				end = pos
				score = 0
			elif score >= cutoff:
				break
	
		kmerleft = len(self.foundkmers) - self.foundkmers.count("") - start - end  + self.foundkmers[:start].count("") + self.foundkmers[(len(self.foundkmers)-end):].count("")


		return kmerleft, start, end
			
		
	
def main(args):
	
	kmertype = args.size
	
	kmerfile = args.kmer
	
	infile = args.input
	
	outfile = args.output
	
	anchorsize = args.anchor
		
	kmers = kmer_database()
	kmers.load(kmerfile)
	
	outputfile = open(outfile, mode ='w' ) 
	excludefile = open(outfile+"_exclude", mode ='w' )

	writefile = outputfile
	header = ""
	name = ""
	seq = ""
	with open(infile, mode = 'r') as f:
		
		for line in f:
			
			if len(line) ==0:
				continue
			if line[0] == ">":
				
				if len(header) and len(seq):
				
					didtrim = "pass"	
					count = kmers.countkmer(name, seq)
					
					if count>= args.cutoff and count >= ( len(seq) - 2 * args.anchor) * args.rcutoff :
						writefile = outputfile
						
					elif count < args.cutoff:
						writefile = excludefile
						kmers.exclude.update(kmers.foundkmers)
						didtrim = "filter"
					else:
						count, start, end  = kmers.trimkmer(seq, args.cutoff, args.rcutoff)
						if count>= args.cutoff :
							writefile = outputfile
							kmers.exclude.update(kmers.foundkmers[:start]+kmers.foundkmers[(len(kmers.foundkmers)-end):])
							header, seq = trim(header, seq, start, end)
							didtrim = "trim"
						else:
							writefile = excludefile
							kmers.exclude.update(kmers.foundkmers)
							didtrim = "filter"
					
					writefile.write(header+"\n"+seq+"\n")
					
					print((name+"\t" + str(len(seq))+ "\t" +str(count) + "\t"+ didtrim))
				
				header = line.strip()
				name = header.split()[0][1:]
				seq = ""
			
			else:
				seq += line.strip()

		if len(header) and len(seq):
			
			didtrim = "pass"		
			count = kmers.countkmer(name, seq)
			
			if count>= args.cutoff and count >= ( len(seq) - 2 * args.anchor) * args.rcutoff :
				writefile = outputfile
						
			elif count < args.cutoff:
				writefile = excludefile
				kmers.exclude.update(kmers.foundkmers)
				didtrim = "filter"		
			else:
				count, start, end  = kmers.trimkmer(seq, args.cutoff, args.rcutoff)
				if count>= args.cutoff :
					writefile = outputfile	
					kmers.exclude.update(kmers.foundkmers[:start]+kmers.foundkmers[(len(kmers.foundkmers)-end):])	
					header, seq = trim(header, seq, start, end)	
					didtrim = "trim"			
				else:
					writefile = excludefile
					kmers.exclude.update(kmers.foundkmers)
					didtrim = "filter"		
						
			writefile.write(header+"\n"+seq+"\n")
					
			print((name+"\t" + str(len(seq))+ "\t" +str(count) + "\t"+ didtrim))
			
		
	outputfile.close()
	excludefile.close()
	
	
	kmeroutfile = args.kmerlist
	if len(kmeroutfile) == 0:
		kmeroutfile = outfile+"_kmer.list"
	with open(kmeroutfile, mode ='w' ) as f:
		
		f.write(">\n"+"\n>\n".join(kmers.kmer_db - kmers.exclude))
		
		
		
		
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-k", "--kmer", help="path to output file", dest="kmer",type=str, required=True)
	parser.add_argument("-s", "--size", help="kmer size", dest="size",type=int, default = 31)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-a", "--anchor", help="path to output file", dest="anchor",type=int, default=5000)
	parser.add_argument("-l", "--list", help="path to kmer output file", dest="kmerlist",type=str, default="")
	parser.add_argument("-c", "--cutoff", help="path to kmer output file", dest="cutoff",type=int, default=500)
	parser.add_argument("-r", "--rcutoff", help="path to kmer output file", dest="rcutoff",type=float, default=0.1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()


