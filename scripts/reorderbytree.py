#!/usr/bin/env python3

import os
import re
import collections as cl
import math
import numpy as np
import scipy 
import pandas as pd
import argparse
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
import kmer_cluster
import kmer_cluster32


def kmerencode(kmer):
	
	kmerint = 0
	for base in kmer:
		
		kmerint *= 4
		if base=='a' or base =='A':
			
			kmerint += 0
			
		elif base=='c' or base =='C':
			
			kmerint += 1
			
		elif base=='g' or base =='G':
			
			kmerint += 2
		elif base=='t' or base =='T':
			
			kmerint += 3
			
	return kmerint

def kmerdecode(kmerint, size = 31):
	
	index = 0
	text = ['A' for a in range(size)]
	
	while kmerint:
		
		text[index] = "ACGT"[kmerint % 4]
		kmerint //= 4
		index += 1
		
	return "".join(text[::-1])


def intencode(value):
	
	code = ""
	
	code = chr( ord('0') + (value % 64) ) + code
	value = value//64
	
	while value:
		
		code = chr( ord('0') + (value % 64) ) + code
		value = value//64
		
	return code

class UPGMANode:
	
	def __init__(self, left=None, right=None, up_dist=0.0, down_dist=0.0, index = None):
		
		self.index = index
		self.left = left
		self.right = right
		self.up_dist = up_dist
		self.down_dist = down_dist
		self.nleaves = 1
		
		if type(left) is type(UPGMANode):
			self.nleaves += left.nleaves + right.nleaves
		else:
			self.nleaves = 1
			
			
	def leaves(self) -> list:
		
		if self is None:
			return []
		elif self.right == None:
			return [self.left]
		
		return self.left.leaves() + self.right.leaves() 
	
	def to_newick(self) -> str:
		
		if self.right == None:
			return self.left + ":" + '{:.7f}'.format(self.up_dist)
		else:
			return ( "(" + ",".join([x.to_newick() for x in [self.left, self.right]]) + "):"+ '{:.6f}'.format(self.up_dist) )

	def allleaves(self):
		
		if type(self.left) is not type(UPGMANode) or  type(self.right) is not type(UPGMANode):
			return self
	
	
		return self.left.allleaves() + self.right.allleaves()

	def reorder(self, dist_matrix, sibling_indexes = None, sign = 0) -> list:
		
		if type(self.left) is not type(UPGMANode) or  type(self.right) is not type(UPGMANode):
			return
		
		leftleaves_indexes = [x.index for x in self.left.allleaves()]
		rightleaves_indexes = [x.index for x in self.right.allleaves()]
		
		self.left.reorder(dist_matrix, rightleaves_indexes, 1)
		self.right.reorder(dist_matrix, leftleaves_indexes, -1)

		if sibling_indexes == None or sign == 0:
			return 
		
		left_distances = [ dist_matrix[i,j] for j in sibling_indexes for i in leftleaves_indexes]
		left_distance_mean = sum(left_distances)/len(left_distances)
		
		right_distances = [ dist_matrix[i,j] for j in sibling_indexes for i in rightleaves_indexes]
		right_distance_mean = sum(right_distances)/len(right_distances)
	
		if sign * right_distance_mean > sign * left_distance_mean:
			
			temp = self.left
			self.left = self.right
			self.right = temp

	
		
class UPGMA:
	
	
	def __init__(self, dist_matrix: np.ndarray, header: list):
		
		self.distances = dist_matrix
		self.header = header
		self.exclude = [0 for x in range(len(header))]+[1 for x in range(len(header))]
		self.build_tree(self.distances, self.header)
		
	def getmindist(self)-> tuple:
		
		min_row, min_col = 0,0
		
		curr_min = np.inf
		for row,values in enumerate(self.work_matrix):
			
			if self.exclude[row]:
				continue
			
			for col, value in enumerate(values):
				
				if self.exclude[col]:
					continue
				
				if value == 0:
					
					return (row, col)
				
				elif value < curr_min:
					
					curr_min = value
					min_row = row
					min_col = col
					
					
		return (min_row, min_col)
	
	
	def build_tree(self, dist_matrix: np.ndarray, header: list) -> UPGMANode:
		
		nodes = [UPGMANode(taxon, index = i) for i, taxon in enumerate(header)]	
	
		headertoindex = {name:i for i,name in enumerate(header)}
		
		self.size = 2*len(header)
		
		self.work_matrix = np.array([row+[np.inf]*len(row) for row in dist_matrix.tolist()]+[[np.inf]*(2*len(row)) for row in dist_matrix.tolist()], dtype=float)
		np.fill_diagonal(self.work_matrix, np.inf)
		
		new_node = None 
		for turn in range(len(nodes)-1):
			
			least_id = self.getmindist()
			least_dist = self.work_matrix[least_id[0], least_id[1]]
			node1, node2 = nodes[least_id[0]], nodes[least_id[1]]
			self.exclude[least_id[0]] = 1
			self.exclude[least_id[1]] = 1
			
			new_node = UPGMANode(node2, node1)
			nodes.append(new_node)
			self.exclude[len(nodes)-1] = 0
			node1.up_dist = least_dist / 2 - node1.down_dist
			node2.up_dist = least_dist / 2 - node2.down_dist
			new_node.down_dist = least_dist / 2
			
			# create new working distance matrix
			self.update_distance( nodes, least_id)
			
		self.tree = new_node 
		
	def update_distance(self, nodes: list, least_id: tuple) -> np.ndarray:
		
		length = len(nodes)
		nleaves1, nleaves2 = nodes[least_id[0]].nleaves, nodes[least_id[1]].nleaves
		nleaves  = nleaves1 + nleaves2 
		
		for i in range(length-1):
			
			if self.exclude[i]:
				continue
			
			
			self.work_matrix[i,length-1] = ( self.work_matrix[i][least_id[0]]*nleaves1 + self.work_matrix[i][least_id[1]]*nleaves2 ) / nleaves
			self.work_matrix[length-1,i] = self.work_matrix[i,length-1]
			
def projectiontree(matrix, header):
	
	if len(header) < 2:
		return header[0] + ";"
	
	matrix = np.matrix(matrix)
	
	dm = np.ones((matrix.shape))
	vars = [matrix[i,i] for i in range(len(matrix))]
	for i in range(len(matrix)):
		
		var1 = max(1.0,vars[i])
		for j in range(i):
			
			var2 = max(1.0,vars[j])
			
			dm [i,j] = 1 - matrix[i,j]/math.sqrt((var1*var2))
			dm [j,i] =  dm [i,j]
			
	tree = UPGMA(dm, header).tree
	tree.reorder(dm)
	
	return tree.to_newick()+";"


def ReorderFile(seqfile, treeorder, name):
	
	with open(seqfile,mode = 'r') as r:
		
		reads = r.read().split(">")[1:]
		
	reordered = [reads[index].splitlines() for index in treeorder]

	if len(name):
		reordered = [name+read[0].replace(" ","\t")+"\n"+"\n".join(read[1:]) for read in reordered]
	else:
		reordered = [read[0].replace(" ","\t")+"\n"+"\n".join(read[1:]) for read in reordered]
	
	with open(seqfile,mode = 'w') as w:
		
		w.write(">"+"\n>".join(reordered)+"\n")
		
	header = [read.splitlines()[0].replace(" ","\t") for read in reordered]
	
	return header

def ReorderMatrix(matrix, order):
	
	
	matrix_rshape = np.reshape(matrix, (len(order), len(order)))
	
	#order_sort = sorted(range(len(order)), key = lambda x: order[x])
	
	newmatrix = []
	for reorderindex, index in enumerate(order):
		
		row = matrix_rshape[index].tolist()
		row_rshape = [row[i] for i in order]
		
		newmatrix.append(row_rshape)
		
	return newmatrix


def MatrixFulltoUpper(matrix, size):
	
	newmatrix = []
	
	for i in range(size):
		
		newmatrix += matrix[(i*size + i) : (i*size + size) ] 
		
	return newmatrix

def ReorderNewickname(text):
	
	locations = [x.span() for x in re.finditer(r'[^(^)^:^,^;]+', text) if text[x.span()[0]-1] != ":"]
	
	last_location = 0
	newtext = ""
	nameindex = 0
	for location in locations:
		
		newtext += text[last_location:location[0]] + str(nameindex)+"_"
		
		last_location = location[1]
		
		nameindex += 1
		
	newtext += text[last_location:]
	
	return newtext

def hashrow(kmerindex , genenum, KmerReader ):
	
	row = tuple(sorted(list(KmerReader.matrix.getkmerrow( kmerindex ) ) ) )
	
	return ( - (abs( len(row) -  genenum/2 - 0.1 ))  , hash(row) )


class KmerData:
	
	def __init__(self, seqfile,  kmersize = 31, msize = 0, leave = []):
		
		self.exclude = set(leave)
	
		self.kmersize = kmersize
		
		self.kmerslist = dict()
		
		self.genekmercounts = []
		
		self.samplesizes = []   
		
		index = 0
		name = ""
		samples = []
		with open(seqfile,mode = 'r') as r:
			for line in r:
				
				line = line.strip()
				
				if len(line):
					if line[0]==">":
						name = line.split()[0][1:]
						if len([x for x in self.exclude if len(x) and x in name]):
							name = ""
							continue
						self.samplesizes.append(0)
						samples.append(name)
						self.genekmercounts.append( 0 )
						index += 1
					else:
						if name != "":
							self.samplesizes[-1] += len(line.strip())
						
						
		self.sampleslist = list(samples)
		
		
		if msize < 20000000:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)
		else:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)
		
	def LoadKmers(self, kmerfile):
		
		kmerindex = 0
		with open(kmerfile, mode = 'r') as f:
			
			for line in f:
				
				if len(line) and line[0] != ">":
					
					self.matrix.addkmer()
					self.kmerslist[kmerencode(line.strip().split()[0])] = kmerindex
					kmerindex += 1
					
		
		
	def AddKmer(self,kmer, index):
		
		kmerindex = self.kmerslist.get(kmer, -1) 
		
		if kmerindex >= 0:
			
			self.genekmercounts[index] += 1
			
			self.matrix.add(kmerindex, index)
			
		return kmerindex
	
	
	def ReadKmers(self, seqfile, outputfile):
		
		intoperator = 4**(self.kmersize-1)
		
		current_k = 0
		reverse_k = 0
		current_size = 0
		
		qname = ""
			
		currentline = ""
		
		with open(outputfile,mode = 'w') as w:
			
			with open(seqfile,mode = 'r') as r:
				
				posi_linestart = 0
				index = -1
				for line in r:
					
					if len(line) ==0:
						continue
					
					if line[0]==">":
						current_k=0
						reverse_k=0
						current_size = 0
						posi_linestart = 0

						qname = line.strip().split()[0][1:]
						if len([x for x in self.exclude if len(x) and x in qname]):
							qname = ""				
							continue
	
						index += 1
						
						w.write(line)
						continue
	
					if qname == "":
						continue

					line_out = bytearray(line.lower().encode())
					for posi,char in enumerate(line.upper()):
						
						if char =="\n":
							continue
						
						if char not in ['A','T','C','G']:
							
							current_k = 0
							reverse_k = 0
							current_size = 0
							
							continue
						
						if current_size >= self.kmersize:
							
							current_k %= intoperator
							current_k <<= 2
							current_k += ['A','C','G','T'].index(char)
							
							reverse_k >>= 2
							reverse_k += (3-['A','C','G','T'].index(char)) * intoperator
							
						else:
							
							
							current_k <<= 2
							current_k += ['A','C','G','T'].index(char)
							reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*current_size))
							
							current_size += 1
							
						if current_size >= self.kmersize :  
							
							kmer = max(current_k, reverse_k)
							
							findindex = self.AddKmer(kmer,index)
							
							if findindex >= 0:
								
								line_out[posi] = line_out[posi] - 32
								
					posi_linestart += len(line)
					
					line_out = line_out.decode()
					
					w.write(line_out)

def annotatefasta(seq, kmer, outputfile, pref = None, leave = []):

	if pref is None:
		name = "_".join(seq.split("/")[-1].split("_")[:2])
	elif len(pref):
		name = pref+"_"
	else:
		name = pref+"_"

	if len(outputfile ) == 0:
		outputfile = seq+"_annotated.fasta"	


	kfilesize = os.stat(kmer).st_size
	
	KmerReader = KmerData(seq, msize = kfilesize, leave = leave)

	KmerReader.sampleslist = [name+x for x in KmerReader.sampleslist]

	KmerReader.LoadKmers(kmer)
	KmerReader.ReadKmers(seq,outputfile)


	
	
	sqmatrix = np.reshape( KmerReader.matrix.SquareMatrix(0), (len(KmerReader.sampleslist), len(KmerReader.sampleslist)) )

	treetext = projectiontree( sqmatrix, KmerReader.sampleslist)
	
	treeorder = [x.split(":")[0].replace("(", "") for x in treetext[:-1].split(",")]
	
	nametoindex = {x:i for i,x in enumerate(KmerReader.sampleslist)}
	
	treeorder = [nametoindex[sample] for sample in treeorder]
	
	headers = ReorderFile(outputfile, treeorder, name)
	
	newnorm = np.array(ReorderMatrix(sqmatrix, treeorder))

	return headers,treetext,newnorm	
	

def main(args):


	headers, treetext, newnorm = annotatefasta(args.seq, args.kmer, args.output, None if args.pref == '>' else args.pref, args.leave.split(","))

	with open(args.output + "_tree.ph", mode = 'w') as f:

		f.write(treetext + "\n")

	np.savetxt(args.output + "_norm.txt", newnorm, delimiter=',', fmt='%1.4f')	
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="seq", type=str, required = True)
	parser.add_argument("-k", "--kmer", help="path to input data file",dest="kmer", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to input data file", dest="output", type=str, default = "")
	parser.add_argument("-p", "--pref", help="path to input data file", dest="pref", type=str, default = '')
	parser.add_argument("-l", "--leave", help="path to input data file", dest="leave", type=str, default = '')
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
