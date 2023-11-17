#!/usr/bin/env python3

import pandas as pd
import os 
import argparse
import sys
import numpy as np
import multiprocessing as mul
import math
import collections as cl


def withinvar(matrix, indexes):
	
	avars = []
	rvars = []
	for index, i in enumerate(indexes):

		size1 = max(1.0 , matrix[i,i])
		#row = dm[ i]
		
		for j in indexes[:index]:
			
			size2 = max(1.0, matrix[j,j])
			match = matrix[i,j]
			dist = 1 - match/math.sqrt((size1 * size2))
			rvars.append(dist)
			avars.append( max( size1 ,  size2 ) -  match )
	
	return rvars,avars	

def withoutvar(matrix, indexes1, indexes2):
	
	avars = []
	rvars = []
	for i in indexes1:

		size1 = max(1.0 , matrix[i,i])		
	
		#row = dm[i,...]
		
		for j in indexes2:
		
			size2 = max(1.0, matrix[j,j])
			match = matrix[i,j]
	
			dist = 1 - match/math.sqrt((size1 * size2)) 
			rvars.append(dist)
			avars.append(  max( size1 ,  size2 ) -  match  )
			
	return rvars, avars
	

class node:
	
	def __init__(self, parent = None, name = "", distance =0.0, index = 0):
		
		self.children = []
		
		self.parent = parent
		
		self.name = name
		
		self.distance = distance
		
		self.index = index
		
		self.annotation = 0
		
		if parent is not None:
			
			parent.children.append(self)
			
	def __str__(self):
		
		if len(self.children) == 0:
			
			return self.name+":"+str(self.distance)
		
		else:
			return "("+",".join([str(x) for x in self.children])+")"+":"+str(self.distance)
	
	def push(self, name = "", distance =0.0, index = 0):
		
		newchild = node(self,name,distance,index)
		
		return newchild
	
	def build(self,text):
		
		if text == "":
			
			return self
		
		elif text[-1] == ";":
			
			text = text[:-1]
			
		current_node = self
		allnames = []
		
		index = 0
		current_name = ""
		for char in text:
			
			if char in [" ", "\'"]:
				
				continue
			
			if char == "(":
				
				current_node = current_node.push()
				
			elif char == ")": 
				
				name,distance = (current_node.name.split(":")+["0.0"])[:2]
				
				name = name.split()[0].split(" ")[0].split("\t")[0]
				
				if len(name)>0:
					
					current_node.name = name
					
					allnames.append(name)
					
					if len(current_node.children) == 0:
						current_node.index = index
						index += 1
					
				else:
					
					allnames.append(current_node.name)
					
				try:
					current_node.distance = float(distance.strip())
					
				except:
					
					current_node.distance = 0.0
					
				current_node = current_node.parent
				
				if len(current_node.name) == 0 :
					
					current_node.name = "bh_"+str(len(allnames))
					
				else:
					print(current_node.name)
				
					
			elif char == "," or char ==";":
				
				name,distance = (current_node.name.split(":")+["0.0"])[:2]
				
				name = name.split()[0].split(" ")[0].split("\t")[0]
				
				if len(name)>0:
					
					current_node.name = name
					
					allnames.append(name)
					
					if len(current_node.children) == 0:
						current_node.index = index
						index += 1
					
				else:
					
					allnames.append(current_node.name)
					
				try:
					current_node.distance = float(distance.strip())
					
				except:
					
					current_node.distance = 0.0
					
					
				if current_node.parent is not None:
					current_node = current_node.parent.push()
					
			else:
				
				current_node.name += char
				
		
		return self
	
	def alloffsprings(self, offset = 0.0):
		
		alloffsprings = [(offset,self)]
		
		for child in self.children:
			
			alloffsprings.extend(child.alloffsprings(offset = child.distance+offset))
			
		return alloffsprings
	
	def allleaves(self):
		return [x for x in self.alloffsprings() if len(x[1].children) == 0]	
	def magnitude(self, exclude = None):
		
		if exclude is None or exclude is self:
			
			alloffsprings = [x[0] for x in self.alloffsprings()]
			
		else:
			
			exclude_offsprings = [x[1] for x in exclude.alloffsprings()]
			
			alloffsprings = [x[0] for x in self.alloffsprings() if x[1] not in exclude_offsprings]
			
			
		return max(alloffsprings+[0.0])


	def findclades(self, dm , cutoff = 0.005, rcutoff = 1.5, acutoff = 200):
		
		if len(self.children) == 0:
			return []
		
		clades = sum([child.findclades(dm, cutoff, rcutoff) for child in self.children], [])

		lnames = [x[1].name for x in self.children[0].allleaves()]
		rnames = [x[1].name for x in self.children[1].allleaves()]
		
		lhaplo = set([tuple(name.split("_")[-3:-1]) for name in lnames])
		rhaplo = set([tuple(name.split("_")[-3:-1]) for name in rnames])

		if lhaplo.intersection(rhaplo):
			rcutoff = 1.5
	
		lindexes = [x[1].index for x in self.children[0].allleaves()]
		rindexes = [x[1].index for x in self.children[1].allleaves()]
				

		lsize = sum([dm[i,i] for i in lindexes])/max(1,len(lindexes))
		rsize = sum([dm[i,i] for i in rindexes])/max(1,len(rindexes))

	
		without,awithout = withoutvar(dm, lindexes, rindexes)

		if len(lindexes) * len(rindexes)  < 10 or min( len(lindexes),  len(rindexes)) < 3:


			if len([x for x,y in zip(without,awithout) if x < 3*cutoff or y < 3*acutoff]) == 0:

				clades += self.children
				
			return clades

		lwithin,awithin = withinvar(dm, lindexes)
		rwithin,awithin = withinvar(dm, rindexes)

		
		within_mean = ( sum(lwithin) + sum(rwithin) ) / max(len(lwithin) + len(rwithin), 1)
		without_mean = sum(without) / len(without)

		if ( min(awithout) >= acutoff  and  without_mean > rcutoff * within_mean ) :
			clades += self.children 
		
		return clades
		
			
	

		
	

def distance_search(root, refs, samples,phyly_groups):
	
	if len(phyly_groups) == 0:
		
		return []
	
	
	refnodes = [root.search(ref) for ref in refs]
	
	results = []
	for group in phyly_groups:
		
		if group in refnodes:
			
			continue
		
		if type(group) == type([]):
			
			refnodes_within_group = [refnode[1] for node in group for refnode in node.get_offsprings(refs)]
			
		else:
			
			refnodes_within_group = [x[1] for x in group.get_offsprings(refs)]
			
		if len(refnodes_within_group) < 1:
			
			continue
		
		if type(group) == type([]):
			
			result = [node for mono_group in group for node in mono_group.distinctive_nodes(refnodes_within_group) if len(node.allhaplotypes(samples))>0]
			
		else:
			
			result = [node for node in group.distinctive_nodes(refnodes_within_group)]
			
		results.extend(result)
		
		
	dis_clusters = []
	for index,  mono in enumerate(results):
		
		mono_cluster = cluster("dis_"+str(index))
		
		mono_cluster.roots.append(mono)
		
		mono_cluster.members = [x[1] for x in mono.alloffsprings()]
		
		mono_cluster.label_members()
		
		dis_clusters.append(mono_cluster)
		
		
	return dis_clusters


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
	
	
	inputfile = args.input
	outputfile = args.output
	
	normfile = args.input + "_norm.txt"
	treefile = args.input + "_tree.ph"

	
	excludenames = set()
	if len(args.exclude):
		with open(args.exclude, mode = 'r') as f:
			excludenames = set([x.strip() for x in f.read().splitlines()])


	normmatrix = np.genfromtxt(normfile, delimiter=',')

	#distmatrix= np.array(normtotdm(normmatrix))

	headers = dict()
	with open(args.input, mode = 'r') as f:
		for line in f:
			if len(line) ==0 or line[0] != ">":
				continue
			lines = line.strip().split()
			headers[lines[0][1:]] = lines[1]
			
	
	

	with open(treefile, mode = 'r') as f:
		
		tree_text = f.read()
	
	tree = node().build(tree_text)
	
	allclades = tree.findclades(normmatrix)


	for i, clade in enumerate(allclades[::-1]):
		
		allleaves = clade.allleaves()
		
		for (dist, leave) in allleaves:
							
			leave.annotation = i + 1
	 
	
	allnames = [x[1].name for x in tree.allleaves()]
	allgroups = [x[1].annotation for x in tree.allleaves()]
	
	newgroups = []
	uniq_group = set([])
	for group in allgroups:
		
		uniq_group.add(group)
		newgroups.append(len(uniq_group))
	
	namegroup = {name: group for name, group in zip (allnames, newgroups) }
	
	groups = [[] for x in range(max(namegroup.values()) + 1)]
	
	for name in allnames:

		groups[namegroup[name]].append(name)
	
	with open(outputfile, mode = 'w') as f:
		
		for group in groups:
			
			f.write( ",".join([x for x in group if headers[x].split(":")[0] not in excludenames]) +"\n" )
	
	
	return 0

def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str, required=True)
	parser.add_argument("-e", "--error", help="path to output file", dest="exclude", type=str, default= "" )
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
