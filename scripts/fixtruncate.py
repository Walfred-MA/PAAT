#!/usr/bin/env python3


import argparse 
import collections as cl
import numpy as np
from statistics import mode

mainchrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']



def combinebreaks(breaks):
	
	all_coordinates = [x for y in breaks for x in y[:2]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	break_groups = []
	break_group = []
	
	break_group_index = []
	break_group_indexs = []
	
	last_coordinate = 0
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if coordinate - last_coordinate < 500:
			
			break_group.append(coordinate)
			break_group_index.append(index)
			
		else:
			
			break_groups.append(break_group)
			break_group_indexs.append(break_group_index)
			
			break_group = [coordinate]
			break_group_index = [index]
			
		last_coordinate = coordinate
		
	break_groups.append(break_group)
	break_group_indexs.append(break_group_index)
	
	return break_groups, break_group_indexs

def getchunks_frombreak(break_groups, break_group_indexs, breaks):
	
	all_queries = [y[2] for y in breaks for x in y[:2]]
	
	chunk_starts = []
	breaks_onquery = cl.defaultdict(list)
	
	chunk_index = 0 
	for (break_group, break_group_index) in zip(break_groups, break_group_indexs):
		
		if break_group == []:
			continue
		
		themode = mode(break_group)
		
		chunk_starts.append(themode)
		
		for coodinate_index in break_group_index:
			
			breaks_onquery[all_queries[coodinate_index]].append(( coodinate_index, chunk_index))
			
		chunk_index += 1
		
	chunks = [ [first, second] for first, second in zip(chunk_starts,chunk_starts[1:])]
	
	chunks_onquery = {query: list(range(breaks[0][1], breaks[1][1])) for query, breaks in breaks_onquery.items()}
	
	
	return chunks, chunks_onquery

def getbreaks(map_locs,labels):
	
	locus_bychr = cl.defaultdict(list)
	for i,locus in enumerate(map_locs):
		
		if ":" not in locus:
			continue
		contig = locus.split(":")[0]
		thebreak = [int(x) for x in locus.split(":")[1][:-1].split('-')]
		locus_bychr[contig].append(thebreak+[i])
		
		
	clusters = cl.defaultdict(list)
	for contig, breaks in locus_bychr.items():
		
		coordinates = [x for y in breaks for x in y[:2]] 
		
		coordinates_sortindex = sorted( range(len(coordinates)), key = lambda x: coordinates[x])
		
		cover = 0
		all_clusters = []
		current_cluster = []
		last_coordi = -100000
		last_label = "XXXX"
		
		for index in coordinates_sortindex:
			
			index2 = index // 2
			current_coordi = coordinates[index]
			current_label = labels[index2]
			if index % 2 == 0:
				cover += 1
				
				if cover == 1 and current_coordi - last_coordi > 10000 and last_label != current_label:
					all_clusters.append(current_cluster)
					current_cluster = [index2]
				else:
					current_cluster.append(index2)
					
			else:
				cover -= 1
				
			last_coordi = current_coordi
			last_label = current_label
			
		all_clusters.append(current_cluster)
		
		for i,cluster in enumerate(all_clusters[1:]):
			
			clusters[(contig,i)] = [breaks[x] for x in cluster]
			
			
	all_chunks_onquery = [[] for x in map_locs]
	for contig, breaks in clusters.items():
		
		alllabels = set([x for y in [labels[x[2]].split(";") for x in breaks] for x in y])
		
		alllabels = ";".join(alllabels)
		
		if alllabels == "NA":
			alllabels = "NA_{}_{}".format(contig[0],contig[1])
			
			
		break_groups, break_group_indexs = combinebreaks(breaks)
		
		chunks,  chunks_onquery = getchunks_frombreak(break_groups, break_group_indexs, breaks)
		
		used_chunkindex = dict()
		for index, chunkindex in chunks_onquery.items():
			
			if len(chunks) < 2:
				chunkindex = -1
				all_chunks_onquery [index] = (alllabels , -1)
			else:
				for chunkindex0 in chunkindex:
					if chunkindex0 not in used_chunkindex:
						used_chunkindex[chunkindex0] = len(used_chunkindex) + 1
						
				all_chunks_onquery [index] = (alllabels , [used_chunkindex[chunkindex0] for chunkindex0 in chunkindex])
				
	return all_chunks_onquery

	
	
def getminmatch(dm, index, contig_indeces, all_loc , types,  group_largestdist):
	
	size = len(all_loc)
	
	work_row = np.asarray(dm[index])
	
	hg38_chunks = []
	mapped_chunks = []
	for contig, indeces in contig_indeces.items():
		
		if len(indeces) < 2 and contig not in mainchrom:
			continue
		
		if contig in mainchrom:
			
			startindex = len(mapped_chunks)
			
		indeces = sorted(indeces, key = lambda x: int(all_loc[x].split(":")[1].split('-')[0]))
		
		thesum = np.zeros(size)
		
		lastdist = 100000000000000
		mapped_chunk = []
		for index2 in indeces:
			row = np.asarray(dm[index2])
			new_sum = thesum + row 
			newdist = np.sum((new_sum - work_row) ** 2)
			rowdist = np.sum((row - work_row) ** 2)
			
			if (newdist <= lastdist and newdist <= rowdist) :
				lastdist = newdist
				thesum = new_sum
				mapped_chunk.append(index2)
			else:  
				mapped_chunks.append([ mapped_chunk, lastdist ])
				if newdist > rowdist:
					mapped_chunk = [index2]
					thesum = row
					lastdist = rowdist
				else:
					thesum = np.zeros(size)
					mapped_chunk = []
					lastdist = 100000000000000	
					
					
		mapped_chunks.append([mapped_chunk,lastdist])
		
		if contig in mainchrom:
			
			hg38_chunks.extend(mapped_chunks[startindex:])
			
	hg38_chunks = [x for x in sorted(hg38_chunks, key = lambda x:x[1])]
	hg38_chunks = hg38_chunks[0] if len(hg38_chunks) else [[],"NA"]
	
	mapped_chunks = [x[0] for x in sorted(mapped_chunks, key = lambda x:(-len(x[0]),x[1] )) if len(x[0]) > 1 and x[1] < sum([group_largestdist[types[i]] for i in x[0]])]
	
	
	return mapped_chunks,hg38_chunks 

def normtodistm(norm, size):
	
	index = 0
	dm = [10000000000 for x in range(size**2)]
	
	for i in range(size):
		
		var1 = norm[i*size + i]
		
		for j in range(i, size):
			
			var2 = norm[j*size + j]
			
			themax = max(var1,var2)
			themin = min(var1,var2)
			match = norm[j*size + i]
			
			dm[j*size+i] = ( match, var1,var2)
			dm[i*size+j] = dm[j*size+i]
			
			
	return dm

def organize_truncate(inputfile, outputfile, normfile):
	
	size = 0
	dmatrix = []
	with open(normfile, mode = 'r') as f:
		
		for line in f:
			
			line = [float(x) for x in line.strip().split(",")]
			dmatrix.append(line)
			size += 1
			
			
			
	dm = np.matrix(dmatrix)
	
	map_locs = []
	all_loc = []
	labels = []
	index = 0
	contig_indeces = cl.defaultdict(list)
	types_indeces = cl.defaultdict(list)
	types = []
	names = []
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.strip().split('\t')
			
			names.append(line[0])
			if "_h1" or "_h2" in line[0]:
				haplo = "_".join(line[0].split('_')[2:4])
				
			types.append(line[1])
			
			types_indeces[line[1]].append(index)
			
			contig_indeces[line[7].split(":")[0]].append(index)
			
			all_loc.append(line[7])
			
			labels.append(line[2])
			
			map_locs.append(line[10])
			index += 1
			
			
	all_chunks_onquery = getbreaks(map_locs, labels)
	
	group_largestdist = dict()
	
	for group, indeces in types_indeces.items():
		
		alldists = []
		
		for index,i in enumerate(indeces):
			
			var1 = dmatrix[i]
			for j in indeces[(index+1):]:
				
				var2 = dmatrix[j]
				
				alldists.append(np.sum( ( np.asarray(var1)  - np.asarray(var2) ) ** 2 ) )
				
		group_largestdist[group] = max(alldists+ [300])
		
	records = []
	
	index = 0
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.strip().split('\t')
			
			tag = line[5]
			
			refgene = line[11] if line[11] != "NA" else line[8]
			
			mapchunks,hg38chunks = getminmatch(dm, index, contig_indeces, all_loc , types,  group_largestdist)
			
			records.append([mapchunks,hg38chunks])
			
			index += 1
			
			
	newrecords = [[] for x in records]
	for i, record in enumerate(records):
		
		if len(record[0]) == 0:
			continue
		
		newrecord = []
		for record0 in record[0]:
			
			if len(record0) < 2:
				
				newrecord.append(record0)
				continue
			
			trace = [x for x in record0]
			
			
			newtraces = sum([],[records[x][0][0] for x in trace if len(records[x][0])])
			
			while len(newtraces) > len(trace):
				
				trace = newtraces
				newtraces = sum([],[records[x][0][0] for x in trace])
				
			newrecord.append(sorted(trace))
			
		newrecords[i] = newrecord
		
	results = []
	index = 0
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.strip().split('\t')
			
			tag = line[5]
			
			refgene = line[11] if line[11] != "NA" else line[8]
			
			mapchunks,hg38chunks = records[index]
			
			maptypes = list(dict.fromkeys([tuple(sorted([names[i] for i in x])) for x in mapchunks]))
			
			chunkinf = all_chunks_onquery[index]
			
			if len(chunkinf) == 0:
				chunkinf = ""
			elif chunkinf[1] == -1:
				chunkinf = chunkinf[0] 
			else:
				chunkinf = "&".join( [ chunkinf[0]+"(part{})".format(str(x)) for x in chunkinf[1] ] )
				
			maptype_str = "&".join([";".join(x) for x in maptypes])
			if maptype_str == "":
				maptype_str = "NA"
				
				
			hg38chunks_str = ";".join([names[x] for x in hg38chunks[0]]) 
			if hg38chunks_str == "":
				hg38chunks_str = "NA"
				
				
			newmaptypes = list(dict.fromkeys([tuple(sorted([names[i] for i in x])) for x in newrecords[index]]))
			newmaptypes_str = "&".join([";".join(x) for x in newmaptypes])
			if newmaptypes_str == "":
				newmaptypes_str = "NA"
				
			newmaptypes2 = list(dict.fromkeys([tuple(sorted([types[i] for i in x])) for x in newrecords[index]]))
			newmaptypes2_str = "&".join([";".join(x) for x in newmaptypes2])
			if newmaptypes2_str == "":
				newmaptypes2_str = "NA"
				
			outline = [line[0], maptype_str , newmaptypes_str ,  newmaptypes2_str ,hg38chunks_str  , str(hg38chunks[1]), chunkinf]
			results.append("\t".join(outline) + "\n")
			
			index += 1
			
			
			
			
	with open(outputfile, mode = 'w') as f:
		for result in results:
			f.write(result)
			
			
			
			
			
def main(args):
	
	organize_truncate(args.input, args.output, args.norm)
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-n", "--norm", help="path to input data file",dest="norm", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

