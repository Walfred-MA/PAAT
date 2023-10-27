#!/usr/bin/env python3


import pandas as pd
import os 
import argparse
import sys
import collections as cl
import re

def makereverse(seq, strd = "-"):
	
	if strd == "+":
		
		return seq
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)


def readfile(querypath):
	
	read = ""
	header = ""
	query = {}
	with open(querypath, mode = 'r') as f:
		for line in f:
			if len(line) == 0:
				continue
			if line[0]==">":
				query[header] = read
				read = ""
				header = line[1:].split()[0]
			else:
				read+=line.strip()
				
	query[header] = read
	read = ""
	
	
	return query

def export_fasta(querypath, infos, outfile_ ):
	
	
	#title = ">{}_{}\t{}:{}-{}{} {}".format(samplename, i, locus[0].split(":")[0], locus[1], locus[2],locus[3],",".join(locus[-1]))
	
	allloci = cl.defaultdict(list)
	for locus in infos:
		allloci[locus[0].split(":")[0]].append(locus)
		
		
	read = ""
	header = ""
	query = {}
	ifstorage = 0
	
	with open(querypath, mode = 'r') as f:
		for line in f:
			if len(line) == 0:
				continue
			if line[0]==">":
				
				if ifstorage:
					
					for locus in allloci[header]:
						
						samplename = "_h".join(locus[0].split("#")[:2])	
						i = locus[-1]
						title = ">{}_{}\t{}:{}-{}{}\t{}".format(samplename, i, locus[0].split(":")[0], locus[1], locus[2],locus[3],",".join(locus[-2]))
						
						seq = read[max(0,locus[1]-1):(locus[2])] if locus[3] == '+' else makereverse(read[max(0,locus[1]-1):(locus[2])])
						
						outfile_.write("{}\n{}\n".format(title, seq ))
						
				
				read = ""
				header = line[1:].split()[0]
				
				ifstorage = 0
				if header in allloci:
					ifstorage = 1
			else:
				if ifstorage:
					read+=line.strip()
	
	if ifstorage:
		
		for locus in allloci[header]:
			
			samplename = "_h".join(locus[0].split("#")[:2])	
			i = locus[-1]
			title = ">{}_{}\t{}:{}-{}{}\t{}".format(samplename, i, locus[0].split(":")[0], locus[1], locus[2],locus[3],",".join(locus[-2]))
			
			seq = read[max(0,locus[1]-1):(locus[2])] if locus[3] == '+' else makereverse(read[max(0,locus[1]-1):(locus[2])])
			
			outfile_.write("{}\n{}\n".format(title, seq ))
		
	read = ""
	
	
	return query
	
	
def fixovermerge(aligns, exonposi, cutoff = 10000, samegene_cutoff = 20000):

	all_coordinates = [x for y in exonposi for x in y[:2]]
	all_names = [y[2].split("_")[-1] for y in aligns]


	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])

	gene_ends = cl.defaultdict(lambda: -samegene_cutoff - 1)
	for index in sort_index:

		coordinate = all_coordinates[index]
		genename = all_names[index//2]
		gene_ends[genename] = coordinate

	lastgene_coordinate = cl.defaultdict(lambda: -samegene_cutoff - 1)
	lastgene_coordinate[""] = -samegene_cutoff - 1
	last_coordinate = -cutoff - 1

	numchunks = 0
	allmidpoints = [ ]

	for index in sort_index:
					
		coordinate = all_coordinates[index]
		genename = all_names[index//2]
	
		if index %2 == 0 :

			
			if numchunks == 0 and coordinate - last_coordinate> cutoff and coordinate - max(list(lastgene_coordinate.values())) > samegene_cutoff:
	
				midpoint =  ( coordinate + last_coordinate ) // 2 

				allmidpoints.append( midpoint )

			numchunks += 1
		
		else:
			
			last_coordinate = coordinate

			if coordinate >= gene_ends[genename]:
				lastgene_coordinate[genename] = -samegene_cutoff - 1
			else:
				lastgene_coordinate[genename] = coordinate 

			numchunks -= 1
					

	all_coordinates = [x for y in aligns for x in y[:2]] + allmidpoints

	lintron = len(aligns) * 2
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])

	allchunks = [ set([]) for x in allmidpoints]


	chunkindex = -1
	for index in sort_index:
		
		if index >= lintron :
			chunkindex += 1		

		else:
			allchunks[chunkindex].add( index//2 )
				
	return [list(x) for x in allchunks], allmidpoints
	

def querycontig_overlap(allaligns, allow_gap = 0, anchor_size = 100):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]

	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group = set([])
	current_group_size = 0 
	current_group_start = 0
	allgroups_sizes = []
	
	number_curr_seq = 0
	coordinate = 0

	last_coordinate = -allow_gap-1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if index %2 == 0 :
			
			number_curr_seq += 1
			
			if number_curr_seq == 1:
				
				current_group_start = coordinate
				
				if coordinate - last_coordinate> allow_gap :
					
					allgroups_sizes.append(current_group_size)
					allgroups.append(current_group)
					
					current_group = set([])
					current_group_size = 0
					
			current_group.add(index//2)
			
		else:
			
			number_curr_seq -= 1
			
			if number_curr_seq == 0:
				
				current_group_size += coordinate - current_group_start  
				
		last_coordinate = coordinate
		
		
	allgroups_sizes.append(current_group_size)
	allgroups.append(current_group)
	
	allgroups_sizes = allgroups_sizes[1:]   
	allgroups = allgroups[1:]
	

	results = []
	for i,group in enumerate(allgroups):
		
		aligns = [allaligns[x] for x in group]

		refgenesizes = [x[4] for x in aligns]
		exonposi = [[x[0], x[1], x[3]] for x in aligns if x[3][:4] == "ENSE"]

		range1 = sorted(exonposi) 

		#print(aligns)
		#print( [x for x in aligns if x[3][:5] == "ENSE0"] )

		if len(exonposi):
			chunks,midpoints = fixovermerge(aligns, exonposi )

		else:
			chunks = [list(range(len(aligns)))]
			midpoints = []
				
			coordinates = [x for y in aligns for x in y[:2]]
			start = max(0,min(coordinates))
			end = max(coordinates)
			genes = list(set([x[3]  for x in aligns if x[3][:5] != "ENSE" ]))

			sizes = [x[1]-x[0] for x in aligns]
			strands = [1 if x[2]=='+' else -1 for x in aligns]

			strandsizes = sum([strd * size for strd, size in zip(strands, sizes)])
			strand = ['-','+'][int(strandsizes >= 0)]

			if end - start > 500:
				results.append([start, end, strand, genes])



			continue

		midpoints.append(100000000000000000)
		lastmidpoint = 0

		#print("mid", midpoints)

		for chunk,midpoint in zip(chunks, midpoints[1:]):
		
			ranges_thischunk = [aligns[x] for x in chunk]
			coordinates = [x for y in ranges_thischunk for x in y[:2]]

			chunk_exonposi = sum([[x[0], x[1]] for x in ranges_thischunk if x[3][:4] == "ENSE"],[])
			chunk_exonposi_range = [min(chunk_exonposi), max(chunk_exonposi) ]


			start = min(coordinates)
			start = max(lastmidpoint, start)

			end = max(coordinates)
			end = min (end, midpoint)

			
			if chunk_exonposi_range[0] < start + anchor_size // 2:
				
				start = max (0, chunk_exonposi_range[0] - anchor_size // 2)

			if chunk_exonposi_range[1] > ( end - anchor_size // 2 ):
			
				end = chunk_exonposi_range[1] + anchor_size // 2

			lastmidpoint  = midpoint
			
			sizes = [x[1]-x[0] for x in ranges_thischunk]
			strands = [x[2] for x in ranges_thischunk]
			genes = list(set([x[3]  for x in ranges_thischunk if x[3][:5] != "ENSE" ]))
			
			if len(genes) == 0:
				genes = ["ENSE"]
				
			#start = start-anchor_size
			#end = end+anchor_size  
				
				
			strandsizes = {"+":0,"-":0}
			for size, strand in zip(sizes, strands):
				
				strandsizes[strand]+=size
				
			if strandsizes["-"] > strandsizes["+"]:
				strand = '-'
			else:
				strand = '+'
		
			#if end - start > 30000 and len(genes) < -1:	
				#print(start, end)
				#print(midpoints)

				#exonposi_s = sorted(exonposi)
				#print(genes)	

				#print(max([abs(y[0]-x[1]) for x,y in zip(exonposi_s, exonposi_s[1:]) ]))
	
			results.append([start, end, strand, genes])

			#print([start, end, strand, genes])

	return results




def refgene_overlap(alloverlapped):
	
	gene_groups = []
	iffind = set([])
	for overlap in alloverlapped:
		
		
		genelist = set(overlap[1])
		
		querylist = overlap[2]
		
		newgroup = [[overlap[0]],genelist, [querylist]]
		
		
		for index,group in enumerate(gene_groups):
			
			if index in iffind:
				
				continue
			
			if len(newgroup[1].intersection(group[1])):		
				
				
				newgroup[2].extend(group[2])
				newgroup[1].update(group[1])
				newgroup[0].extend(group[0])
				
				iffind.add(index)				
				
		gene_groups.append(newgroup)
		
		
	gene_groups = [group for index,group in enumerate(gene_groups) if index not in iffind]
	
	return gene_groups


def main(args):
	
	inputpath = args.input
	outpath = args.output
	querypath = args.query
	refgenepath = args.refgenes
	hg38,chm13 = args.ref.split(",")
	
	refgenes_titles = {}
	refgenes_sizes = {}
	name  = ""
	with open(refgenepath, mode = 'r') as f:
		
		for line in f:
			if len(line.strip()) == 0:
				continue
			if line[0]==">":
				name = line.split()[0].split("_")[0]
				refgenes_titles[name] = line
				refgenes_sizes[name] = 0
			else:
				refgenes_sizes[name] +=len(line.strip())
				
				
	table = pd.read_csv(inputpath, header = None , sep = "\t")
	
	eachcontig = cl.defaultdict(list)
	refcontig = cl.defaultdict(list)
	
	allgenes = {}
	for gene,title in refgenes_titles.items():
		
		title = title.strip()
		name = title.split()[0] + "_" + title.split()[-1]

		contig = title.split()[1].split(":")[0]
		
		therange = title.split()[1].split(":")[1][:-1]
		sign = 1
		if therange[0] == '-':
			therange = therange[1:]	
			sign = -1

		start,end = map(int,therange.split("-")) 
		start *= sign
	

		refcontig[contig].append([start,end,title.split()[1][-1],name,refgenes_sizes[name.split("_")[0]]])
		
		allgenes[gene] = len(allgenes)
		
	allgenes["ENSE0"] = 0
	
	for index, row in table.iterrows():
		
		contig = row.iat[0]
		
		gene = row.iat[5]
		
		if gene not in allgenes:
			allgenes[gene] = len(allgenes)
			
		length = row.iat[6]
		
		eachcontig[contig].append(index)
		
		
	allqueryloci = []
	
	"""
	refcontig = {}
	for contig, usecols in refcontig.items():

		maxgenesize = max([refgenes_sizes[row[3].strip()] for row in usecols]) 

		loci = querycontig_overlap(usecols, max(30000, maxgenesize))

		loci = [[contig]+x for x in loci]

		allqueryloci.extend(loci)
	"""
	
	for contig, indexes in eachcontig.items():
		
		subdata = table.iloc[indexes]
		
		maxgenesize = max(list(subdata[6]))
		mingenesize = min(list(subdata[6]))		
		
		usecols = subdata[[2,3,4,5,6]].values.tolist()

		loci = querycontig_overlap(usecols, 2*args.anchor, 100)

		loci = [[contig]+x for x in loci]
		
		allqueryloci.extend(loci)


	strand_arg = {"+":"", "-":"--reverse-complement"}
	total_loci = len(allqueryloci)
	output_regions = cl.defaultdict(list)

	with open(querypath, mode = 'r') as f:
		
		querypaths = {x.split()[0]:x.split()[1] for x in f.read().splitlines()}
	

	querypaths["CHM13"] = chm13
	querypaths["HG38"] = hg38
	
	strand_arg = {"+":"", "-":"--reverse-complement"}
	
	try:
		os.remove(outfile)
	except:
		pass
		
	try:
		os.remove(outfile+"_info.txt")
	except:
		pass
	
	allqueryloci = sorted(allqueryloci)
	
	allqueryloci_exd = []
	for i, locus in enumerate(allqueryloci):
				
		size = abs(locus[1] - locus[2])
			
		if size < 5*args.anchor:
			
			extendsize = 5*args.anchor - size
			
			locus[1] = max(locus[1] - extendsize//2, 0)
			
			locus[2] += extendsize//2
		
		allqueryloci_exd.append(locus)

	lastcontig = ""
	lastend = -10000
	for i, locus in enumerate(allqueryloci_exd):
		
		if lastcontig == locus[0] and locus[1] < lastend:
		
			midpoint = min( max( (locus[1]+lastend)//2 , allqueryloci[i-1][2] ) , allqueryloci[i][1] )

			allqueryloci_exd[i-1][2] = midpoint
			locus[1] = midpoint


		lastcontig = locus[0]
		lastend = locus[2]
		
			
	
	outfile = outpath
	for i, locus in enumerate(allqueryloci_exd):
		
		samplename = "_h".join(locus[0].split("#")[:2]) 
		
		queryfile = querypaths.get(samplename,hg38)
		
		if queryfile == hg38 and "NC_0609" in samplename:
			queryfile = querypaths["CHM13"]
			locus[0] = locus[0].split("#")[-1]
		samplename = samplename.split(":")[0]   
		
		
			
		title = ">{}_{}\t{}:{}-{}{}\t{}".format(samplename, i, locus[0].split(":")[0], locus[1], locus[2],locus[3],",".join(locus[-1]))
		
		if total_loci < 1000:
			cm = "echo \"{}\" >> {}".format(title, outfile)
			os.system(cm)
			
			cm = "samtools faidx {} {}:{}-{} {}  | grep -v '>' >> {}".format(queryfile,locus[0].split(":")[0], max(0,locus[1]), locus[2], strand_arg[locus[3]], outfile)
			os.system(cm)
			
			cm = "echo $\'\\n\' >> {}".format(outfile)
			os.system(cm)
		else:
			output_regions[queryfile].append(locus+[i])
				

	if len(output_regions):
		outfile_ = open(outfile, mode = "w")
	
		for query, infos in output_regions.items():
		
			export_fasta(query, infos, outfile_)
	
		outfile_.close()
	

	
		
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str,required=True)
	parser.add_argument("-q", "--query", help="query path", dest="query",type=str,required=True)
	parser.add_argument("-g", "--refgenes", help="refgene file path", dest="refgenes", type=str,required=True)
	parser.add_argument("-r", "--ref", help="refgene file path", dest="ref", type=str,required=True)
	parser.add_argument("-a", "--anchor", help="anchor size", dest="anchor", type=int,default=3000)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
