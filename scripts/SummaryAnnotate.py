#!/usr/bin/env python3


import subprocess
import collections as cl
import argparse
import os
from StretcherReader import findglobalalign
from maskstretcher import maskfile
from CompareTrans import Comparegenes
import re

class GenodeDB:
	
	def __init__(self, path):
		
		self.pattern = re.compile(r'(?:gene_name=)([^;]+)')
		self.path = path
		self.gene = cl.defaultdict(list)
		
		
		with open(self.path, mode = 'r') as r:
			
			for line in r:
				if len(line) == 0 or line[0]=="#":
					continue
				
				row = line.split("\t")
				
				
				if row[2] != "gene":
					continue
				
				gene_name = re.findall( self.pattern ,row[8])
				
				if len(gene_name) == 1:
					
					geneinfo = [int(row[3]), int(row[4])] + gene_name
					
					self.gene[row[0]].append(geneinfo)
					
	def OverlapChr(self, chr , queries):
		
		results = [[] for x in queries]
		
		refs = self.gene[chr]
		
		length = len(refs) + len(queries)
		lgenes = len(refs)
		
		coordinates = [x for y in refs + queries for x in y[:2]] 
		names = [x[2] for x in refs]
		
		coordinates_sortindex = sorted( range(length * 2), key = lambda x: coordinates[x])
		current_refs = []
		current_queries = []
		for index in coordinates_sortindex:
			
			index2 = index // 2
			if index % 2 == 0:
				
				if index2 < lgenes:
					
					current_refs.append(index2)
					
					for query in current_queries:
						
						results[query].append(index2)
						
				else:
					
					current_queries.append(index2 - lgenes)
					results[index2 - lgenes].extend(current_refs)
					
			else:
				
				if index2 < lgenes:
					
					current_refs.remove(index2)
					
				else:
					
					current_queries.remove(index2 - lgenes)
					
		outputs = cl.defaultdict(list)
		for i,result in enumerate(results):
			
			queryname = queries[i][2]
			
			outputs[queryname].extend([refs[refindex][2] for refindex in result])
			
		return outputs
	
	def Overlap(self, locations):
		
		results = cl.defaultdict(list)
		for chr, locations_bychr in locations.items():
			
			result_chr = self.OverlapChr(chr, locations_bychr)
			
			for queryname,genes in result_chr.items():
				
				results[queryname].extend(genes)
				
		return results
	
	
	def getgenename(self, location, pref):
		
		if len(location):
			
			locations_dic = {location.split(":")[0]:[list(map(int,location.split(":")[1][:-1].split("-")))+[""]]}
			
			type_locs = list(set(sum ( [loc for name, loc in self.Overlap(locations_dic).items()],[])))
			type_locs = [x for x in list(type_locs) if x.startswith(pref)]
			
		shortnames = set(type_locs)
		newtype_locs = []
		for gene in type_locs:
			if '-' in gene and gene.split('-')[0] in shortnames:
				continue
			newtype_locs.append(gene)
			
		type_locs=newtype_locs
		type_locs = ";".join(type_locs)
		
		if len(type_locs) == 0:
			type_locs = "NA"
			
		return type_locs
	
	
def analyzesv(svs, exon_q, exon_r):
	
	exon_q = [list(map(int,x.split("_"))) for x in exon_q.split(";") if len(x)]
	exon_r = [list(map(int,x.split("_"))) for x in exon_r.split(";") if len(x)]
	
	
	for sv in svs:
		
		sv = [x for x in sv]
		sv[1] = int(sv[1])
		sv[2] = int(sv[2])
		
		if "Insert" in sv[0] or "Duplication" in sv[0]:
			
			for exon in exon_q:
				
				if sv[2] - sv[1] + exon[1] - exon[0] - (max(sv[2], exon[1]) - min(sv[1], exon[0]) ) >0:
					
					return "SVonExon"
		else:
			for exon in exon_r:
				
				if sv[2] - sv[1] + exon[1] - exon[0] - (max(sv[2], exon[1]) - min(sv[1], exon[0]) ) >0:
					
					return "SVonExon"
				
	return ""

def  analyzecv(refpath, gffpath, geneliftover, transliftover):
	
	
	if geneliftover == "NA" or len(transliftover) == 0:
		return ""
	
	geneliftover = geneliftover.split(";")
	
	for gene in geneliftover:
		
		for (gene2, tran) in transliftover:
			
			if gene2 == gene or Comparegenes(refpath, gffpath, gene, tran):
				
				return ""
			
			
	return "Convertion"



def mergeloci(loci):
	
	
	allcontigs = cl.defaultdict(list)
	for name, (locus, simi) in loci.items():
		
		if locus != "NA":
			contig,therange = locus.split(":")
			therange = [int(x) for x in therange[:-1].split('-')]
			allcontigs[contig].append([name, therange])
			
	allclusters = []
	for contig, coordis in allcontigs.items():
		
		allnames = [x[0] for x in coordis]
		allcoordis = sum([x[1] for x in coordis], [])
		allcoordis_sortindex = sorted(range(len(allcoordis)), key = lambda x: allcoordis[x])
		
		clusters = [[]]
		coverage = 0
		for sortindex in allcoordis_sortindex:
			
			if sortindex % 2 == 0:
				coverage+=1
				clusters[-1].append(allnames[sortindex//2])
			else:
				coverage-=1
				if coverage == 0:
					clusters.append([])
		allclusters.extend(clusters[:-1])
		
		
	nametoloci = dict()
	for index, cluster in enumerate(allclusters):
		
		allhaplos = cl.defaultdict(list)
		for name in cluster:
			if "_h" not in name:
				haplo = "CHM13" if "NC_0609" in name else "HG38" if "v" not in "_".join(name.split("_")[2:]) and name.startswith("HLA") == False else "_".join(name.split("_")[2:4])
			else:
				haplo = "_".join(name.split("_")[2:4])
				
			allhaplos[haplo].append((loci[name][1],loci[name][0],name))
			
			
			for haplo, allnames in allhaplos.items():
				
				allnames = sorted(allnames)
				
				for j,name in enumerate(allnames):
					
					nametoloci [name[-1]] = str(index) + "_" + str(j)
					
					
	return nametoloci


def cleananchors(cigar_text, anchorsize = 100, ifreverse = 0):
	
	cigars = re.findall(r'\d+[A-Za-z]+',cigar_text)
	
	if ifreverse:
		cigars  = cigars[::-1]
		
	lastcigar= ''
	leftsize = anchorsize 
	index = 0
	
	cigars_new = []
	for index,cigar in enumerate(cigars):
		
		size_str = re.findall(r'\d+',cigar)[0]
		size_str_len = len(size_str)
		
		csize = int(size_str)
		ctype = cigar[size_str_len]
		theseq = cigar[(size_str_len+1):]
		
		if ctype in ['M','X','=','D']:
			
			if csize > leftsize:
				
				newsize = csize - leftsize
				qseqstart = len(theseq) - csize 
				if ifreverse:
					newseq = theseq[qseqstart:(qseqstart+newsize)]
					if ctype == 'X':
						newseq = theseq[:newsize] + newseq
				else:
					
					newseq = theseq[( len(theseq)- newsize):]
					if ctype == 'X':
						newseq = theseq[(qseqstart-newsize):qseqstart] + newseq
						
				lastcigar = str(newsize)+ctype+newseq
				
			leftsize = max(0, leftsize- csize)
			
			if leftsize == 0:
				break
		else:
			cigars_new.append(cigar)
			
	if ifreverse:
		cigars = "".join(cigars[(index+1):][::-1]) +  lastcigar + "".join(cigars_new)
	else:   
		cigars = "".join(cigars_new)+lastcigar+"".join(cigars[(index+1):])
		
		
	return cigars


def main(args):
	
	genedb= GenodeDB(args.database)
	refdb = args.ref
	nametogenenames= dict()
	nameandloci=dict()
	nametocopies = cl.defaultdict(str)
	nametograph = dict()
	nametoanno = dict()
	nametotype = dict()
	nametocigar = cl.defaultdict(str)
	nametosvs = cl.defaultdict(list)
	nametoexons =  cl.defaultdict(str)
	nametoliftover = cl.defaultdict(list)
	
	assignref = cl.defaultdict(lambda : ["NA","",10000.0])
	reftypes = set()
	
	usekmernum = cl.defaultdict(lambda : "3000")
	if len(args.kmerinfo):
		with open(args.kmerinfo, mode = 'r') as f:
			for line in f:
				line = line.strip().split()
				usekmernum[line[0]] = line[2]
				
				
	typeindex = -1
	with open(args.typeannotate, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				
				typeindex += 1
				line = line.split()
				
				typename = "_".join(line[-1].split('_')[:2])+ "_"+str(typeindex)
				allalleles = line[-1].split(',')
				
				for allele in allalleles:
					nametotype[allele] = typename
					
	with open(args.refmatch, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				
				line = line.split()
				
				if line[-1] == "Ref":
					line[3] = line[2]
					nametogenenames[line[0]] = genedb.getgenename(line[2], line[0].split("_")[0][:3])
					reftypes.add(nametotype[line[0]])
					line[1] = line[0]
					
				nametoanno[line[0]] = [line[-1], line[4], line[2], line[1], line[3], line[5]]
				nameandloci[line[0]] = (line[3] , float(line[5]))
				nametoliftover[line[0]] = line
				
				if float(line[5]) < min(0.1, assignref[line[0]][2]):
					assignref[line[0]] = (line[1] , line[3],float(line[5]))
					
					
	with open(args.liftover1, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				
				line = line.split()
				nametoanno[line[0]] = [line[-1], line[4], line[2], line[1], line[3], line[5]]
				nameandloci[line[0]] = (line[3], float(line[5]))
				nametoliftover[line[0]] = line
				
				if assignref[line[0]][2] == 10000.0:
					assignref[line[0]] = (line[1],line[3],0.0)
					
					
	additional_liftover = []
	with open(args.liftover2, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				
				line = line.split()
				nametoliftover[line[0]] = line
				additional_liftover.append(line[0])
				
	for name in additional_liftover:
		
		info = nametoliftover[name]
		liftto = info[1]
		nextlift = info[1]
		if nextlift != 'NA':
			liftto = nextlift
			nextlift = nametoliftover[nextlift][1]
		nametoliftover[name][1] = nextlift
		
		nametoanno[name][4] = nametoliftover[name][3]
		
		if nextlift != 'NA':
			nameandloci[name] = (nametoliftover[nextlift][3],float( nametoliftover[name][5]))
			
			
	allloci = mergeloci(nameandloci)
	with open(args.graphfile, mode = 'r') as f:
		for line in f:
			if len(line.strip()) and line[0] == 'L':
				line = line.split()
				nametograph[line[1]] = "\t".join([line[3], line[5]])
				
				
				
	nametoproteincoding = cl.defaultdict(list)
	with open(args.genecount, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				
				line = line.split()
				if  float(line[6]) > 90.0 and line[4] == "protein_coding":
					
					nametocopies[line[0]] += "{}:{}:{};".format(line[1].replace(":",""), line[3].split('-')[0], round(float(line[6]), 2) )
					
					nametoproteincoding[line[0]].append([ line[3].split('-')[0], line[1].replace(":","")]) 
					
					
					
	with open(args.genecount, mode = 'r') as f:
		
		for line in f:
			if len(line.strip()):
				line = line.split()
				if len(nametocopies[line[0]]) == 0:
					
					nametocopies[line[0]] += "{}:{}:{};".format( line[1].replace(":",""),line[3].split('-')[0], round(float(line[6]), 2) )
					
	with open(args.exon, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				line = line.split('\t')
				nametoexons[line[0]] = line[1].strip()
				
	with open(args.svfile, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				line = line.strip().split('\t')
				if len(line) > 1:
					
					if "Insert" in line[2] or "Duplication" in line[2]:
						sv_text = [line[2],line[3], line[4], line[-1]]
					else:
						sv_text = [line[2],line[5], line[6], line[-1]]
						
					nametosvs[line[0]].append(sv_text)
					
	with open(args.svfile2, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()):
				line = line.strip().split('\t')
				if len(line) > 1:
					if "Insert" in line[2] or "Duplication" in line[2]:
						sv_text = [line[2],line[3], line[4], line[-1]]
					else:
						sv_text = [line[2],line[5], line[6], line[-1]]
						
					nametosvs[line[0]].append(sv_text)
					
	name_ordered = []
	with open(args.input, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()) and line[0] == '>':
				name_ordered.append(line.split()[0][1:])
				
	folder = args.folder + "/"
	for name in name_ordered:
		
		thefile = folder + name + "_o.txt"
		
		reversefile = thefile + "_reverse"
		if os.path.isfile(reversefile):
			thefile = reversefile
			
		cigar = ""
		if os.path.isfile(thefile ):
			maskfile(thefile)
			score, cigar = findglobalalign(thefile)
			#cigar = cigar.translate(str.maketrans('TCGNatcgn','AAAAAAAAA')).replace('A','') 
			
			cigarrsize =  sum([ int(x[:-1]) for x in re.findall(r'\d+[M=DX]',cigar)])
			location = nameandloci[name][0].split(":")[1][:-1].split('-')
			
			location = [int(location[0]), int(location[1])]
			leftanchor = min(100, location[0])
			rightanchor = min(100,  location[1] - (leftanchor + cigarrsize) ) 
			
			cigar = cleananchors(cleananchors(cigar, leftanchor),rightanchor,1)
			
		nametocigar[name] = cigar
		
	"""
	typefix = cl.defaultdict(int)
	typeallcounts = cl.defaultdict(set)
	for name in name_ordered:

		thetype = nametotype[name]

		#genecounts = ";".join([x.split(":")[1] for x in  nametocopies[name].split(";")[:-1]])

		genecounts = ";".join(nametoproteincoding[name])

		typeallcounts[thetype].add(genecounts)
		typefix[(thetype, genecounts)] += 1

	pair_sort = sorted(list(typefix.keys()), key = lambda x: typefix[x], reverse = 1)

	typecounts = cl.defaultdict(int)
	for (thetype, genecounts) in pair_sort:
		typefix[(thetype, genecounts)] = typecounts[thetype]
		typecounts[thetype] += 1
	"""
		
	with open(args.output, mode = 'w') as f:
		
		for name in name_ordered:
			
			
			thetype = nametotype[name]
			#genecounts = ";".join([x.split(":")[1] for x in  nametocopies[name].split(";")[:-1]])
			"""
			genecounts = ";".join(nametoproteincoding[name])
			if typefix[(thetype, genecounts)] :

				thetype += "_sub"+str(typefix[(thetype, genecounts)])   
			"""
			
			refgenename = nametogenenames.get(nametoanno[name][3], genedb.getgenename(  nametoanno[name][4] , name.split("_")[0][:3]) if ":" in nametoanno[name][4] else "NA"  )
			note1 = analyzesv(nametosvs[name], nametoexons[name], nametoexons[nametoanno[name][3]])
			note2 = analyzecv(args.ref, args.gff,refgenename, nametoproteincoding[name])
			
			note = ";".join([ x for x in [note1, note2] if len(x)])
			if len(note) == 0:
				note = "NA"
				
			svs = nametosvs[name]   
			if len(svs) == 0:
				sv_text = "NA"
			else:
				sv_text =  ";".join(["_".join(sv) for sv in svs])
				
			if nametoanno[name][0] == 'Pri':
				nametoanno[name][0] = 'Ref' if thetype in reftypes  else 'Alt'
				
			anno = nametoanno[name]
			
			if len(nametocopies[name]) == 0:
				nametocopies[name] = "NA"
				
				
			if anno[1] == "1" or nametocopies[name] != "NA":
				groupname = name.split('_')[0]
				
				
				if not [x for x in nametocopies[name].split(";") if ":" in x and x.split(":")[1].startswith(groupname)]:
					
					thetag = "Decoy"
					
				else:
					thetag = "Exon"
					
			else:
				thetag = "Intron"
				
				
				
				
				
			anno = [anno[0],thetag,anno[2],anno[3],anno[5], anno[4]]
			
			exonloc = nametoexons[name]
			if len(exonloc) == 0:
				exonloc = "NA"
				
				
				
			output = [name,thetype, refgenename ,"Para_" +str(allloci.get(name, thetype)), nametocopies[name], "\t".join(anno),  nametoliftover[name][1]+":"+ nametoliftover[name][3], usekmernum["_".join(name.split("_")[2:])] , note ,sv_text, exonloc,nametograph[name],  assignref[name][0]+":"+assignref[name][1]+":"+nametocigar[name]]
			
			f.write("\t".join(output)+"\n")
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-t", "--type", help="path to input data file",dest="typeannotate", type=str, required=True)
	parser.add_argument("-r", "--refmatch", help="path to input data file",dest="refmatch", type=str, required=True)
	parser.add_argument("-g", "--graphfile", help="path to input data file",dest="graphfile", type=str, required=True)
	parser.add_argument("-c", "--genecount", help="path to input data file",dest="genecount", type=str, required=True)
	parser.add_argument("-l1", "--liftover1", help="path to input data file",dest="liftover1", type=str, required=True)
	parser.add_argument("-l2", "--liftover2", help="path to input data file",dest="liftover2", type=str, required=True)
	parser.add_argument("-s", "--svfile", help="path to input data file",dest="svfile", type=str, required=True)
	parser.add_argument("-f", "--folder", help="path to input data file",dest="folder", type=str, required=True)
	parser.add_argument("-s2", "--svfile2", help="path to input data file",dest="svfile2", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
	parser.add_argument("-d", "--database", help="path to input data file",dest="database", type=str, required=True)
	parser.add_argument("-k", "--kmerinfo", help="path to input data file",dest="kmerinfo", type=str, required=True)
	parser.add_argument("-e", "--exon", help="path to input data file",dest="exon", type=str, required=True)
	parser.add_argument("-ref", "--ref", help="path to input data file",dest="ref", type=str, required=True)
	parser.add_argument("-gff", "--gff", help="path to input data file",dest="gff", type=str, required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
