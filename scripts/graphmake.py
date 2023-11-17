#!/usr/bin/env python3

import argparse
import os
import subprocess
import datetime
from graphfilesplit import filesplit
from graphinserts import insertdistact
from graphcigar import graphcigar
from graphcleanrepeats import cleanrepeats
from graphtolinear import graphtolinear

def runblastn(ref, query, folder, output,nthreads):
	
	refname = ref.split("/")[-1]
	
	dbpath = folder+refname+"_db"
	
	cmd = "blastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity 95 -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( query,dbpath , output,nthreads)
	
	os.system(cmd)
	
def alignnewseq(folder,graphfile, thefile,nthreads):
	
	refname = graphfile.split("/")[-1]
	
	dbpath = folder+refname+"_db"
	
	output = thefile + "_align.out"
	
	cmd = "blastn -task megablast  -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30  -perc_identity 95 -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( thefile,dbpath , output, nthreads)
	
	os.system(cmd)
	
	
def addnewseq(folder, graphfile, addfile,nthreads):
	
	#cmd = "cat {} | grep \"^S\" | awk '{{print \">\"$2\"\\n\"$3}}' > {}".format( graphfile ,  graphfasta )
	#os.system(cmd)
	
	filealign = addfile + "_blastn.out"
	fileinserts = addfile + "_inserts.fasta"
	
	try:
		os.remove(fileinserts)
	except:
		pass
		
		
	#cmd = "minimap2 -c -x map-hifi --dual=yes {} {} > {}".format(graphfile,addfile,filealign)
	#os.system(cmd)
		
	runblastn(graphfile, addfile, folder, filealign,nthreads)

	insertfiles = insertdistact( filealign, addfile, fileinserts, 1)
	
	if len(insertfiles):
		os.system("cat {} >> {}".format(fileinserts, graphfile))
		
	return insertfiles	

def runcleanrepeats(folder, graphfile, newgraphfile,nthreads):

	#try:
		#os.remove(newgraphfile)
	#except:
		#pass

	cleanrepeats(graphfile, newgraphfile, nthreads)

	#filesplit(graphfile, folder, 0)

	#files = [folder+thefile for thefile in  os.listdir(folder) if thefile[-3:]==".fa" ]
	
	#for afile in files:

		#outfile = afile+"_norep.fasta"

		#cleanrepeats(afile, outfile, nthreads)

		#os.system("cat {} >> {}".format( outfile , newgraphfile ))	
	
	return newgraphfile

def editname(thefile, graphfile):

	size = 0
	with open(thefile, mode = 'r') as r:
		for line in r:
			if len(line) and line[0] != '>':
				size += len(line.strip())			


	with open(thefile, mode = 'r') as r:
		with open(graphfile, mode = 'w') as w:
			for line in r:
				if len(line) and line[0] == '>':
					
					header = line.strip().split()
					header[0] += "_{}_{}".format(0, size)

					if len(header[0]) > 50 -2:
						header[0] = ">"+"_".join(header[0].split("_")[2:])
					w.write("\t".join(header)+"\n")
				else:
					w.write(line)


def creategraph(folder, graphfile,nthreads ):
	
	
	allfiles = [x[1] for x in sorted([( 10e19 + os.path.getsize(folder+"/"+file) if "NC_" in file else 10e20 + os.path.getsize(folder+"/"+file) if  ("_chr" in file and 'v' not in file ) else  10e18 + os.path.getsize(folder+"/"+file) if  ( 'v' in file )  else os.path.getsize(folder+"/"+file) , folder+"/"+file ) for file in os.listdir(folder) if file[-3:]==".fa"], reverse = 1)]
	
	if len(allfiles) == 0:
		return

	editname(allfiles[0],graphfile)		
	
	graphfile_filename = graphfile.split("/")[-1]
	
	dbpath = folder+graphfile_filename+"_db"
	
	dbcmd = "makeblastdb -in {}  -dbtype nucl -parse_seqids -out {}".format(graphfile, dbpath)
	
	os.system(dbcmd)
	
	for addfile in allfiles[1:]:
		
		insertfiles = addnewseq(folder, graphfile, addfile,nthreads)
		
		if len(insertfiles):
			os.system(dbcmd)
			
def aligngraph(folder, graphfile, graphalign,nthreads):
	
	try:
		os.remove(graphalign)
	except:
		pass	


	graphfile_filename = graphfile.split("/")[-1]

	dbpath = folder+graphfile_filename+"_db"

	dbcmd = "makeblastdb -in {}  -dbtype nucl -parse_seqids -out {}".format(graphfile, dbpath)

	os.system(dbcmd)
	
	#allfiles = [x[1] for x in sorted([( 10e19 + os.path.getsize(folder+"/"+file) if "NC_" in file else 10e20 + os.path.getsize(folder+"/"+file) if  ("_chr" in file and 'v' not in file ) else os.path.getsize(folder+"/"+file) , folder+"/"+file ) for file in os.listdir(folder) if file[-3:]==".fa"], reverse = 1)]
	
	allfiles = [x[1] for x in sorted([( 10e19 + os.path.getsize(folder+"/"+file) if "NC_" in file else 10e20 + os.path.getsize(folder+"/"+file) if  ("_chr" in file and 'v' not in file ) else  10e18 + os.path.getsize(folder+"/"+file) if  ( 'v' in file )  else os.path.getsize(folder+"/"+file) , folder+"/"+file ) for file in os.listdir(folder) if file[-3:]==".fa"], reverse = 1)]

	for thefile in allfiles:
		
		output = thefile + "_align.out"
		
		cmd = "blastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30 -perc_identity 95 -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( thefile,dbpath , output,nthreads)
		
		os.system(cmd)
		
		queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, thefile , output)
		
		with open(graphalign, mode = 'a') as w:
			w.write("{}\t{}\t{}\t{}\t{}\n".format(queryname, fullpath, fullcigar, pathranges,qranges))
					
					
					
def main(args):
	
	folder = args.input +  datetime.datetime.now().strftime("%y%m%d_%H%M%S/")
	#folder = "./tempx/"

	folder = folder.replace("%","_")
	folder2 = folder + "/graphtemp/"

	try:
		os.system("rm -rf {}".format(folder))
		os.mkdir(folder)
	except:
		pass

	try:
		os.system("rm -rf {}".format(folder2))
		os.mkdir(folder2)
	except:
		pass
		
	inputfile_name = args.input.split("/")[-1]

	graphfile_raw = folder + inputfile_name+"_graphwithrep.FA"

	graphfile = args.input +"_graph.FA"
	graphalign = args.input + "_allgraphalign.out"
	graphlinear = args.input + "_lineargraph.gaf"
	
	filesplit(args.input, folder, 0)
	
	if args.create:
		creategraph(folder, graphfile_raw,args.thread)

	if args.fine:

		runcleanrepeats(folder2, graphfile_raw, graphfile,args.thread)

	if args.align:

		aligngraph(folder, graphfile, graphalign,args.thread)

	if args.linear:
				
		graphtolinear(graphalign, args.input, graphfile, graphlinear, args.thread, folder+"stretcherouts/")	
	
	try:
		os.system("rm -rf {}".format(folder))
		pass
	except:
		pass	



def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-c", "--create", help="if generate fine graph", dest="create", type=int,default = 1)
	parser.add_argument("-f", "--fine", help="if generate fine graph", dest="fine", type=int,default = 1)
	parser.add_argument("-a", "--align", help="if align to graph", dest="align", type=int,default = 1)
	parser.add_argument("-l", "--linear", help="if align to graph", dest="linear", type=int,default = 1)
	parser.add_argument("-t", "--thread", help="if align to graph", dest="thread", type=int,default = 1)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
