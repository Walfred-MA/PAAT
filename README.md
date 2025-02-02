# Pangenome Allele Annotation Toolkits (PATs)
The toolkits and pipeline for Pangenome Allele annotation from pangenomes are organized into four components:

AlleleSearch: A pipeline to search and extract pangenome alleles.

Filtration: A pipeline to partition initial gene groups, extract exclusive k-mers for each gene group, and filter out sequences of no interest, such as tandems.

Compile: A pipeline to construct k-mer matrices and phylogenetic trees then compile them into database files used by Ctyper.

Annotation: A pipeline to annotate all pangenome-alleles to be used to interpret the results
How to run this pipeline:

1. First, you will need find the interested sequences from pangenome assemblies. This can be achived by the snakefile in AlleleSearch folder. You will need to provide the path of pangenome assemblies and the prefixes or names of genes or gene families you interested in. Then run the snakemake file and it will generate a list of a .fasta file for each of genes or gene families you interested in. 
2. Next, you will need to distract informative k-mers from each gene or gene families of interest. This can be achieved by the snakefile in Kmers folder. You will need to compile a custom c++ script in src/kmer_haplotyping and put into the script folder. After running this pipeline, in the partitions/ folder, it will output one or multiple matrices (if this gene family can be partited into several smaller groups), with their kmer list included.
3. Then, you will need to compile the sequences and k-mers into databases used for genotyping, which we called matrices files. This can be achived by the snakefile in Compile folder. For each matrix from the previous step, it will output a _matrix.txt file in the same folder, which is the database file you need.
4. Final, you will need to concatenate all  _matrix.txt file together to a single database file and index it for ctyper to run. 
5. (optional) Additional, to annotate all PAs in the database, you may run snakefile at Annotate/ folder. It will generate a list of annotation results, concatenate all *_annotatesummaryfix.txt to generate annotation summary file to be used.


If you are interested in using any ideas or methods presented here in your own projects, please contact us at wangfeim@usc.edu and mchaisso@usc.edu. We will respond to requests regarding the annotation of new pangenome assemblies.
