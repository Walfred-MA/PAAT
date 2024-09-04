# Pangenome Allele Annotation Toolkits (PATs)
The toolkits and pipeline for Pangenome Allele annotation from pangenomes are organized into four components:

AlleleSearch: A pipeline to search and extract pangenome alleles.

Filtration: A pipeline to partition initial gene groups, extract exclusive k-mers for each gene group, and filter out sequences of no interest, such as tandems.

Compile: A pipeline to construct k-mer matrices, phylogenetic trees, and pan-gene graphs, then compile them into database files used by Ctyper.

Annotation: A pipeline to annotate PAs, which aids in interpreting Ctyper's genotyping results.

We provide limited maintenance for the code in this repository. While we do not guarantee compatibility with all environments, we offer assistance upon request to help run the pipeline for annotating new pangenome assemblies. For other purposes or individual use of specific scripts, maintenance and assistance may be provided on a case-by-case basis.

If you are interested in using any ideas or methods presented here in your own projects, please contact us at wangfeim@usc.edu and mchaisso@usc.edu. We will respond to requests regarding the annotation of new pangenome assemblies.
