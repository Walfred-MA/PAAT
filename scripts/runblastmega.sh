#!/bin/sh

sizecut=$4
blastn -task megablast  -query $1 -db $2 -outfmt 17 -word_size 50 -num_threads $3 -evalue 1e-300 -dust yes -lcase_masking  -perc_identity 99 -qcov_hsp_perc 0.1  -max_target_seqs 100 -out - | awk -v cutoff=$sizecut '{for(i=1;i<=NF;i++){if($i ~ /AS:i:/){split($i,a,":"); if(a[3]>cutoff){print $0}}}}' | awk '{for(i=1;i<=NF;i++){if($i ~ /PI:f:/){split($i,a,":"); if(a[3]>95){print $0}}}}' 
