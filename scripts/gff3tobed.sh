#!/bin/bash

input=$1
type=$2
output=$3

cat $input  | awk -v var=$type '{if ($3 == var) print $1,$4,$5,$7,$9}' | awk -v var=$type"_id=[^;]+" -v var2="gene_id=[^;]+"  'BEGIN{OFS="\t"}; match($5,var) { match1 = match($5,var); start1 = RSTART; len1 = RLENGTH ; match2 = match($5,var2);start2 = RSTART; len2 = RLENGTH  ; print $1,$2,$3,$4, substr($5,start1+8,len1-8),substr($5,start2+8,len2-8) } '  > $output
