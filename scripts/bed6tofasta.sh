#!/bin/bash

theref=$1
input=$2
anchor=$3
output=$4

cat $input | while read chrom start end strd name gene; 
do 
	start2=$(( $start > $anchor  ?  ( $start - $anchor ) : 0  )) ;
	start=$(( $start - $anchor )) ; end=$(( $end + $anchor ))  ;
	echo ">$name $chrom:$start-$end$strd $gene"  ;
	if [[ $strd == "+" ]] ;
	then		
		samtools faidx $theref "$chrom:$start2-$end" | grep -v "^>" ;
	else
		samtools faidx  $theref "$chrom:$start2-$end" -i  | grep -v "^>" ;
	fi ;
done > $output ; 
