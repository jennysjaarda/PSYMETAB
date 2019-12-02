#!/bin/bash
inp=$1
cat $inp | awk '{ print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6}'| awk -F'\t' -vOFS='\t' '{ gsub("0", "X", $5) ; gsub("0", "X", $6) ; gsub("I", "X", $6);  print }' > $(basename $inp .bim).snp
