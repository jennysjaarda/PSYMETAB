#!/bin/bash
inp=${1%%.bed}
plink --bfile $inp --recode A-transpose --out tmp
cat tmp.traw | cut -f 1-6 --complement | tail -n +2 | sed 's%\t%%g'| sed 's/NA/9/g' > `basename $inp.geno`
