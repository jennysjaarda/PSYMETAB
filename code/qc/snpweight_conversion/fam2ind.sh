#!/bin/bash
inp=$1
cat $inp | awk '{ print $2" "$5" "$1}' | sed -e 's% 2 % F %' -e 's% 1 % M %' | tr ' ' '\t' > $(basename $inp .fam).ind
