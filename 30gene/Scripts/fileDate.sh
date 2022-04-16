#! /bin/bash

# This script generates a tsv file with the names of the vcf and its fileDates

for i in $(find . -name '*.vcf')
do
	j=$(grep '##fileDate' $i)
	k=${i%/*}
	l=$(du $i)
	printf "${k##*/}\t$(basename $i)\t${j/#*=}\t${l%\t}\n" %b
done

