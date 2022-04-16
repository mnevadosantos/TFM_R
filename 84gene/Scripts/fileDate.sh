#! /bin/bash

# Este script genera un fichero tsv con los nombres de los vcf y sus fileDates

for i in $(find . -name '*.vcf')
do
	j=$(grep '##fileDate' $i)
	k=${i%/*}
	l=$(du $i)
	printf "${k##*/}\t$(basename $i)\t${j/#*=}\t${l%\t}\n" %b
done

