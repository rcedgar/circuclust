#!/bin/bash -e

if [ x$circuclust == x ] ; then
	if [ -s ../bin/circuclust ] ; then
		circuclust=../bin/circuclust
	else
		circuclust=`which circuclust`
	fi
fi

mkdir -p ../test_output
cd ../test_output

# test1.fa is an arbitrary nt sequence
# test.fa is a set of rotated and/or revcomp'd versions of test1
# 	they should all cluster together
#	when rotated and oritented they should all be identical to test1

$circuclust \
  -cluster ../test_data/test.fa \
  -id 0.9 \
  -tsvout test_cluster.tsv \
  -alnout test_cluster.aln

$circuclust \
  -align ../test_data/test.fa \
  -db ../test_data/test1.fa \
  -rotated test_rotated.fa \
  -tsvout test_align.tsv \
  -alnout test_align.aln
