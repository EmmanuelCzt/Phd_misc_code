#!/bin/bash

OUTFILE=$1

for i in `ls *.gz`
do
	md5sum $i >> ${OUTFILE}.txt
done