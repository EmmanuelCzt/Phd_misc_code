#!/bin/bash

input=$1
mkdir alfred_parsed

zcat ${input} | grep "^ME" | cut -f 2- | sed 's/ /\t/g' > alfred_parsed/`basename ${input} .tsv.gz`_ME.tsv
zcat ${input} | grep "^MQ" | cut -f 2- | sed 's/ /\t/g'> alfred_parsed/`basename ${input} .tsv.gz`_MQ.tsv
zcat ${input} | grep "^CM" | cut -f 2- | sed 's/ /\t/g'> alfred_parsed/`basename ${input} .tsv.gz`_CM.tsv
zcat ${input} | grep "^IZ" | cut -f 2- | sed 's/ /\t/g'> alfred_parsed/`basename ${input} .tsv.gz`_IZ.tsv
zcat ${input} | grep "^TC" | cut -f 2- | sed 's/ /\t/g'> alfred_parsed/`basename ${input} .tsv.gz`_TC.tsv