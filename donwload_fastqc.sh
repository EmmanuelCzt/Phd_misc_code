#!/bin/bash

export PATH=$PATH:/home/emmanuel/software/fastqc_v0.11.8/FastQC

url=$1

wget --spider --force-html -r -l2 ${url} 2>&1 | grep '^--' | awk '{ print $3 }'| grep -v '\.\(css\|js\|png\|gif\|jpg\)$' > url.list

for i in `cat url.list`
do
  wget $i |& tee -a ftp_fastqc.log
done

mkdir fastqc

for i in `ls`
do
  fastqc -o fastqc/ $i |& tee -a ftp_fastqc.log
done
# Scrap all urls of a webpage
