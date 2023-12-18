#!/bin/bash

url=$1

wget --spider --force-html -r -l2 ${url} 2>&1 | grep '^--' | awk '{ print $3 }'| grep -v '\.\(css\|js\|png\|gif\|jpg\)$' > url.list

for i in `cat url.list`
do
  wget $i | tee -a ftg.log
done

# Scrap all urls of a webpage
# --spider check the availability of a webpage 
# --force-html treat the downladed page as html url
# -r recursive retrieving 
# -l2 depth of recursivity (2)