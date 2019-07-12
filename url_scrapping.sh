#!/bin/bash

url=$1
url_list=$2

wget --spider --force-html -r -l2 ${url} 2>&1 | grep '^--' | awk '{ print $3 }'| grep -v '\.\(css\|js\|png\|gif\|jpg\)$' > ${url_list}.list

# Scrap all urls of a webpage
