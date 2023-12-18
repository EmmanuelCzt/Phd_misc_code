#!/bin/bash

GENE_TAB=$1
PREFIX=$2

#All gene set
cut -f 1 ${GENE_TAB} | tail -n +2 > ${PREFIX}\_all.txt
#Lineage Acquired
awk '$7=="LNG" || $7=="LNG-PAR" {print $1}' ${GENE_TAB} > ${PREFIX}\_lineageAcq.txt
#Species Acquired
awk '$7=="INDE" {print $1}' ${GENE_TAB} > ${PREFIX}\_SpeciesSpe.txt

#Ampliconic only
awk '$5=="Ampliconic" {print $1}' ${GENE_TAB} > ${PREFIX}\_ampliconic.txt
#Ampliconic Lineage acquired
awk '$5=="Ampliconic" && ($7=="LNG" || $7=="LNG-PAR") {print $1}' ${GENE_TAB} > ${PREFIX}\_ampliconic_lineageAcq.txt
#Ampliconic Species acquired
awk '$5=="Ampliconic" && $7=="INDE" {print $1}' ${GENE_TAB} > ${PREFIX}\_ampliconic_SpeciesSpe.txt


#Non-ampliconic
awk '$5=="N.A." {print $1}' ${GENE_TAB} > ${PREFIX}\_NonAmp.txt
#Non-Ampliconic Lineage acquired
awk '$5=="N.A." && ($7=="LNG" || $7=="LNG-PAR") {print $1}' ${GENE_TAB} > ${PREFIX}\_NonAmp_lineageAcq.txt
#Non-Ampliconic Species acquired
awk '$5=="N.A." && $7=="INDE" {print $1}' ${GENE_TAB} > ${PREFIX}\_NonAmp_SpeciesSpe.txt

