
#!/bin/bash

INFILE=$1
CHAINS=$2
CHROMSIZES=$3


/home/emmanuel/software/bigWigToBedGraph ${INFILE}.bw ${INFILE}.bedGraph #Strip infile name

CHRFORMAT=$(cut -f 1 ${INFILE}.bedGraph | head -n +1)

echo ${CHRFORMAT}

if (( $CHRFORMAT >= 1 && ${CHRFORMAT} <= 22 ))
then
	echo "ensembl chr format .... converting to UCSC"
	sed -E 's/^(([0-9]|[XY]|M)([0-9]|T|.))/chr\1/' ${INFILE}.bedGraph > ${INFILE}_UCSC.bedGraph

	/home/emmanuel/software/liftOver ${INFILE}_UCSC.bedGraph ${CHAINS} ${INFILE}_UCSC_hg38.bedGraph ${INFILE}.unmapped

	sort -k1,1 -k2,2n ${INFILE}_UCSC_hg38.bedGraph | bedtools merge -c 4 -o mean -d "-1" > ${INFILE}_UCSC_hg38_sorted_mergedMean.bedGraph

	/home/emmanuel/software/bedGraphToBigWig ${INFILE}_UCSC_hg38_sorted_mergedMean.bedGraph ${CHROMSIZES} ${INFILE}_UCSC_hg38_sorted_mergedMean.bw
else
	echo "UCSC chr format"

	/home/emmanuel/software/liftOver ${INFILE}.bedGraph ${CHAINS} ${INFILE}_hg38.bedGraph ${INFILE}.unmapped

	sort -k1,1 -k2,2n ${INFILE}_hg38.bedGraph | bedtools merge -c 4 -o mean -d "-1" > ${INFILE}_hg38_sorted_mergedMean.bedGraph

	/home/emmanuel/software/bedGraphToBigWig ${INFILE}_hg38_sorted_mergedMean.bedGraph ${CHROMSIZES} ${INFILE}_hg38_sorted_mergedMean.bw
fi

rm ${INFILE}*.bedGraph

#| bedtools merge -c 4 -o sum 