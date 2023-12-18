configfile: "config.yaml"

#SAMPLES = ["SRR5048133", "SRR5048134"]

rule all:
    input:
        expand("mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}_SNP_q10.vcf", sample=config["samples"])


#rule fastqc:
#    input:
#        "fastq/{sample}.R1.fastq.gz"
#    output:
#        html="qc/fastqc/{sample}.html",
#        zip="qc/fastqc/{sample}.zip"
#    params: ""
#    log:
#        "logs/fastqc/{sample}.log"
#    wrapper:
#        "0.30.0/bio/fastqc"



#rule hisat2:
#    input:
#        fastq1="fastq/{sample}_1.fastq.gz",
#        fastq2="fastq/{sample}_2.fastq.gz",
#    output:
#        temp("star_{sample}_{ref_genome_name}.sam")
#    params:
#        hisat_index="~/Documents/HISAT_index/{ref_genome_name}/{ref_genome_name}"
#    threads: 32
#    log:
#        "logs/{sample}_{ref_genome_name}_hisat2.log"
#    shell:
#        "(hisat2 -p {threads} --dta \
#        -x {params.hisat_index} \
#        -1 {input.fastq1} \
#        -2 {input.fastq2} \
#        -S {output}) 2> {log}"


#rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = ["fastq/{sample}.R1.fastq.gz"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = ["fastq/{sample}.R2.fastq.gz"] #optional
#    output:
        # see STAR manual for additional output files
        "star/pe/{sample}/Aligned.out.bam"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index="/media/louis/dd2To/Genomes/{ref_genome_name}_fasta/STAR_index/",
        # optional parameters
        extra=""
    threads: 8
    wrapper:
        "0.31.0/bio/star/align"


rule samtools_sort:
    input:
        "~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/{sample}.cram"
    output:
        temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted.bam")
    threads: 12
    log:
        "logs/{sample}_{ref_genome_name}_samtools_sort.log"
    shell:
        "(samtools sort -@ {threads} \
            -o {output} \
            {input}) 2> {log}"


rule uniq_reads:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted.bam"
    output:
        temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq.bam")
    params: qual=10
    log:
        "logs/{sample}_{ref_genome_name}_uniq.log"
    shell:
        "(samtools view -q {params.qual} -b {input.bam} > {output}) \
            2> {log}"

rule samtools_index:
    input: 
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq.bam"
    output: 
        temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq.bam.bai")
    params:
        "" # optional params string
    shell:
        "(samtools index {input.bam} {output})"

#rule stringtie:
    input:
        bam="bam/{sample}_{ref_genome_name}_uniq.bam",
        bai="bam/{sample}_{ref_genome_name}_uniq.bam.bai"
#    output: protected("stringtie/{sample}_{ref_genome_name}.gtf")
    threads: 16
    params: "{sample}"
    log: "logs/{sample}_{ref_genome_name}_stringtie.log"
    shell:
        "(stringtie -p {threads} -o {output} \
         -l {params} {input.bam}) 2> {log}"

#rule bigwig_forward:
    input:
        bam="bam/{sample}_{ref_genome_name}_uniq.bam"
#    output: protected("bigwig_STAR/{sample}_{ref_genome_name}_BPM_all_forward.bw")
    params:
        filterRNAstrand="forward",
        normalizeUsing="BPM",
        binSize=20,
        smoothLength=40
    log: "logs/{sample}_{ref_genome_name}_BPM_all_forward.bw"
    shell:
        "(bamCoverage -b {input.bam} --filterRNAstrand {params.filterRNAstrand} \
        --normalizeUsing {params.normalizeUsing}  --binSize {params.binSize} \
        --smoothLength {params.smoothLength} -o {output}) 2> {log}"

#rule bigwig_reverse:
    input:
        bam="bam/{sample}_{ref_genome_name}_uniq.bam"
    output: protected("bigwig_STAR/{sample}_{ref_genome_name}_BPM_all_reverse.bw")
    params:
        filterRNAstrand="reverse",
        normalizeUsing="BPM",
        binSize=20,
        smoothLength=40
    log: "logs/{sample}_{ref_genome_name}_BPM_all_reverse.bw"
    shell:
        "(bamCoverage -b {input.bam} --filterRNAstrand {params.filterRNAstrand} \
        --normalizeUsing {params.normalizeUsing}  --binSize {params.binSize} \
        --smoothLength {params.smoothLength} -o {output}) 2> {log}"

#rule htseq:
    input:
        bam="bam/{sample}_{ref_genome_name}_uniq.bam",
        string="stringtie/{sample}_{ref_genome_name}.gtf",
        bai="bam/{sample}_{ref_genome_name}_uniq.bam.bai",
        #bw_forward="bigwig_STAR/{sample}_{ref_genome_name}_BPM_all_forward.bw",
        #bw_reverse="bigwig_STAR/{sample}_{ref_genome_name}_BPM_all_reverse.bw"
    output:
#        protected("htseq_curated_rmsk_class_id_chrX_whXACTT113/{sample}_{ref_genome_name}.txt")
    params:
        r="name",
        format="bam",
        stranded="reverse",
        type="exon",
        gtf="/media/louis/2dd2To/annotation/annotation/GENCODE_v29/gencode.v29.annotation_XACT_T113_wh_officialXACT.gtf",
        idattr="class_id"
    log: "logs/{sample}_{ref_genome_name}_repeatmasker_htseq.log"
    shell:
        "(htseq-count -r {params.r} --format={params.format} \
        --stranded={params.stranded} --type={params.type} \
        --idattr={params.idattr} {input.bam} {params.gtf} > {output}) 2> {log}"

#rule AddOrReplaceReadGroups:
    input:
        bam="bam/{sample}_{ref_genome_name}_uniq.bam",
        bai="bam/{sample}_{ref_genome_name}_uniq.bam.bai",
        htseq="htseq_curated_rmsk_class_id/{sample}_{ref_genome_name}.txt"
#    output: "{sample}_{ref_genome_name}_rg_added_sorted.bam"
    params:
        SO="coordinate",
        RGID="id",
        RGLB="library",
        RGPL="platform",
        RGPU="machine",
        RGSM="sample"
    log:
        "logs/{sample}_{ref_genome_name}_AddOrReplaceReadGroups.log"
    shell:
        "(java -jar /media/louis/2dd2To/bin/picard.jar  AddOrReplaceReadGroups \
            I={input.bam} \
            O={output} \
            SO={params.SO} \
            RGID={params.RGID} \
            RGLB={params.RGLB} \
            RGPL={params.RGPL} \
            RGPU={params.RGPU} \
            RGSM={params.RGSM}) 2> {log}"

rule MarkDuplicates:
    input: "~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq.bam"
    output: temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam")
    params:
            CREATE_INDEX="true",
            VALIDATION_STRINGENCY="SILENT",
            REMOVE_DUPLICATES="true",
            ASSUME_SORTED="true",
            M="output.metrics"
    log:
        "logs/{sample}_{ref_genome_name}_MarkDuplicates.log"
    shell:
        "(java -jar /home/emmanuel/software/picard/picard.jar MarkDuplicates \
            I={input} \
            O={output} \
            CREATE_INDEX={params.CREATE_INDEX} \
            VALIDATION_STRINGENCY={params.VALIDATION_STRINGENCY} \
            REMOVE_DUPLICATES={params.REMOVE_DUPLICATES} \
            ASSUME_SORTED={params.ASSUME_SORTED} \
            M={params.M}) 2> {log}"

rule samtools_index_csi:
    input: "~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam"
    output: temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam.csi")
    params:
        "" # optional params string
    shell:
        "(samtools index -c {input})"

rule samtools_index_bai:
    input: "~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam"
    output: temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam.bai")
    params:
        "" # optional params string
    shell:
        "(samtools index {input})"

rule SplitNCigarReads:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam"
    output:
        temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split.bam")
    params:
        ref="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa"
    log:
        "logs/{sample}_{ref_genome_name}_SplitNCigarReads.log"
    shell:
        "(java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar SplitNCigarReads \
            -R {params.ref} \
            -I {input.bam} \
            -O {output}) 2> {log}"

rule BaseRecalibrator:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split.bam"
    output: temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split_recal.table")
    params:
        ref="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa",
        #knownIndels="/media/louis/2dd2To/bin/gatk-4.1.2.0/Mills_and_1000G_gold_standard.indels.{ref_genome_name}.vcf.gz",
        dbsnp146="/run/user/1000/gvfs/smb-share:server=172.27.25.120,share=data/Louis/dbSNP/dbsnp_146.{ref_genome_name}.vcf.gz"
    log:
        "logs/{sample}_{ref_genome_name}_BaseRecalibrator.log"
    shell:
        "(java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar BaseRecalibrator -I {input.bam} \
            -R {params.ref} \
            --known-sites {params.knownIndels} \
            --known-sites {params.dbsnp146} \
            -O {output}) 2> {log}"

rule ApplyBQSR:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split.bam",
        bqsr="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split_recal.table"
    output:
        temp("~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split_recal.bam")
    params:
        ref="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa"
    log:
        "logs/{sample}_{ref_genome_name}_ApplyBQSR.log"
    shell:
        "(java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar ApplyBQSR \
            -R {params.ref} \
            -I {input.bam} \
            -bqsr {input.bqsr} \
            -O {output}) 2> {log}"

rule HaplotypeCaller:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp_split_recal.bam"
    params:
        ref="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa",
        stand_call_conf=10.0
    output:
        temp("mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}.vcf")
    log:
        "logs/{sample}_{ref_genome_name}_HaplotypeCaller.log"
    shell:
        "(java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar HaplotypeCaller \
            -R {params.ref} \
            -I {input.bam} \
            --dont-use-soft-clipped-bases \
            -stand-call-conf {params.stand_call_conf} \
            -O {output}) 2> {log}"

rule IndexFeatureFile:
    input: "mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}.vcf"
    output: temp("mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}.vcf.idx")
    shell:
        "java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar IndexFeatureFile \
            -F {input}"

rule VariantFiltration:
    input:
        vcf="mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}.vcf",
        index="mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}.vcf.idx"
    output:
        "mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}_filtered.vcf"
    params:
        ref="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa",
    log:
        "logs/{sample}_{ref_genome_name}_VariantFiltration.log"
    shell:
        "(java -jar /home/emmanuel/software/gatk-package-4.1.5.0/gatk-package-4.1.5.0-local.jar VariantFiltration \
        -R {params.ref} \
        -V {input.vcf} \
        -window 50 \
        -cluster 3 \
        --filter-name FS \
        -filter 'FS > 30.0' \
        --filter-name QD \
        -filter 'QD < 2.0' \
        -O {output}) 2> {log}"

rule mpileup:
    input:
        bam="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam",
        bai="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam.bai",
        csi="~/Documents/Omics_dataset/DNASeq/1kGProject/Brits_EUAnc/bam/{sample}_sorted_uniq_dedupp.bam.csi"
    output: "mpileup/1kGP_GBEU/{sample}_1kGEU_{ref_genome_name}_SNP_q10.vcf"
    params:
        threads=12,
        adjust=0,
        fasta="/home/emmanuel/Documents/Genomes/{ref_genome_name}/{ref_genome_name}.fa",
        filter="\'QUAL<10 || DP<10\'",
        region="chrX"
    log: "logs/{sample}_{ref_genome_name}_mpileup.log"
    shell:
        "(bcftools mpileup --threads {params.threads} \
        -I -C {params.adjust} \
        -r {params.region} \
        -f {params.fasta} {input.bam} | \
        bcftools call --threads {params.threads} -m | \
        bcftools filter --threads {params.threads} \
        -e {params.filter} > {output}) 2> {log}"

#rule bedtools_intersect:
    input:
        a="mpileup/{sample}_{ref_genome_name}_chrX_q10.vcf"
#    output: "mpileup/{sample}_{ref_genome_name}_informativechrX_q10.vcf"
    params:
        b="SNP/bed_informativeSNPH9_{ref_genome_name}.bed"
    log: "logs/{sample}_{ref_genome_name}_bedtools_intersect.log"
    shell:
        "(bedtools intersect -a {input.a} \
        -b {params.b} > {output}) \
        2> {log}"


#rule CoverageBamFile10reads:
#    input: "bam/{sample}_{ref_genome_name}_uniq.bam"
#    output: "bedCov/{sample}_{ref_genome_name}_10reads.bed"
#    log:
#        "logs/{sample}_{ref_genome_name}_CoverageCamFile10read.log"
#    shell:
#        "bedtools genomecov -ibam {input} \
#            -bg | awk -F '\t' '($4 >= 10) {{print $0}}' \
#            > {output} 2> {log}"


#rule intersectVcfBed:
#    input:
#        vcfdbsnp="~/ngs_bin/gatk-4.0.8.1/dbsnp_146.{ref_genome_name}.vcf",
#        bed="bedCov/{sample}_{ref_genome_name}_10reads.bed"
#    output: "bedCov/{sample}_{ref_genome_name}_dbsnp146.vcf"
#    log:
#        "logs/{sample}_{ref_genome_name}_intersectVcfBed.log"
#    shell:
#        "bedtools intersect -a {input.vcfdbsnp} \
#            -b {input.bed} \
#            > {output} 2> {log}"
