
#!/bin/bash

stringtie_list=$1
scallop_list=$2
annot_file=$3
assembly=$4

mkdir meta_annotation


stringtie --merge -o meta_annotation/stringtiemerge_stringtie\_${assembly}_noref.gtf -f 0.05 ${stringtie_list}
stringtie --merge -G ${annot_file} -o meta_annotation/stringtiemerge_stringtie\_${assembly}_ref.gtf -f 0.05 ${stringtie_list}
stringtie --merge -o meta_annotation/stringtiemerge_scallop\_${assembly}_noref.gtf -f 0.05 ${scallop_list}
stringtie --merge -G ${annot_file} -o meta_annotation/stringtiemerge_scallop\_${assembly}_ref.gtf -f 0.05 ${scallop_list}

/home/emmanuel/software/taco-v0.7.3.Linux_x86_64/taco_run -v -o meta_annotation/taco_scallop --gtf-expr-attr cov --filter-min-expr 0.0 --isoform-frac 0.05 ${scallop_list}
/home/emmanuel/software/taco-v0.7.3.Linux_x86_64/taco_run -v -o meta_annotation/taco_stringtie --gtf-expr-attr cov --filter-min-expr 0.0 --isoform-frac 0.05 ${stringtie_list}
