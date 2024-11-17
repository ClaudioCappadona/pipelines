#!/bin/bash

module load conda
source activate bcftools

grouping=/home/path/rare_var_analysis/grouping
exomes_info=/home/path/exomes_info/new
annotation=/home/path/annotation/vep/annotations

group=$1
test=skat
suffix=_snponly
freq=$2
set=$3
feature=$4

#freq="1%"
#freq="5%"

if [ "$freq" == "1%" ]; then
	vcf_file=$grouping/allChr_95_95_rare.vcf.gz
fi

if [ "$freq" == "5%" ]; then
	vcf_file=$annotation/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier_annotated_canonic.vcf.gz
fi


if [ "$feature" == "vars" ]; then
	feat="1,3"
fi

if [ "$feature" == "genes" ]; then
	feat="3"
fi



groups=( "lof_high" "lof_high_missense_del" "lof_high_missense")
freqs=( "1%" "5%")
sets=("training" "test")
features=("vars" "genes")




header () {
#get samples ID from VCF file header

bcftools query \
	-S $exomes_info/ML/exomes_process_info_filter_noNA_ML_${set}.covar.samples \
	-f '%ID [%GT] \n' \
	-H \
	$vcf_file | \
	head -n 1 | tr ":" " " | \
	sed 's/GT//g' | tr " " "\n" | \
	cut -f 2 -d "]" | \
	sed '/^$/d' | \
	paste -s -d ' ' | tr -s " " | tr " " "\t"
}

#echo "a"

body () {
#get only genotypes information from body of VCF file
bcftools query \
	-R $grouping/$group/$test/docker_fork_skatoh_ML_training_${freq}/sign_single_variants \
	-f '%ID [ %GT]\n' \
	-S $exomes_info/ML/exomes_process_info_filter_noNA_ML_${set}.covar.samples \
	$vcf_file | \
	grep -f $grouping/${group}/$test/docker_fork_skatoh_ML_training_${freq}/sign_single_variants.vcfids | \
	tr -s " " | cut -f 2- -d " " | tr "\t" " " | tr " " "\t" 

	#| python $grouping/$group/$test/docker_fork_skatoh/feature_select.py
}

varid_gene () {
bcftools query \
	-R $grouping/$group/$test/docker_fork_skatoh_ML_training_${freq}/sign_single_variants \
	-f '%ID %INFO/CSQ \n' \
	-S $exomes_info/ML/exomes_process_info_filter_noNA_ML_${set}.covar.samples \
	$vcf_file | \
	grep -f $grouping/${group}/$test/docker_fork_skatoh_ML_training_${freq}/sign_single_variants.vcfids | \
	tr -s " " | cut -f 1,4 -d "|" | tr "|" " " | cut -f $feat -d " " | tr " " "_"
	
}


#if [ -z "$group" ]; then
if [[ ! " ${groups[@]} " =~ " ${group} " ]]; then
	 echo "The script requires 4 arguments:"
	 echo "- first argument (variants group) can be:" ${groups[@]}
	 echo "- second argument (freq) can be:" ${freqs[@]}
	 echo "- third argument (set) can be:" ${sets[@]}
	 echo "- fourth argument (features) can be:" ${features[@]}
	 exit
else
	if [ -z "$set" ]; then
		echo "second argument (set) can be either 'training' or 'test'"
	else
		header
		#body
		#paste varid_gene body
		paste <(varid_gene) <(body) --delimiters '\t' | tr "\n\n" "\n"
	fi
fi
