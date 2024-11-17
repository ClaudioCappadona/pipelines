#!/bin/bash


pwd; hostname; date

module load conda/anaconda3
source activate python3.7

annotation=/home/path/annotation/vep/annotations
grouping=/home/path/rare_var_analysis/grouping
vcf_file=$grouping/allChr_95_95_rare.vcf.gz

cd $grouping


################# Group file per gene divided by vars classes #######################################

ALL_single_variants () {

	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	#gzip -cd $vcf_file | grep -v "#" | cut -f 3 > ${FUNCNAME[0]}/ids &
	gzip -cd allChr_95_95_rare.vcf.gz | grep -v "#" | cut -f 3  | sed 's/_/:/' | sed 's/_/\//g' | sed 's/\//_/' > ${FUNCNAME[0]}/ids &
	gzip -cd $vcf_file | grep -v "#" | cut -f 8 | cut -f 2-4 -d "|" > ${FUNCNAME[0]}/var_info

	wait
	paste ${FUNCNAME[0]}/ids ${FUNCNAME[0]}/var_info -d "|" | paste - ${FUNCNAME[0]}/ids > ${FUNCNAME[0]}/${FUNCNAME[0]}.grp

	wait 
	if [ -s ${FUNCNAME[0]}/ALL_single_variants.grp ]; then
		rm ${FUNCNAME[0]}/ids &
		rm ${FUNCNAME[0]}/var_info
	fi
	wait
}
#ALL_single_variants &


ALL_single_rare_variants () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $vcf_file \
	-o ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--snponly
}
#ALL_single_rare_variants


lof_high () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
		
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $grouping/$vcf_file \
	-o $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "IMPACT[=]HIGH" \
	--excludevep "Consequence[=]stop lost" \
	--excludevep "Consequence[=]start_lost" \
	--excludevep "Consequence[=]transcript_amplification" \
	--excludevep "Consequence[=]transcript_ablation" \
	--snponly
}
#lof_high &

lof_high () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
		
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $grouping/$vcf_file \
	-o $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snp_indel.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "IMPACT[=]HIGH" \
	--excludevep "Consequence[=]stop lost" \
	--excludevep "Consequence[=]start_lost" \
	--excludevep "Consequence[=]transcript_amplification" \
	--excludevep "Consequence[=]transcript_ablation"
	#--snponly
}
#lof_high &

lof_high_single_vars () {

cd $grouping
if [ -s $grouping/lof_high/lof_high_snponly.grp ]; then
	grpfile=$grouping/lof_high/lof_high_snponly.grp
	touch $grouping/lof_high/${FUNCNAME[0]}_snponly.grp
	
	while read -r line; do
		vars=( $line )
			for (( i=1 ; i<${#vars[@]} ; i++)); do
				echo ${vars[0]}_${vars[$i]} $'\t' ${vars[$i]} >> $grouping/lof_high/${FUNCNAME[0]}_snponly.grp
			done
	done < $grpfile
fi	                                
}
#lof_high_single_vars

missense_deleterious () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $vcf_file \
	-o ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "REVEL_score[>]0.5" \
	--includevep "Consequence[=]missense_variant" \
	--snponly
}
#missense_deleterious &

missense_deleterious () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $vcf_file \
	-o ${FUNCNAME[0]}/${FUNCNAME[0]}_snp_indel.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "REVEL_score[>]0.5" \
	--includevep "Consequence[=]missense_variant"
	#--snponly
}
#missense_deleterious &

missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
	
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $vcf_file \
	-o ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "Consequence[=]missense_variant" \
	--snponly
}
#missense &

missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
	
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $vcf_file \
	-o ${FUNCNAME[0]}/${FUNCNAME[0]}_snp_indel.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includeinfo "AF[<]0.01" \
	--includevep "Consequence[=]missense_variant"
	#--snponly
}
#missense 

wait

########### merging grp files ###############

lof_high_missense_del () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/merge_snp_file_CC.py \
	--snpfiles lof_high/lof_high_snponly.grp,missense_deleterious/missense_deleterious_snponly.grp \
	--outfile ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp
	
}
#lof_high_missense_del

lof_high_missense_del () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/merge_snp_file_CC.py \
	--snpfiles lof_high/lof_high_snp_indel.grp,missense_deleterious/missense_deleterious_snp_indel.grp \
	--outfile ${FUNCNAME[0]}/${FUNCNAME[0]}_snp_indel.grp
	
}
#lof_high_missense_del


lof_high_missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/merge_snp_file_CC.py \
	--snpfiles lof_high/lof_high_snponly.grp,missense/missense_snponly.grp \
	--outfile ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly.grp
	
}
#lof_high_missense

lof_high_missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/merge_snp_file_CC.py \
	--snpfiles lof_high/lof_high_snp_indel.grp,missense/missense_snp_indel.grp \
	--outfile ${FUNCNAME[0]}/${FUNCNAME[0]}_snp_indel.grp
	
}
#lof_high_missense




######################## grp files for variants with MAF < 5% ###############


annotations=/home/path/annotation/vep/annotations
vcf_file_annot=allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier_annotated_canonic.vcf.gz
exomes_info=/home/path/exomes_info/new


lof_high () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
		
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $annotations/$vcf_file_annot \
	-o $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snponly_5%.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includesamples $exomes_info/exomes_process_info_filter_noNA.covar.samples \
	--includeinfo "AF[<]0.05" \
	--includevep "IMPACT[=]HIGH" \
	--excludevep "Consequence[=]stop lost" \
	--excludevep "Consequence[=]start_lost" \
	--excludevep "Consequence[=]transcript_amplification" \
	--excludevep "Consequence[=]transcript_ablation" \
	--snponly
}
#lof_high &


missense_deleterious () {
        if [ ! -d ${FUNCNAME[0]} ]; then
                mkdir ${FUNCNAME[0]}
        fi

        python /home/path/software/TRAPD/make_snp_file_CC.py \
        -v $annotations/$vcf_file_annot \
        -o $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snponly_5%.grp \
        --vep  \
        --snpformat "EPACTS" \
        --genecolname SYMBOL \
        --includesamples $exomes_info/exomes_process_info_filter_noNA.covar.samples \
        --includeinfo "AF[<]0.05" \
        --includevep "REVEL_score[>]0.5" \
        --includevep "Consequence[=]missense_variant" \
        --snponly
}
#missense_deleterious &


missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi
	
	python /home/path/software/TRAPD/make_snp_file_CC.py \
	-v $annotations/$vcf_file_annot \
	-o $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snponly_5%.grp \
	--vep  \
	--snpformat "EPACTS" \
	--genecolname SYMBOL \
	--includesamples $exomes_info/exomes_process_info_filter_noNA.covar.samples \
	--includeinfo "AF[<]0.05" \
	--includevep "Consequence[=]missense_variant" \
	--snponly
}
#missense


########### merging grp files ###############

lof_high_missense_del () {
        if [ ! -d ${FUNCNAME[0]} ]; then
                mkdir ${FUNCNAME[0]}
        fi

        python /home/path/software/TRAPD/merge_snp_file_CC.py \
        --snpfiles $grouping/lof_high/lof_high_snponly_5%.grp,$grouping/missense_deleterious/missense_deleterious_snponly_5%.grp \
        --outfile $grouping/${FUNCNAME[0]}/${FUNCNAME[0]}_snponly_5%.grp
        
}
lof_high_missense_del


lof_high_missense () {
	if [ ! -d ${FUNCNAME[0]} ]; then
		mkdir ${FUNCNAME[0]}
	fi

	python /home/path/software/TRAPD/merge_snp_file_CC.py \
	--snpfiles $grouping/lof_high/lof_high_snponly_5%.grp,$grouping/missense/missense_snponly_5%.grp \
	--outfile ${FUNCNAME[0]}/${FUNCNAME[0]}_snponly_5%.grp
	
}
lof_high_missense


wait


############################ Pathways_grp_file ########################################################

#pathways="IFN hemostatic_genes our_PRS_genes our_procoagulants our_anticoagulants Renieri lectin surfactant"
#pathways="IFN complement lectin surfactant Renieri hemostatic_genes our_PRS_genes"
pathways="custom_IFN"

groups="lof_high lof_high_missense_del lof_high_missense"
#groups="ALL_single_variants"

pathways_grp_file_snponly () {

for pathway in $pathways; do
	echo -n "" > $pathway/${pathway}_snponly.grp

	for group in $groups; do
		echo ${pathway}_${group} > $pathway/path_name.txt
		grep -w -f $pathway/${pathway}.txt $group/${group}_snponly.grp | \
		tr "\t" "\n" | \
		grep -v -w -f $pathway/${pathway}.txt | \
		tr "\n" "\t" | \
		paste $pathway/path_name.txt - >> $pathway/${pathway}_snponly.grp
	done
done
}
#pathways_grp_file_snponly


pathways_grp_file_snp_indel () {

for pathway in $pathways; do
	echo -n "" > $pathway/${pathway}_snp_indel.grp

	for group in $groups; do
		echo ${pathway}_${group} > $pathway/path_name.txt
		grep -w -f $pathway/${pathway}.txt $group/${group}_snp_indel.grp | \
		tr "\t" "\n" | \
		grep -v -w -f $pathway/${pathway}.txt | \
		tr "\n" "\t" | \
		paste $pathway/path_name.txt - >> $pathway/${pathway}_snp_indel.grp
	done
done
}
#pathways_grp_file_snp_indel

############################ Hallmarks_grp_file ########################################################
wait
groups="lof_high lof_high_missense_del lof_high_missense" 


hallmarks_GSEA_grp_file_snponly () {

hallmarks=$(cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | cut -f 1 | cut -f 2- -d "_")

for hallmark in $hallmarks; do

	if [ ! -d $grouping/hallmarks_GSEA/$hallmark ]; then
		mkdir -p $grouping/hallmarks_GSEA/$hallmark
	fi

	echo -n "" > $grouping/hallmarks_GSEA/$hallmark/${hallmark}_snponly.grp

	cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | grep $hallmark | cut -f 3- | tr "\t" "\n" > $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt

	for group in $groups; do
		echo ${hallmark}_${group} > $grouping/hallmarks_GSEA/$hallmark/tmp.txt
		grep -w -f $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt $grouping/$group/${group}_snponly.grp | \
		tr "\t" "\n" | \
		grep -v -w -f $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt | \
		tr "\n" "\t" | \
		paste $grouping/hallmarks_GSEA/$hallmark/tmp.txt - >> $grouping/hallmarks_GSEA/$hallmark/${hallmark}_snponly.grp
	done
done
}
#hallmarks_GSEA_grp_file_snponly

hallmarks_GSEA_grp_file_snp_indel () {

hallmarks=$(cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | cut -f 1 | cut -f 2- -d "_")

for hallmark in $hallmarks; do

	if [ ! -d $grouping/hallmarks_GSEA/$hallmark ]; then
		mkdir -p $grouping/hallmarks_GSEA/$hallmark
	fi

	echo -n "" > $grouping/hallmarks_GSEA/$hallmark/${hallmark}_snp_indel.grp

	cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | grep $hallmark | cut -f 3- | tr "\t" "\n" > $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt

	for group in $groups; do
		echo ${hallmark}_${group} > $grouping/hallmarks_GSEA/$hallmark/tmp.txt
		grep -w -f $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt $grouping/$group/${group}_snp_indel.grp | \
		tr "\t" "\n" | \
		grep -v -w -f $grouping/hallmarks_GSEA/$hallmark/${hallmark}.txt | \
		tr "\n" "\t" | \
		paste $grouping/hallmarks_GSEA/$hallmark/tmp.txt - >> $grouping/hallmarks_GSEA/$hallmark/${hallmark}_snp_indel.grp
	done
done
}
#hallmarks_GSEA_grp_file_snp_indel

wait
#}




exit
