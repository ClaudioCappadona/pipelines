#!/bin/bash 

pwd; hostname; date

grouping=/mnt/path/rare_var_analysis/grouping
exomes_info=/mnt/path/exomes_info/new
vcf_file=$grouping/allChr_95_95_rare.vcf.gz
ref=/mnt/path/esomi_project_backedup/reference

#module load conda
#source activate python3.7


	
####################### calculate_burder_per_hallmark_and_group #######################################

hallmarks=$(cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | cut -f 1 | cut -f 2- -d "_")

#suffix="_snp_indel"
suffix="_snponly"
#add=""
add="docker_fork_skatoh"
#testset="skat b.collapse b.wcnt"
testset="skat"


burder_per_hallmark_and_group () {

#cd $HOME/software/EPACTS_SKATOh
#docker build -t fastapi --force-rm -f ./dockerfile .

cd $grouping/hallmarks_GSEA

#making folders if they are absent
for hallmark in $hallmarks; do
	for test in $testset; do
		if [ ! -d $grouping/hallmarks_GSEA/$hallmark/$test/$add ]; then
			mkdir -p $grouping/hallmarks_GSEA/$hallmark/$test/$add
		fi

		if [ ! -s $grouping/hallmarks_GSEA/$hallmark/$test/$add/${hallmark}${suffix}.${test} ]; then
	
	docker run  --rm \
	-v /home/$PI:/home/$PI:Z \
	-v /mnt/path/$user/ \
	-w /EPACTS_SKATOh \
	-p 8034:8034 \
	epacts-docker_dockerfile_debug:latest  /usr/local/bin/epacts group \
	--vcf $vcf_file \
	--groupf $grouping/hallmarks_GSEA/$hallmark/${hallmark}${suffix}.grp \
	--out $grouping/hallmarks_GSEA/$hallmark/$test/$add/${hallmark}${suffix}.${test} \
	--max-maf 0.01 \
	--ped $exomes_info/exomes_process_info_filter_noNA.covar.ped \
	--pheno PHENO \
	--missing ./. \
	--no-plot \
	--cov SEX \
	--cov age_at_2020 \
	--cov age_sex \
	--cov age_age \
	--cov PC1 \
	--cov PC2 \
	--cov PC3 \
	--cov PC4 \
	--cov PC5 \
	--cov PC6 \
	--cov PC7 \
	--cov PC8 \
	--cov PC9 \
	--cov PC10 \
	--run 4 \
	--skat-o \
	--skato-h \
	--test $test
#--ped $exomes_info/exomes_process_info_filter.covar.ped \	

	fi

	done
done

wait 
}
#burder_per_hallmark_and_group 

wait 

####################### hallmarks burden results recap  #######################################

hallmarks_burden_results_recap () {

hallmarks=($(cat $grouping/hallmarks_GSEA/h.all.v2022.1.Hs.symbols.gmt | cut -f 1 | cut -f 2- -d "_"))

#testset="skat b.collapse b.wcnt"
testset="skat"

cat $grouping/hallmarks_GSEA/${hallmarks[0]}/skat/$add/${hallmarks[0]}${suffix}.skat.epacts | head -n 1 > $grouping/hallmarks_GSEA/skat_${add}${suffix}.recap.tmp
cd $grouping/hallmarks_GSEA 

for hallmark in ${hallmarks[@]}; do
	for test in $testset; do
		cat $grouping/hallmarks_GSEA/$hallmark/$test/$add/${hallmark}${suffix}.${test}.epacts | grep -v "#" >> $grouping/hallmarks_GSEA/skat_${add}${suffix}.recap.tmp
	done
done 

cat $grouping/hallmarks_GSEA/skat_${add}${suffix}.recap.tmp | sort -gk 10 |  sed 's/#//g' > $grouping/hallmarks_GSEA/skat_${add}${suffix}.recap
rm $grouping/hallmarks_GSEA/skat_${add}${suffix}.recap.tmp
}
hallmarks_burden_results_recap            















  

######################################### cases allele count per pathway #########################################################
#pathways="IFN hemostatic_genes our_PRS_genes our_procoagulants our_anticoagulants"
pathways="IFN"

pathway_cases_count () {
for pathway in $pathways; do
#groups="lof_high"
#groups="lof_high lof_high_missense_del lof_high_missense ALL_single_variants"

	python /home/path/software/TRAPD/count_cases_CC.py \
	-v $vcf_file \
	--snpfile $grouping/$pathway/${pathway}${suffix}.grp \
	-o $grouping/$pathway/${pathway}${suffix}_controlcount.txt \
	--samplefile $exomes_info/controls_sample_id.txt \
	--snpformat "EPACTS" &

	python /home/path/software/TRAPD/count_cases_CC.py \
	-v $vcf_file \
	--snpfile $grouping/$pathway/${pathway}${suffix}.grp \
	-o $grouping/$pathway/${pathway}${suffix}_casecount.txt \
	--samplefile $exomes_info/cases_sample_id.txt \
	--snpformat "EPACTS" &
	 
done

wait 

for pathway in $pathways; do
paste $grouping/$pathway/${pathway}${suffix}_controlcount.txt $grouping/$pathway/${pathway}${suffix}_casecount.txt > $grouping/$pathway/${pathway}${suffix}_controlcase_count &
done

}
#pathway_cases_count


merge_pathway_epacts_and_counts () {
cat $grouping/$pathway/$test/${pathway}${add}${suffix}.${test}.epacts | cut -f 1-9 | paste - $grouping/$pathway/${pathway}${add}${suffix}_casecount.txt > $grouping/$pathway/$test/tempfile
wait
cat $grouping/$pathway/$test/${pathway}${add}${suffix}.${test}.epacts | cut -f 10-13 | paste $grouping/$pathway/$test/tempfile - > $grouping/$pathway/$test/${pathway}${add}${suffix}.${test}.epacts_info
}
#merge_pathway_epacts_and_counts



exit
