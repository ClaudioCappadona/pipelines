#!/bin/bash 


pwd; hostname; date

grouping=/mnt/path/rare_var_analysis/grouping
exomes_info=/mnt/path/exomes_info/new
vcf_file=$grouping/allChr_95_95_rare.vcf.gz
ref=/mnt/path/reference

module load conda
source activate python3.7

	
####################### calculate_burder_per_pathway_and_group #######################################

pathways="IFN"
#"hemostatic_genes our_PRS_genes our_procoagulants our_anticoagulants complement lectin surfactant Renieri"
#pathways="complement lectin surfactant Renieri hemostatic_genes our_PRS_genes"
#pathways="hemostatic_genes our_PRS_genes our_procoagulants our_anticoagulants complement lectin surfactant Renieri"
#suffix="_snp_indel"
suffix="_snponly"
#add=""
add="docker_fork_skatoh_user"
#testset="skat b.collapse b.wcnt"
testset="skat"


burder_per_pathway_and_group () {

#cd $HOME/software/EPACTS_SKATOh
#docker build -t fastapi --force-rm -f ./dockerfile .

cd $grouping

#get dockerfile
#wget https://raw.githubusercontent.com/ClaudioCappadona/EPACTS_SKATOh/SKATOh_test/dockerfile

for pathway in $pathways; do
	for test in $testset; do
		if [ ! -d $grouping/$pathway/$test/$add ]; then
			mkdir -p $grouping/$pathway/$test/$add
		fi


	docker run  --rm --user 101407 \
	-v /home/$PI/$USER:/home/$PI/$USER:Z \
	-v /mnt/path/$PI/$USER/ \
	-w /EPACTS_SKATOh \
	-p 8034:8034 \
	epacts-docker_dockerfile_debug:latest  /usr/local/bin/epacts group \
	--vcf $vcf_file \
	--groupf $grouping/$pathway/${pathway}${suffix}.grp \
	--out $grouping/$pathway/$test/$add/${pathway}${suffix}.${test} \
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

	done
done
}
burder_per_pathway_and_group 

wait 

####################### top variants to gene annotation ##########################
#top_vars_to_gene () {
#for pathway in $pathways; do
#	for test in $testset; do
		#top_vars =( $( cat $grouping/$pathway/$test/$add/${pathway}${suffix}.${test} | awk '$10 < 0.05' | cut -f 7 ) )

#	done
#done
#}
#top_vars_to_gene
####################### pathways burden results recap  #############################


pathways_burden_results_recap () {

pathways=(IFN complement lectin surfactant Renieri hemostatic_genes our_PRS_genes)


#suffix="_debug"
#add=""
add="docker_fork"
#testset="skat b.collapse b.wcnt"
testset="skat"

cat $grouping/${pathways[0]}/skat/$add/${pathways[0]}${suffix}.skat.epacts | head -n 1 > $grouping/pathways_skat${suffix}.recap
cd $grouping 

for pathway in ${pathways[@]}; do
	for test in $testset; do
		cat $grouping/$pathway/$test/$add/${pathway}${suffix}.${test}.epacts | grep -v "#" >> $grouping/pathways_skat${suffix}.recap
	done
done 

}
#pathways_burden_results_recap








                                
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
