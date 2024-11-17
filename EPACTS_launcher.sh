#!/bin/bash 

grouping=/mnt/path/rare_var_analysis/grouping
exomes_info=/mnt/path/exomes_info/new

ref=/mnt/path/reference
scripts=/mnt/path/rare_var_analysis/scripts

annotation=/home/path/annotation/vep/annotations

#freq="_1%"
freq="_5%"

if [ "$freq" == "_1%" ]; then
	vcf_file=$grouping/allChr_95_95_rare.vcf.gz
	nfreq=0.01
fi

if [ "$freq" == "_5%" ]; then
	vcf_file=$annotation/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier_annotated_canonic.vcf.gz
	nfreq=0.05
fi



####################### test GENE burder_per vars group #######################################
#groups="lof_high lof_high_missense_del lof_high_missense"
groups="lof_high_missense"
#groups="  ALL_single_variants missense_deleterious missense"

testset="skat"
#suffix="_debug"
suffix=_snponly${freq}
#testset="skat b.collapse b.wcnt"

#add="docker_fork_skatoh"
add=docker_fork_skatoh${freq}



burder_per_group () {

#cd $HOME/software/EPACTS_SKATOh
#docker build -t fastapi --force-rm -f ./dockerfile .

cd $grouping

for group in $groups; do
 for test in $testset; do

 
if [ ! -d "$grouping/$group/$test/$add" ]; then
        mkdir -p $grouping/$group/$test/$add
fi

#if [ ! -s "$grouping/$group/$test/${group}.*epacts" ]; then

        docker run  --rm --user 101407 \
        -v /home/$PI:/home/$PI:Z \
        -v /mnt/path/$user/ \
        -w /EPACTS_SKATOh \
        -p 8033:8033 \
        epacts-docker_dockerfile_debug:latest  /usr/local/bin/epacts group \
        --vcf $vcf_file \
        --groupf $grouping/$group/${group}${suffix}${freq}.grp \
        --out $grouping/$group/${test}/$add/${group}${suffix}.${test} \
        --ped $exomes_info/exomes_process_info_filter_noNA.covar.ped \
        --pheno PHENO \
        --missing ./. \
        --max-maf $nfreq \
        --no-plot \
        --cov SEX \
        --cov age_at_2020 \
        --cov age_sex \
        --cov age_age \
        --unit 10 \
        --run 10 \
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
        --test $test \
        --skat-o \
        --skato-h
#--unit 50 \
#--min-mac 2 \    
#
#--run 10 \  
#--ped $exomes_info/exomes_process_info_filter.covar.ped \ 
	done
done

#--groupf $grouping/$group/${group}${add}${suffix}.grp \

#fi
#statgen/epacts:dev group \  \
#--ped $exomes_info/exomes_process_info_filter_noNA.covar.ped \
#--ped $exomes_info/exomes_process_info_filter.covar.ped \
}
#burder_per_group


#loading burden_per_group_single_vars function
#. $scripts/burden_per_group_single_vars.sh


####################### gene burden results filter and sort  #############################

gene_burden_results_filter_sort () {

for group in $groups; do
	for test in $testset; do
		cd $grouping/$group/${test}/$add/
		
		sed -i 's/#//g' $grouping/$group/${test}/$add/${group}${suffix}.${test}.epacts

		Rscript --vanilla $scripts/adjust_pvalue.R $grouping/$group/${test}/$add/${group}${suffix}.${test}.epacts

	done
done
}
#gene_burden_results_filter_sort




get_samples_per_var () {
cd $grouping

for group in $groups; do

	bash  $scripts/make_single_variants_grp.sh $grouping/$group/${group}${suffix}.grp  > $grouping/$group/${group}${suffix}_single_vars.grp
											
	for test in $testset; do
		if [ -s $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples ]; then
			rm $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples
			touch $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples
		fi
		while read -r line; do
		echo $line | cut -f 2- -d " " | tr "\t" "\n" | tr ":" "\t" | tr "_" "\t" | awk -F "\t" '{ {print $1"\t" $2"\t"$2}}' > $grouping/$group/${group}${suffix}_single_vars_coord.tmp
		varid=$(echo $line | cut -f 2 -d " ")
		varname=$(echo $line | cut -f 1 -d " ")
		perl $scripts/count_var_group_tabix.pl $vcf_file $grouping/$group/${group}${suffix}_single_vars_coord.tmp $varid $varname>> $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples

		done < $grouping/$group/${group}${suffix}_single_vars.grp
	done
done
}
get_samples_per_var

samples_occurrencies () {

cd $grouping

for group in $groups; do

	bash  $scripts/make_single_variants_grp.sh $grouping/$group/${group}${suffix}.grp  > $grouping/$group/${group}${suffix}_single_vars.grp
	for test in $testset; do
		if [ -s $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples ]; then
			echo "ok"
		fi
	done
done

}
#samples_occurrencies

wait 

exit
