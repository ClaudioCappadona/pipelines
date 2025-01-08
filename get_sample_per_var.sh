#!/bin/bash
#SBATCH --job-name=get_sample_per_var                                         # Job name
##SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=claudio.cappadona@humanitasresearch.it       # Where to send mail
#SBATCH --ntasks=1
#SBATCH --nodelist=node8                                            # Run on a N CPU
#SBATCH --mem=30gb                                          # Job memory request
##SBATCH --time=24:00:00                                # Time limit hrs:min:sec
#SBATCH --output=get_sample_per_var%j.log                                # Standard output and error log

pwd; hostname; date

grouping=/mnt/backup/asselta/ccappadona/esomi_project_backedup/rare_var_analysis/grouping
exomes_info=/mnt/backup/asselta/ccappadona/esomi_project_backedup/exomes_info/new
vcf_file=$grouping/allChr_95_95_rare.vcf.gz
ref=/mnt/backup/asselta/ccappadona/esomi_project_backedup/reference
scripts=/mnt/backup/asselta/ccappadona/esomi_project_backedup/rare_var_analysis/scripts

#module load conda
#source activate python3.7


####################### SETTINGS  #######################################
groups="lof_high lof_high_missense_del lof_high_missense"
#groups="lof_high_missense_del"

#groups="  ALL_single_variants missense_deleterious missense"

testset="skat"
#suffix="_debug"

#freq="_1%"
freq="_5%"

#testset="skat b.collapse b.wcnt"
#add=""
add=docker_fork_skatoh_ML_training${freq}
#add="docker_fork_skatoh_5%"



suffix=_snponly${freq}


#######################################################################

for group in $groups; do
for test in $testset; do

        cd $grouping


        make_single_variants_grp () {

                grpfile=$grouping/$group/${group}_ML_training${suffix}.grp

                while read -r line; do

                        vars=( $line )
                        for (( i=1 ; i<${#vars[@]} ; i++)); do
                        echo ${vars[0]}_${vars[$i]} $'\t' ${vars[$i]}
                        done
                done < $grpfile
        }
        #make_single_variants_grp > $grouping/$group/${group}_ML_training${suffix}_single_vars.grp


        get_samples_per_var () {

                if [ -s $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples ]; then
                        rm $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples
                        touch $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples
                fi
                while read -r line; do
                echo $line | cut -f 2- -d " " | tr "\t" "\n" | tr ":" "\t" | tr "_" "\t" | awk -F "\t" '{ {print $1"\t" $2"\t"$2}}' > $grouping/$group/${test}/$add/${group}_single_vars_coord.tmp
                varid=$(echo $line | cut -f 2 -d " ")
                varname=$(echo $line | cut -f 1 -d " ")
                perl $scripts/count_var_group_tabix.pl $vcf_file $grouping/$group/${test}/$add/${group}_single_vars_coord.tmp $varid $varname>> $grouping/$group/${test}/$add/${group}${suffix}.${test}.samples

                done < $grouping/$group/${group}_ML_training${suffix}_single_vars.grp

        }
        #get_samples_per_var

        sign_single_variants () {

                cat $grouping/$group/${test}/$add/${group}${suffix}.$test.epacts.sort.padj | cut -f 4 | cut -f 2 -d "_" | awk '{print "^"$0"_"}' | \
                grep -f - $grouping/$group/${test}/$add/${group}${suffix}.$test.samples | tr "\t" "\n" | grep "/" | tee $grouping/$group/${test}/$add/sign_single_variants.txt | \
                cut -f 2- -d "_" | tr ":" "_" | tr "/" "_" | tee $grouping/$group/${test}/$add/sign_single_variants.vcfids | \
                cut -f 1,2 -d "_" | awk -F "_" '{ print $1"\t"$2"\t"$2 }' > $grouping/$group/${test}/$add/sign_single_variants

        }
        sign_single_variants

done
done
