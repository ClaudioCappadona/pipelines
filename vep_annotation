#!/bin/bash
#SBATCH --job-name=vep                                          # Job name
##SBATCH --mail-type=END,FAIL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --partition=cpu
##SBATCH --mail-user=claudio.cappadona@humanitasresearch.it       # Where to send mail
#SBATCH --ntasks=5
#SBATCH --nodelist=node5                                            # Run on a N CPU
#SBATCH --mem=50gb                                          # Job memory request
##SBATCH --time=24:00:00                                # Time limit hrs:min:sec
#SBATCH --output=vep%j.log                                # Standard output and error log

pwd; hostname; date

vqsr="/home/path/ccappadona/esomi_project_backedup/joint_call_2023/vqsr"


vep_annotation () {
docker run --rm -v /home/path:/home/path:Z -v /mnt/backup/path/ccappadona/:/mnt/backup/path/ccappadona/ -w /data -p 8034:8034   labduga/vep104_root \
/opt/vep/src/ensembl-vep/vep --cache \
--refseq \
-i ${vqsr}/allChr_95_95recalibrated_demultiplex_max_mc_id_proj19.vcf.gz \
-o /home/path/ccappadona/esomi_project_backedup/annotation/vep/annotations_2023/allChr_95_95recalibrated_demultiplex_max_mc_id_proj19_annotated_canonic.vcf.gz \
--force_overwrite \
--use_transcript_ref \
--dir_cache /opt/vep/.vep \
--dir_plugins /opt/vep/.vep/Plugins \
--species homo_sapiens \
--use_transcript_ref \
--assembly GRCh38       \
--no_stats \
--fork 18 \
--offline \
--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--buffer_size 500 \
--vcf \
--compress_output bgzip \
--biotype \
--hgvs \
--symbol \
--canonical \
--regulatory \
--gene_phenotype \
--variant_class \
--plugin DisGeNET,file=/home/path/share/vep_data/DisGeNET/all_variant_disease_associations_final.tsv.gz \
--plugin CADD,/home/path/share/vep_data/CADD/whole_genome_SNVs.tsv.gz \
--plugin Mastermind,/home/path/share/vep_data/Mastermind/mastermind_cited_variants_reference-2021.08.03-grch38.vcf.gz,0,0,1 \
--plugin dbNSFP,/home/path/share/vep_data/dbNSFP/dbNSFP4.2a_grch38.gz,REVEL_score,REVEL_rankscore,SIFT_pred,FATHMM_pred,PROVEAN_pred,MutationAssessor_pred,VEST4_score,LRT_score,LRT_pred,MutationTaster_pred,Polyphen2_HVAR_pred \
--plugin dbscSNV,/home/path/share/vep_data/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
--plugin SpliceAI,snv=/home/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/home/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
--plugin MaxEntScan,/home/path/share/vep_data/MaxEntScan/fordownload,SWA
}
#vep_annotation

vep_annotation () {
docker run --rm -v /home/path:/home/path:Z -v /mnt/backup/path/ccappadona/:/mnt/backup/path/ccappadona/ -w /data -p 8034:8034   labduga/vep104_root \
/opt/vep/src/ensembl-vep/vep --cache \
--refseq \
-i ${vqsr}/allChr_95_95recalibrated_demultiplex_max_mc_id_proj19.vcf.gz \
-o /home/path/ccappadona/esomi_project_backedup/annotation/vep/annotations_2023/allChr_95_95recalibrated_demultiplex_max_mc_id_proj19_annotated_canonic.vcf.gz \
--force_overwrite \
--use_transcript_ref \
--dir_cache /opt/vep/.vep \
--dir_plugins /opt/vep/.vep/Plugins \
--species homo_sapiens \
--use_transcript_ref \
--assembly GRCh38       \
--no_stats \
--fork 18 \
--offline \
--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--buffer_size 500 \
--vcf \
--compress_output bgzip \
--biotype \
--af_1kg \
--af_gnomad \
--hgvs \
--symbol \
--canonical \
--regulatory \
--gene_phenotype \
--variant_class \
--plugin DisGeNET,file=/home/path/share/vep_data/DisGeNET/all_variant_disease_associations_final.tsv.gz \
--plugin CADD,/home/path/share/vep_data/CADD/whole_genome_SNVs.tsv.gz \
--plugin Mastermind,/home/path/share/vep_data/Mastermind/mastermind_cited_variants_reference-2021.08.03-grch38.vcf.gz,0,0,1 \
--plugin dbNSFP,/home/path/share/vep_data/dbNSFP/dbNSFP4.2a_grch38.gz,REVEL_score,REVEL_rankscore,SIFT_pred,FATHMM_pred,PROVEAN_pred,MutationAssessor_pred,VEST4_score,LRT_score,LRT_pred,MutationTaster_pred,Polyphen2_HVAR_pred \
--plugin dbscSNV,/home/path/share/vep_data/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
--plugin SpliceAI,snv=/home/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/home/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
--plugin MaxEntScan,/home/path/share/vep_data/MaxEntScan/fordownload,SWA
}
vep_annotation

#### notes ####
#--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
#fork 4 is the advised value for multithreding here: https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#faster
#buffer size 5000 is the advised value

#-entrypoint /usr/bin/perl
#--fasta  /home/path/ccappadona/reference/vep_references_refseq/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \

#-o /home/path/ccappadona/esomi_project_backedup/annotation/vep/annotations/CCRL2_annotated.vcf \



exit
