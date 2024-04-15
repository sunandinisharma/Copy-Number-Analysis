#!/bin/sh
#SBATCH --job-name=cnvkit
#SBATCH --mem=8G
#SBATCH --time=28:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sunandini@unmc.edu
#SBATCH --licenses=common
#SBATCH --chdir=/lustre/work/javeediqbal/shared/weiwei/P01

module load anaconda/4.3

source activate cnvkit


pfx=${1}

normalBam=/common/javeediqbal/sunandini/CopywriteR/Tonsil1_Group7_23/Alignment/Tonsil1_Group7_23.recalibrated.final.bam

bed=/common/javeediqbal/sunandini/Sequencing_Pipeline/SureSelectAllExonV7_S31285117_hs_hg38/for.cnvKITSureSelectAllExonV7_hg38_S31285117_Covered.bed

ref="/common/javeediqbal/sunandnini/reference_GRCh38_tophat/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/genome.fa"

bam=${pfx}/Alignment/${pfx}.recalibrated.final.bam

access=/common/javeediqbal/sunandini/PROGRAMS/cnvkit/data/access-100kb.hg38.bed

mkdir -p ${pfx}/cnvkit_vTonsil1

cnvkit.py batch $bam --normal $normalBam -t $bed --access $access -f $ref --output-reference ${pfx}/cnvkitSubmit_100kb_vTonsil1/20kb_reference.cnn -d ${pfx}/cnvkitSubmit_100kb_vTonsil1

cnvkit.py scatter ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cnr -s ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cns -o ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final-scatter.pdf

cnvkit.py diagram ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cnr -s ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cns -o ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final-diagram.pdf

cnvkit.py call ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cns -y -v ${pfx}/Variants/Mutect2_GDC_FLAGGEDdbSNPfiltered.PON/${pfx}.single_sample.filtered.vcf -o ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.mutect2.call.cns
cnvkit.py call ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.cns -y -v ${pfx}/Variants/${pfx}.varscan.snp.vcf -o ${pfx}/cnvkitSubmit_100kb_vTonsil1/${pfx}.recalibrated.final.varscan.call.cns

conda deactivate
