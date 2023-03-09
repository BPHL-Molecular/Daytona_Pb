#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=daytonapb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --output=daytonapb.%j.out
#SBATCH --error=daytonapb.%j.err
#SBATCH --time=3-00


#module load singularity
module load apptainer

#demultiplex hifi reads by lima
hifireads=$(realpath ./hifireads/*.bam)
m13barcodes=$(realpath ./M13barcodes/*.fasta)
probes=$(realpath ./enrichment_probes/*.fasta)
biosample=$(realpath ./biosample_sheet/*.csv)
dem=$(realpath ./hifireads_dem)
topdir=$(realpath ./)
lima $hifireads $m13barcodes dem.bam --hifi-preset ASYMMETRIC --biosample-csv $biosample --store-unbarcoded --dump-clips --log-level INFO --split-named
mv ${topdir}/dem* $dem

nextflow run daytonapb.nf -params-file params.yaml

rm ${dem}/dem*
sort ./output/*/report.txt | uniq > ./output/sum_report.txt
sed -i '/sampleID\treference/d' ./output/sum_report.txt
sed -i '1i sampleID\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag\tpangolin_version\tlineage\tSOTC' ./output/sum_report.txt

cat ./output/assemblies/*.fasta > ./output/assemblies.fasta
singularity exec /apps/staphb-toolkit/containers/nextclade_2021-03-15.sif nextclade --input-fasta ./output/assemblies.fasta --output-csv ./output/nextclade_report.csv