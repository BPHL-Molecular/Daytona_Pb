# Daytona_Pb
A pipeline is used to analyze SARS-CoV-2 sequencing data from PacBio equipment.

## What to do
The pipeline can analyze SARS-CoV-2 sequencing data from PacBio equipment. The sample's sequencing quality, mapping/alignment with reference, coverage, SNPs, annotation, consensus sequence, lineage, phylogeny, etc can be analyzed.   

## Prerequisites
Nextflow should be installed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed.

Singularity/Apptainer is also needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.

PacBio SMRTLINK stand-alone tools (lima, mimux, pbindex, pbmm2, VCFCons) are needed. About how to install them, please see the file "How_to_install_smrtlink_tools.txt".

## Before running
1. Add the reads file in bam format (*.bam) from PacBio equipment to the folder "hifireads"
2. Add the file of M13 barcodes in fasta (*.fasta) format to the folder "M13barcodes".
3. Add enrichment probes file in fasta format (*.fasta) to the folder "enrichment_probes". 
4. Add biosample sheet file in csv format (*.csv) to the folder "biosample_sheet". Only two columns are allowed in the sheet. One column is barcode name, the other column is sample name. A sample sheet can be found in the folder "biosample_sheet".
         
## How to run
### 
1. open file "parames.yaml", set the full paths to the folders of input, probe, output, and reference. 
2. get into the top of the pipeline directory, then run 
```bash
./daytonapb.sh
```
If SLURM is used, you can run
```bash
sbatch ./daytonapb.sh
```

## Results
All results can be found in the directory /output.
