process pbsmrtlink {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${params.output}/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #echo ${params.input}/${x}_1.fastq.gz >> xfile.txt
   
   mkdir -p ${params.output}/assemblies
   mkdir -p ${params.output}/variants
   mkdir -p ${params.output}/vadr_error_reports
   mkdir -p ${params.output}/${x}
   cp ${params.input}/${x}.bam ${params.output}/${x}
   
   #SMTLINK standalone tools
   #trim molecular loop adapters and tag reads with UMIs
   mimux --probes ${params.probe}/*.fasta --probe-report ${x}.mimus_probe_report.tsv --log-level INFO ${params.output}/${x}/${x}.bam ${params.output}/${x}/${x}.mimux_trimmed.bam
   pbindex ${params.output}/${x}/${x}.mimux_trimmed.bam
   
   #align trimmed reads to Wuhan reference
   pbmm2 align --sort --preset HiFi --log-level INFO ${params.reference}/*.fasta ${params.output}/${x}/${x}.mimux_trimmed.bam ${params.output}/${x}/${x}.mapped.bam
   
   #calculate statistics
   singularity exec docker://staphb/samtools:1.12 samtools coverage ${params.output}/${x}/${x}.mapped.bam -o ${params.output}/${x}/${x}.coverage.txt
   
   #samtools mpileup --min-BQ 1 -f ${params.reference}/*.fasta -s -o ${params.output}/${x}/${x}.mapped.bam.mpileup ${params.output}/${x}/${x}.mapped.bam
   #samtools depth -q 0 -Q 0 -o ${params.output}/${x}/${x}.mapped.bam.depth ${params.output}/${x}/${x}.mapped.bam
   singularity exec docker://staphb/samtools:1.12 samtools mpileup --min-BQ 1 -f ${params.reference}/*.fasta -s -o ${params.output}/${x}/${x}.mapped.bam.mpileup ${params.output}/${x}/${x}.mapped.bam
   singularity exec docker://staphb/samtools:1.12 samtools depth -q 0 -Q 0 -o ${params.output}/${x}/${x}.mapped.bam.depth ${params.output}/${x}/${x}.mapped.bam
   
   #call variants
   singularity exec docker://staphb/bcftools:1.16 bcftools mpileup --open-prob 25 --indel-size 450 --gap-frac 0.01 --ext-prob 1 --min-ireads 3 --max-depth 1000 --max-idepth 5000 --seed 1984 -h 500 -B -a FORMAT/AD -f ${params.reference}/*.fasta ${params.output}/${x}/${x}.mapped.bam | singularity exec docker://staphb/bcftools:1.16 bcftools call -mv -Ov | singularity exec docker://staphb/bcftools:1.16 bcftools norm -f ${params.reference}/*.fasta - | singularity exec docker://staphb/bcftools:1.16 bcftools filter -e 'QUAL < 20' - > ${params.output}/${x}/${x}.variants_bcftools.vcf

   #bcftools mpileup --open-prob 25 --indel-size 450 --gap-frac 0.01 --ext-prob 1 --min-ireads 3 --max-depth 1000 --max-idepth 5000 --seed 1984 -h 500 -B -a FORMAT/AD -f ${params.reference}/*.fasta ${params.output}/${x}/${x}.mapped.bam | bcftools call -mv -Ov | bcftools norm -f ${params.reference}/*.fasta - | bcftools filter -e 'QUAL < 20' - > ${params.output}/${x}/${x}.variants_bcftools.vcf

   #generate viral consensus sequence based on major allele in VCF
   VCFCons ${params.reference}/*.fasta ${params.output}/${x}/${x} --min_coverage 2 --min_alt_freq 0.5 --vcf_type bcftools --input_depth ${params.output}/${x}/${x}.mapped.bam.depth --input_vcf ${params.output}/${x}/${x}.variants_bcftools.vcf

   #align viral consensus sequence to Wuhan reference
   pbmm2 align --sort --preset HiFi --log-level INFO ${params.reference}/*.fasta ${params.output}/${x}/${x}.vcfcons.frag.fasta ${params.output}/${x}/${x}.consensus_mapped.bam

   """
}
