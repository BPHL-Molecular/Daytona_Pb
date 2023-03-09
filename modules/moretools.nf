process moretools {
    input:
        val mypath
    output:
        stdout
        
    $/
    #!/usr/bin/env python3
    import subprocess
    
    items = "${mypath}".strip().split("/")
    #print(items[-1])
    
    filepath1 = "${mypath}"+"/"+items[-1]+".coverage.txt"
    #print(filepath1)
    with open(filepath1, 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        #print(ref_name)
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        #print(reads_mapped)
        mapq = stats[8]
        
   
    #### get the total number of reads of a BAM file
    proc_cx = subprocess.run('singularity exec docker://staphb/samtools:1.12 samtools view -c ' + "${mypath}/" + items[-1] + '.bam', shell=True, capture_output=True, text=True, check=True)
    wc_out = proc_cx.stdout.rstrip()
    clean_reads = int(wc_out)
    
    #Get percentage of mapped reads/clean reads
    percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
    #print(percent_map)
    
    #Gather QC metrics for consensus assembly
    filepath2 = "${mypath}"+"/"+items[-1]+".vcfcons.fasta"
    with open(filepath2, 'r') as assem:
        header = assem.readline()
        header = header.rstrip()
        bases = assem.readline()
        bases = bases.rstrip()
        num_bases = len(bases)
        ns = bases.count('N')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)
        #print(called)
        #print(end)
        #print(pg)
    #Rename header in fasta to just sample name
    subprocess.run("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2, shell=True, check=True)
    #print("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2)
    
    #QC flag
    pg_flag = ''
    dp_flag = ''
    qc_flag = ''
    pangolin = "v4.1.2_pdata-v1.13"
    vadr_flag = ''
    lineage = ''
    sotc_out = ''
    if float(pg) < 79.5:
        pg_flag = 'FAIL: Percent genome < 80%'
        qc_flag = qc_flag + pg_flag
    else:
        if float(depth) < 100:
            dp_flag = 'FAIL: Mean read depth < 100x'
            qc_flag = qc_flag + dp_flag
        if qc_flag == '':
            qc_flag = qc_flag + 'PASS'
    #print(qc_flag)
    
    if qc_flag == 'PASS':
        subprocess.run("cp "+filepath2+" "+"${params.output}"+"/assemblies/", shell=True, check=True)   
        subprocess.run('cp ' + "${mypath}" + '/' + items[-1] + '.variants_bcftools.vcf ' + "${params.output}"+'/variants/', shell=True, check=True)
    
        #Run VADR
        out_log = open("${mypath}/"+items[-1]+'.out', 'w')
        err_log = open("${mypath}/"+items[-1]+'.err', 'w')
        
        subprocess.run("singularity exec -B "+"${mypath}"+":/data docker://staphb/vadr:1.5.1 /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 /data/" + items[-1]+".vcfcons.fasta > " + "${mypath}"+"/"+items[-1]+".vadrtrimmed.fasta", shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+":/data docker://staphb/vadr:1.5.1 v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ /data/"+items[-1]+".vadrtrimmed.fasta -f /data/"+"vadr_results --noseqnamemax", shell=True, stdout=out_log, stderr=err_log, check=True)
        
        #Parse through VADR outputs to get PASS or REVIEW flag
        
        with open("${mypath}"+"/vadr_results/vadr_results.vadr.pass.list", 'r') as p_list:
            result = p_list.readline()
            result = result.rstrip()
            if result == items[-1]:
                vadr_flag = 'PASS'

        with open("${mypath}"+"/vadr_results/vadr_results.vadr.fail.list", 'r') as f_list:
            result = f_list.readline()
            result = result.rstrip()
            if result == items[-1]:
                vadr_flag = 'REVIEW'

        #Copy VADR error report to main analysis folder for easier review
        if vadr_flag == 'REVIEW':
            subprocess.run("cp " + "${mypath}"+"/vadr_results/vadr_results.vadr.alt.list "+"${params.output}"+"/vadr_error_reports/"+items[-1]+".vadr.alt.list", shell=True, check=True)

        #Run pangolin
        
        
        subprocess.run('singularity exec -B '+"${mypath}"+":/data docker://staphb/pangolin:4.1.2-pdata-1.13" + ' pangolin -o /data /data/'+items[-1]+".vcfcons.fasta", shell=True, check=True)

        #Get lineage
        proc = subprocess.run("tail -n 1 "+"${mypath}"+"/lineage_report.csv | cut -d \',\' -f 2", shell=True, check=True, capture_output=True, text=True)
        lineage = proc.stdout.rstrip()
        print(lineage)

        #Run nextclade
        subprocess.run("singularity exec -B "+"${mypath}"+":/data docker://nextstrain/nextclade:latest nextclade dataset get --name sars-cov-2 --output-dir /data/dataset", shell=True, check=True)
        subprocess.run("singularity exec -B "+"${mypath}"+":/data docker://nextstrain/nextclade:latest nextclade run --input-dataset /data/dataset --output-csv /data/nextclade_report.csv /data/"+items[-1]+".vcfcons.fasta", shell=True, check=True)

        #Parse nextclade output and screen for sotc
        sotc_v = "${params.sotc}".split(',')
        with open("${mypath}"+"/nextclade_report.csv", 'r') as nc:
                header = nc.readline()
                c_results = nc.readline()
                c_results = c_results.rstrip()
                data = c_results.split(',')
                #print(data)
                sotc = []
                for v in sotc_v:
                    print(v)
                    if v in data:
                        sotc.append(v)
                sotc_out = (',').join(sotc)

        out_log.close()
        err_log.close()
        #subprocess.run("rm -r " + "${mypath}"+"/dataset", shell=True, check=True)
    else:
        vadr_flag = 'NA'
        lineage = 'NA'
        sotc_out = 'NA'
        
    with open("${mypath}"+"/report.txt", 'w') as report:
        header = ['sampleID', 'reference', 'start', 'end', 'num_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag', 'pangolin_version', 'lineage', 'SOTC']
        report.write('\t'.join(map(str,header)) + '\n')
        results = [items[-1], ref_name, start, end, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag, pangolin, lineage, sotc_out]
        report.write('\t'.join(map(str,results)) + '\n')
    /$
}