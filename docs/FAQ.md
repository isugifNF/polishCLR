---
title: "FAQ"
sidebar: toc
maxdepth: 1
sort: 1
---

### What does "[Pipeline error] Parameter(s) [...] is(are) not valid in the pipeline!" mean?

If you receive an error similar to:

```
N E X T F L O W  ~  version 22.04.5
Launching `https://github.com/isugifNF/polishCLR` [focused_sax] DSL2 - revision: 98b864ca3a [check_valid_params]
[Pipeline error] Parameter(s) [other, mito] is(are) not valid in the pipeline!
```

You may have mispelled one of the parameters. By default, Nextflow quietly ignores any misspelled or extra parameters. Since polishCLR sometimes requires quite a number of files/parameters, it was necessary to add an “unrecognized parameter” catch-all kind of error.

Right now the valid parameters are hardcoded in a `parameters_valid` set within main.nf:

```
 def parameters_valid = ['help','monochrome_logs','outdir', 
 'primary_assembly','alternate_assembly','paternal_assembly','maternal_assembly','mitochondrial_assembly', 
 'illumina_reads','pacbio_reads','k','species','meryldb','falcon_unzip','same_specimen','steptwo','step', 
 'queueSize','account','threads','clusterOptions', 
 'queue-size', 'cluster-options', 
 'parallel_app','bzcat_app','pigz_app','meryl_app','merqury_sh','pbmm2_app','minimap2_app','samtools_app', 
 'gcpp_app','bwamem2_app','freebayes_app','bcftools_app','merfin_app','pbcstat_app','hist_plot_py', 
 'calcuts_app','split_fa_app','purge_dups_app','get_seqs_app','gzip_app','busco_app','busco_lineage', 
 'parallel_params','pbmm2_params','minimap2_params','gcpp_params','bwamem2_params','freebayes_params', 
 'purge_dups_params','busco_params','check_software','profile'] as Set 
```

You can also trigger this error by using empty input and extra flags such as `--other` and `--mito`. 

```
mkdir data 
cd data
touch pri.fa alt.fa mit.fa ill_R1.fastq.bz2 ill_R2.fastq.bz2 pac.subreads.bam
cd ..

nextflow run isugifNF/polishCLR \
  --primary_assembly data/pri.fa \
  --alternate_assembly data/alt.fa \
  --mitochondrial_assembly data/mit.fa \
  --illumina_reads "data/ill_{R1,R2}.fastq.bz2" \
  --pacbio_reads data/pac.subreads.bam \
  --species "BugName" \
  --k "21" \
  --falcon_unzip true \
  --step 1 \
  -stub-run \
  -profile local \
  --other hi --mito "doesn't_exist"
```