---
title: "FAQ"
sidebar: toc
maxdepth: 1
sort: 4
---

### How to I check that software dependencies are available?

Use the `--check_software` flag to print out paths to software dependencies and their versions:

```
nextflow run isugifNF/polishCLR -r main \
  --check_software \
  -profile local
```

Which should print out something like:

```
N E X T F L O W  ~  version 22.04.5
Pulling isugifNF/polishCLR ...
 downloaded from https://github.com/isugifNF/polishCLR.git
Launching `https://github.com/isugifNF/polishCLR` [amazing_lamarck] DSL2 - revision: e0325a97a5 [main]
executor >  local (1)
[38/26cb53] process > check_software [100%] 1 of 1 ✔
===== Dependencies check =====
parallel   .... good.  GNU parallel 20160622
bzcat      .... good.  bzip2, a block-sorting file compressor. Version 1.0.8, 13-Jul-2019.
pigz       .... good.  pigz 2.6
meryl      .... good.  meryl 1.3
pbmm2      .... good.  pbmm2 1.9.0 Using: pbmm2 : 1.9.0 (commit v1.9.0) pbbam : 2.1.0 (commit v2.1.0) pbcopper : 2.0.0 (commit v2.0.0-54-g1ce9870) boost : 1.77 htslib : 1.15 minimap2 : 2.15 zlib : 1.2.11
minimap2   .... good.  2.24-r1122
samtools   .... good.  samtools 1.15.1
gcpp       .... good.  gcpp 2.0.2-2.0.2
bwa-mem2   .... good.  bwa-mem v
freebayes  .... good.  version: v1.3.6
bcftools   .... good.  bcftools 1.15.1
merfin     .... good.  merfin 1.0
pbcstat    .... good.  Usage: aa_pb [options] <PAF_FILE> ...
calcuts    .... good.  Usage: calcuts [<options>] <STAT> ...
split_fa   .... good.  Usage: split_fa [<options>] <STAT> ...
purge_dups .... good.  Version: 1.2.5
get_seqs   .... good.  Usage: get_seqs [<options>] <DUPs.BED> <FASTA>
gzip       .... good.  gzip 1.12
busco      .... good.  BUSCO 5.4.2
```

or some combination of `software ... need to install.`


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