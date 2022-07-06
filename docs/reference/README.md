---
title: "Reference"
sidebar: toc
maxdepth: 1
sort: 4
---

# References

**2009** bwa, update in 2013, and again in 2019 to bwa-mem2

* Li, H. and Durbin, R., 2009. Fast and accurate short read alignment with Burrows–Wheeler transform. bioinformatics, 25(14), pp.1754-1760.
* Li, H., 2013. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997.
* Vasimuddin, M., Misra, S., Li, H. and Aluru, S., 2019, May. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS) (pp. 314-324). IEEE.

**2012** FreeBayes

* Garrison, E. and Marth, G., **2012**. [Haplotype-based variant detection from short-read sequencing](https://api.semanticscholar.org/CorpusID:15153602). arXiv preprint arXiv:1207.3907.

**2013** FALCON, FALCON-unzip, FALCON-Phase

* Chin, C.S., Alexander, D.H., Marks, P., Klammer, A.A., Drake, J., Heiner, C., Clum, A., Copeland, A., Huddleston, J., Eichler, E.E. and Turner, S.W., **2013**. [Nonhybrid, finished microbial genome assemblies from long-read SMRT sequencing data](https://api.semanticscholar.org/CorpusID:205421576). Nature methods, 10(6), pp.563-569.
* Chin, C.S., Peluso, P., Sedlazeck, F.J., Nattestad, M., Concepcion, G.T., Clum, A., Dunn, C., O'Malley, R., Figueroa-Balderas, R., Morales-Cruz, A. and Cramer, G.R., **2016**. [Phased diploid genome assembly with single-molecule real-time sequencing](https://api.semanticscholar.org/CorpusID:11465102). Nature methods, 13(12), pp.1050-1054.
  * p-contigs=primary contigs; a-contigs=alternative contigs (local alt seqs); falcon-unzip (haplotype phasing) generates p-contigs, and final haplotig set (h-contigs which is more specific than a-contigs). 
* Kronenberg, Z.N., **Rhie, A.**, Koren, S., Concepcion, G.T., Peluso, P., Munson, K.M., Porubsky, D., Kuhn, K., Mueller, K.A., Low, W.Y. and Hiendleder, S., **2021**. [Extended haplotype-phasing of long-read de novo genome assemblies using Hi-C. Nature communications](https://www.nature.com/articles/s41467-020-20536-y), 12(1), pp.1-10.
  * "Thus, we suggest the following genome assembly workflow: (1) partially phased long-read assembly, (2) FALCON-Phase on primary contigs and haplotigs, (3) scaffolding with HI-C data, and (3) FALCON-Phase on scaffolds.

**2014** GRAAL, instaGRAAL (update in 2020, utilizes GPUs)

* Marie-Nelly, H., Marbouty, M., Cournac, A., Flot, J.F., Liti, G., Parodi, D.P., Syan, S., Guillén, N., Margeot, A., Zimmer, C. and Koszul, R., 2014. [High-quality genome (re) assembly using chromosomal contact data](https://api.semanticscholar.org/CorpusID:457083). Nature communications, 5(1), pp.1-10.
* Baudry, L., Guiglielmoni, N., Marie-Nelly, H., Cormier, A., Marbouty, M., Avia, K., Mie, Y.L., Godfroy, O., Sterck, L., Cock, J.M. and Zimmer, C., 2020. [instaGRAAL: chromosome-level quality scaffolding of genomes using a proximity ligation-based scaffolder](https://api.semanticscholar.org/CorpusID:219730049). Genome biology, 21(1), pp.1-22.

**2015** Longranger, BUSCO (updates in 2021)

* Bishara, A., Liu, Y., Weng, Z., Kashef-Haghighi, D., Newburger, D.E., West, R., Sidow, A. and Batzoglou, S., 2015. [Read clouds uncover variation in complex regions of the human genome](https://genome.cshlp.org/content/25/10/1570.short). Genome research, 25(10), pp.1570-1580.
* Manni, M., Berkeley, M.R., Seppey, M., Simao, F.A. and Zdobnov, E.M., 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. arXiv preprint arXiv:2106.11799.

**2016** minimap2, gEVAL

* Li, H., 2016. [Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences](https://api.semanticscholar.org/CorpusID:3618555). Bioinformatics, 32(14), pp.2103-2110.
* Li, H., 2018. [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778). Bioinformatics, 34(18), pp.3094-3100.
* Chow, W., Brugger, K., Caccamo, M., Sealy, I., Torrance, J. and Howe, K., 2016. [gEVAL—a web-based browser for evaluating genome assemblies](https://api.semanticscholar.org/CorpusID:18844893). Bioinformatics, 32(16), pp.2508-2510.
* Howe, K., Chow, W., Collins, J., Pelan, S., Pointon, D.L., Sims, Y., Torrance, J., Tracey, A. and Wood, J., 2021. [Significantly improving the quality of genome assemblies through curation](https://api.semanticscholar.org/CorpusID:231304018). Gigascience, 10(1), p.giaa153.
  * gEVAL is a browser based method for evaluating quality of genome assemblies
  * "This is especially timely in the context of emerging projects that aim to assemble the genomes of very large numbers of species to highest quality possible, including the Vertebrate Genomes Project (VGP), the Darwin Tree of Life Project (DToL, darwintreeoflife.org), and the overarching Earth Biogenome Project [1, 10]."
  * "Before being loaded into gEVAL, all assemblies are run through a nextflow [38] pipeline that performs contamination detection and separation or removal as described in Table 1, combined with removal of trailing Ns [38]." 

**2016** Juicer, Juicebox

* Durand, N.C., Shamim, M.S., Machol, I., Rao, S.S., Huntley, M.H., Lander, E.S. and Aiden, E.L., 2016. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. Cell systems, 3(1), pp.95-98.
* Durand, N.C., Robinson, J.T., Shamim, M.S., Machol, I., Mesirov, J.P., Lander, E.S. and Aiden, E.L., 2016. Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom. Cell systems, 3(1), pp.99-101.

**2017** SALSA, Canu (TrioCanu)

* Ghurye, J., Pop, M., Koren, S., Bickhart, D. and Chin, C.S., 2017. [Scaffolding of long read assemblies using long range contact information](https://api.semanticscholar.org/CorpusID:4029050). BMC genomics, 18(1), pp.1-11.
* Koren, S., Walenz, B.P., Berlin, K., Miller, J.R., Bergman, N.H. and Phillippy, A.M., 2017. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome research, 27(5), pp.722-736.
* Koren, S., Rhie, A., Walenz, B.P., Dilthey, A.T., Bickhart, D.M., Kingan, S.B., Hiendleder, S., Williams, J.L., Smith, T.P. and Phillippy, A.M., 2018. De novo assembly of haplotype-resolved genomes with trio binning. Nature biotechnology, 36(12), pp.1174-1182.

**2018** purge_haplotigs, purge_dups

* Roach, M.J., Schmidt, S.A. and Borneman, A.R., 2018. [Purge Haplotigs: allelic contig reassignment for third-gen diploid genome assemblies.](https://api.semanticscholar.org/CorpusID:54039252) BMC bioinformatics, 19(1), pp.1-10.
  * purge_haplotigs
* Guan, D., McCarthy, S.A., Wood, J., Howe, K., Wang, Y. and Durbin, R., 2020. [Identifying and removing haplotypic duplication in primary genome assemblies](https://api.semanticscholar.org/CorpusID:202030660). Bioinformatics, 36(9), pp.2896-2898.
  * C source code at https://github.com/dfguan/purge_dups
  * Pipeline outline: (1) minimap2 (li, 2016), (2) create windows by contigs and self align, (3) remove haplotigs, (4) chain overlaps.. something about the shorter contig. (more detail in [Supplementary Material](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/36/9/10.1093_bioinformatics_btaa025/1/btaa025_supplementary_data.pdf?Expires=1632508562&Signature=cDc1UMOPwGVGmfhQstXvuQ-QEWbFPyi9kNOxCRmH4uYEcBDYk3uuDXGQNcLo~uXA75MF5V9o4S-Rx2XYmBm8dLtHb0HEZ1l9VkfcsHTw-V8xRMwFEoLabq-P0kmRoLPGVoxXVpkZ4V9FJ76KGJ4N9DV0TtYJ7nz1To4FuTGK2nPXWDRWcM-7NegCCt9JfJmCZE-r8oDYdr07ogMwhc-Qei25ornGvAHs2U7rb4guiWLOsOgJ1XNdx3w3Fapn2NOlixxBDGnNvmUaP~p1Y2xyncKg-6M4TaNq0LlIxEaNnjEsZdm8JTi4Bio22mAjZUFPH30E-9w955YuJ4aaudVDow__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)).
  * "Following this [Scaff10x] with a round of polishing with Arrow closed a number of gaps, reducing contig number further and increasing contig N50" Wait... arrow merges contigs? or maybe it's Scaff10x.
  * "To our knowledge, scaffolders that use long-range information, such as Scaff10X with linked reads or SALSA with Hi-C data, do not handle heterozygous overlaps. We therefore recommend applying purge_dups directly after initial assembly, prior to scaffolding."
  * "In conclusion, purge_dups can significantly improve genome assemblies by removing overlaps and haplotigs caused by sequence divergence in heterozygous regions." ... removes false dups, while retaining assembly completeness, improves scaffolding
  * Supplemental
  
  ```
  # === input/output variables
  pfs=*.pfs                # raw Pacbio read alignment PAF files
  asm=all_p_ctg.fasta      # primary assembly..um do I include mito and haplo here?
  
  # === Purge dups commands
  pbcstat $pfs       # will generate PB.base.cov and PB.stat
  calcuts PB.stat > cutoffs 2> calcults.log
  split_fa $asm > $asm.split.fa
  minimap2 -xasm5 -DP $asm.split.fa $asm.split.fa > $asm.split.self.paf
  purge_dups -2 -T cutoffs -c PB.base.cov $asm.split.self.paf > dups.bed 2> purge_dups.log
  get_seqs dups.bed $asm > purged.fa 2> hap.fa        # so it separates here..haplotigs sent to stderr?
  ```


**2020** Merqury

* **Rhie, A.**, Walenz, B.P., Koren, S. and Phillippy, A.M., **2020**. [Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies](https://api.semanticscholar.org/CorpusID:214727187). Genome biology, 21(1), pp.1-27.
  * 
  > For a comparison of multiple assemblies, we demonstrate Merqury on haplotype-resolved (TrioCanu [10]), pseudo-haplotype (FALCON-Unzip [9]), and mixed-haplotype (Canu [33]) assemblies of this hybrid genome.

**2021** merfin, mitoVGP, VGP assembly pipeline

* **Formenti, G.**, **Rhie, A.**, Walenz, B.P., Thibaud-Nissen, F., Shafin, K., Koren, S., Myers, E.W., Jarvis, E.D. and Phillippy, A.M., 2021. [Merfin: improved variant filtering and polishing via k-mer validation](https://api.semanticscholar.org/CorpusID:236145752). bioRxiv.
* **Formenti, G.**, **Rhie, A.**, Balacco, J., Haase, B., Mountcastle, J., Fedrigo, O., Brown, S., Capodiferro, M.R., Al-Ajli, F.O., Ambrosini, R. and Houde, P., 2021. [Complete vertebrate mitogenomes reveal widespread repeats and gene duplications](https://api.semanticscholar.org/CorpusID:233433542). Genome biology, 22(1), pp.1-22.
* **Rhie, A.**, McCarthy, S.A., Fedrigo, O., Damas, J., **Formenti, G.**, Koren, S., Uliano-Silva, M., Chow, W., Fungtammasan, A., Kim, J. and Lee, C., 2021. [Towards complete and error-free genome assemblies of all vertebrate species](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8081667/). Nature, 592(7856), pp.737-746.
  * "Genome heterozygosity posed additional problems, because homologous haplotypes in a diploid or polyploid genome are forced together into a single consensus by standard assemblers, sometimes creating false gene duplications."
  * Website: https://vertebrategenomesproject.org
  * "To our knowledge, this was the first systematic analysis of many sequence technologies, assembly algorithms, and assembly parameters applied on the same individual" heh, that would be fun
  * "After fixing a function in the PacBio FALCON software that caused artificial breaks in contigs between stretches of highly homozygous and heterozygous haplotype sequences (Supplementary Note 1, Table 2), ..." did we fix this as well?
  * VGP assembly pipeline (v1.0): haplotype-separated CLR contigs, scaffolding with linked reads, optical maps and Hi-C, gap filling, base call polishing, manual curation (extended data **Figs 2a** (polishing after scaffolding), 3a).
  * VGP assembly flowchart (Extended Data Fig 3): purge dups -> scaffold -> polish {arrow, longranger+FreeBayes, longranger+FreeBayes} "with binned reads" means reads by contig?

<details><summary>Expandable notes</summary>

  * 
  > FALCON and FALCON-Unzip were run with default parameters, except for computing the overlaps. Raw read overlaps were computed with DALIGNER parameters -k14 -e0.75 -s100 -l2500 -h240 -w8 to better reflect the higher error rate in early PacBio sequel I and II. Pread (preassembled read) overlaps were computed with DALIGNER parameters -k24 -e.90 -s100 -l1000 -h600 intending to collapse haplotypes for the FALCON step to better unzip genomes with high heterozygosity rate. FALCON-Unzip outputs both a pseudo-haplotype and a set of alternate haplotigs that represent the secondary alleles. We refer to these outputs as the primary contig set (c1) and alternate contig set (c2).
  * 
  >  To reduce these false duplications, we ran Purge_Haplotigs13, first during curation (VGP v1.0 pipeline) and then later after contig formation (VGP v1.5 pipeline). To do the former, Purge_Haplotigs was run on the primary contigs (c1), and identified haplotigs were mapped to the scaffolded primary assembly with MashMap286 for removal. In the latter, identified haplotigs were moved from the primary contigs (c1) to the alternate haplotig set (p2). The remaining primary contigs were referred to as p1; p2 combined with c2 was referred to as q2. Later, in the VGP v1.6 pipeline, we replaced Purge_Haplotigs with Purge_Dups14, a new program developed by several of the authors in response to Purge_Haplotigs not removing partial false duplication at contig boundaries. Purging also removes excessive low-coverage (junk) and high-coverage (repeats) contigs. To calculate the presence and overall success of purging false duplications, we used a k-mer approach (Supplementary Methods, Supplementary Fig. 6).
  *
  > To polish bases in both haplotypes with minimal alignment bias, we concatenated the alternate haplotig set (c2 in v1.0 or q2 in v1.5–1.6) to the scaffolded primary set (s3) and the assembled mitochondrial genome (mitoVGP in v1.6). We then performed another round of polishing with Arrow (smrtanalysis 5.1.0.26412) using PacBio CLR reads, aligning with pbalign --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 –concordant --algorithm=blasr --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1 and consensus polishing with variantCaller --skipUnrecognizedContigs haploid -x 5 -q 20 -X120 –v --algorithm=arrow. While this round of polishing resulted in higher QV for all genomes herein considered, we noticed that it was particularly sensitive to the coverage cutoff parameter (-x). This is because Arrow generates a de novo consensus from the mapped reads without explicitly considering the reference sequence. Later, we found that the second round of Arrow polishing sometimes reduced the QV accuracy for some species. Upon investigation, this issue was traced back to option -x 5, which requires at least 5 reads to call consensus. Such low minimum requirements can lead to uneven polishing in low coverage regions. To avoid this behaviour, we suggest to increase the -x close to the half sequence coverage (for example, 30× when 60× was used for assembly) and check QV before moving forward.

</details>

**2021** ag100pest update

* **Childers, A.K.**, **Geib, S.M.**, **Sim, S.B.**, Poelchau, M.F., **Coates, B.S.**, Simmonds, T.J., **Scully, E.D.**, Smith, T.P., Childers, C.P., Corpuz, R.L. and Hackett, K., 2021. [The USDA-ARS Ag100Pest Initiative: High-Quality Genome Assemblies for Agricultural Pest Arthropod Research](https://api.semanticscholar.org/CorpusID:235896348). Insects, 12(7), p.626.
  * Figure 1: general workflow
  * Bioproject: https://www.ncbi.nlm.nih.gov/bioproject/555319
  * "Ag100Pest began by using continuous long reads (CLRs) for assembly (details not presented herein) as the improved HiFi procedure [33] had not yet been developed"


## Online Videos

* [2018 Nov 8 - Best Practices for Rapid Reference Quality Genome Assembly - Webinar](https://youtu.be/HTJnhuFpM30)
