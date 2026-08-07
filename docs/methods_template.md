# Methods — NGS Data Processing and Analysis

*Template for publication-grade methods sections. Adapt sections based on modality (RNAseq, WES, WGS). Parameters reflect pipeline defaults; update values if modified for your study.*

---

## Data Processing and Bioinformatics

All raw sequencing data were processed using a custom Snakemake workflow (Mölder et al., 2021) implementing GATK Best Practices (Van der Auwera et al., 2013) with modality-specific adaptations for RNA-seq, whole-exome sequencing (WES), and whole-genome sequencing (WGS). The human reference genome GRCh38 primary assembly (Schneider et al., 2017) was used throughout. Gene annotation was obtained from GENCODE v48 (Frankish et al., 2021).

### Quality Control and Read Preprocessing

Raw paired-end (or single-end) FASTQ files were assessed for sequencing quality using FastQC v0.12.1 (Andrews, 2010) before and after adapter removal. Adapter sequences were detected and trimmed, and low-quality bases were removed using fastp v1.0.1 (Chen et al., 2018) with a minimum Phred quality score of 20 and a minimum read length of 50 bp. For paired-end libraries, automatic adapter detection was enabled (`--detect_adapter_for_pe`). Aggregate quality metrics across all samples were summarized using MultiQC v1.20 (Ewels et al., 2016).

### Sequence Alignment

**[Use for RNAseq]** Trimmed reads were aligned to the GRCh38 reference genome using the splice-aware aligner HISAT2 v2.2.1 (Kim et al., 2019) with a GENCODE v48-derived splice site-aware index. Alignments were performed with the `--dta` flag enabled to facilitate downstream transcript quantification. Full read group information (sample identity, library, platform, and platform unit) was embedded in each alignment file.

**[Use for WES/WGS]** Trimmed reads were aligned to the GRCh38 reference genome using bwa-mem2 (Vasimuddin et al., 2019), an AVX-accelerated implementation of the BWA-MEM algorithm (Li, 2013). The `-Y` flag was used to enable soft-clipping of supplementary alignments. Full read group information (sample identity, library, sequencing platform, and platform unit) was embedded in each alignment file.

Resulting alignments were coordinate-sorted using SAMtools v1.21 (Li et al., 2009). PCR and optical duplicates were flagged (but not removed) using GATK v4.6.2.0 MarkDuplicates (McKenna et al., 2010) with `VALIDATION_STRINGENCY=LENIENT`. Alignment quality metrics, including total reads, mapping rate, and properly paired read fraction, were computed using SAMtools flagstat.

### Base Quality Score Recalibration (DNA only)

**[Use for WES/WGS only]** Base quality scores were recalibrated following GATK Best Practices using BaseRecalibrator and ApplyBQSR from GATK v4.6.2.0. Known variant sites used for recalibration included dbSNP v138 (Sherry et al., 2001), the Mills and 1000 Genomes gold standard indel set, and the Broad Institute known indels resource (all lifted to GRCh38 coordinates). For WES samples, recalibration was restricted to capture kit target regions.

### Variant Calling (DNA only)

**[Use for WES/WGS only]** Germline variant calling followed the GATK joint genotyping workflow. Per-sample genomic variant call files (gVCFs) were generated using GATK HaplotypeCaller v4.6.2.0 in `-ERC GVCF` mode, which emits reference confidence blocks alongside variant calls. For WES samples, calling was restricted to capture kit intervals. Per-sample gVCFs were consolidated into a GenomicsDB workspace using GenomicsDBImport, and joint genotyping across the cohort was performed using GenotypeGVCFs, producing a multi-sample VCF containing SNPs and indels.

### Gene Expression Quantification (RNA only)

**[Use for RNAseq only]** Gene-level read counts were quantified from deduplicated BAM files using featureCounts v2.1.1 (Liao et al., 2014) from the Subread package. Reads were assigned to exonic features (`-t exon`) and summarized at the gene level (`-g gene_id`) using GENCODE v48 annotation. Only primary alignments were counted (`--primary`). For paired-end libraries, fragments (not individual reads) were counted using `--countReadPairs`. Strandedness was set to unstranded (`-s 0`) [*adjust to `-s 1` (forward) or `-s 2` (reverse) as appropriate for library preparation protocol*]. Per-sample count files were merged into a single gene-by-sample count matrix for downstream differential expression analysis.

### Genetic Ancestry Estimation

Genetic ancestry was estimated for all samples using a containerized pipeline incorporating GATK 3.8, PLINK v1.9 (Purcell et al., 2007), and ADMIXTURE v1.3.0 (Alexander et al., 2009). Deduplicated BAM files were processed through genotype likelihood estimation at approximately 10,000 ancestry-informative marker (AIM) sites, followed by merging with the 1000 Genomes Project Phase 3 reference panel (The 1000 Genomes Project Consortium, 2015), which comprises 2,504 individuals from 26 populations across five continental super-populations (AFR, AMR, EAS, EUR, SAS). Supervised ADMIXTURE analysis was performed at K=5 (continental super-populations) and optionally at K=23 (sub-populations). Principal component analysis (PCA) was used to project study samples onto the reference coordinate space. Ancestry proportions and PCA coordinates were visualized using R v4.4.2 with ggplot2.

---

## References

1. Alexander DH, Novembre J, Lange K. Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*. 2009;19(9):1655–1664. doi:10.1101/gr.094052.109

2. Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. Available online at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

3. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*. 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560

4. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016;32(19):3047–3048. doi:10.1093/bioinformatics/btw354

5. Frankish A, Diekhans M, Jungreis I, et al. GENCODE 2021. *Nucleic Acids Research*. 2021;49(D1):D916–D923. doi:10.1093/nar/gkaa1087

6. Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*. 2019;37(8):907–915. doi:10.1038/s41587-019-0201-4

7. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*. 2013;1303.3997. doi:10.48550/arXiv.1303.3997

8. Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 2009;25(16):2078–2079. doi:10.1093/bioinformatics/btp352

9. Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*. 2014;30(7):923–930. doi:10.1093/bioinformatics/btt656

10. McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*. 2010;20(9):1297–1303. doi:10.1101/gr.107524.110

11. Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33. doi:10.12688/f1000research.29032.2

12. Purcell S, Neale B, Todd-Brown K, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics*. 2007;81(3):559–575. doi:10.1086/519795

13. Schneider VA, Graves-Lindsay T, Howe K, et al. Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly. *Genome Research*. 2017;27(5):849–864. doi:10.1101/gr.213611.116

14. Sherry ST, Ward MH, Kholodov M, et al. dbSNP: the NCBI database of genetic variation. *Nucleic Acids Research*. 2001;29(1):308–311. doi:10.1093/nar/29.1.308

15. The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015;526(7571):68–74. doi:10.1038/nature15393

16. Van der Auwera GA, Carneiro MO, Hartl C, et al. From FastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Current Protocols in Bioinformatics*. 2013;43:11.10.1–11.10.33. doi:10.1002/0471250953.bi1110s43

17. Vasimuddin Md, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In: *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*. 2019:314–324. doi:10.1109/IPDPS.2019.00041
