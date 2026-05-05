# HAWEI
Haplotype-Resolved Alignment-Based Window Expression Inference (HAWEI)

**IMPORTANT!!** **HAWEI is currently not a fully developed software package. It was developed primarily for personal investigation**.  Others are very welcome to use, inspect, adapt, or extend HAWEI if the workflow is useful for related problems.

## What is HAWEI?
HAWEI is a haplotype-aware RNA-seq analysis workflow designed for cases where a target gene is represented by multiple highly similar haplotypes and conventional transcript quantification can lose resolution because of extensive multi-mapping.

HAWEI was developed around gene-level targeted analyses rather than broad transcriptome-wide deployment.
Users should inspect intermediate outputs and verify whether the selected panels are biologically and technically sensible.
If you use HAWEI in your own analyses, independent validation is strongly recommended.

## Repository status
HAWEI is currently maintained as a **working research repository**.  
It is most appropriate for users who are comfortable reading scripts, checking intermediate files, and modifying commands when needed.

## Basic usage
```bash
bash HAWEI.sh -i <aligned_haplotype_fasta> -o <route_dir> --fq-manifest <samples.tsv>

-i, --input
An aligned full-length haplotype FASTA. All haplotype sequences must be aligned on the same coordinate system.

-o, --outdir
Output directory for the analysis route.

--fq or --fq-manifest
RNA-seq input reads, provided either as:
a directory of paired-end FASTQ files, or
a manifest table with columns: sample, r1, r2

## Example
bash HAWEI.sh \
  -i Elovl5b_ref_align.ABCD.fasta \
  -o route_elovl5b_VO_liver \
  --fq-manifest manifests/VO_liver.tsv \
  --nt 8
Main options
bash HAWEI.sh -h

Commonly used options include:
--len
Window lengths to build, default: 300,400,800,1200
--nt
Number of threads
--stop-after
Stop after a specific stage: build, cws, select, or ecboot
--panels
Restrict downstream analysis to selected panels
--samples
Restrict downstream analysis to selected samples
--min-ungapped
Minimum ungapped length required to retain a haplotype window

