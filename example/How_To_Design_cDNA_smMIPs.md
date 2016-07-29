Overview
========

This document describes how to design cDNA-smMIPs for specific set of genes
using the *design_mips_for_all_transcripts.py* script.

The script performs the following steps:

1. Load Ensembl gene model and transcript sequences
2. Exhaustively generate all possible cDNA-smMIPs against all transcripts for the target gene(s)
3. Determine copy number of the individual probes and compute MIPS predicted performance scores based on the MIPGEN 2.0 logistic model by Boyle et al. See Note below on restrictions for commercial usage.
4. Generate Fastq file with probe and target sequence for cDNA-smMIPs with performance score above threshold.

The user can then map the fastq file to the reference genome with a splicing-aware mapper (e.g. STAR). We then recommend the use of a genome browser such as [IGV](https://www.broadinstitute.org/igv/) to select the cDNA-smMIPs that cover the target transcripts as desired.

Requirements
============

1. To compute smMIP performance scores (based on the MIPGEN 2.0 model by Boyle et al, Bioinformatics 2014), it is necessary
to have the read mapper [bwa](http://bio-bwa.sourceforge.net/) installed. Using the *--bwa* option the path to the bwa binary can be provided.
In this case also a whole-genome reference sequence fasta indexed with SAMtools and bwa must be provided.
2. The Python script relies on the [pysam](http://pysam.readthedocs.io/en/latest/installation.html) package for access to VCF files. 
3. The Python packages *scipy*, *numpy*, *pandas*, *htseq*, and *pybedtools* must be installed. These packages can all be installed with *pip* ('pip install <package name>').
4. A C++ compiler (such as *gcc* on Linux) as the script uses scipy.weave to compile inline C++ code. 
5. SAMTools and tabix can be found [here](http://htslib.org).
6. To avoid issues with gene names mapping to multiple Ensembl identifiers, the script requires the user to specify the Ensembl gene identifiers for the target genes. 

Download Ensembl cDNA and gene description files
================================================

To design cDNA-smMIPs for any gene, you first need to download the Ensembl files with the cDNA sequence for all transcripts, and the file with the gene model (transcripts, exons, order of exons and genomic coordinates): 

```
 wget ftp://ftp.ensembl.org/pub/grch37/release-84/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
 wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
```

Provide bgzipped VCF file indexed with Tabix with variants to be excluded
=========================================================================

There are two VCF files that can be specified:

1. A VCF file which contains variants that should be excluded from the extension and ligation probes of the smMIPs. This is usually a VCF file with common variants (SNPs and indels) from 1000 Genomes or Exac.
2. A VCF file which contains variants that should be targeted by the smMIP, for instance, for allele-specific expression analysis.

To bgzip a VCF file *variants.vcf* and index it with tabix, do
```
    bgzip variants.vcf
    tabix -p vcf variants.vcf.gz
```
The second step creates a file called *variants.vcf.gz.tbi*; this is the index file that allows random access.

For this example, the VCF file *ExAC.r0.3.sites.vep.vcf.gz.AC_1000_SNPs.vcf.gz.PNN.vcf.gz* is provided, which contains common SNPs only for the PNN gene. In the example we use it as both the source of variants that should not be contained in the probes (*--exclude_variants_vcf*) and the source of target variants that should be contained in the target sequence of cDNA-smMIPs (*--target_variants_vcf*).

Optional and required arguments
===============================

```
    required named arguments:
      --ensembl_gff ENSEMBL_GFF
                            Gff3 file with Ensembl gene, transcript and exon
                            definitions. May be gzipped (default: None)
      --ensembl_cDNA ENSEMBL_CDNA
                            Fasta file with cDNA of all Ensembl transcripts. May
                            be gzipped (default: None)
      --output_file_prefix OUTPUT_FILE_PREFIX
                            Prefix of names output files. Output files will be
                            gzipped. (default: None)
      --target_ensembl_file TARGET_ENSEMBL_FILE
                            File with Ensembl gene IDs (ENSG..) of genes for which
                            cDNA-smMIPs should be designed (default: None)
      --length_UMI_extension LENGTH_UMI_EXTENSION
                            Number of Ns (length of UMI) on extension probe.
                            Choose value between 0 and 10. Total length of UMI
                            sequence should not exceed 10 nt. (default: None)
      --length_UMI_ligation LENGTH_UMI_LIGATION
                            Number of Ns (length of UMI) on ligation probe. Choose
                            value between 0 and 10. Total length of UMI sequence
                            should not exceed 10 nt. (default: None)

    optional arguments:
      -h, --help            show this help message and exit
      --remove_transcript_biotypes REMOVE_TRANSCRIPT_BIOTYPES
                            Comma-separated list of biotypes. No smMIPs will be
                            designed for these transcripts (default: unprocessed_p
                            seudogene,retained_intron,nonsense_mediated_decay)
      --exclude_variants_vcf EXCLUDE_VARIANTS_VCF
                            VCF file (bgzip'ed and indexed with tabix) with
                            variants that should not be present in the extension
                            or ligation probe (default: None)
      --target_variants_vcf TARGET_VARIANTS_VCF
                            VCF file (bgzip'ed and indexed with tabix) with
                            variants that should be covered by the target sequence
                            of a cDNA-smMIP (default: None)
      --target_seq_length TARGET_SEQ_LENGTH
                            Length of cDNA-smMIPs target sequence (sequence
                            between extension and ligation probe) (default: 100)
      --calculate_performance_scores
                            Compute performance scores for MIPS using MIPGEN 2.0
                            logistic model (copyright Evan Boyle, used by
                            permission.) Commercial use of this option is not
                            allowed. (default: False)
      --ref_genome REF_GENOME
                            Fasta reference sequence (indexed by Samtools and
                            BWA). Required for calculating MIP performance score.
                            (default: None)
      --bwa BWA             Path to BWA binary. (default: bwa)
      --output_candidate_mips_to_fastq
                            Outputs MIP probe and target sequence to Fastq file
                            for mapping with a splicing-aware mapper. Can be used
                            to select cDNA-smMIPs. (default: False)
      --fastq_minimum_performance_score FASTQ_MINIMUM_PERFORMANCE_SCORE
                            If outputting candidate cDNA-smMIPs to Fastq, only
                            output MIPs with performance score above this
                            threshold. Only effective if
                            calculate_performance_scores is also specified
                            (default: 0.9)

```

Unique Molecular Identifiers (UMIs)
===================================

The script has two options to specify the length of the unique molecular identifier (UMI). 
The UMI can be put in two positions: 1) between the sequence primer and the extension probe and 2) between the sequence primer and ligation probe.
It is possible to have both. 
The length of the extension and ligation UMI can be specified with the options *--length_UMI_extension* and *--length_UMI_ligation*. 
By default, the extension UMI is of length 9 and the ligation UMI is of length 0. For all smMIPs in the manuscript these are the value that were used.
It is not recommended to specify UMI sequences that together are longer than 10 nt.

**NOTE** The values for *--length_UMI_extension* and *--length_UMI_ligation* must be provided to the script that counts molecules for each smMIP 
as described in the tutorial *How_To_Estimate_Relative_and_Differential_expression*.


Example
=======

In this example we will design cDNA-smMIPs for the PNN gene, one of the genes included in the ASE experiment in the manuscript.
We have created a gff3 file with records for just this gene, Homo_sapiens.GRCh37.82.gff3.gz.PNN.gz. However, it also works with the full *Homo_sapiens.GRCh37.82.gff3.gz* file, only slower.

First, we have to create a file called *ensembl_genes.txt* which contains all the Ensembl identifiers of the genes for which we would like to design cDNA-smMIPs. In this case, the file contains the Ensembl ID for the *PNN* gene:
 
    ENSG00000100941

Next, design cDNA-smMIPs for this gene with the following command:

```bash
 python design_mips_for_all_transcripts.py --ensembl_gff files/Homo_sapiens.GRCh37.82.gff3.gz.PNN.gff3.gz --ensembl_cDNA Homo_sapiens.GRCh37.cdna.all.fa.gz \
         --output_file_prefix PNN.mips --target_ensembl_file files/ensembl_genes.txt \
         --length_UMI_extension 9 --length_UMI_ligation 0 \
         --exclude_variants_vcf files/ExAC.r0.3.sites.vep.vcf.gz.AC_1000_SNPs.vcf.gz.PNN.vcf.gz \
         --target_variants_vcf files/ExAC.r0.3.sites.vep.vcf.gz.AC_1000_SNPs.vcf.gz.PNN.vcf.gz \
         --calculate_performance_scores --output_candidate_mips_to_fastq --ref_genome /data1/genomes/hg19.fa
```
Here, the prefix of all the output files is specified with the required parameter *--output_file-prefix*.

*NOTE* The file '/data1/genomes/hg19.fa' should be replaced by the one on your system. The fie *Homo_sapiens.GRCh37.82.gff3.gz.PNN.gff3.gz* is included in the example, but *Homo_sapiens.GRCh37.cdna.all.fa.gz* must be downloaded as described above. 


This command should run in about 10 seconds, which includes the mapping of the probes with *bwa*. The following output files should have been created:

```
PNN.mips.designed_genes.txt  PNN.mips.filtered_by_score.fastq.gz  PNN.mips.without_scores.txt.gz  PNN.mips.with_scores.txt.gz
```

*PNN.mips.with_scores.txt.gz* is a table containing all the smMIPs designed. *PNN.mips.filtered_by_score.fastq.gz* is the corresponding Fastq file with the probe and target sequences for each cDNA-smMIP whose predicted performance score is above the threshold.

An example of a cDNA-smMIP that was produced is:

```
    ctcaacctcagcctcagCTTCAGCTTCCCGATATCCGACGGTAGTGTNNNNNNNNNttgagccagataaagaatgtaaa
```
The smMIP-backbone (which contains the sequence primer binding sites and the UMI sequence) is indicated by upper case letter.
Here the the lower-case sequence represent the ligation probe sequence (left of backbone) and the extension sequence (right of the backbone).
In this example there is only an UMI sequence between the extension probe sequence and the backbone  (the 9 N's).


Next, the Fastq file with candidate smMIPs can be mapped to the genome with a splicing-aware mapper such as STAR:

```
 STAR --runThreadN 4 --genomeDir /data1/genomes/STAR_2.5.1b_modified/hg19/hg19/ \
      --readFilesIn PNN.mips.filtered_by_score.fastq.gz --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate  --outSAMunmapped Within \
      --outFileNamePrefix PNN.mips.filtered_by_score.fastq.gz.STAR.
```
*NOTE* You need to replace */data1/genomes/STAR_2.5.1b_modified/hg19/hg19/* with the actual directory of the STAR genome index on your system. It needs to be generated with STAR first (see STAR manual).

The resulting BAM file (PNN.mips.filtered_by_score.fastq.gz.STAR.Aligned.sortedByCoord.out.bam) needs to be indexed with SAMtools (samtools index PNN.mips.filtered_by_score.fastq.gz.STAR.Aligned.sortedByCoord.out.bam). Finally, the BAM file can be viewed in the IGV genome browser in order to select the cDNA-smMIPs.


In this example, inpection of the BAM with the IGV browser will show that only cDNA-smMIPs have been designed for the first exons of the transcripts. For this particular example, this is because the common variants overlap with the probe sequences of smMIPs that cover the exons towards the 3' end.


Limitations
===========

At the moment, the design script removes all candidate smMIPs for which the probe sequences overlap with a variant specified in the input VCF. 
Unlike the MIPGEN pipeline for DNA resequencing, in such cases the script does *not* design smMIPs where the probe sequence is specific to one of the SNP alleles.
This feature is not included because the downstream differential expression analysis cannot yet handle allele-specific probes (note that this is a different issue from measuring allele-specific expression; in that case the variant overlaps with the target sequence of the smMIP but not one of the smMIP extension/ligation probes). 
If no suitable cDNA-smMIPs can be found that cover particular exons, please try running the script with different values for the option *--target_seq_length*.

At the moment, the script also does not automatic tiling of smMIPs across transcripts; it is up to the user to select the cDNA-smMIPs that best cover the relevant transcripts in a particular experimental system. 
In contrast with DNA-resequencing, where the same DNA is shared by all cells (in principle), transcript abundance is highly cell-type specific. 
We therefore leave it to the user to choose which are the relevant transcripts for which smMIPs should be designed. 
A reference RNA-Seq data set can be helpful in deciding which are relevant transcripts.  



Note on commercial usage 
========================

The code to compute MIP performance scores is based on the MIPGEN 2.0 package to design MIPS for DNA, by Boyle et al, [Bioinformatics, 18(30):2670-2, 2014](http://www.ncbi.nlm.nih.gov/pubmed/24867941).
The copyright is held by Evan Boyle. Commercial use of anything that makes use of the MIP performance scores is not allowed. 

For details on commercial usage, please see the [MIPGEN website](http://shendurelab.github.io/MIPGEN/).
