Overview
========

This document describes how to estimate allelic ratios for specific 
coding variants using cDNA-smMIPs.

There are three steps:

1. Remove UMI sequence from the read sequence and put them in the read ID, produces a new pair of Fastq files (CDNA_smMIPS_analysis.jar).
2. Map the trimmed reads from new Fastq file to the genome using an RNA-Seq mapper (STAR).
3. Obtain allele-specific molecule counts for all targeted variants (count_umis_at_mips_targets_ase.py).

In this approach, reads are mapped to the genome in order to count the SNP alleles at positions specified in terms of a genome reference.


Step 1: Produce UMI-trimmed Fastq files
=======================================

For this example, the UMIs are removed from the reads and added to the read IDs with the following command:

```
java -Xmx2g -jar CDNA_smMIPS_analysis.jar program=REMOVEUMIS \
                 fastq1=fastq/JB_K562HEK_3-1_R1_09-03-2016.fastq.gz fastq2=fastq/JB_K562HEK_3-1_R2_09-03-2016.fastq.gz \
                 outputfastq1=fastq/stripped.JB_K562HEK_3-1_R1_09-03-2016.fastq.gz \
                 outputfastq2=fastq/stripped.JB_K562HEK_3-1_R2_09-03-2016.fastq.gz \
                 ligationUmiLength=0 extensionUmiLength=9
```

Again, the *ligationUMILength* and *extensionUMILength* should be specified by the user based on the parameters that were used for the cDNA-smMIP design.

The first record in the unstripped forward and reverse Fastq files are:

Forward, original
```
 fastq/JB_K562HEK_3-1_R1_09-03-2016.fastq.gz
 @NB501280:19:H75YHAFXX:1:11101:4542:1249 1:N:0:GTCGGTAA
 GGCCTCTTTTCTTCCCCTTCCCCATCTCCAAGGCCTCCCCGCGCATGACACTACCGTCGGATCGTGCGTGTGTCGGTAAA
 +
 AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEE
```
Forward, UMI stripped (read length unchanged, but read ID is changed)
```
 fastq/stripped.JB_K562HEK_3-1_R1_09-03-2016.fastq.gz
 @NB501280:19:H75YHAFXX:1:11101:4542:1249::CATGCGCGG 1:N:0:GTCGGTAA
 GGCCTCTTTTCTTCCCCTTCCCCATCTCCAAGGCCTCCCCGCGCATGACACTACCGTCGGATCGTGCGTGTGTCGGTAAA
 +
 AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEE
```
Reverse, original
```
 fastq/JB_K562HEK_3-1_R2_09-03-2016.fastq.gz
 @NB501280:19:H75YHAFXX:1:11101:4542:1249 2:N:0:GTCGGTAA
 CATGCGCGGGGAGGCCTTGGAGATGGGGAAGGGGAAGAAAAGAGGCCCTTCAGCTTCCCGATTACGGATCTCGTATGTGT
 +
 A/AAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
```
Reverse, UMI stripped (read length changed, ID changed)
``` 
 fastq/stripped.JB_K562HEK_3-1_R2_09-03-2016.fastq.gz
 @NB501280:19:H75YHAFXX:1:11101:4542:1249::CATGCGCGG 2:N:0:GTCGGTAA
 GGAGGCCTTGGAGATGGGGAAGGGGAAGAAAAGAGGCCCTTCAGCTTCCCGATTACGGATCTCGTATGTGT
 +
 EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
```
The read sequence of the forward read is unchanged after the stripping procedure, as
in this example only the extension probe (the second read in sequencing) contains the UMI.
However, the ID is changed from *@NB501280:19:H75YHAFXX:1:11101:4542:1249* to 
*@NB501280:19:H75YHAFXX:1:11101:4542:1249::CATGCGCGG*, where CATGCGCGG is the UMI sequence 
taken from the reverse read (first nine nucleotides of the reverse read).
The reverse read does change in length as this contains the UMI sequence in the first nine nucleotides.

## Output

This step creates two new fastq files with stripped UMIs. In this examples, the two output files are:
```
 fastq/stripped.JB_K562HEK_3-1_R1_09-03-2016.fastq.gz 
 fastq/stripped.JB_K562HEK_3-1_R2_09-03-2016.fastq.gz
```


Step 2: Map trimmed reads to reference genome
=============================================

The next step is to align the trimmed reads to a reference genome. It is essential that the read IDs contain the 
original UMI sequence, as step 3 requires the UMIs to infer the individual molecules from the mapped reads. As
long as the read IDs containing the UMI sequence are preserved, any dedicated RNA-Seq mapper can be used in principle.
However, for allele-specific analysis one should always keep mapping bias in mind.

We used STAR [Dobin et al, 2013] to map reads to the genome with the following command:

```
 STAR --runThreadN 4 --genomeDir /data1/genomes/STAR_2.5.1b_modified/hg19/hg19/ \
      --readFilesIn fastq/stripped.JB_K562HEK_3-1_R1_09-03-2016.fastq.gz fastq/stripped.JB_K562HEK_3-1_R2_09-03-2016.fastq.gz \
      --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --outSAMunmapped Within --outFileNamePrefix STAR.JB_K562HEK_3-1. \
```

**Note*** that the directory */data1/genomes/STAR_2.5.1b_modified/hg19/hg19/* refers to the STAR-index of the reference sequence hg19. 
On your system you need to first generate this index with STAR (see STAR manual) and replace the directory with the appropriate one on your system. The loading of the STAR index when STAR starts can take a few minutes. Also, you may need to install STAR on your system.

## Output

Among other files created by STAR, the main output file is 
```
 STAR.JB_K562HEK_3-1.Aligned.sortedByCoord.out.bam
```

This bam files is sorted by coordinate. It must still be indexed with Samtools, so that the script in Step 3 can have random access to specific genomic regions in the BAM file.
Use the following command, which assumes SAMtools [htslib.org](http://htslib.org) is installed: 
```
 samtools index STAR.JB_K562HEK_3-1.Aligned.sortedByCoord.out.bam
```

Step 3: Get allele counts from BAM file
=======================================

The third step determines for each target variants what the molecule counts are for all cDNA-smMIPs that cover it.

## Input requirements

The target variants should be specified in the *ase_mip_design_file*. It should have the same columns as the file with cDNA-smMIPs for expression analysis, and *in addition* a column that describes the target variant. It is of the following format (CHROM,POS,REF,ALT joined by an underscore):
```
 chr_{chromomosome id}_{position}_{reference allele}_{non-reference allele}
```
For example:
```
 chr_1_236749649_T_C
```

To obtain the counts, run the following command:
```
 python count_umis_at_mips_targets_ase.py --ase_mip_design_file files/cdna_mips_ase-K562-HEK293.txt \
        --bamfile STAR.JB_K562HEK_3-1.Aligned.sortedByCoord.out.bam \
        --output_file counts.JB_K562HEK_3-1.txt
```

## Output

The output file (in this example, counts.JB_K562HEK_3-1.txt) is a tab-separated file that lists allele counts (molecule counts based on the UMIs) for each target variant and each cDNA-smMIP. 
Thus, if a certain variant is covered by multiple smMIPs, there will be a line for each probe. This makes it possible to check consistency of estimates between smMIPs.

Limitations
==========

At the moment, the script can only count alleles for SNPs, but not small indels.

The read mapping may introduce a bias in the allele counts. In principle it should be easy to use a more unbiased approach and use a de novo
assembler to infer haplotypes from the cDNA-smMIPs reads for a given smMIP, and then do variant calling based on that.  

