# cDNA-smMIPs

## Overview

Single-molecule inversion probes (smMIPs) can be used to resequence small DNA fragments (100-112 nt) 
and allow single-molecule counting through unique molecular identifiers. 
They were developed by the group of [Jay Shendure](http://krishna.gs.washington.edu/), who used it to identify disease mutations in large cohorts. 
See [PMID:23160955](http://www.ncbi.nlm.nih.gov/pubmed/23160955) and [PMID:23382536](http://www.ncbi.nlm.nih.gov/pubmed/23382536) for details.

While the efficiency of smMIPs is lower compared to that of regular primers and variable from smMIP to smMIP,
it turns out that performance of a single smMIP is quite robust.
As a result, when applied to cDNA (reverse-transcribed RNA), cDNA-smMIPs can be used to accurately quantify differential gene expression 
(by comparing the same probe between conditions) and allele-specific expression.
With multiple smMIPs per transcript it is also possible to quantify relative gene expression.

A manuscript describing cDNA-smMIPs has been submitted.

## What is here?

This repository contains software and instructions for designing cDNA-smMIps and statistical analysis of the resulting data.

- [A Python script to design cDNA-smMIPs](src/python/design_mips_for_all_transcripts.py)
- [A Python script to obtain molecule counts from Fastq files](src/python/get_molecule_counts_from_fastq.py)
- [A Python script to estimate relative and differential expression](src/python/estimate_expression.py)
- [A Python script to obtain allele counts for allele-specific expression analysis](src/python/count_umis_at_mips_targets_ase.py)

## Documentation

Examples and documentation on how to use the above scripts are in the [example](example/) folder:

- [How to design cDNA-smMIPs](example/How_To_Design_cDNA_smMIPs.md)
- [How to estimate relative and differential expression](example/How_To_Estimate_Relative_and_Differential_expression.md)
- [How to estimate allelic ratios](example/How_To_Estimate_Allelic_Ratios.md)

Note that the script to count the molecules requires [this Java program](src/java/CDNA_smMIPS_analysis.jar). 
It is also required to preprocess Fastq files for the ASE analysis, as described in the How To's in the [example](example/) directory.

## Designing cDNA-smMIPs

This repository contains Python scripts to design cDNA-smMIPs for a set of genes, using cDNA transcript sequences from Ensembl to design the probes against.
The directory [example](example/) contains an example detailing how to design cDNA-smMIPs.
A VCF file with genomic variants that should be excluded from the smMIP probes may be specified.
It is also possible to provide a VCF file with variants that should be targeted by the cDNA-smMIPs (ie, covered by the sequence between the two probes) for ASE analysis.

## Examples with real data

The directory [example](example/) contains a tutorial based on a real data set, from Fastq-file to analysis of differential expression with
the statistical model that we developed for cDNA-smMIPs.
There is also a tutorial for allele-specific expression analysis with cDNA-smMIPs.

## Installation of software

Java and Python is required to run the software. Running on a non-Unix OS is not supported, as the scripts require access to the GCC compiler for compilation
of certain inline functions through Scipy.weave. 

The Python scripts can be run stand-alone provided the correct libraries have been installed. 
We recommend the use of [Anaconda](https://www.continuum.io/why-anaconda) to manage these.
Apart from standard scientific Python libraries such as Pandas, Numpy and Scipy, a few dedicated bioinformatics
libraries are required, including [Pysam](http://pysam.readthedocs.io/en/latest/api.html), [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html) and 
[pybedtools](https://pythonhosted.org/pybedtools/). 
These can be installed with the Python package installer *pip* (the pip version of Pysam is more recent than the Anaconda one) or through conda itself.

## Questions

For questions or comments about the software, contact Kees.
