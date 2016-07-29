Overview
========

This document illustrates the workflow for estimating differential
expression using cDNA-smMIPs.

For this analysis it is assumed that the file with cDNA-smMIP probe sequences and 
corresponding Fastq files from the cDNA-smMIPs sequencing experiment
are available.

Also, the user must provide a study design file that lists the condition (string identifier)
for each experiment (=Fastq file).

In this example it is assumed that there is one pair (forward, reverse) 
of Fastq files for each sample/experiment.


The key steps are:

1. Create study design file;
2. Obtain molecule counts for each experiment and merge into a single file containing counts for all experiments;
3. Estimate normalized log-expression values and differential expression between conditions.


Step 1: Study design file
=================

The study design file is a file with four tab-separated columns:

1. condition
2. experiment
3. forward_fastq
4. reverse_fastq

The first row of the file should contain these exact column labels in this order.

Each row describes a cDNA-smMIP capture experiment.  The *condition* column should contain a string-identifier for the condition.
The *experiment* column should be a unique identifier for the experiment (ie, a row).
Note that if one has generated replicates (biological or technical), multiple experiments will share
the same condition name.

See *files/pbmc_study_design.txt* for an example of a study design file.


Step 2: Obtain molecule counts for each experiment and merge
====================================================

The next step is to obtain the molecule counts from the Fastq files. This can be achieved with the script called *get_molecule_counts_from_fastq.py*.
To see what input it needs, type 
```
python get_molecule_counts_from_fastq.py -h
```
at the Unix command line.

## Input

It requires the following input parameters to be specified:

1. --study_design_file The study design file as described above.
2. --output_file_prefix The prefix of the various output file.
3. --mips_design_file The cDNA-MIPs design file. It should be a tab-separated file with the following columns: *ext_probe_sequence*,*lig_probe_sequence*,*unique_probe_id*.
4. --extension_umi_length The length of the UMI (number of Ns) between the sequence primer and the extension probe.
5. --ligation_umi_length The length of the UMI (number of Ns) between the sequence primer and the ligation probe.

***NOTE*** The length of the UMIs should have been specified when the cDNA-smMIPs were designed. 

***NOTE*** By default, the script will assume that *java* is present, and that the Java-package *CDNA_smMIPS_analysis.jar* (used to obtain the counts from Fastq files) is present in the current path. Alternative locations can be specified with the options *--java* and *--jar_count*.

## Output

The *get_molecule_counts_from_fastq.py* script produces the file *{output_file_prefix}.merged_counts.txt*, where *{output_file_prefix}* is the value that was specified at the command line for this option.

It also produces an output file for each experiment, these are files ending with *molecule_count.log*. These contain additional statistics, such as the fraction of read pairs matched to a smMIP, the number of extnsion-ligation dimers detected, etc.

## Example

To perform this step for the PBMC data in this example, run the following command at the command line:

```
python get_molecule_counts_from_fastq.py --study_design_file files/pbmc_study_design.txt \
       --output_file_prefix pbmc_counts --mips_design_file files/cdna_mips_geuvadis.txt \ 
       --extension_umi_length 9 --ligation_umi_length 0 --fastq_dir fastq/
```

This should produce the file *pbmc_counts.merged_counts.txt*, which should have 381 lines.
This file can be imported into Excel if necessary.
The columns from the *mip_design_file* are listed first, followed by columns giving various statistics. The most relevant ones are:

1. numObservedUMIs, which is the number of unique UMIs observed for reads matched to a given probe.
2. numObservedUMIsCorrected, which is the number of UMIs observed but corrected for sequencing errors (see manuscript for details).
3. condition, which is the condition specified in the study design file.
4. experiment, which is the experiment specified in the study design file.

Step 3: Estimate relative and differential expression
=====================================================

The third step consists of running the Bayesian model on the count data to estimate relative smMIP expression (similar to log(RPKM) expression values)
and to estimate differential expression between conditions.

The two main reasons for using the Bayesian model are:

1. To systematically integrate counts from different replicates for the same condition into a condition-level estimate of expression;
2. To quantify differential expression (log2 fold change) as well as the uncertainty (standard deviation) in the estimates of DE, while correcting for certain smMIP artifacts.

If one is not interested in uncertainty, one may simply estimate expression values for normalizing the numObservedUMIsCorrected counts to log(numObservedUMIsCorrected per million). These may be sufficiently accurate for your purpose.

## Estimating relative expression

To estimate relative expression of the cDNA-smMIPs for the different conditions, use the following command:

```
python estimate_expression.py --merged_counts_file pbmc_counts.merged_counts.txt --output_file_prefix bayesian_model
```

This will run the Bayesian model to estimate normalized log-expression values. These will be output to the file *bayesian_model.bayesian_estimate_log_expression.txt*. Here, the *.bayesian_estimate_log_expression.txt* part is added to the prefix specified with the program option * --output_file_prefix*.

This output file will have a column for each of the conditions that are detected in the merged counts file *pbmc_counts.merged_counts.txt*. In this example, the two conditions are *RPMI* and *Candida*, as specified in the study design file in Step 1.

**NOTE** This step assumes that the merged_counts_file *pbmc_counts.merged_counts.txt* was produced by the script *get_molecule_counts_from_fastq.py*. 
This script creates a *condition* in the merged_counts_file that the script *estimate_expression.py* uses to determine which samples/experiments belong to
same condition and are treated as replicates. 

## Estimating differential expression

To estimate differential expression, you can simply specify which pair of conditions to compare using the *--compare_conditions* program option:

```
python estimate_expression.py --merged_counts_file pbmc_counts.merged_counts.txt \
                              --load_model_from_pickle --output_file_prefix bayesian_model \ 
                              --compare_conditions RPMI/Candida
```
The two conditions that should be compared have to be separated by a '/'. Multiple pairs of conditions can be specified by separating them with a comma (',').Although it is not particularly useful in this example, as there only two conditions, this can be done like so:

```
python estimate_expression.py --merged_counts_file pbmc_counts.merged_counts.txt \
                              --load_model_from_pickle --output_file_prefix bayesian_model \
                              --compare_conditions RPMI/Candida,RPMI/Candida
```
This will cause the RPMI/Candida comparison to be done twice.

### Output of differential expression analysis

The differential expression analysis performed by the option *--compare_conditions* creates an output file *output_file_prefix.bayesian_estimate_log2_differential_expression.txt* with two columns for each pair of conditions. For this example, the top rows are:

|                         | mean_log2_DE_RPMI_Candida | stdev_log2_DE_RPMI_Candida |
| ----------------------- | ------------------------- | -------------------------- |
| ENSG00000011304.1       | 0.124927987015            | 0.342951043889             | 
| ENSG00000011304.10      | -0.338230710403           | 0.415207159972             |

The third column records the standard deviation of the *posterior* distribution for the DE estimate. Thus, to get an approximate 95% confidence interval,
one can use (mean_log2_DE_RPMI_Candida - 2*stdev_log2_DE_RPMI_Candida, mean_log2_DE_RPMI_Candida + 2*stdev_log2_DE_RPMI_Candida).


## Additional options

The *estimate_expression.py* script has a number of additional options:

1. *--excel* will save output files also in Excel format;
2. *--save_model_as_pickle* will save the results from the MCMC inference to a file, so that the MCMC inference does not have to be re-run when one wants to estimate differential expression between different pairs of conditions.
3. *--load_model_from_pickle* If one has previously run *estimate_expression.py* with the *--save_model_as_pickle* flag, then specifying this flag will not run MCMC again but instead load it from the file. ***NOTE*** that requires that the *same* *--output_file_prefix* is used for load and save as the name of the file with MCMC results (either to save or load from) is constructed from the *--output_file_prefix* value.

For example, the save model feature may be used as follows:
```
# fit Bayesian model and estimate relative expression 
python estimate_expression.py --merged_counts_file pbmc_counts.merged_counts.txt \
                              --save_model_as_pickle --output_file_prefix bayesian_model
# Estimate differential expression between specified conditions
python estimate_expression.py --merged_counts_file pbmc_counts.merged_counts.txt \
                              --load_model_from_pickle --output_file_prefix bayesian_model \
                              --compare_conditions RPMI/Candida
```
Here, the first command will run the MCMC model and save the MCMC results to a file. The second command will not re-run the MCMC inference for the Bayesian model, but instead load it from the file *bayesian_model.stan_model.pkl.stan_cdna_smmips_mt.pkl*, and then estimate differential expression for probe between the RPMI and Candida condition.

Limitations
===========

At the moment, the *estimate_expression.py* allows differential expression only to be estimated between pairs of conditions. 



