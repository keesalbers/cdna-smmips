import os, sys
import pandas as pd
import gzip
import argparse
import subprocess


def load_study_design(fn):
    x = pd.read_table(fn, sep = "\t")
    assert x.columns[0] == "condition"
    assert x.columns[1] == "experiment"
    assert x.columns[2] == "forward_fastq"
    assert x.columns[3] == "reverse_fastq"
    assert x.shape[0] == len(set(x['experiment'])), "Experiment names must be unique"
    return x

def get_counts(design, args):
    design = design.copy()
    all_counts = []
    for i in design.index:
        prefix_experiment = args.output_file_prefix + "." + design.loc[i,"condition"] + "." + design.loc[i,"experiment"] 
        design.loc[i,"prefix_count"] = prefix_experiment
        if args.fastq_dir != "":
            fastq1 = args.fastq_dir + "/"
        else:
            fastq1 = ""
        fastq1 += design.loc[i, "forward_fastq"]
        if args.fastq_dir != "":
            fastq2 = args.fastq_dir + "/"
        else:
            fastq2 = ""
        fastq2 += design.loc[i, "reverse_fastq"]
        if not os.path.exists(fastq1):
            sys.stderr.write("Fastq file %s does not exist for condition %s and experiment %s.\n" % (fastq1, design.loc[i,"condition"], design.loc[i,"experiment"]))
            sys.exit(1)
        
        if not os.path.exists(fastq2):
            sys.stderr.write("Fastq file %s does not exist for condition %s and experiment %s.\n" % (fastq2, design.loc[i,"condition"], design.loc[i,"experiment"]))
            sys.exit(1)

        cmd = "%s -Xmx4g -jar %s PROGRAM=COUNTUMIS BLOCKCOUNTUMIS=0 Fastq1=%s Fastq2=%s inputDesignFile=%s outputFile=%s ligationUmiLength=%d extensionUmiLength=%d seedSequenceHashLength=5 blockSize=100000000" % (args.java, args.jar_count, fastq1, fastq2, args.mips_design_file, prefix_experiment, int(args.ligation_umi_length), int(args.extension_umi_length))

        print "\tProcessing experiment %s for condition %s:\n\t %s\n" % (design.loc[i,"experiment"], design.loc[i,"condition"],cmd)

        try:
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            with open("%s.molecule_count.log" % prefix_experiment, 'w') as flog: 
                for line in output:
                    flog.write(line)
        except Exception,e:
            sys.stderr.write("There was an error:\n%s\n" % str(e))
            sys.exit(2)

        # READ COUNTS OUTPUT FILE
        counts_file = "%s.merged.statistics.txt" % prefix_experiment
        counts = pd.read_table(counts_file, sep = "\t")
        counts["condition"] = design.loc[i,"condition"]
        counts["experiment"] = design.loc[i,"experiment"]
        all_counts.append(counts)

    return design, all_counts

# run get_molecule_counts_from_fastq.py --study_design_file pbmc_study_design.txt --output_file_prefix pbmc_counts--mips_design_file cdna_mips_geuvadis.txt --extension_umi_length 9 --ligation_umi_length 0 --fastq fastq/ --jar_count scripts/CDNA_smMIPS_analysis.jar
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--java", help = "Path to java binary.", default = "java")
    parser.add_argument("--fastq_dir", help = "Directory with fastq files (prepended to fastq filenames in study design file)", default = "")
    parser.add_argument("--jar_count", help = "Path to JAR with counting program", default = "CDNA_smMIPS_analysis.jar")


    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--study_design_file", help="Tab-separated file with four columns condition,experiment,forward_fastq,reverse_fastq", required = True)
    requiredNamed.add_argument("--mips_design_file", help="File with probe sequences cDNA-smMIPs", required = True)
    requiredNamed.add_argument("--output_file_prefix", help ="Prefix of output files. ", required = True) 

    requiredNamed.add_argument("--extension_umi_length", help = "Length of UMI on the side of the extension probe. Specify zero if not present.", required = True)
    requiredNamed.add_argument("--ligation_umi_length", help = "Length of UMI on the side of the ligation probe. Specify zero if not present.", required = True)
    

    args = parser.parse_args()
    
    sys.stdout.write( "Checking study design...\n")
    sys.stdout.flush()
    design = load_study_design(args.study_design_file)
    sys.stdout.write("Done.\n")

    sys.stdout.write( "Counting molecules...\n")
    sys.stdout.flush()
    ann_design, all_counts = get_counts(design, args)
    sys.stdout.write("Done.\n")
    
    sys.stdout.write( "Merging counts...\n")
    merged_counts = pd.concat(all_counts, ignore_index=True)

    output_file = args.output_file_prefix + ".merged_counts.txt"
    sys.stdout.write("Wrote merged counts for all experiments to %s\n" % output_file)
    merged_counts.to_csv(output_file, sep = "\t", index = False)

    sys.stdout.write("Finished\n")
if __name__ == "__main__":
    main()
