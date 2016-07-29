import os, sys
import pandas as pd
import itertools
import pybedtools as pbt
import numpy as np
import pysam
import HTSeq
import scipy.weave
from scipy.weave import converters
import gzip
import argparse


RCBASE = { 'A':'T',
           'C':'G',
           'G':'C',
           'T':'A' }

def kprint(msg):
    sys.stdout.write(msg)
    sys.stdout.flush()

class Shendure_scores:
    def __init__(self, args):
        self.args = args
        self.ext_fastq_file = self.args['output_prefix'] + ".ext_probe.fastq.gz"
        self.lig_fastq_file = self.args['output_prefix'] + ".lig_probe.fastq.gz"
        self.ext_sam = self.args['output_prefix'] + ".ext_probe.fastq.gz.sam.gz"
        self.lig_sam = self.args['output_prefix'] + ".lig_probe.fastq.gz.sam.gz"




    def get_shendure_score(self, row):
        scan_seq = row['target_seq'] #row['scan_target_sequence']
        ext_seq = row['ext_probe_seq']
        lig_seq = row['lig_probe_seq']
        scan_start = 0 # row['cmip_mip_transcript_start']
        ext_count = row['ext_probe_copy_count']
        lig_count = row['lig_probe_copy_count']
        chrom = "chr_mip" # "%s" % row['cmip_transcript']
        mip = shendure_get_mip_score.PlusMIP(chrom = chrom, scan_start = scan_start, scan_target_sequence = scan_seq, ext_probe_sequence = ext_seq, lig_probe_sequence = lig_seq, ext_probe_copy = ext_count, lig_probe_copy = lig_count)
        return mip.get_score()

    def get_shendure_score_weave(self, row):
        scan_seq = row['target_seq'] #row['scan_target_sequence']
        ext_seq = row['ext_probe_seq']
        lig_seq = row['lig_probe_seq']
        scan_start = 0 # row['cmip_mip_transcript_start']
        ext_count = row['ext_probe_copy_count']
        lig_count = row['lig_probe_copy_count']


        code = """
        /*
        Code is taken from MIPGEN 2.0.
        MIPgen: optimized modeling and design of molecular inversion probes for targeted resequencing.
        Boyle et al, Bioinformatics 18(30):2670-2, 2014
        COPYRIGHT EVAN BOYLE, 2014.
        COMMERCIAL USE IS NOT PERMITTED.
        */
        using namespace std;
        std::string ext_probe_sequence = ext_seq;
        std::string lig_probe_sequence = lig_seq;
        std::string scan_target_sequence = scan_seq;
        double ligation_arm_length = double(lig_probe_sequence.size());
        double extension_arm_length = double(ext_probe_sequence.size());
        double scan_size = double(scan_seq.size());
        double ext_probe_copy = double(ext_count);
        double lig_probe_copy = double(lig_count);

        std::string ligation_junction = lig_probe_sequence.substr(0,2);

        map<string, double> junction_scores;
        junction_scores["AA"] = 0.0;
        junction_scores["AC"] = 0.35;
        junction_scores["AG"] = 0.046;
        junction_scores["AT"] = 0.079;
        junction_scores["CA"] = 0.34;
        junction_scores["CC"] = 0.22;
        junction_scores["CG"] = 0.55;
        junction_scores["CT"] = -0.071;
        junction_scores["GA"] = 0.35;
        junction_scores["GC"] = 0.92;
        junction_scores["GG"] = 0.24;
        junction_scores["GT"] = 0.48;
        junction_scores["TA"] = -0.46;
        junction_scores["TC"] = -0.35;
        junction_scores["TG"] = -0.25;
        junction_scores["TT"] = -0.98;

        if (ext_probe_sequence.find("N") < (unsigned int) extension_arm_length || lig_probe_sequence.find("N") < (unsigned int) ligation_arm_length) {
            return_val = -1000.0;
        } else {

            string last_base = scan_target_sequence.substr(0, 1);
            double run_count = 0;
            for (int i = 1; i< int(scan_size); i++)
            {
                string current_base = scan_target_sequence.substr(i, 1);
                if (current_base == "G" || current_base == "C")
                {
                    if (last_base == "G" || last_base == "C") {}
                    else
                    {
                        run_count++;
                        last_base = current_base;
                    }
                }
                else
                {
                        if ("A" == last_base || "T" == last_base) {}
                        else {
                        run_count++;
                        last_base = current_base;
                    }
                }
            }
            run_count++;
            double bases_per_switch = scan_size / run_count;
            double ext_g_count = (double) std::count(ext_probe_sequence.begin(), ext_probe_sequence.end(), 'G');
            double lig_g_count = (double) std::count(lig_probe_sequence.begin(), lig_probe_sequence.end(), 'G');
            double target_g_count = (double) std::count(scan_target_sequence.begin(), scan_target_sequence.end(), 'G');

            double ext_gc_count = (double) std::count(ext_probe_sequence.begin(), ext_probe_sequence.end(), 'C') + ext_g_count;
            double lig_gc_count = (double) std::count(lig_probe_sequence.begin(), lig_probe_sequence.end(), 'C') + lig_g_count;
            double target_gc_count = (double) std::count(scan_target_sequence.begin(), scan_target_sequence.end(), 'C') + target_g_count;

            double ext_a_count = (double) std::count(ext_probe_sequence.begin(), ext_probe_sequence.end(), 'A');
            double lig_a_count = (double) std::count(lig_probe_sequence.begin(), lig_probe_sequence.end(), 'A');
            double target_a_count = (double) std::count(scan_target_sequence.begin(), scan_target_sequence.end(), 'A');

            double ext_length = extension_arm_length;
            double lig_length = ligation_arm_length;
            double target_length = scan_size > 250 ? 250 : scan_size;

            double ext_gc_content = ext_gc_count/ext_length;
            double lig_gc_content = lig_gc_count/lig_length;
            double target_gc_content = target_gc_count/scan_size;

            double ext_g_content = ext_g_count/ext_length;
            double lig_g_content = lig_g_count/lig_length;
            double target_g_content = target_g_count/scan_size;

            double ext_a_content = ext_a_count/ext_length;
            double lig_a_content = lig_a_count/lig_length;
            double target_a_content = target_a_count/scan_size;

            double junction_score = junction_scores[ligation_junction];

            double log_ext_copy = ext_probe_copy > 100 ? 2 : log10(ext_probe_copy);
            double log_lig_copy = lig_probe_copy > 100 ? 2 : log10(lig_probe_copy);
            double exponent =
                -35.0464 - 4 +
                -1.974282*bases_per_switch +
                2.63667*bases_per_switch*ext_gc_content +
                2.540741*bases_per_switch*lig_gc_content +
                -0.006488*bases_per_switch*target_length +
                -0.018137*ext_a_content +
                7.795877*pow(ext_a_content, 2) +
                -5.576753*ext_a_content*ext_g_content +
                -0.274062*ext_a_content*ext_length +
                8.695568*ext_a_content*lig_g_content +
                -2.014126*ext_a_content*log_ext_copy +
                3.163087*ext_a_content*log_lig_copy +
                -19.900678*pow(ext_gc_content, 2) +
                35.747084*ext_gc_content +
                -2.082136*ext_g_content +
                7.204324*pow(ext_g_content, 2) +
                11.73888*ext_g_content*lig_g_content +
                -2.173235*ext_g_content*log_lig_copy +
                -7.123214*ext_g_content*target_gc_content +
                1.068617*ext_length +
                -0.008666*pow(ext_length, 2) +
                -0.555599*ext_length*ext_gc_content +
                -0.289857*ext_length*lig_g_content +
                -0.009621*ext_length*lig_length +
                0.119863*ext_length*log_ext_copy +
                2.405833*junction_score +
                -1.764289*junction_score*lig_gc_content +
                2.112564*junction_score*lig_g_content +
                0.656183*junction_score*log_ext_copy +
                -3.099451*junction_score*target_a_content +
                -2.097335*junction_score*target_gc_content +
                6.542827*pow(lig_a_content, 2) +
                -8.885956*lig_a_content +
                5.297562*lig_a_content*ext_g_content +
                13.042436*lig_a_content*lig_gc_content +
                -12.333361*lig_a_content*lig_g_content +
                19.866468*lig_gc_content +
                -13.698276*pow(lig_gc_content, 2) +
                -0.517301*lig_gc_content*lig_length +
                -2.846777*lig_gc_content*log_lig_copy +
                13.64081*lig_gc_content*target_a_content +
                13.614709*lig_g_content +
                -4.759165*lig_g_content*ext_gc_content +
                -7.48883*lig_g_content*lig_gc_content +
                -2.308594*lig_g_content*log_ext_copy +
                -12.640154*lig_g_content*target_a_content +
                1.164626*lig_length +
                -0.010326*pow(lig_length, 2) +
                -0.354197*lig_length*ext_gc_content +
                0.111448*lig_length*log_ext_copy +
                0.238095*lig_length*target_a_content +
                -4.161632*log_ext_copy +
                -2.728864*log_ext_copy*ext_gc_content +
                0.641717*log_ext_copy*log_lig_copy +
                3.738798*log_ext_copy*target_gc_content +
                -1.98457*log_lig_copy +
                2.362253*log_lig_copy*ext_gc_content +
                3.467229*log_lig_copy*target_gc_content +
                -18.443242*target_a_content +
                20.89245*pow(target_a_content, 2) +
                -0.048679*target_a_content*target_length +
                -50.249451*pow(target_gc_content, 2) +
                27.132716*target_gc_content +
                -0.050633*target_gc_content*target_length +
                -20.772366*target_g_content +
                -60.796481*pow(target_g_content, 2) +
                26.630245*target_g_content*target_a_content +
                87.162648*target_g_content*target_gc_content +
                0.030256*target_g_content*target_length +
                0.032811*target_length;
            return_val = double( pow(2.71828, exponent)/(1 + pow(2.71828, exponent)) );
        }
        """
        score = scipy.weave.inline(code, ['scan_seq','ext_seq','lig_seq','ext_count','lig_count'], compiler = 'gcc', headers = ["<vector>","<map>","<cmath>","<string>","<algorithm>","<map>" ], extra_compile_args = ['-O3'])
        return score


    def map_bwa_se(self, ref_file, fastq_file, sam_file, skip_mapping = False):

        if not skip_mapping:
            cmd1 = "%s aln %s %s > %s.sai" % (self.args['bwa'], ref_file, fastq_file, fastq_file)
            kprint("Executing "+cmd1)
            res1 = os.system(cmd1)
            assert res1 == 0, "Error in bwa alignment: " + cmd1

            cmd2 = "%s samse -n 100 %s %s.sai %s | gzip > %s" % (self.args['bwa'], ref_file, fastq_file, fastq_file, sam_file)
            kprint("Executing "+cmd2)
            res2 = os.system(cmd2)
            assert res2 == 0, "Error in bwa alignment: " + cmd2

            os.remove("%s.sai" % fastq_file)

    def make_fastq_file(self, seqs, fn):
        assert fn.endswith('.fastq.gz')
        f = gzip.open(fn,'w')
        for i, seq in enumerate(seqs):
            f.write("@seq%d.%s\n" % (i,seq))
            f.write(seq + "\n+\n")
            f.write("".join(['A']*len(seq))+"\n")
        f.close()


    def get_counts_dict(self, seqs, samfile):
        copy_count = np.ones( (len(seqs)), dtype='int')*-1    
        cd = {}
        f = gzip.open(samfile, 'r')
        for line in f:
            if line[0] == '@':
                continue
            line = line.rstrip("\n")
            dat = line.split('\t')
            rid = dat[0]
            assert rid.startswith('seq')
            srid = rid.split('.')
            read_no = int(srid[0][3:])
            probe_seq = srid[1]
            if dat[2] == '*':
                count = 100 # unmapped
            else:
                count = 1 # the actual location 

                if len(dat) >= 20:
                    assert dat[19].startswith('XA')
                    xa = dat[19].split(';')
                    for ml in xa:
                        if ml != '':
                            mls = ml.split(',')
                            nm = int(mls[-1])
                            if nm == 0:
                                count += 1
                else:
                    # bwa omits XA tag if there are more than 100 mapping locations or only 1 location
                    if int(dat[4]) == 0: # mapping quality has to be zero for 100 copies, oth
                        count = 100
            copy_count[read_no] = count
            cd[ seqs[read_no] ] = count
        f.close()
        assert (copy_count == -1).sum() == 0
        return cd, copy_count


    def map_bwa_se_probes(self, ext_probe_seqs, lig_probe_seqs):
        if True:
            print "Making extension probe fastqs"
            self.make_fastq_file(ext_probe_seqs, self.ext_fastq_file)
            print "Mapping extension probes"
            self.map_bwa_se(self.args['ref_genome_fasta'], self.ext_fastq_file, self.ext_sam)

            print "Making ligation probe fastqs"
            self.make_fastq_file(lig_probe_seqs, self.lig_fastq_file)
            print "Mapping ligations probes"
            self.map_bwa_se(self.args['ref_genome_fasta'], self.lig_fastq_file, self.lig_sam)

        return self.get_counts_dict(ext_probe_seqs, self.ext_sam)[0], self.get_counts_dict(lig_probe_seqs, self.lig_sam)[0], []


    def compute_scores(self):
        tv = pd.read_table(self.args['candidate_mip_file'], sep="\t", compression = "gzip")
        
        ext_probes = list(set(tv['ext_probe_seq']))
        lig_probes = list(set(tv['lig_probe_seq']))

        print "Computing BWA copy count numbers for %d extension probes and %d ligation probes\n" % (len(ext_probes), len(lig_probes))

        # to get copynumber counts, do single-end BWA mapping
        ext_probe_counts, lig_probe_counts, pair_mapq = self.map_bwa_se_probes(ext_probes, lig_probes)
        otv = pd.DataFrame(tv)
        otv['index'] = range(otv.shape[0])
        otv = otv.set_index('index')
        otv['ext_probe_copy_count'] = -1
        otv['lig_probe_copy_count'] = -1
        for i in otv.index:
            ext_probe_seq, lig_probe_seq = otv.loc[i,"ext_probe_seq"], otv.loc[i,"lig_probe_seq"]
            otv.loc[i,'ext_probe_copy_count'] = ext_probe_counts[ext_probe_seq]
            otv.loc[i,'lig_probe_copy_count'] = lig_probe_counts[lig_probe_seq]
            if ext_probe_counts[ext_probe_seq] >= 100 and  lig_probe_counts[lig_probe_seq] >= 100:
                otv.loc[i,'shendure_logistic_score'] = 0.0
            else:
                otv.loc[i,'shendure_logistic_score'] = self.get_shendure_score_weave(otv.loc[i])
            if i % 5000 == 1:
                print "Finished",i
        otv['ext_probe_copy_count'] = pd.Series(otv['ext_probe_copy_count'], dtype='int')
        otv['lig_probe_copy_count'] = pd.Series(otv['lig_probe_copy_count'], dtype='int')
        otv.to_csv(self.args['output_file'], sep = "\t", index = None, compression = 'gzip')
        if self.args['remove_fastq_sam']:
            os.remove(self.ext_sam)
            os.remove(self.lig_sam)
            os.remove(self.ext_fastq_file)
            os.remove(self.lig_fastq_file)
            

class Exon:
    def __init__(self, gff, transcript):
        # reference to transcript instance
        self.transcript = transcript
        # copy exon feature from gff 
        self.gff = pbt.interval_constructor(gff)
    def length(self):
        l = self.gff.end - self.gff.start
        assert l >= 0
        return l
    def getid(self):
        return self.gff.attrs['exon_id']
    def set_transcript_start(self, tstart):
        self.tstart = tstart

class Transcript:
    def __init__(self, gff, gene):
        # chromosome
        self.chrom = gff.chrom 
        self.gff = pbt.interval_constructor(gff)
        # hash map from integer position to exon instance
        self.exonid_to_exon = {}
        self.number_to_exon = {}
        self.exonid_to_number = {}
        self.transcript_id = gff.attrs['transcript_id']
        self.seq = None

    def add_exon(self, chrom, firstbase, lastbase, exon_id):
        if not self.firstbase_to_exon.has_key(firstbase):
            self.firstbase_to_exon = {}
            self.firstbase_to_exon[firstbase] = {}
        else:
            if self.firstbase_to_exon[firstbase].has_key(lastbase):
                exon  = self.firstbase_to_exon[firstbase][lastbase]
                sys.stderr.write("Exon with coordinates already exists: %s\n" % (exon.get_id()))
                return exon

        exon = Exon(firstbase, lastbase, exon_id, self)
        self.firstbase_to_exon[firstbase][lastbase] = exon
        assert self.exonid_to_exon.has_key(exon_id) == False
        self.exonid_to_exon[exon_id] = exon
        return exon

    def length(self):
        l = 0
        for i, e in self.number_to_exon.iteritems():
            l += e.gff.end - e.gff.start
        assert l >= 0
        return l
    def add_exon_gff(self, gff):
        exonid = gff.attrs['exon_id']
        if not self.exonid_to_exon.has_key(exonid):
            exon = Exon(gff, self)
            number = int(gff.attrs['rank'])
            self.exonid_to_exon[exonid] = exon
            self.number_to_exon[number] = exon
            self.exonid_to_number[exonid] = number
        else:
            assert 1 == 0
        return exon
    def number_keys_sorted(self):
        return sorted(self.number_to_exon.keys())
   
    def get_sequence_exonid(self, exonid):
        tpos = self.exonid_to_exon[exonid].tstart
        tend = tpos + self.exonid_to_exon[exonid].length()
        return self.seq[tpos:tend]


    def kprint(self):
        numbers = sorted(self.number_to_exon.keys())
        for num in numbers:
            v = self.number_to_exon[num]
            exonid = v.gff.attrs['exon_id']
            print "\t\tExon",exonid, v.gff.chrom, v.gff.start, v.gff.end
    def set_exon_transcript_coordinates(self):
        numbers = sorted(self.number_to_exon.keys())
        tpos = 0
        for num in numbers:
            v = self.number_to_exon[num]
            v.set_transcript_start(tpos)
            tpos += v.length()
        return tpos 

class Gene:
    def __init__(self, gff, genecollection):
        self.id_to_transcripts = {}
        self.exonid_to_exon = {}
        self.gff = pbt.interval_constructor(gff)
        self.genecollection = genecollection
        self.transcripts_cdna_status = {}

    def add_transcript(self, gff):
        transcript_id = gff.attrs['transcript_id']
        gene_id = gff.attrs['Parent'][gff.attrs['Parent'].index('ENSG'):]
        assert gene_id == self.gff.attrs['gene_id']
        assert transcript_id not in self.id_to_transcripts
        if not self.id_to_transcripts.has_key(transcript_id):
            self.id_to_transcripts[transcript_id] = Transcript(gff, self)
 
    def add_transcript_exon_gff(self, gff):
        # NOTE that this is a transcript-specific exon. An exon with the same exon-id can be present in multiple transcripts at different positions.
        transcript_id = gff.attrs['Parent'][gff.attrs['Parent'].index('ENST'):]
        transcript = self.id_to_transcripts[transcript_id]
        exon = transcript.add_exon_gff(gff)
        if not self.exonid_to_exon.has_key(gff.attrs['exon_id']):
            self.exonid_to_exon[gff.attrs['exon_id']] = []
        self.exonid_to_exon[gff.attrs['exon_id']].append(exon)
        return exon

    def kprint(self):
        for tr, v in self.id_to_transcripts.iteritems():
            print "\tTranscript", tr, v.gff.attrs['gene_biotype']
            v.kprint()
    
    # ADD TRANSCRIPT SEQUENCES FROM FASTA FROM ENSEMBL
    # checks if the transcript length defined by the class is equal to the transcript sequence length from ENSEMBL
    def add_transcript_sequences(self, cdna):
        added = 0
        not_added = 0
        for tid, t in self.id_to_transcripts.iteritems():
            if not tid in cdna:
                self.transcripts_cdna_status[tid] = 0
                not_added += 1
            else:
                tlen = t.set_exon_transcript_coordinates()
                seq = cdna[tid]
                # check if the transcriptclass total exon length is equal to the fasta sequence length for this transcript
                assert tlen == len(seq)
                t.seq = seq
                self.transcripts_cdna_status[tid] = 1
                added += 1
        return added, not_added

    def design_mips_for_transcript(self, transcript_id, design_options = {}):
        # target variant is a row in dataframe created by get_transcripts_for_target_SNPs 
        if self.transcripts_cdna_status[transcript_id] == 0:
            print "Skipping candidates for transcript %s, because no cDNA was provided for this transcript." % transcript_id
            return None
        transcript = self.id_to_transcripts[transcript_id]
        assert transcript.gff.strand == self.gff.strand
        transcript_seq = transcript.seq.seq

        variants = []
        exons = []

        # Approach:
        # 1. seq: get full transcript sequence
        # 2. list: annotate exclude variants wrt transcript, check ref seq
        # 3. list: annotate target variants wrt transcript, check ref seq
        # 4. c-code to make list of all possible smMIPS
        # 5. returns a list of target variants that are overlapping with target sequence

        transcript_exclude_variants_list = []
        transcript_target_variants_list = []
        exon_starts = []


        for exnum in sorted(transcript.number_to_exon.keys()):
            exon = transcript.number_to_exon[exnum]
            eid = exon.getid()

            tstart = exon.tstart
            
            # QC
            if exnum > 1:
                if self.gff.strand == '+':
                    assert exon.gff.start >= transcript.number_to_exon[exnum-1].gff.end
                else:
                    assert exon.gff.end <= transcript.number_to_exon[exnum-1].gff.start
            
            # set exon boundaries
            if exnum > 1:
                exon_starts.append(exon.tstart)
                

            # annotate exon variants wrt to transcript
            if eid in self.genecollection.exonid_to_exclude_variants:
                for var in self.genecollection.exonid_to_exclude_variants[eid]:
                    var_idx = var[0]
                    var_transcript_start = int(exon.tstart + var[1])
                    var_transcript_end = int(exon.tstart + var[2])
                    
                    # check if reference seq of variant is consistent with transcript sequence 
                    if var_transcript_start == var_transcript_end and self.gff.strand == '+':
                        assert transcript_seq[var_transcript_start] == self.genecollection.exclude_variants['variants'][var_idx].ref


                    transcript_exclude_variants_list.append(var_idx)
                    transcript_exclude_variants_list.append(var_transcript_start)
                    transcript_exclude_variants_list.append(var_transcript_end)
            
            if eid in self.genecollection.exonid_to_target_variants: 
                for var in self.genecollection.exonid_to_target_variants[eid]:
                    var_idx = var[0]
                    var_transcript_start = int(exon.tstart + var[1])
                    var_transcript_end = int(exon.tstart + var[2])
                    
                    if var_transcript_start == var_transcript_end and self.gff.strand == '+':
                        assert transcript_seq[var_transcript_start] == self.genecollection.target_variants['variants'][var_idx].ref

                    transcript_target_variants_list.append(var_idx)
                    transcript_target_variants_list.append(var_transcript_start)
                    transcript_target_variants_list.append(var_transcript_end)

        num_exclude_variants = len(transcript_exclude_variants_list)/3
        num_target_variants = len(transcript_target_variants_list)/3
        len_mip_target_seq = design_options['len_mip_target_seq']
        num_exon_starts = len(exon_starts)

        code = """
        std::string tseq = transcript_seq;
        int tlen = (int) tseq.size();
        int mip_len = int(len_mip_target_seq) + 40;
        std::map<int, std::map<int,int> > position_to_exclude_variant_index;
        std::map<int, std::map<int,int> > position_to_target_variant_index;                        
        std::map<int, int> exon_start_pos;
        // std::cout << "tlen:" << int(tlen) << " mip_len: " << mip_len << std::endl;
       
        py::list ret; // overall result list

        for (int i = 0; i < num_exon_starts; i++) exon_start_pos[exon_starts[i]] = 1;

        for (int i = 0; i < num_exclude_variants; i++) {
           int var_idx = transcript_exclude_variants_list[i*3]; 
           int var_tstart = transcript_exclude_variants_list[i*3+1];  
           int var_tend = transcript_exclude_variants_list[i*3+2];
           for (int tpos = var_tstart; tpos <= var_tend; tpos++) {
                position_to_exclude_variant_index[tpos][var_idx] = 1;
           }  
        }
        for (int i = 0; i < num_target_variants; i++) {
           int var_idx = transcript_target_variants_list[i*3]; 
           int var_tstart = transcript_target_variants_list[i*3+1];  
           int var_tend = transcript_target_variants_list[i*3+2];
           for (int tpos = var_tstart; tpos <= var_tend; tpos++) {
                position_to_target_variant_index[tpos][var_idx] = 1;
           }  
        }
        std::map<int,int> ext_probe_exclude_variants;
        std::map<int,int> lig_probe_exclude_variants;
        std::map<int,int> covered_target_variants;
        for (int mip_start = 0; mip_start < tlen; mip_start++) {
            // std::cout << " " << int(mip_start) << std::endl;
            if (mip_start + mip_len >= tlen) break;
            for (int ext_probe_len = 16; ext_probe_len<25; ext_probe_len++) {
                int ext_probe_spans_boundary = 0;
                int lig_probe_spans_boundary = 0;
                int target_seq_spans_boundary = 0;
                // check ext probe exclude variants
                for (int probe_pos = 0 ; probe_pos < ext_probe_len; probe_pos++) {
                    int pos = mip_start + probe_pos;
                    if (pos>0 && exon_start_pos.find(pos) != exon_start_pos.end())
                        ext_probe_spans_boundary = 1;
                    std::map<int, std::map<int,int> >::const_iterator it = position_to_exclude_variant_index.find(pos);
                    if (it != position_to_exclude_variant_index.end()) {
                        for (std::map<int,int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                            ext_probe_exclude_variants[it2->first]=1;
                        }
                    }
                }

                int lig_probe_len = 40 - ext_probe_len;
                int transcript_target_start = mip_start + ext_probe_len;
                int transcript_target_end = transcript_target_start+int(len_mip_target_seq) - 1;
                // std::cout << " target_start: " << transcript_target_start << " " << transcript_target_end << std::endl; 
                
                // check ligation probe exclude variants
                for (int probe_pos = 0; probe_pos < lig_probe_len; probe_pos++) {
                    int pos = transcript_target_end+1 + probe_pos;
                    if (pos>0 && exon_start_pos.find(pos) != exon_start_pos.end())
                        lig_probe_spans_boundary = 1;
                    
                    std::map<int, std::map<int,int> >::const_iterator it = position_to_exclude_variant_index.find(pos);
                    if (it != position_to_exclude_variant_index.end()) {
                        for (std::map<int,int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                            lig_probe_exclude_variants[it2->first]=1;
                        }
                    }
                }
                // check target variants that are covered
                for (int pos = transcript_target_start ; pos <= transcript_target_end; pos++) {
                    if (pos>0 && exon_start_pos.find(pos) != exon_start_pos.end()) target_seq_spans_boundary=1;
                    
                    std::map<int, std::map<int,int> >::const_iterator it = position_to_target_variant_index.find(pos);
                    if (it != position_to_target_variant_index.end()) {
                       // std::cout << "  found target variant @ pos " << pos << std::endl;
                       for (std::map<int,int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                            // check if this variant does not also overlap the extension probe    
                            if (pos == transcript_target_start && position_to_target_variant_index[pos-1].find(it2->first) != position_to_target_variant_index[pos-1].end()) continue;
                            
                            if (pos == transcript_target_end && position_to_target_variant_index[pos+1].find(it2->first) != position_to_target_variant_index[pos+1].end()) continue;

                            covered_target_variants[it2->first] = 1;
                        }
                    }
                }
                // std::cout << "   covered_target_variants.size(): " << covered_target_variants.size() << " " << ext_probe_exclude_variants.size() << " " << lig_probe_exclude_variants.size()  <<  std::endl;
                
                if (ext_probe_exclude_variants.size() == 0 && lig_probe_exclude_variants.size() == 0) {
                    py::list mip;
                    mip.append(mip_start); // start of mip and extension probe
                    mip.append(transcript_target_start); // start of target sequence
                    mip.append(transcript_target_end+1); // start of ligation probe
                    mip.append(mip_start+mip_len);
                    mip.append(ext_probe_spans_boundary);
                    mip.append(lig_probe_spans_boundary);
                    mip.append(target_seq_spans_boundary);
                    py::list l_target_vars;
                    for (std::map<int,int>::const_iterator it = covered_target_variants.begin(); it != covered_target_variants.end(); it++) l_target_vars.append(it->first);
                    mip.append(l_target_vars);

                    ret.append(mip);
                }

                /*
                For the moment we exclude mips for which the probes overlap one of the exclude variants
                py::list l_ext_vars, l_lig_vars, l_target_vars; 
                
                for (std::map<int,int>::const_iterator it = ext_probe_exclude_variants.begin(); it != ext_probe_exclude_variants.end(); it++) l_ext_vars.append(int(it->first));
                for (std::map<int,int>::const_iterator it = lig_probe_exclude_variants.begin(); it != lig_probe_exclude_variants.end(); it++) l_lig_vars.append(it->first);


                mip.append(l_ext_vars);
                mip.append(l_lig_vars);
                */
            }
            return_val = ret;

        }
        // std::cout << std::endl;
        """
        res = scipy.weave.inline(code, ['transcript_exclude_variants_list','transcript_target_variants_list', 'num_exclude_variants','num_target_variants','transcript_seq','len_mip_target_seq','exon_starts','num_exon_starts'], compiler = 'gcc', headers = ["<vector>","<map>"], extra_compile_args = ['-O3'])
        mips = []
        if res is not None:
            for mip in res:
                #if len(mip[7])!=0:
                #    1/0
                nmip = (transcript_seq[mip[0]:mip[1]] +"/"+transcript_seq[mip[2]:mip[3]], mip[4],mip[5],mip[6], mip[7], transcript_seq[mip[1]:mip[2]], mip[1], mip[2])
                mips.append(nmip)
            return mips
    def design_mips_for_all_transcripts(self, design_options = {}):
        all_mips = {}
        
        # first design mips for all transcripts
        for transcript_id in self.id_to_transcripts.keys():
            mips = self.design_mips_for_transcript(transcript_id, design_options = design_options)
            if mips is None:
                continue
            for mip in mips:
                mip_ext_lig = mip[0]
                if mip_ext_lig not in all_mips:
                    all_mips[mip_ext_lig] = {}
                target_vars = ".vars/" + "/".join([str(var_idx) for var_idx in mip[4]])
                mip_description = "%s.%d.%d.%d%s" % (mip[5], mip[1],mip[2],mip[3], target_vars)
    
                if mip_description not in all_mips[mip_ext_lig]:
                    all_mips[mip_ext_lig][mip_description] = {}
                store_id = "%s:%d-%d" % (transcript_id, mip[6], mip[7])
                all_mips[mip_ext_lig][mip_description][store_id] = 1
        return all_mips 


def get_mip_oligosequence(ext, lig, design_options):
    umi_ext = 'N' * design_options['length_UMI_extension']
    umi_lig = 'N' * design_options['length_UMI_ligation']

    mip_seq = lig.lower() + umi_lig + design_options['backbone_sequence'].upper() + umi_ext + ext.lower()
    return mip_seq


def vcf_contains_indel(v):
    ref = v.ref
    min_dlen = 0
    max_dlen = 0
    for alt in v.alts:
        dlen = len(alt) - len(ref)
        if dlen < 0:
            if dlen < min_dlen:
                min_dlen = dlen
        elif dlen > 0:
            if dlen > max_dlen: 
                max_dlen = dlen
    return (min_dlen != 0 | max_dlen !=0), min_dlen, max_dlen

def vcf_get_unique_id(v):
    return ":".join([str(v.chrom),str(v.id), str(v.pos), str(v.ref), str(",".join(v.alts))])

class GeneCollection:
    def __init__(self):
        self.geneid_to_gene = {}
        self.transcriptid_to_geneid = {}
        self.exonid_to_exclude_variants = {}
        self.exonid_to_target_variants = {}
        self.geneid_to_name = {}
        self.name_to_geneid = {}
        self.exclude_variants = {}
        self.target_variants = {}

        self.exclude_variants['variantid_to_idx'] = {}
        self.exclude_variants['variants'] = []
        self.exclude_variants['variant_ids'] = []
 
        self.target_variants['variantid_to_idx'] = {}
        self.target_variants['variants'] = []
        self.target_variants['variant_ids'] = []
 
  
    def get_all_ensembl_genes(self):
        return self.geneid_to_gene.keys()[:]

    def get_all_gene_names(self):
        return self.name_to_geneid.keys()[:]

    def get_ensembl_ids(self, geneNames):
        ensembl_ids = []
        for name in geneNames:
            try:
                ensembl_id = self.name_to_geneid[name]
                ensembl_ids.append(ensembl_id)
            except KeyError:
                sys.stderr.write("Could not find ensembl ID corresponding to gene %s\n" % name)
        return ensembl_ids
        
    def size(self):
        return len(self.geneid_to_gene)
    
    def add_gene(self, gff):
        assert gff.fields[2] == 'gene'
        gene_id  = gff.attrs['gene_id']
        assert gene_id.startswith('ENSG')
        assert gene_id not in self.geneid_to_gene
        gene = Gene(gff, self)
        self.geneid_to_gene[gene_id] = gene

        # update name lists
        gene_name = gff.attrs['Name']
        self.geneid_to_name[gene_id] = gene_name
        if gene_name in self.name_to_geneid:
            sys.stderr.write("Warning: %s has multiple Ensembl IDs: %s %s.\nUsing the first one.\n" % (gene_name, self.name_to_geneid[gene_name], gene_id))
        else:
            self.name_to_geneid[gene_name] = gene_id

    def add_transcript(self, gff):
        assert gff.fields[2] == 'transcript'
        transcript_id = gff.attrs['transcript_id']
        gene_id = gff.attrs['Parent'][gff.attrs['Parent'].index('ENSG'):]
        if transcript_id in self.transcriptid_to_geneid:
            assert self.transcriptid_to_geneid[transcript_id] == gene_id
        else:
            self.transcriptid_to_geneid[transcript_id] = gene_id
        if gene_id in  self.geneid_to_gene:
            gene = self.geneid_to_gene[gene_id]
            gene.add_transcript(gff)
            return 0
        else:
            return 1

 
    def add_exon_from_transcript(self, gff):
        assert gff.fields[2] == 'exon'
        transcript_id = gff.attrs['Parent'][gff.attrs['Parent'].index('ENST'):]
        if transcript_id in self.transcriptid_to_geneid:
            gene_id = self.transcriptid_to_geneid[transcript_id]
        else:
            # transcript is not in dictionary, for instance because it was a processed_transcript
            return 1
        if gene_id in self.geneid_to_gene:
            gene = self.geneid_to_gene[gene_id]
            gene.add_transcript_exon_gff(gff)
            return 0
        else:
            # in this case the gene was not in the collection.
            # This may happen when the gene is a pseudogene that was not therefore not added
            return 2
    
    def kprint(self):
        for g,v in self.geneid_to_gene.iteritems():
            print "Gene",g, v.gff.attrs['gene_name'], v.gff.strand
            v.kprint()
    def iteritems(self):
        return self.geneid_to_gene.iteritems()

    def add_transcript_sequences(self, cdna):
        added, not_added = 0,0
        for gid, g in self.iteritems():
            stats = g.add_transcript_sequences(cdna)
            added += stats[0]
            not_added += stats[1]
        return added, not_added

    def overlap_variants(self, vcf_file, exonid_to_variants, variants):
        # note that an exon may be present in multiple transcript in different ranks.
        vcf = pysam.VariantFile(vcf_file,'r')
        assert type(exonid_to_variants) == dict
        assert type(variants) == dict
        tot_overlapping_variants = 0
        for gid, g in self.iteritems():
            for eid,exlist in g.exonid_to_exon.iteritems():
                ex = exlist[0] # We only need the genomic coordinates of the exon. 
                               # The exon may be present in multiple transcripts, in that case len(exlist)>1
                nrec = 0
                assert ex.gff.start < ex.gff.end
                for vcfrec in vcf.fetch(ex.gff.chrom, ex.gff.start, ex.gff.end):
                                        
                    variantid = vcf_get_unique_id(vcfrec)
                    vidx = -1
                    if variantid not in variants['variantid_to_idx']:
                        vidx = len(variants['variants'])
                        variants['variantid_to_idx'][variantid] = vidx
                        variants['variants'].append(vcfrec)
                        variants['variant_ids'].append(variantid)
                    else:
                        vidx = variants['variantid_to_idx'][variantid]

                    # compute variant positions relative to start of exon (5')

                    has_indel, len_del, len_ins = vcf_contains_indel(vcfrec)
                    var_start = vcfrec.pos - ex.gff.start - 1
                    var_end = vcfrec.pos + abs(len_del) - ex.gff.start - 1
                    
                    #assert var_start >= 0
                    #assert var_end <= ex.gff.end - ex.gff.start                    
                    # in some cases genomic variant may partially overlap the exon.
                    if var_start < 0:
                        var_start = 0
                    if var_end >= ex.gff.end - ex.gff.start:
                        var_end = ex.gff.end - ex.gff.start
 
                    if nrec == 0:
                        exonid_to_variants[eid] = []
                    exonid_to_variants[eid].append([vidx, var_start, var_end])
                    nrec += 1
                tot_overlapping_variants += nrec
        return tot_overlapping_variants

    def overlap_variants_exclude_from_probes(self, vcf_file):
        stats = self.overlap_variants(vcf_file, self.exonid_to_exclude_variants, self.exclude_variants)
        return stats

    def overlap_target_variants(self, vcf_file):
        stats = self.overlap_variants(vcf_file, self.exonid_to_target_variants, self.target_variants)
        return stats
 
    def design_mips_for_genes(self, output_file, ensembl_genes = [], design_options = {}):
        with gzip.open(output_file, 'w') as f:
            f.write("unique_probe_id\text_probe_seq\tlig_probe_seq\ttarget_seq\tmip_sequence\text_probe_overlaps_exon_boundary\tlig_probe_overlaps_exon_boundary\ttarget_seq_overlaps_exon_boundary\ttarget_variants\tensembl_gene\tgene_name\tensembl_transcripts\n")
            num_genes = 0
            for geneid in ensembl_genes:
                gene_mip_count = 0
                try:
                    gene = self.geneid_to_gene[geneid]
                except KeyError:
                    print "Could not find gene",geneid
                    continue
                gene_name = self.geneid_to_name[geneid]
                if num_genes % 10 == 1:
                    print "\tdesigned for %d genes" % num_genes
                mips = gene.design_mips_for_all_transcripts(design_options = design_options)
                for mip_seq in mips:
                    ext, lig = mip_seq.split('/')
                    for description in mips[mip_seq]:
                        #des_trs.append(description+";" + ",".join(mips[mip_seq][description].keys()))
                        des_split = description.split('.')
                        target_seq = des_split[0]
                        overlap_exon = "\t".join(des_split[1:4])
                        variants = des_split[4].split('/')
                        if des_split[4]!='vars/':
                            var_str = ";".join([self.target_variants['variant_ids'][int(var_idx)] for var_idx in variants[1:]])
                        else:
                            var_str = "NoVariants"

                        transcripts = ",".join(mips[mip_seq][description].keys())
                        mip_oligo_seq = get_mip_oligosequence(ext, lig, design_options)
                        gene_mip_count += 1
                        unique_probe_id = "smmip_%s_%d" % (geneid, gene_mip_count)
                        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (unique_probe_id, ext, lig, target_seq, mip_oligo_seq, overlap_exon, var_str, geneid, gene_name, transcripts))
                num_genes += 1
            print "\tdesigned for %d genes. Done." % num_genes
def read_file_with_ids(filename):
    id_dict = {}
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if line == '':
                break
            x = line.split()
            for xid in x:
                id_dict[xid] = 1
    return id_dict.keys()

def write_file_with_ids(filename, list_ids):
    with open(filename, 'w') as f:
        for wid in list_ids:
            f.write("%s\n" % str(wid))

def output_candidate_mips_to_fastq(candidate_mips_file, fastq_file, filter_performance_score = False, min_score = -1.0):
    assert candidate_mips_file.endswith('.gz')
    mips = pd.read_table(candidate_mips_file, compression="gzip")
    do_filter = False
    if filter_performance_score:
        if 'shendure_logistic_score' in mips.columns:
            do_filter = True
        else:
            sys.stderr.write("WARNING: no MIPS performance score found in candidate mips file\n")

    with gzip.open(fastq_file,'w') as fq:
        for i in mips.index:
            if do_filter and float(mips.loc[i,"shendure_logistic_score"]) < min_score:
                continue
            seq = mips.loc[i,"ext_probe_seq"] + mips.loc[i, "target_seq"] + mips.loc[i,"lig_probe_seq"]
            if do_filter:
                fq.write("@%s_score_%1.2f\n" % (mips.loc[i, "unique_probe_id"], float(mips.loc[i,"shendure_logistic_score"])))
            else:
                fq.write("@%s\n" % (mips.loc[i, "unique_probe_id"]))
            fq.write(seq + "\n")
            fq.write('+\n')
            fq.write('A'*len(seq)+"\n")




def main():
    # TODO add transcript start, end.
    # make file that says which smMIPS cover a variants from the target variant file.
    #  run design_mips_for_all_transcripts.py --ensembl_gff gff3.gz --ensembl_cDNA cdna.fa.gz --output_file mips.txt --target_ensembl_file genes.txt --exclude_variants_vcf test/ExAC.r0.3.sites.vep.vcf.gz.AC_1000_SNPs.vcf.gz --target_variants_vcf test/ExAC.r0.3.sites.vep.vcf.gz.AC_1000_SNPs.vcf.gz
    print """\nThis program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"""

    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--remove_transcript_biotypes", help = "Comma-separated list of biotypes. No smMIPs will be designed for these transcripts", default = 'unprocessed_pseudogene,retained_intron,nonsense_mediated_decay')
    parser.add_argument("--exclude_variants_vcf", help = "VCF file (bgzip'ed and indexed with tabix) with variants that should not be present in the extension or ligation probe")
    parser.add_argument("--target_variants_vcf", help = "VCF file (bgzip'ed and indexed with tabix) with variants that should be covered by the target sequence of a cDNA-smMIP")         
    parser.add_argument("--target_seq_length", help = "Length of cDNA-smMIPs target sequence (sequence between extension and ligation probe)", default = 100)
    parser.add_argument("--calculate_performance_scores", help = "Compute performance scores for MIPS using MIPGEN 2.0 logistic model (copyright Evan Boyle, used by permission.) Commercial use of this option is not allowed. ", action = "store_true")
    parser.add_argument("--ref_genome", help="Fasta reference sequence (indexed by Samtools and BWA). Required for calculating MIP performance score.")
    parser.add_argument("--bwa", help = "Path to BWA binary.", default = "bwa")
    parser.add_argument("--output_candidate_mips_to_fastq", help = "Outputs MIP probe and target sequence to Fastq file for mapping with a splicing-aware mapper. Can be used to select cDNA-smMIPs.", action = "store_true")
    parser.add_argument("--fastq_minimum_performance_score", help ="If outputting candidate cDNA-smMIPs to Fastq, only output MIPs with performance score above this threshold. Only effective if calculate_performance_scores is also specified", default = 0.9)


    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--ensembl_gff", help="Gff3 file with Ensembl gene, transcript and exon definitions. May be gzipped", required = True)
    requiredNamed.add_argument("--ensembl_cDNA", help="Fasta file with cDNA of all Ensembl transcripts. May be gzipped", required = True)
    requiredNamed.add_argument("--output_file_prefix", help ="Prefix of names output files. Output files will be gzipped.", required = True) 
    requiredNamed.add_argument("--target_ensembl_file", help = "File with Ensembl gene IDs (ENSG..) of genes for which cDNA-smMIPs should be designed", required = True)
    requiredNamed.add_argument("--length_UMI_extension", help = "Number of Ns (length of UMI) on extension probe. Choose value between 0 and 10. Total length of UMI sequence should not exceed 10 nt.", required = True)
    requiredNamed.add_argument("--length_UMI_ligation", help = "Number of Ns (length of UMI) on ligation probe. Choose value between 0 and 10. Total length of UMI sequence should not exceed 10 nt.", required = True)
 
   

    args = parser.parse_args()

    # check conditional arguments
    if args.calculate_performance_scores:
        if args.ref_genome is None:
            sys.stderr.write("Need to specify a genome reference file (Fasta) for calculation of smMIP performance scores.\n")
            sys.exit(1)

    if args.calculate_performance_scores:
        print """\nThe code to compute MIP performance scores is based on
        MIPGEN2.0 (Boyle et al, Bioinformatics, 18(30):2670-2, 2014).
        Copyright Evan Boyle.
        COMMERCIAL USE IS NOT PERMITTED.
        Please visit http://shendurelab.github.io/MIPGEN/ for details.
        \n"""
    
    candidate_mip_file = args.output_file_prefix+".without_scores.txt.gz"
    candidate_mip_file_with_scores = args.output_file_prefix + ".with_scores.txt.gz"
    fastq_file_without_scores = args.output_file_prefix+".without_scores.fastq.gz"
    fastq_file_with_scores = args.output_file_prefix+".filtered_by_score.fastq.gz"
    output_file_designed_genes = args.output_file_prefix + ".designed_genes.txt"

    design_options = {'len_mip_target_seq':int(args.target_seq_length),
                      'backbone_sequence' : 'CTTCAGCTTCCCGATATCCGACGGTAGTGT', # DO NOT CHANGE unless you really know what you are doing.
                      'length_UMI_extension' : int(args.length_UMI_extension),
                      'length_UMI_ligation' : int(args.length_UMI_ligation) }
    assert design_options['length_UMI_extension'] >= 0 and design_options['length_UMI_extension'] < 16
    assert design_options['length_UMI_ligation'] >= 0 and design_options['length_UMI_ligation'] < 16

    if False: 
        data_dir = "./" # "/Volumes/Kees_SD_1/MacbookAir/data/GenomeAnnotation/"
        fname_ensembl_gff = data_dir + "ftp.ensembl.org.pub.grch37.release-84.gff3.homo_sapiens.Homo_sapiens.GRCh37.82.gff3.gz"
        fname_cdna = data_dir + "ftp.ensembl.org.pub.grch37.release-84.fasta.homo_sapiens.cdna.Homo_sapiens.GRCh37.cdna.all.fa.gz"
        #fname_prefix_output_files = "mips.txt.gz"
        #fname_candidate_mips = "%s.target_variants.transcripts.candidate_mips.txt" % (fname_prefix_output_files)


        remove_transcript_biotypes = ['unprocessed_pseudogene','retained_intron','nonsense_mediated_decay'] 

        ref_genome_fasta = '/home/kalbers/ref/hg19.fa'
        candidate_mip_file = "mips.txt.gz"
        candidate_mip_file_with_scores = candidate_mip_file + ".with_scores.txt.gz"
    else:
        remove_transcript_biotypes = args.remove_transcript_biotypes.split(',')
        print "Transcript with these biotypes will be excluded:"," ".join(remove_transcript_biotypes)

 
    test_transcript = 'ENST00000508384'
    #fname_cdna += "." + test_transcript + ".fa"
    if True:
        sys.stderr.write("Loading transcripts from %s\n" % args.ensembl_gff)
        genes_gff = pbt.BedTool(args.ensembl_gff)
 
        # add genes and transcripts
        print "Reading gene, transcript and exon definitions from", args.ensembl_gff
        gc = GeneCollection()
        stats = {'transcript':{ 0: 0, 1: 0}}
        for gff in genes_gff:
            if gff.fields[2] == "gene":
                if gc.size() % 100 == 1:
                    print "\tRead %d genes" % gc.size()
                #if gc.size() == 10:
                #    break 
                gc.add_gene(gff)
                   
            elif gff.fields[2] == "transcript" and gff.attrs['biotype'] not in remove_transcript_biotypes:
                res = gc.add_transcript(gff)
                stats['transcript'][res] += 1 
            elif gff.fields[2] == "exon":
                gc.add_exon_from_transcript(gff)    
        # adding variants from VCF
        if args.exclude_variants_vcf is not None:
            print "Overlapping exons with variants to be excluded from probes from", args.exclude_variants_vcf
            num = gc.overlap_variants_exclude_from_probes(args.exclude_variants_vcf)
            print "\t%d variants overlapped"% num
        else:
            print "No VCF file with variants to exclude from probes specified. Are you sure?"

        if args.target_variants_vcf is not None:
            print "Overlapping exons with target variants from", args.target_variants_vcf
            num = gc.overlap_target_variants(args.target_variants_vcf)
            print "\t%d variants overlapped"% num
        else:
            print "No VCF file with target variants specified."


        # LOAD ENSEMBL cDNA for ALL TRANSCRIPTS
        print "LOADING cDNA DATA from",args.ensembl_cDNA
        cdna = dict( (s.name.split('.')[0], s) for s in HTSeq.FastaReader(args.ensembl_cDNA) )
        print "Adding to GeneCollection"
        stats = gc.add_transcript_sequences(cdna)
        print "Added %d transcripts, could not find cDNA for %d transcripts" % (stats[0],stats[1]) 
        

        if args.target_ensembl_file is None:
            ensembl_genes = gc.get_all_ensembl_genes()
        else:
            ensembl_genes = []
            if args.target_ensembl_file is not None:
                ensembl_genes.extend(read_file_with_ids(args.target_ensembl_file))
                print "Read %d ensembl gene ids from %s" % (len(ensembl_genes), args.target_ensembl_file)
            
            #if args.target_genes_file is not None:
            #    gene_names = read_file_with_ids(args.target_genes_file)
            #    ensembl_genes_2 = gc.get_ensembl_ids(gene_names)
            #    print "Read %d gene names, resulting in %d ensembl gene ids, from %s" % (len(gene_names), len(ensembl_genes_2), args.target_genes_file)
            #    ensembl_genes.extend(ensembl_genes_2)

        print "Designing cDNA-smMIPs for %d genes" % len(ensembl_genes)
        gc.design_mips_for_genes(candidate_mip_file, ensembl_genes = ensembl_genes, design_options = design_options)
        write_file_with_ids(output_file_designed_genes, ["%s\t%s" % (e,gc.geneid_to_name[e]) for e in ensembl_genes])

    if args.calculate_performance_scores:
        score_args = {'candidate_mip_file' : candidate_mip_file,
                      'ref_genome_fasta' :  args.ref_genome,
                      'output_file' : candidate_mip_file_with_scores,
                      'output_prefix' : candidate_mip_file_with_scores + ".temporary_file",
                      'remove_fastq_sam' : True,
                      'bwa' : args.bwa}
        scores = Shendure_scores(score_args)
        scores.compute_scores()
        filter_performance_score = True
        candidate_mips_file_for_fastq = candidate_mip_file_with_scores
        fastq_file = fastq_file_with_scores 
    else:
        candidate_mips_file_for_fastq = candidate_mip_file
        filter_performance_score = False
        fastq_file = fastq_file_without_scores

    if args.output_candidate_mips_to_fastq:
        print "Writing cDNA-smMIPs to Fastq file",fastq_file
        if filter_performance_score:
            print "\trequiring MIPS performance score > ",args.fastq_minimum_performance_score
            print "\tMIPS with extension and ligation copy count >= 100 are given score of zero."
        output_candidate_mips_to_fastq(candidate_mips_file_for_fastq, fastq_file, filter_performance_score = filter_performance_score, min_score = args.fastq_minimum_performance_score)

if __name__ == "__main__":
    main()
