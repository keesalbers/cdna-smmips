import pysam
import pandas as pd
import numpy as np
import argparse
import scipy.weave
import sys
from scipy.weave import converters
import matplotlib.pyplot as plt

def reverse_complement(kmer):
    rc = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    res = ""
    l = len(kmer)
    for i in range(l):
        res = res + rc[kmer[l-i-1]]
    return res

class Fasta(object):
    def __init__(self, fname, is_cdna = False):
       self.fname = fname
       self.ff = pysam.FastaFile(fname)
       if is_cdna:
           self.set_gene_to_transcript_dict()
    
    def set_gene_to_transcript_dict(self):
       print "Making gene to transcript dictionary"
       self.transcript_to_annot = {}
       self.ensembl_gene_to_transcripts = {}
       f = open(self.fname,'r')
       for line in f:
            if line.startswith('>'):
                d = line[1:].split()
                transcript = d[0]
                annot = d[2]
                if not "protein_coding" in d[5]:
                    continue
                assert d[3].startswith('gene:')
                ensembl_gene = d[3].split(':')[1]
                assert transcript not in self.transcript_to_annot
                self.transcript_to_annot[transcript] = annot

                if not ensembl_gene in self.ensembl_gene_to_transcripts:
                    self.ensembl_gene_to_transcripts[ensembl_gene] = [transcript]
                else:
                    self.ensembl_gene_to_transcripts[ensembl_gene].append(transcript)
       f.close()

class MIPS(object):
    def __init__(self, mip_design_file = ""):
        self.mips_table = pd.read_table(mip_design_file, sep = "\t")
        self.variant_column = "variant"
        
        self.required = ["ext_probe_sequence", "lig_probe_sequence","unique_probe_id",self.variant_column]
        for req in self.required:
            if not req in self.mips_table.columns:
                sys.stderr.write("Required columns in MIPS design file: " + " ".join(self.required) + "\n")
                raise "ERROR"

        self.unique_probe_ids = list(self.mips_table['unique_probe_id'])
        assert len(self.unique_probe_ids) == len(set(self.unique_probe_ids)), "Probe IDs are not unique!"
        self.setup_probe_seq()
        self.setup_ase_variants()

    def setup_probe_seq(self):
        self.ext_probe_sequence =list(self.mips_table["ext_probe_sequence"])
        self.lig_probe_sequence =list(self.mips_table["lig_probe_sequence"])
        self.lig_probe_sequence_rc = []
        self.ext_probe_sequence_rc = []
        for seq in self.ext_probe_sequence:
            self.ext_probe_sequence_rc.append(reverse_complement(seq))
        for seq in self.lig_probe_sequence:
            self.lig_probe_sequence_rc.append(reverse_complement(seq))

    def setup_ase_variants(self):
        vdict = {}
        self.variants_ase = []
        for i in self.mips_table.index:
            varstr = self.mips_table[self.variant_column].loc[i]
            if not varstr in vdict:
                var = varstr.split('_')
                variant = {'CHROM': "chr" + var[1], 'START':int(var[2]), 'END':int(var[2]), 'REF':var[3], 'ALT':var[4] }
                vdict[varstr] = 1
                self.variants_ase.append(variant)
   
    def get_ase_variants(self):
        return self.variants_ase[:]
            
    def match_read_pairs_to_probe_dict(self, list_forward_reads, list_reverse_reads):
        """
        Returns dictionary of index in probe_list -> index in forward/reverse read list
        probe_index = self.match_read_pairs_to_probe_index(forward_reads, reverse_reads)
        NOTE: reverses read sequence to original as in the FASTQ for the matching
        """
        L = 76
        assert len(list_forward_reads) == len(list_reverse_reads)
        # get ligation probe sequence from reads
        list_forward_reads_seq = []
        for i, read in enumerate(list_forward_reads):
            assert read.is_read1
            if read.is_reverse:
                list_forward_reads_seq.append(reverse_complement(read.seq[-L:]))
            else:
                list_forward_reads_seq.append(read.seq[:L])

        list_reverse_reads_seq = []
        for i, read in enumerate(list_reverse_reads):
            assert read.is_read2
            if read.is_reverse:
                list_reverse_reads_seq.append(reverse_complement(read.seq[-L:]))
            else:
                list_reverse_reads_seq.append((read.seq[:L]))


        match_probe = self.match_read_pairs_to_probe_index_check_dimer(list_forward_reads_seq, list_reverse_reads_seq)
        assert len(match_probe) == len(list_forward_reads_seq), "Internal error"
        sort_by_probe = {}
        for idx, probe_index in enumerate(match_probe):
            try:
                sort_by_probe[probe_index].append(idx)
            except KeyError:
                sort_by_probe[probe_index] = [idx]

        if True:
            for k,v in sort_by_probe.iteritems():
                print "Probe_index %d: num matched reads: %d" % (k, len(v))

        return sort_by_probe

    def match_read_pairs_to_probe_index(self, forward_reads, reverse_reads):
        """
        Performs exhaustive matching against all probes.
        Returns numpy array with probe index.
        """
        ext_probe_sequence = self.ext_probe_sequence
        lig_probe_sequence_rc = self.lig_probe_sequence_rc
        assert len(forward_reads) == len(reverse_reads), "Internal error"
        assert len(ext_probe_sequence) == len(lig_probe_sequence_rc)
        # numpy array to store the index of the mip probe
        match_probe = np.zeros( (len(forward_reads)), dtype = 'int') - 1
        max_num_mismatches = 1
 

        code = """
        std::vector<std::string> lig_probe_seq, ext_probe_seq;
        for (int p = 0; p < ext_probe_sequence.length(); p++) {
            std::string lseq = lig_probe_sequence_rc[p];
            std::string eseq = ext_probe_sequence[p];
            lig_probe_seq.push_back(lseq);
            ext_probe_seq.push_back(eseq);
        }
        for (int r = 0; r < forward_reads.length(); r++) {
            std::string fseq = forward_reads[r];
            std::string rseq = reverse_reads[r];
            
            int best_match_probe = -1;
            int best_match = 1000;
            for (int p = 0; p < lig_probe_seq.size(); p++) {
                int nm_lig = 0; // mismatches ligation probe
                int nm_ext = 0; // mismatched extension probe
                std::string & lig = lig_probe_seq[p];
                std::string & ext = ext_probe_seq[p];
                for (size_t b = 0; b < lig.size(); b++) {
                    if (fseq[b] != lig[b]) nm_lig++;
                }
                int L = (int) rseq.size();
                for (size_t b = 0; b < ext.size(); b++) {
                    if (rseq[b] != ext[b]) nm_ext++; 
                }
                if (nm_lig <= max_num_mismatches && nm_ext <= max_num_mismatches) {
                    // std::cout << "r: " << r << " p: " << p << " nm: " << nm_lig << "\t" << nm_ext << std::endl;

                    if (nm_ext + nm_lig < best_match) {
                        best_match = nm_ext + nm_lig;
                        best_match_probe = p;
                    }
                }
            }
            match_probe[r] = best_match_probe;
        }
        """
        scipy.weave.inline(code, ['ext_probe_sequence','lig_probe_sequence_rc','max_num_mismatches','forward_reads','reverse_reads', 'match_probe'], compiler = 'gcc', headers = ["<vector>"])
        return match_probe

    def match_read_pairs_to_probe_index_check_dimer(self, forward_reads, reverse_reads):
        """
        Performs exhaustive matching against all probes.
        Returns numpy array with probe index.
        """
        ext_probe_sequence = self.ext_probe_sequence
        lig_probe_sequence_rc = self.lig_probe_sequence_rc
        ext_probe_sequence_rc = self.ext_probe_sequence_rc
        lig_probe_sequence = self.lig_probe_sequence
        assert len(forward_reads) == len(reverse_reads), "Internal error"
        assert len(ext_probe_sequence) == len(lig_probe_sequence_rc)
        # numpy array to store the index of the mip probe
        match_probe = np.zeros( (len(forward_reads)), dtype = 'int') - 1
        max_num_mismatches = 1
        num_dimers = -1; 

        code = """
        num_dimers = 0;
        std::vector<std::string> lig_probe_seq, ext_probe_seq, lig_probe_seq_rc, ext_probe_seq_rc;
        for (int p = 0; p < ext_probe_sequence.length(); p++) {
            std::string lseq = lig_probe_sequence[p];
            std::string eseq = ext_probe_sequence[p];
            lig_probe_seq.push_back(lseq);
            ext_probe_seq.push_back(eseq);

            std::string lseq2 = lig_probe_sequence_rc[p];
            std::string eseq2 = ext_probe_sequence_rc[p];
            lig_probe_seq_rc.push_back(lseq2);
            ext_probe_seq_rc.push_back(eseq2);
        }

        //for (int r = 0; r < 100; r++) { // forward_reads.length(); r++) {
        for (int r = 0; r < forward_reads.length(); r++) {
            std::string fseq = forward_reads[r];
            std::string rseq = reverse_reads[r];
            
            int best_match_probe = -1;
            int best_match = 1000;
            for (int p = 0; p < lig_probe_seq.size(); p++) {
                int nm_lig = 0; // mismatches ligation probe
                int nm_ext = 0; // mismatched extension probe
                std::string & lig = lig_probe_seq_rc[p];
                std::string & ext = ext_probe_seq[p];
                for (size_t b = 0; b < lig.size(); b++) {
                    if (fseq[b] != lig[b]) nm_lig++;
                }
                int L = (int) rseq.size();
                for (size_t b = 0; b < ext.size(); b++) {
                    if (rseq[b] != ext[b]) nm_ext++; 
                }
                if (nm_lig <= max_num_mismatches && nm_ext <= max_num_mismatches) {
                    // std::cout << "r: " << r << " p: " << p << " nm: " << nm_lig << "\t" << nm_ext << std::endl;

                    if (nm_ext + nm_lig < best_match) {
                        best_match = nm_ext + nm_lig;
                        best_match_probe = p;
                    }
                }
            }
            match_probe[r] = best_match_probe;
            
            if (best_match_probe != -1) { 
                bool hit = false;
                int kmer = 15; // look for substrings of 10 probe nt in read and calculate distance
                int ext_seq_index_in_ligation_read = 1000000; // will contain first match.
                int lig_seq_index_in_extension_read = 1000000; // will contain first match
                
                int lenUmiLigation = 0;    // in BAM file UMIs have been trimmed, so these sequences will not contain any UMIs.
                int lenUmiExtension = 0;   // in BAM file UMIs have been trimmed.

                const std::string & ligationRead = fseq;
                const std::string & extensionRead = rseq;
                const std::string & _lig_probe_seq = lig_probe_seq[best_match_probe];
                const std::string & _ext_probe_seq = ext_probe_seq[best_match_probe];
                const std::string & _lig_probe_seq_rc = lig_probe_seq_rc[best_match_probe];
                const std::string & _ext_probe_seq_rc = ext_probe_seq_rc[best_match_probe];
                // std::cout << std::endl << _ext_probe_seq.length()  << " " <<  _lig_probe_seq.length() << " " << _ext_probe_seq_rc.length() << " " << _lig_probe_seq_rc.length() << std::endl;
                // std::cout << fseq << " " << rseq << std::endl;
                for (int ext_probe_start = 0; ext_probe_start < _ext_probe_seq.length() - kmer; ext_probe_start++) {
                    
                    for (int index = 0; index < ligationRead.length() - _ext_probe_seq_rc.length(); index++) {
                        int nm = 0;
                        // countMismatches(_ext_probe_seq_rc, ligationRead.substring(index, index + kmer));
                        for (int i = 0; i < (int) _ext_probe_seq_rc.size(); i++)  {
                            if (_ext_probe_seq_rc[ext_probe_start+i] != ligationRead[index+i]) nm++;
                        }

                        if (nm <= 1) {
                            hit = true;
                            ext_seq_index_in_ligation_read = index - (_ext_probe_seq.length()-kmer-ext_probe_start) - lenUmiLigation - _lig_probe_seq.length();
                            break;                        
                        }
                        // std::cout << "nm: " << nm << std::endl;
                    }
                }
                
                // check ligation probe in extension read        
                for (int lig_probe_start = 0; lig_probe_start < _lig_probe_seq.length() - kmer; lig_probe_start++) {
                    //String lig_probe_seq = this.ligationProbeSequence.substring(lig_probe_start, lig_probe_start + kmer);
                    
                    for (int index = 0; index < extensionRead.length() - kmer; index++) {
                        //int nm = SeqUtil.countMismatches(lig_probe_seq, extensionRead.substring(index, index + kmer));
                        int nm = 0;
                        for (int i = 0; i < (int) _lig_probe_seq.size(); ++i) {
                            if (_lig_probe_seq[lig_probe_start+i] != extensionRead[index+i]) nm++;
                        }
                        if (nm <= 1) {
                            hit = true;
                            lig_seq_index_in_extension_read = index - lig_probe_start - lenUmiExtension - _ext_probe_seq.length();
                            break;
                        }
                    }
                }
                // std::cout << "ext_seq_index_in_ligation_read: " << ext_seq_index_in_ligation_read << " lig_seq_index_in_extension_read: " << lig_seq_index_in_extension_read << std::endl;
                //score[0] = ext_seq_index_in_ligation_read;
                //score[1] = lig_seq_index_in_extension_read;
                if (hit) {
                    match_probe[r] = -2;
                    num_dimers++;
                }
            } // check dimer
        }
        std::cout << "Number of dimers detected: " << num_dimers << std::endl;

        """
        scipy.weave.inline(code, ['ext_probe_sequence','lig_probe_sequence','ext_probe_sequence_rc','lig_probe_sequence_rc', 'max_num_mismatches','forward_reads','reverse_reads', 'match_probe','num_dimers'], compiler = 'gcc', extra_compile_args = ['-O3'],headers = ["<vector>"])
        # print "Number of dimers detected:",num_dimers
        return match_probe

def get_mip_data(mip_design_file):
    mips = pd.read_table(mip_design_file, sep = "\t")
    
    for r in mips.index:
        ext_seq = mips.loc[r,"ext_probe_sequence"]
        lig_seq = mips.loc[r,"lig_probe_sequence"]
        probe_id = mips.loc[r,"probe_id"]
        scan_seq = mips.loc[r,"scan_target_seq"]
        1/0
    return mips_ercc

def count_num_cigar_match(read):
    count = 0
    for tup in read.cigartuples:
        if tup[0] == 0:
            count += tup[1]
    return count

def get_umi_histogram(read_forward, read_reverse, match_options):
    umi_hist = {}
    for qname in read_forward:
        if read_reverse.has_key(qname):
            umi = qname.split(':')[-1]
            if umi.count('N') == 0:
                umi_hist[umi] = umi_hist.get(umi,0) + 1
    return umi_hist

def make_umi_to_seq_dictionary(list_bamreads_forward = [], list_bamreads_reverse = []):
    """
    Makes umi from extension and ligation read.
    Only includes read pairs where the UMI does not contain an N.
    Returns matched dictionaries from umi to read sequence.
    """
    res_forward, res_reverse = {},{}
    for read_forward, read_reverse in zip(list_bamreads_forward, list_bamreads_reverse):
        seq_id = read_forward.qname.split(':')
        assert read_forward.qname == read_reverse.qname
        umi = seq_id[-2] + seq_id[-1]
        if umi.count('N') == 0:
            try:
                res_forward[umi].append(read_forward.seq)
                res_reverse[umi].append(read_reverse.seq)
            except KeyError:
                res_forward[umi] = [read_forward.seq]
                res_reverse[umi] = [read_reverse.seq]
    return res_forward, res_reverse

def get_umi_from_qname(qname):
    seq_id = qname.split(':')
    umi = seq_id[-2] + seq_id[-1]
    return umi

def get_dict_consensus_sequence(dict_seq = {}, consensus_options = {}):
    """
    Uses UMI sequence in sequence identifier to build consensus sequence for forward and reverse separately.
    NOTE: assumes that duplicate UMIs belong to the same original molecule.
          This will not be true if UMIs from different smMIPS are combined.
    """
    consensus_seq = {}
    for seq_id,list_seq in dict_seq.iteritems():
        consensus_seq[seq_id] = seqan_get_consensus_sequence_with_realignment(list_seq)

    return consensus_seq

def filter_read(read, target, match_options):
    """
    Returns true if read should be filtered out based on overlap with target sequence and 
    match_options parameters
    """
     # check if minimum read length is reached
    if read.is_read1 and read.query_length < match_options.min_query_length_forward:
        return True
    elif read.is_read2 and read.query_length < match_options.min_query_length_reverse:
        return True

    if target != []:
        fraction_overlap = float(read.get_overlap(target[1], target[2])) / float(read.query_length)
        if fraction_overlap < match_options.min_fraction_read_overlap_target:
            return True

        num_matched_nt = count_num_cigar_match(read)
        fraction_matched_nt = float(num_matched_nt) / float(read.query_length)

        if fraction_matched_nt < match_options.min_fraction_matched_nt:
            return True

    return False

def get_pair_mapped_as_list_seq(forward_reads, reverse_reads):
    """
    Takes dictionary of read_id -> bam_read and converts it to matched list
    of forward read seq, reverse read seq, and read_id.
    NOTE Removes reads with unmapped mate.
    """
    list_forward_reads_seq, list_reverse_reads_seq, list_read_id = [],[],[]
    for read_id, read in forward_reads.iteritems():
        if reverse_reads.has_key(read_id):
            list_forward_reads_seq.append(read.seq)
            list_reverse_reads_seq.append(reverse_reads[read_id].seq)
            list_read_id.append(read_id)
    return list_forward_reads_seq, list_reverse_reads_seq, list_read_id

def get_pair_mapped_as_list(forward_reads, reverse_reads):
    """
    Takes dictionary of read_id -> bam_read and converts it to matched list
    of forward read, reverse read.
    NOTE Removes reads with unmapped mate.
    """
    list_forward_reads, list_reverse_reads = [],[]
    for read_id, read in forward_reads.iteritems():
        if reverse_reads.has_key(read_id):
            list_forward_reads.append(read)
            list_reverse_reads.append(reverse_reads[read_id])
    return list_forward_reads, list_reverse_reads

def get_reads_for_target(target, bam, match_options, do_filter = False):
    """
    Returns dictionary of id -> bam_read for forward reads and reverse reads.
    Performs filtering using filter parameters in match_options
    """
    read_iter = bam.fetch(target[0], target[1], target[2])
    n = 0
    nall = 0
    forward_reads = {}
    reverse_reads = {}
    for read in read_iter:
        nall += 1
        if do_filter:
            if filter_read(read, target, match_options):
                continue
        if read.is_read1:
            forward_reads[read.qname] = read
        elif read.is_read2:
            reverse_reads[read.qname] = read
        else:
            1/0

        n += 1
    return forward_reads, reverse_reads

def get_all_reads_from_bamfile(bam, match_options):
    """
    Returns dictionary of id -> bam_read for forward reads and reverse reads.
    Performs filtering using filter parameters in match_options
    """
    n = 0
    nall = 0
    forward_reads = {}
    reverse_reads = {}
    for read in bam:
        nall += 1
        if filter_read(read, [], match_options):
            continue
        if read.is_read1:
            forward_reads[read.qname] = read
        elif read.is_read2:
            reverse_reads[read.qname] = read
        else:
            1/0

        n += 1
    return forward_reads, reverse_reads

def get_majority_vote(nuc_dict = {}):
    m = -1
    basecall = 'nocall'
    for nuc, count in nuc_dict.iteritems():
        if count > m:
            m = count
            basecall = nuc
    return basecall

def determine_readcount_threshold(readcount_to_numumis, threshold_fraction_data):
    if len(readcount_to_numumis) == 0:
        return 0
    readcounts = np.array(sorted(readcount_to_numumis.keys()))
    numumis = np.zeros( readcounts.shape )
    
    for idx, rc in enumerate(readcounts):
        numumis[idx] = readcount_to_numumis[rc]

    a = readcounts*numumis
    num_read_threshold = readcounts[a.shape[0]-1-np.where((a[::-1]).cumsum()/(1.0*a.sum())>threshold_fraction_data)[0][0]]
    return num_read_threshold


def check_umi_errors(pdict, base_error_rate = 0.006):
    umis = []
    counts = []
    for umi in pdict:
        umis.append(umi)
        c = pdict[umi]
        count = sum(c.values())
        counts.append(count)
    counts = np.array( counts )
    idx = np.argsort(counts)[::-1]
    counts = list(counts[idx])
    umis = [ umis[i] for i in idx]
    
    plt.hist(counts);


    A = max(counts)
    L = len(umis[0]) # length of UMI
    p_umi_error = 1.0 - scipy.stats.binom.cdf(0, n = L, p = base_error_rate) # prob at least one seq error in the UMI

    binom_table = np.zeros( (A+1, A+1) )
    for i in range(A+1):
        for j in range(A+1):
            binom_table[i,j] = 1-scipy.stats.binom.cdf(i,n = j, p = p_umi_error)


    code = """
    std::map < std::string, int> umi_to_index;
    for (size_t i = 0; i < umis.size(); i++) {
        umi_to_index[ umis[i] ] =(int) i;
    }
    const char nuc[] = "ACGT";
    int tot_match = 0;
    int single_match = 0;
    for (int i = (int) counts.size()-1; i >= 0; i--) {
        std::string umi = umis[i];
        int num_match = 0;
        std::vector<int> single_matches, double_matches;
        int single_count = 0;

        for (size_t p = 0; p < umi.size(); ++p) {
            const char nu = umi[p];
            for (int m = 0; m < 4; ++m) {
                const char nm = nuc[m];
                if (nm != nu) {
                    std::string mumi(umi);
                    mumi[p] = nm;
                    if (umi_to_index.find(mumi) != umi_to_index.end()) {
                        if (i != umi_to_index[mumi]) {
                            single_matches.push_back(umi_to_index[mumi]);
                            single_count += int(counts[ umi_to_index[mumi] ]);
                        }
                    }
                    for (size_t p2 = 0; p2 < umi.size(); ++p2) {

                        if (p2 != p) {
                            const char nu2 = umi[p2];
                            for (int m2 = 0; m2 < 4; ++m2) {
                                const char nm2 = nuc[m2];
                                if (nm2 != nu2) {
                                    std::string mumi2(mumi);
                                    mumi2[p2] = nm2;
                                    if (umi_to_index.find(mumi2) != umi_to_index.end()) {
                                        if (i != umi_to_index[mumi2]) {
                                            double_matches.push_back(umi_to_index[mumi2]);
                                        }
                                    }
    
                                }
                            }
                        }
                   }

                }
            
            }
        }
        num_match = single_matches.size() + double_matches.size();
        single_match += single_matches.size();
        tot_match += single_matches.size() ;
        if (single_matches.size() > 0 ) {
            std::cout << "umi[" << i << "]: " << umi << " num matches: " << single_matches.size() << " num_reads_single_mut: " << single_count << " double_matches: " << double_matches.size();
            std::cout << " counts[" << i << "]: " << int(counts[i]) << " counts single matches: ";
            for (size_t j = 0; j < single_matches.size(); ++j) {
                int c = counts[single_matches[j]];
                std::cout << " " << c;
            }
            std::cout << " idx single matches: ";
            for (size_t j = 0; j < single_matches.size(); ++j) {
                int c = single_matches[j];
                std::cout << " " << c;
            }

            std::cout << std::endl;
        }
    }
    std::cout << "Total number of matches: " << tot_match << std::endl;
    std::cout << "Single matches: " << single_match << std::endl;
    """
    scipy.weave.inline(code, ['umis','counts'], compiler = 'gcc',  type_converters = converters.blitz, extra_compile_args = ['-O3'],headers = ["<vector>","<map>"])


    1/0


def count_umis(mips, bamfile, match_options):
    bamfiles = [f for f in bamfile.split(',') if f not in [',','']]
    bam = {}
    for bf in set(bamfiles):
        bam[bf] = pysam.AlignmentFile(bf, "rb")


    # output dataframe
    columns = ['CHROM','POS','ALT','REF','bamfile','unique_probe_id','A','C','G','T','N','umi_hist','umi_readcount_threshold','no_umi_A','no_umi_C','no_umi_G','no_umi_T','no_umi_N' ]
    res = pd.DataFrame(columns = columns)
    num_rows = 0

    for variant in mips.get_ase_variants(): # [{'START': 160113872, 'ALT': 'C', 'CHROM': 'chr6', 'REF': 'T', 'END': 160113872}]:#mips.get_ase_variants():
        print "Variant:", variant
        assert variant['START'] == variant['END'] # only SNPs for now..
        # target is a window around the target ASE variant such that also the mate read is fetched.
        # what if the variant 
        target = [ variant['CHROM'], max(variant['START'] - match_options.window_around_ase_variant, 0), variant['END'] + match_options.window_around_ase_variant]
        for bamfilename in bam:
            print "\tBAMFILE",bamfilename 
            dict_forward_reads, dict_reverse_reads = get_reads_for_target(target, bam[bamfilename], match_options, do_filter = False)
            print "\tNumber of reads for",target,":",len(dict_forward_reads)
        
            # identify mapped pairs and create lists
            list_forward_reads, list_reverse_reads = get_pair_mapped_as_list(dict_forward_reads, dict_reverse_reads)
            print "\t",len(list_forward_reads),len(dict_forward_reads)

            # make sure that the mate reads is caught for the majority of reads
            if not float(len(list_forward_reads)) > 0.95 * float(len(dict_forward_reads)):
                sys.stderr.write("Warning. float(len(list_forward_reads)) > 0.95 * float(len(dict_forward_reads)) for %s:%d in %s\n" % (variant['CHROM'],variant['START'], bamfilename))
            
            print "\tNumber of read pairs:",len(list_forward_reads)
        
            # match each read pair to a MIP, so that we can compute consensus sequence for reads with same UMI
            dict_match_probe_idx = mips.match_read_pairs_to_probe_dict(list_forward_reads, list_reverse_reads)
            dict_qname_to_probe_idx = {}
            probe_umi_hist = {}
            for probe_idx in dict_match_probe_idx:
                probe_umi_hist[probe_idx] = {} # use below for consensus calling
                for read_idx in dict_match_probe_idx[probe_idx]:
                    dict_qname_to_probe_idx[ list_forward_reads[read_idx].qname ] = probe_idx
            num_qual_skipped = 0
            qnames_seen = {}
            for pileupcolumn in bam[bamfilename].pileup(variant['CHROM'], variant['START']-1, variant['END']+1,**{'max_depth':1000000} ):
                if pileupcolumn.pos == variant['START']-1: #pileup is zero-based
                    # print ("\tcoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    unmatched_read = 0
                    # print "\tpc:",pileupcolumn.pos
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.
                            # print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position])) 
                            qname = pileupread.alignment.query_name
                            qnames_seen[qname] = 1
                            qbase = pileupread.alignment.query_sequence[pileupread.query_position]
                            qqual = pileupread.alignment.query_qualities[pileupread.query_position]
                            if qqual < match_options.min_base_qual:
                                num_qual_skipped += 1
                                continue
                            umi = get_umi_from_qname( qname )
                            if umi.count('N') == 0:
                                if dict_qname_to_probe_idx.has_key(qname):
                                    probe_idx = dict_qname_to_probe_idx[ qname ]
                                    pdict = probe_umi_hist[probe_idx]
                                    if not pdict.has_key(umi):
                                        pdict[umi] = {'A':0,'C':0,'G':0,'T':0,'N':0}
                                    probe_umi_hist[probe_idx][umi][qbase] += 1
                                else:
                                    unmatched_read += 1
            print "\tNumber of unpaired reads:",unmatched_read
            print "\tNumber of low-quality nucleotides:",num_qual_skipped
            print "\tNumber of read-pairs seen in pileup:",len(qnames_seen.keys())
            # now output for each variant a table for allele counts in each probe.
            for probe_idx in probe_umi_hist:
                probe_id = "unmatched"
                if probe_idx >= 0:
                    probe_id = mips.unique_probe_ids[probe_idx]
                rd = {'CHROM':variant['CHROM'],
                      'POS':variant['START'],
                      'ALT':variant['ALT'],
                      'REF':variant['REF'],
                      'bamfile' : bamfilename, 
                      'unique_probe_id' : probe_id,
                      'A' : 0,
                      'C' : 0,
                      'G' : 0,
                      'T' : 0,
                      'N' : 0,
                      'no_umi_A' : 0,
                      'no_umi_C' : 0,
                      'no_umi_G' : 0,
                      'no_umi_T' : 0,
                      'no_umi_N' : 0,
                      'umi_hist' : 'NA',
                      'umi_readcount_threshold':'NA'}
                pdict = probe_umi_hist[probe_idx]
                
                # update readcount histogram
                readcount_to_numumis = {}
                for umi in pdict:
                    readcount = sum(pdict[umi].values())
                    if not readcount_to_numumis.has_key(readcount):
                        readcount_to_numumis[readcount] = 1
                    else:
                        readcount_to_numumis[readcount] += 1

                min_num_reads = determine_readcount_threshold(readcount_to_numumis, match_options.threshold_fraction_data)
                print "\tprobe_idx[%d] read_count_threshold: %d" % (probe_idx, min_num_reads)
                
                # determine umi readcount threshold to avoid counting UMIs 
                # due to sequencing errors
                for umi in pdict:
                    readcount = sum(pdict[umi].values())
                    if readcount >= min_num_reads:
                        basecall = get_majority_vote(pdict[umi])
                        rd[basecall] += 1
                    # this counts alleles without taking into account the umi's, for purpose of accuracy comparison.
                    for basecall in ['A','C','G','T','N']:
                        rd['no_umi_'+basecall] += pdict[umi][basecall]

                rd['umi_hist'] = ",".join(["%d:%d" % (k, readcount_to_numumis[k]) for k in sorted(readcount_to_numumis.keys())])
                rd['umi_readcount_threshold'] = min_num_reads
                
                res.loc[num_rows] = pd.Series(rd)
                num_rows += 1
            res.to_csv(match_options.output_file, sep = "\t", index = None)



    for bamfilename in bam: 
        bam[bamfilename].close()
    
    return res

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ase_mip_design_file", help="mip design file", required = True)
    parser.add_argument("--bamfile", help="indexed BAM file", required = True)
    parser.add_argument("--output_file", help="tab-separated output file with counts in columns numObserved and numObservedUMIs", required = True)

    parser.add_argument("--min_query_length_forward", type = int, default = 80, help="minimum length of forward read")
    parser.add_argument("--min_query_length_reverse", type = int, default = 71, help="minimum length of reverse read")
    parser.add_argument("--min_fraction_read_overlap_target", type = float, default = 0.90, help="minimum fraction of read that should overlap MIP target sequence") 
    parser.add_argument("--min_fraction_matched_nt", type = float, default = 0.95, help = "minimum fraction of bases in reads that should have cigar string MATCH")
    parser.add_argument("--window_around_ase_variant", type = int, default = 100000)
    parser.add_argument("--min_base_qual", type = int, default = 0)
    parser.add_argument("--threshold_fraction_data", type = float, default = 0.95, help = "minimum read count for UMI to be included is value above which threshold_fraction_data fraction of reads is accounted for")
    args = parser.parse_args()
    
    mips = MIPS(args.ase_mip_design_file)
    
    res_df = count_umis(mips, args.bamfile, args)
    res_df.to_csv(args.output_file, sep = "\t", index = None)
    return res_df

# run count_umis_at_mips_targets_ase.py --ase_mip_design_file ase_design_jan_2016/cdna_mips_ase-K562-HEK293.txt --bamfile /Volumes/Kees_SD_1/Manuscripts/cDNA-smMIPS/data/for_paper/experiment_ase-K562-HEK293-1/original_fastq/JB_K562HEK_8-2_R1_09-03-2016.stripped.fastq.gz.star.Aligned.sortedByCoord.out.bam

if __name__ == "__main__":
    res_df = main()


    



