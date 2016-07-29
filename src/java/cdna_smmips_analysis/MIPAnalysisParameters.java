/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.util.HashMap;

/**
 *
 * @author Z009157
 */
public class MIPAnalysisParameters {
    
    int thresholdExtLigDimerDistance; 
    int dimerKMERsize;
    int maxProbeMismatchesDimerDetection;
    int max_ext_mismatch, max_lig_mismatch; // maximum number of mismatches between extension and ligation probe
    int read_num_blocks;
    boolean createMIPKmerHash; // tells whether hash of kmer counts for each target should be contructed
                                      // -1 : don't do anything, 0 : create it, 1 : use it.
    int mipTargetKMERSize;
    String forwardFileMatched,  reverseFileMatched;
    String forwardFileUnmatched,  reverseFileUnmatched;
    boolean estimateFragmentLength;
    boolean blockCountUMIs;
    boolean outputUnmatched; // output unmatched reads to fastq for mapping.
    boolean outputMatched;
    int blockSize;
    
    double fractionReadsExplained;
    
    HashMap<MIPTarget, MIPTargetSequenceStatistics> statsWithKMERFreqsForFiltering;
    boolean filterReadsUsingKMERs;
    
    String kmerHistogramOutputFile;
    
    // remove UMIs in output to Fastq
    boolean removeUMIs;
    
    public MIPAnalysisParameters()
    {
        thresholdExtLigDimerDistance = 20;
        dimerKMERsize = 15;
        maxProbeMismatchesDimerDetection = 1;
        max_ext_mismatch = 2; // Ma
        max_lig_mismatch = 2; 
        
        read_num_blocks = -1;
        mipTargetKMERSize = 10;
        outputMatched = false;
        outputUnmatched = false;
        blockSize = 10000;
        blockCountUMIs = true;
        createMIPKmerHash = false;
        fractionReadsExplained = 0.95;
        
        statsWithKMERFreqsForFiltering = null;
        filterReadsUsingKMERs = false;
        kmerHistogramOutputFile = "";
        
        removeUMIs = false;
    }

    
}
