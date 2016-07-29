/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.ArrayList;
import java.io.FileWriter;
import java.io.BufferedWriter;

/**
 *
 * @author Z009157
 */

public class MIPKmerCounts {
    HashMap<String, Integer> kmerFreqs;
    int kmerLength;
    int maxFreq;
    public MIPKmerCounts(int kmerLength) {
        this.kmerLength = kmerLength;
        kmerFreqs = new HashMap<String, Integer>();
    }
    /***
     * Adds all the kmers in the ligation and extension read. Ligation read is RC'd.
     * @param ligationRead
     * @param extensionRead 
     */
    void addReadPair(String ligationRead, String extensionRead)   {
        String ligationReadRC = SeqUtil.reverseComplementRandomN(ligationRead);
        this.addKmers(ligationReadRC);
        this.addKmers(extensionRead);
    }
    int getNumObs() {
        int n = 0;
        for (Integer v : this.kmerFreqs.values()) {
            n += v;
        }
        return n;
    }
    
    
    
    void addKmers(String seq) {
        for (int i = 0; i < seq.length() - this.kmerLength; i++) {
            String kmer = seq.substring(i,i+this.kmerLength);
            Integer count = this.kmerFreqs.get(kmer);
            if (count == null) {
                this.kmerFreqs.put(kmer, 1);
                if (this.maxFreq == 0) maxFreq = 1;
            } else {
                this.kmerFreqs.put(kmer, count+1);
                if (count+1>this.maxFreq) {
                    this.maxFreq = count + 1;
                }
            }
        }
    }
    
    int countNumberOfKMERSBelowThreshold(String seq, int maxFreq) {
        
        int count = 0;
        for (int i = 0; i < seq.length() - this.kmerLength; i++) {
            String kmer = seq.substring(i,i+this.kmerLength);
            int f = 0;
            Integer freq  = this.kmerFreqs.get(kmer);
            if (freq != null) f = (int)(freq);
            if (f<=maxFreq) count++;
        }
        return count;
    }
    
    boolean removeReadPair(String ligationRead, String extensionRead, int rule) {
        String ligationReadRC = SeqUtil.reverseComplementRandomN(ligationRead);
        int numKmers = ligationRead.length() + extensionRead.length() - 2*this.kmerLength;
        if (rule == 0) {
            // count number of singleton kmers in ligationRC and extension read. 
            // If more than 50% of kmers are singletons then remove the read pair 
            int freqThreshold = (int) (0.90*((double)this.maxFreq));
            int numSingletons = this.countNumberOfKMERSBelowThreshold(ligationReadRC, freqThreshold) + this.countNumberOfKMERSBelowThreshold(extensionRead, freqThreshold);
            double fraction = ((double) numSingletons) / ((double) numKmers);
            if (fraction>0.8) return true; else return false;
        }
        return false;
    }
    
    void printSummary(String id, BufferedWriter output) throws Exception {
        TreeMap<Integer, ArrayList<String> > ordered = new TreeMap<Integer, ArrayList<String> >();
        for (String kmer : this.kmerFreqs.keySet()) {
            Integer count = this.kmerFreqs.get(kmer);
            if (!ordered.containsKey(count)) {
                ordered.put(count, new ArrayList<String>());
            }
            ordered.get(count).add(kmer);
        }
        
        for (Integer count  : ordered.descendingKeySet()) {
            for (String kmer : ordered.get(count)) {
                String out = String.format("%s\t%s\t%d\n",id, kmer, count);
                output.write(out);
            }
        }
  
    }
}
