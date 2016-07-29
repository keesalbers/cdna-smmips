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
public class MakeFragment {
    int overlap, relpos, fragmentLength;
    double confidence;
    
    String read1, read2;
    HashMap<String, Integer> kmerToPos1, kmerToPos2;
    static int KMER = 5;
    static int GAP = 1;
    public MakeFragment(String read1, String read2) {
        this.read1 = read1;
        this.read2 = read2;
        init();
        match();        
    }
    
    final void init() {
        overlap = -1;
        fragmentLength = -1;
        confidence = -1.0;
        kmerToPos1 = new HashMap<String, Integer>();
        kmerToPos2 = new HashMap<String, Integer>();
        buildHash(read1, kmerToPos1);
        buildHash(read2, kmerToPos2);
    }
    final void buildHash(String read, HashMap<String, Integer> kmerToPos) {
        int rlen = read.length();
        for (int i = 0; i+KMER < rlen; i += GAP) {
            kmerToPos.put(read.substring(i,i+KMER), i);
        }
    }
    final void match() {
        HashMap<Integer, Integer> posToCount = new HashMap<Integer, Integer>();
        for (String kmer : kmerToPos1.keySet()) {
            Integer pos1 = kmerToPos1.get(kmer);
            Integer pos2 = kmerToPos2.get(kmer);
            if (pos2 != null) {
                int startRead2 = pos1-pos2;
                Integer k = posToCount.get(startRead2);
                if (k == null) {
                    posToCount.put(startRead2, 1);
                } else {
                    posToCount.put(startRead2, k+1);
                }
            }
        }
        int max = -1;
        int bestpos = -1;
        for (int rp : posToCount.keySet()) {
            int count = posToCount.get(rp);
            if (count > max) {
                bestpos = rp;
                max = count;
            }
        }
        if (max == -1) {
            this.overlap = -1;
            this.fragmentLength = -1;
            this.confidence = 1.0;
        } else {
            if (bestpos+read2.length()>0 && bestpos<read1.length()) {
                int left = this.max(bestpos, 0);
                int right = this.min(read1.length(), bestpos+read2.length());
                overlap = right-left;
                fragmentLength = this.max(read1.length(), bestpos+read2.length())-this.min(bestpos, 0);
                this.confidence = ((double) max)/((double) this.overlap-KMER);
                /*
                System.out.printf("L1:%d L2:%d max:%d bestpos:%d overlap: %d fragmentLength: %d confidence: %1.2f \n", read1.length(), read2.length(), max, bestpos, this.overlap, this.fragmentLength, this.confidence);
                if (bestpos>=0) {
                    System.out.printf("%s\n",read1); for (int i = 0; i < bestpos;i++) System.out.printf(" "); 
                    System.out.printf("%s\n",read2);
                }
                */
            } else {
                overlap = -1;
                fragmentLength = -1;
                this.confidence = 1.0;
            }
        }
        //System.out.printf("L1:%d L2:%d max:%d bestpos:%d overlap: %d fragmentLength: %d confidence: %1.2f \n", read1.length(), read2.length(), max, bestpos, this.overlap, this.fragmentLength, this.confidence);
    }
    int min(int a, int b) {
        if (a<b) {
            return a;
        } else {
            return b;
        }
    }
   int max(int a, int b) {
        if (a>b) {
            return a;
        } else {
            return b;
        }
    }
}
