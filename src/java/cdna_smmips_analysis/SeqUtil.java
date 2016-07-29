/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;
import java.math.*;
import java.util.*;
/**
 *
 * @author Z009157
 */
public class SeqUtil {
    public static final char[] nucleotides = {'A','C','G','T','N' };
    static final Random rn = new Random();
    static public char complement(char ch) throws Exception {
            if (ch == 'a' || ch == 'A') {
               return 'T';
            } else if (ch == 'c' || ch == 'C') {
               return 'G';
            } else if (ch == 'g' || ch == 'G') {
               return 'C';
            } else if (ch == 't' || ch == 'T') {
               return 'A';
            } else throw new Exception("Unknown base");
    }
    
    static public char complementRandomN(char ch) {
            if (ch == 'a' || ch == 'A') {
               return 'T';
            } else if (ch == 'c' || ch == 'C') {
               return 'G';
            } else if (ch == 'g' || ch == 'G') {
               return 'C';
            } else if (ch == 't' || ch == 'T') {
               return 'A';
            } else  {
                return SeqUtil.nucleotides[SeqUtil.rn.nextInt(4)];
            }
    }
    
    static public String reverseComplement(String seq) throws Exception {
        int L = seq.length();        
        char[] nseq = new char[L];
        for (int i = 0; i<L; ++i) nseq[L-i-1] = complement(seq.charAt(i));        
        return String.copyValueOf(nseq);
    }
    
    static public String reverseComplementRandomN(String seq)  {
        int L = seq.length();        
        char[] nseq = new char[L];
        for (int i = 0; i<L; ++i) nseq[L-i-1] = complementRandomN(seq.charAt(i));        
        return String.copyValueOf(nseq);
    }
    
    static public int countMismatches(String a, String b) {
        if (a.length() != b.length()) return a.length() + b.length();
        int nm = 0;
        for (int i = 0; i < a.length(); i++) {
            if (a.charAt(i)!=b.charAt(i)) {
                nm++;
            }
        }
        return nm;
    }
    
}
