/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.util.HashMap;
/**
 * Note. mip properties from the pipeline marked by //X were removed to put fewer restrictions on input file.
 */


/**
 *
 * @author Z009157
 */
public class MIPTarget {
    String uniqueID;
    String extensionProbeSequence, ligationProbeSequence, mipSequence, ligationProbeSequenceRC, extensionProbeSequenceRC;
    //X String chrom;
    //X String targetSequence;
    //X Integer extensionProbeStart, extensionProbeEnd, ligationProbeStart, ligationProbeEnd;
    //X String designStrand; // was probe designed for target or reverse-complement of target?
    //X String targetStrand; // was the target sequence itself on the forward strand or the reverse strand?
    
    int indexInMIPTargetArray;
    
    MIPTargetRNAProperties rnaProperties; // null if it is not an RNA smMIP
    
    public MIPTarget(Row row, HashMap <String, String> columnNames) throws Exception {
        this.init(row, columnNames);
    }
   
    
    final void init(Row row, HashMap <String, String> columnNames) throws Exception {
   
       
        
        this.uniqueID = row.get(columnNames.get("uniqueID")).toString();
        this.extensionProbeSequence = row.get(columnNames.get("extensionProbeSequence")).toString().toUpperCase();
        this.ligationProbeSequence = row.get(columnNames.get("ligationProbeSequence")).toString().toUpperCase();        
        this.mipSequence = row.get(columnNames.get("mipSequence")).toString().toUpperCase();
        
        //X this.targetSequence = row.get(columnNames.get("targetSequence")).toString().toUpperCase();
        
        //X this.chrom = row.get(columnNames.get("chrom")).toString();
        
        //X this.extensionProbeStart = Integer.decode(row.get(columnNames.get("extensionProbeStart")).toString());
        //X this.extensionProbeEnd = Integer.decode(row.get(columnNames.get("extensionProbeEnd")).toString());
        //X this.ligationProbeStart = Integer.decode(row.get(columnNames.get("ligationProbeStart")).toString());
        //X this.ligationProbeEnd = Integer.decode(row.get(columnNames.get("ligationProbeEnd")).toString());
        //X this.designStrand = row.get(columnNames.get("designStrand")).toString();
        
        //X if (row.hasColumn(columnNames.get("targetStrand"))) {
        //X    this.targetStrand = row.get(columnNames.get("targetStrand")).toString();
        //X }
       
        //X if (row.hasColumn(columnNames.get("spansEJB"))) {
        //X     rnaProperties = new MIPTargetRNAProperties(row, columnNames);
        //X }
        
        // Reverse complement ligation probe sequence
        this.ligationProbeSequenceRC = SeqUtil.reverseComplement(ligationProbeSequence);
        this.extensionProbeSequenceRC = SeqUtil.reverseComplement(extensionProbeSequence);
    }
    public boolean isOK() {
        boolean ok1 = true;
        boolean ok2 = (this.extensionProbeSequence.length() + this.ligationProbeSequence.length() >= 40);
        boolean ok3 = true; //(this.mipSequence.length() >= 70);
        boolean ok4 = true; //(this.targetSequence.length() > 100);
        boolean ok5 = true; //X (this.extensionProbeStart>0);
        boolean ok6 = true; //X (this.extensionProbeEnd>0);        
        boolean ok7 = true; //X (this.ligationProbeStart>0);        
        boolean ok8 = true; //X (this.ligationProbeEnd>0);        
        
        boolean ok = ok1 & ok2 & ok3 & ok4 & ok5 & ok6 & ok7 & ok8;
        if (!ok) {
            System.out.printf("Not ok: %b %b %b %b %b %b %b %b\n",ok1,ok2,ok3,ok4,ok5,ok6,ok7,ok8);
        }
        return ok;
    }
    /***
     * Calculate distance of extension probe in ligation read from beginning of read;
     * calculate distance of ligation probe in extension read from beginning of read;
     * @param ligationRead
     * @param extensionRead
     * @param lenUmiLigation
     * @param lenUmiExtension
     * @return 
     */
    public boolean checkExtensionLigationDimerFullRead(String ligationRead, String extensionRead, int lenUmiLigation, int lenUmiExtension, int max_num_mismatch, int[] score) throws Exception {
        boolean hit = false;
        int kmer = 15; // look for substrings of 10 probe nt in read and calculate distance
        int ext_seq_index_in_ligation_read = 1000000; // will contain first match.
        int lig_seq_index_in_extension_read = 1000000; // will contain first match
        
        // check extension probe in ligation read        
        for (int ext_probe_start = 0; ext_probe_start < this.extensionProbeSequence.length() - kmer; ext_probe_start++) {
            String ext_probe_seq_rc = SeqUtil.reverseComplement(this.extensionProbeSequence.substring(ext_probe_start, ext_probe_start + kmer));
            
            for (int index = 0; index < ligationRead.length() - ext_probe_seq_rc.length(); index++) {
                int nm = SeqUtil.countMismatches(ext_probe_seq_rc, ligationRead.substring(index, index + kmer));
                if (nm <= max_num_mismatch) {
                    hit = true;
                    ext_seq_index_in_ligation_read = index - (this.extensionProbeSequence.length()-kmer-ext_probe_start) - lenUmiLigation - this.ligationProbeSequence.length();
                    /*
                            if (ext_seq_index_in_ligation_read<0) {
                        System.out.printf("\nligationRead: %s\n", ligationRead);
                        System.out.printf("ext_probe_seq_rc: %s\n", ext_probe_seq_rc);                    
                        System.out.printf("index: %d ext_probe_start: %d lenUmiLigation: %d ligationProbeSequence.length(): %d\n", index,ext_probe_start,lenUmiLigation,this.ligationProbeSequence.length());
                        System.out.printf("ext_seq_index_in_ligation_read: %d\n", ext_seq_index_in_ligation_read);
                    //}
                    */
                    break;                        
                }
            }
        }
        
        // check ligation probe in extension read        
        for (int lig_probe_start = 0; lig_probe_start < this.ligationProbeSequence.length() - kmer; lig_probe_start++) {
            String lig_probe_seq = this.ligationProbeSequence.substring(lig_probe_start, lig_probe_start + kmer);
            
            for (int index = 0; index < extensionRead.length() - kmer; index++) {
                int nm = SeqUtil.countMismatches(lig_probe_seq, extensionRead.substring(index, index + kmer));
                if (nm <= 1) {
                    hit = true;
                    lig_seq_index_in_extension_read = index - lig_probe_start - lenUmiExtension - this.extensionProbeSequence.length();
                    break;
                }
            }
        }
        score[0] = ext_seq_index_in_ligation_read;
        score[1] = lig_seq_index_in_extension_read;
        return hit;
    }
    
    
    public void checkExtensionLigationDimer(String ligationRead, String extensionRead, int lenUmiLigation, int lenUmiExtension, int[] dimerScore) {
        
        // get potential probe sequences from ligation read. 
        
        // The following is supposed to be the RC of the extension probe sequence IF it is present in the ligationRead 
        String ligReadExtensionProbeRC = ligationRead.substring(lenUmiLigation + this.ligationProbeSequenceRC.length(), lenUmiLigation + this.ligationProbeSequenceRC.length() + this.extensionProbeSequence.length() );
        int nmLigReadExtensionProbeRC = SeqUtil.countMismatches(ligReadExtensionProbeRC, this.extensionProbeSequenceRC);
        // System.out.printf("D: %d %d\n",ligReadExtensionProbeRC.length(), this.extensionProbeSequenceRC.length());
        String extReadLigationProbe = extensionRead.substring(lenUmiExtension + this.extensionProbeSequence.length(), lenUmiExtension + this.extensionProbeSequence.length() + this.ligationProbeSequence.length());     
        int nmExtReadLigationProbe = SeqUtil.countMismatches(extReadLigationProbe, this.ligationProbeSequence);
       
        dimerScore[0] = nmLigReadExtensionProbeRC;
        dimerScore[1] = nmExtReadLigationProbe;
    } 
            
    
    public void print() {
        System.out.printf("mip %d %s %s %s\n", this.indexInMIPTargetArray, this.uniqueID, this.extensionProbeSequence, this.ligationProbeSequence);
    }
}
