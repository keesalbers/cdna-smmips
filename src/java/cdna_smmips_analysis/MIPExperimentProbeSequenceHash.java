/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Set;
import java.util.Collection;
/**
 *
 * @author Z009157
 */
public class MIPExperimentProbeSequenceHash {
    // maps probe sequences to MIP targets
    protected HashMap <String, HashMap<Integer, MIPTarget> > extensionProbeSequenceToTarget;
    protected HashMap <String, HashMap<Integer, MIPTarget> > ligationRCProbeSequenceToTarget;
    
    // same but for the longest minimal probe sequence
    protected HashMap <String, HashMap<Integer, MIPTarget> > minExtensionProbeSequenceToTarget;
    protected HashMap <String, HashMap<Integer, MIPTarget> > minLigationRCProbeSequenceToTarget;
    
    protected ArrayList<MIPTarget> mipTargets;
    int extensionProbeMinLength; // minimum length across all MIP probes in experiment
    int ligationProbeMinLength; // minimum length across all MIP probes in experiment
    int seedSequenceHashLength;
    boolean doExhaustiveMatch;
    
    public class LookupResult {
        public MIPTarget mip;
        public String umiExtension, umiLigation;
        int result, diff;
        int[] dimerDistanceScore;
        // int[] dimerScore;
        
        public MIPTarget getMip() {
            return mip;
        }

        public String getUmiExtension() {
            return umiExtension;
        }

        public String getUmiLigation() {
            return umiLigation;
        }

        public int getResult() {
            return result;
        }
        
        public LookupResult() {
            mip = null;
            result = -1;
            umiExtension = "";
            umiLigation = "";
            dimerDistanceScore = new int[2];
            //dimerScore = new int[2];
        }
        public LookupResult(MIPTarget _mip, String _umiExtension,  String _umiLigation, int _result) {
            mip = _mip;
            umiExtension = _umiExtension;
            umiLigation = _umiLigation;
            result = _result;
            dimerDistanceScore = new int[2];
            //dimerScore = new int[2];
        } 
    };
    
    public MIPExperimentProbeSequenceHash(ArrayList<MIPTarget> _mipTargets, int seedSequenceHashLength) throws Exception {
        this.mipTargets = _mipTargets; //TODO maybe shouldn't store this.
        if (seedSequenceHashLength<0) {
            this.doExhaustiveMatch = true;
            this.seedSequenceHashLength = 8;
        } else {
            this.seedSequenceHashLength = seedSequenceHashLength;
            this.doExhaustiveMatch = false;
        }
        this.setProbeSequenceHashTables(this.mipTargets);
    }
    
    private final int setMinSequenceHash(HashMap <String, HashMap<Integer, MIPTarget> > probeSequenceToTarget, HashMap <String, HashMap<Integer, MIPTarget> > minProbeSequenceToTarget) {
        int minLength = 1000000;
        for (String seq : probeSequenceToTarget.keySet()) {
            if (seq.length() < minLength) {
                minLength = seq.length();
            }
        }
        for (String seq : probeSequenceToTarget.keySet()) {
            String kseq = seq.substring(0, minLength);
            if (!minProbeSequenceToTarget.containsKey(kseq)) {
                minProbeSequenceToTarget.put(kseq, new HashMap<Integer, MIPTarget>());
            }             
            for (Integer i : probeSequenceToTarget.get(seq).keySet()) {
                MIPTarget mip = probeSequenceToTarget.get(seq).get(i);
                minProbeSequenceToTarget.get(kseq).put(i, mip);
            }
        }
        return minLength;
    }
    
    final void setProbeSequenceHashTables(ArrayList<MIPTarget> mipTargets) throws Exception {
        extensionProbeSequenceToTarget = new HashMap <String, HashMap<Integer, MIPTarget> > ();
        ligationRCProbeSequenceToTarget = new HashMap <String, HashMap<Integer, MIPTarget> > ();
        
        for (int i = 0; i < mipTargets.size(); ++i) {
            // ligation probe
            MIPTarget mip = mipTargets.get(i);
            String probeSeq = mip.extensionProbeSequence.substring(0, seedSequenceHashLength);
            if (!extensionProbeSequenceToTarget.containsKey(probeSeq)) {
                extensionProbeSequenceToTarget.put(probeSeq, new HashMap<Integer, MIPTarget>());
            } 
            HashMap<Integer, MIPTarget> map = extensionProbeSequenceToTarget.get(probeSeq);
            if (map.containsKey(mip.indexInMIPTargetArray)) throw new Exception("Error.");
            map.put(mip.indexInMIPTargetArray, mip);
            
            // ligation probe. Use RC of probe sequence.
            probeSeq = mip.ligationProbeSequenceRC.substring(0, seedSequenceHashLength);
            if (!ligationRCProbeSequenceToTarget.containsKey(probeSeq)) {
                ligationRCProbeSequenceToTarget.put(probeSeq, new HashMap<Integer, MIPTarget>());
            } 
            
            map = ligationRCProbeSequenceToTarget.get(probeSeq);
            if (map.containsKey(mip.indexInMIPTargetArray)) throw new Exception("Error.");
            map.put(mip.indexInMIPTargetArray, mip);                        
        }
        
        minExtensionProbeSequenceToTarget = new HashMap <String, HashMap<Integer, MIPTarget> > ();       
        minLigationRCProbeSequenceToTarget = new HashMap <String, HashMap<Integer, MIPTarget> > ();
        
        this.extensionProbeMinLength = setMinSequenceHash(this.extensionProbeSequenceToTarget, this.minExtensionProbeSequenceToTarget);
        this.ligationProbeMinLength = setMinSequenceHash(this.ligationRCProbeSequenceToTarget, this.minLigationRCProbeSequenceToTarget);
        
        if (this.extensionProbeMinLength>seedSequenceHashLength)
            this.extensionProbeMinLength = seedSequenceHashLength;
        if (this.ligationProbeMinLength>seedSequenceHashLength)
            this.ligationProbeMinLength = seedSequenceHashLength;
        if (this.doExhaustiveMatch) {
            System.out.printf("Matching probe sequences exhaustively.\n");
        } else {
            System.out.printf("extensionProbeMinLength: %d ligationProbeMinLength: %d (used for finding candidate matches)\n",extensionProbeMinLength, ligationProbeMinLength);
        }
    }
    
    void findProbe(String ligationRead, String extensionRead) {
    // This finds most of the reads have extension probe at position 9 (length UMI) in second read.
        System.out.printf("\n");
        for (String k : minExtensionProbeSequenceToTarget.keySet()) {
            int f1 = ligationRead.indexOf(k,0);
            int f2 = extensionRead.indexOf(k,0);
            if (f1 != -1 || f2 != -1) {
                System.out.printf("\textensionRead  k: %s f1: %d f2:%d", k, f1, f2);
                for (Integer i : minExtensionProbeSequenceToTarget.get(k).keySet())
                    System.out.printf(" [%d ligRC:%s ext:%s]", i, minExtensionProbeSequenceToTarget.get(k).get(i).ligationProbeSequenceRC, minExtensionProbeSequenceToTarget.get(k).get(i).extensionProbeSequence);
                System.out.printf(" extension read: %s \n", extensionRead);
            }
        }
        for (String k : minLigationRCProbeSequenceToTarget.keySet()) {
            int f1 = ligationRead.indexOf(k,0);
            int f2 = extensionRead.indexOf(k,0);
            if (f1 != -1 || f2 != -1) {
                System.out.printf("\tligationRead k: %s f1: %d f2:%d", k, f1, f2);
                for (Integer i : minLigationRCProbeSequenceToTarget.get(k).keySet())
                    System.out.printf(" [%d lig:%s ext:%s]", i,minLigationRCProbeSequenceToTarget.get(k).get(i).ligationProbeSequence, minLigationRCProbeSequenceToTarget.get(k).get(i).extensionProbeSequence);
                System.out.printf(" ligation read: %s \n", ligationRead);
            }
        }
        System.out.printf("\n");
    }
    
    
    
    /***
     * 
     * @param readPair
     * @param umi
     * @return 
     */
    LookupResult getMIPTarget(String ligationRead, String extensionRead, int lenUmiLigation, int lenUmiExtension, MIPAnalysisParameters mipAnalysisParameters) throws Exception {
        HashMap<Integer, MIPTarget> mapLigation, mapExtension;
        LookupResult res = new LookupResult();
        HashMap<MIPTarget,Integer> candidates = new HashMap<MIPTarget,Integer>();
        
        res.result = 0;
        res.umiExtension = extensionRead.substring(0, lenUmiExtension);
        res.umiLigation = ligationRead.substring(0,lenUmiLigation);
    
       
        if (!this.doExhaustiveMatch) {
            String rcMinProbeSeq = ligationRead.substring(lenUmiLigation, lenUmiLigation + ligationProbeMinLength);
            mapLigation = this.minLigationRCProbeSequenceToTarget.get(rcMinProbeSeq);
            if (mapLigation != null) {
                for (Integer i : mapLigation.keySet()) {
                    candidates.put(mapLigation.get(i),1);
                }
            }

            String minProbeSeq = extensionRead.substring(lenUmiExtension, lenUmiExtension + extensionProbeMinLength);
            mapExtension = this.minExtensionProbeSequenceToTarget.get(minProbeSeq);
            if (mapExtension!=null) {
                // System.out.print("mapExtension.length(): " + mapExtension.size() + "\n");
                for (Integer i : mapExtension.keySet()) {
                    int count = 0;
                    MIPTarget mip = mapExtension.get(i);
                    if (candidates.containsKey(mip)) 
                        count = candidates.get(mip);
                    count += 1;
                    candidates.put(mip,count);
                }
            }
        }
        
        
        
        
        
        if (candidates.size() == 0 || this.doExhaustiveMatch) {
            for (int i = 0; i < this.mipTargets.size();i++) {
                candidates.put(this.mipTargets.get(i),1);
            }
            // System.out.print("No candidates\n");
        }
        
        HashMap<MIPTarget,Integer> selected = new HashMap<MIPTarget,Integer>();
        for (MIPTarget mip : candidates.keySet()) {
             String ligReadProbe = ligationRead.substring(lenUmiLigation, lenUmiLigation + mip.ligationProbeSequenceRC.length());
             String extReadProbe = extensionRead.substring(lenUmiExtension, lenUmiExtension + mip.extensionProbeSequence.length());
             int numLigMismatch = SeqUtil.countMismatches(ligReadProbe, mip.ligationProbeSequenceRC);
             int numExtMismatch = SeqUtil.countMismatches(extReadProbe, mip.extensionProbeSequence);
             int score = numLigMismatch + numExtMismatch;
             if (numLigMismatch <= mipAnalysisParameters.max_lig_mismatch && numExtMismatch <= mipAnalysisParameters.max_ext_mismatch) {
                 selected.put(mip, score);
             }
        }
        int bestScore = 10000;
        int diff = 10000;
        MIPTarget bestMIP = null;
        res.result = -1;
        for (MIPTarget mip : selected.keySet()) {
            int score = selected.get(mip);
            if (score < bestScore) {
                bestMIP = mip;
                diff = bestScore - score;
                bestScore = score;
                res.result = score;
                res.diff = diff;
            }
        }
        
        /*
        if (bestMIP == null) {
            // do exhaustive search
            this.findProbe(ligationRead, extensionRead);
        }
        */
       
        
        res.mip = bestMIP;
        if (res.mip != null) {
            // res.mip.checkExtensionLigationDimer(ligationRead, extensionRead, lenUmiLigation, lenUmiExtension, res.dimerScore);
            res.mip.checkExtensionLigationDimerFullRead(ligationRead, extensionRead, lenUmiLigation, lenUmiExtension, mipAnalysisParameters.maxProbeMismatchesDimerDetection, res.dimerDistanceScore);
        }
        
        return res;
    }
}
