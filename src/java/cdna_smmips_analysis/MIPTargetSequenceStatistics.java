/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;
import java.util.TreeSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Map;

/**
 *
 * @author Z009157
 */
public class MIPTargetSequenceStatistics {
    
    
    
    public class UMIHistogram {
        TreeMap<Integer, Integer> readCountToNumUMIs;
        public UMIHistogram(TreeMap<String, Integer> umis) {
            readCountToNumUMIs = new TreeMap<Integer, Integer>();
            for (Map.Entry<String, Integer> e : umis.entrySet()) {
                int numReads = e.getValue();
                if (readCountToNumUMIs.containsKey(numReads)) {
                    int count = readCountToNumUMIs.get(numReads);
                    readCountToNumUMIs.put(numReads, count + 1);
                } else {
                    readCountToNumUMIs.put(numReads, 1);
                }            
            }
        }

        public int getReadCountThreshold(double fractionReadsExplained) {
            int totReads = 0;
            for (Integer readcount : this.readCountToNumUMIs.descendingKeySet()) {
                int numUMIs = this.readCountToNumUMIs.get(readcount);
                totReads += numUMIs * readcount;
            }
            
            int cumReads = 0;
            for (Integer readcount : this.readCountToNumUMIs.descendingKeySet()) {
                int numUMIs = this.readCountToNumUMIs.get(readcount);
                cumReads += numUMIs * readcount;
                double fractionExplained = ((double)(cumReads))/((double)(totReads));
                if (fractionExplained >= fractionReadsExplained) {
                    return readcount;
                }
            }
            return 1;
        }
        
        public int getCorrectedUMICount(int readCountThreshold) {
            int numUMIs = 0;
            for (Integer readcount : this.readCountToNumUMIs.descendingKeySet()) {
                if (readcount >= readCountThreshold) {
                    numUMIs += this.readCountToNumUMIs.get(readcount);
                }
            }
            return numUMIs;
        }
        
        public String getString() {
            String s = "";
            int i = 0;
            for (Integer count : readCountToNumUMIs.keySet()) {
                if (i>0) s+= ",";
                s += String.format("%s:%d", count,readCountToNumUMIs.get(count));
                i++;
            }
            if (s.equals("")) {
                s = "NA";
            }
            return s;
        }
        
    };
    
    int numCount;
    int numUMIsWithNs;
    TreeMap<String, Integer> umis;
    TreeMap<String, Integer> umiToDimers;
    TreeMap<Integer, Integer> fragLengthHistogram;
    boolean outputUmiHistogram;
        
    MIPKmerCounts mipKmerCounts; // Keeps kmer frequency distribution to screen out crappy reads that are not ext-lig dimers.
    
    // Adds observations of ms to this
    public void addStatistics(MIPTargetSequenceStatistics ms) {
        this.numCount += ms.numCount;
        // add umi observations
        for (String umi : ms.umis.keySet()) {
            this.umis.put(umi, ms.umis.get(umi));
        }
        for (String umi : ms.umiToDimers.keySet()) {
            this.umiToDimers.put(umi, ms.umiToDimers.get(umi));
        }
        // add fraglength observations
        for (Integer fl : ms.fragLengthHistogram.keySet()) {
            this.addFragLengthObservation(fl, ms.fragLengthHistogram.get(fl));
        }
    }
    public MIPTargetSequenceStatistics(MIPAnalysisParameters mipAnalysisParameters) {
        numCount = 0;
        umis = new TreeMap<String,Integer>();
        fragLengthHistogram = new TreeMap<Integer, Integer>();
        umiToDimers = new TreeMap<String, Integer>();
        fragLengthHistogram = new TreeMap<Integer, Integer>();
        outputUmiHistogram = false;
        numUMIsWithNs = 0;
        mipKmerCounts = new MIPKmerCounts(mipAnalysisParameters.mipTargetKMERSize);
    }
    public int incNumCount() {
        numCount++;
        return numCount;
    }
    
    public int incNumUMIsWithNs() {
        this.numUMIsWithNs++;
        return this.numUMIsWithNs;
    }
    
    /***
     * Store observation for umi. Only gets stored as count if isDimer is false.
     * @param umi
     * @param inc
     * @param isDimer 
     */
    public void addUMI(String umi, int inc, boolean isDimer) {
        Integer cnt = umis.get(umi);
        if (!isDimer) {
            if (cnt == null) {
                umis.put(umi, inc);
            } else {
                umis.put(umi, cnt+inc);
            }
        } else {
            cnt = this.umiToDimers.get(umi);
            if (cnt == null) {
                umiToDimers.put(umi, inc);
            } else {
                umiToDimers.put(umi, cnt+inc);
            }
        }
    }
    
    public int getTotalNumberOfDimers()
    {
        int sum = 0;
        for (int value : this.umiToDimers.values()) {
            sum += value;
        }
        return sum;
    }
    
    public int addFragLengthObservation(int fragLength, int inc) {
        Integer c = fragLengthHistogram.get(fragLength);
        if (c == null) {
            fragLengthHistogram.put(fragLength, inc);
            return 1;
        } else {
            fragLengthHistogram.put(fragLength, c+inc);
            return c+1;
        }        
    }
    public double getMeanFragmentLength() {
        double count = 0;
        double weighted = 0.0;
        for (int v : this.fragLengthHistogram.keySet()) {
            double c = (double) this.fragLengthHistogram.get(v);
            weighted += c * ((double)v);
            count += c;
        }
        return weighted/count;
    }
    
    public String getUMIHistogramString() {
        if (!outputUmiHistogram) return "NA";
        TreeMap<Integer, TreeMap<String, Integer> > x = new TreeMap<Integer, TreeMap<String, Integer> >();
        for (Map.Entry<String, Integer> e : this.umis.entrySet()) {
            TreeMap<String, Integer> m = x.get(e.getValue());
            if (m == null) {
                m = new TreeMap<String, Integer>();
                x.put(e.getValue(),m);
            }
            m.put(e.getKey(),1);            
        }
        
        String s = "";
        for (Integer count : x.descendingKeySet()) {
            for (Map.Entry<String, Integer> e : x.get(count).entrySet()) {
                s += String.format("[%s:%d]", e.getKey(),count);
            }
        }
        if (s.equals("")) {
            s = "NA";
        }
        return s;
    }
    

    
    public String getFragmentHistogramString() {
        String s = "";
        int count = 0;
        for (int v : this.fragLengthHistogram.keySet()) {
            int c =  this.fragLengthHistogram.get(v);
            count += c;
        }
        for (int v : this.fragLengthHistogram.keySet()) {
            int c =  this.fragLengthHistogram.get(v);     
            int pct = (int) (((double) c)/((double) count)*100.0);
            s += String.format("[%d:%d]", v,pct);
        }
        if (s.equals("")) {
            s = "NA";
        }

        return s;
    }
    
    public void updateKMERHash(String ligationRead, String extensionRead) {
        this.mipKmerCounts.addReadPair(ligationRead, extensionRead);
    }
    /***
     * Returns header string as columns separated by tabs
     * @return 
     */
    static public String getHeaderString() {
        return "numObserved\tnumObservedUMIs\tnumObservedDimer\tnumUMIsWithDimer\tmeanFragmentLength\tfragmentLengthHistogram\tnumReadsPerUMIHistogram\tUMIReadCountThreshold\tnumObservedUMIsCorrected";
    }
    /***
     * Returns data as string, separated by tabs, in order given by getHeaderString
     * @return 
     */
    public String getRow(MIPAnalysisParameters mipAnalysisParameters) {
        UMIHistogram umiHist = new UMIHistogram(this.umis);
        int readCountThreshold = umiHist.getReadCountThreshold(mipAnalysisParameters.fractionReadsExplained);
        int correctedUMICount = umiHist.getCorrectedUMICount(readCountThreshold);
        return String.format("%d\t%d\t%d\t%d\t%1.2f\t%s\t%s\t%d\t%d", numCount, umis.size(), this.getTotalNumberOfDimers(), this.umiToDimers.size(), getMeanFragmentLength(),this.getFragmentHistogramString(), umiHist.getString(), readCountThreshold, correctedUMICount); 
    }
}
