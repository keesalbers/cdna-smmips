/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;
import java.lang.StringBuilder;
import java.util.HashMap;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.util.zip.GZIPOutputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
/**
 *
 * @author Z009157
 */
public class MIPProbeCoverageAnalysisBlock {
    // References to instances
    MIPExperiment mipExperiment;
    ProcessFastqPairs pfq;
    BufferedWriter[] matchedReadsOutputFiles;
    BufferedWriter[] unmatchedReadsOutputFiles;
    // updated by this class
    ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>> targetStatistics;
    HashMap<MIPTarget, MIPTargetSequenceStatistics> targetStatisticsForKmerCounting;
    // minimum distance between extension and ligation probe sequence in the SAME read
    // Below this threshold a read (and consequently the pair) will be classified as extension/ligation dimer
    
    MIPAnalysisParameters mipAnalysisParameters;
    
   
    
    
    public MIPProbeCoverageAnalysisBlock(MIPExperiment _mipExperiment, ProcessFastqPairs _pfq, MIPAnalysisParameters mipAnalysisParameters) throws Exception {
        this.mipAnalysisParameters = mipAnalysisParameters;
        this.mipExperiment = _mipExperiment;
        this.pfq = _pfq;
        this.targetStatistics = null;
        this.matchedReadsOutputFiles = null;
        this.targetStatistics = new ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics> > ();
        this.targetStatisticsForKmerCounting = this.getNewStatisticsBlock();    
        this.setFilteredReadFileNames();
    }
    
    /*
    public void setOutputUnmatched()
    {
        this.outputUnmatched = true;
    }
    */
    
    HashMap<MIPTarget, MIPTargetSequenceStatistics> getNewStatisticsBlock() {
        HashMap<MIPTarget, MIPTargetSequenceStatistics> t_stat = new HashMap<MIPTarget, MIPTargetSequenceStatistics>();
        for (int i = 0; i < this.mipExperiment.mipTargets.size(); ++i) {
            MIPTarget mip = this.mipExperiment.mipTargets.get(i);
            t_stat.put(mip, new MIPTargetSequenceStatistics(this.mipAnalysisParameters));
        }
        return t_stat;
    }
    
    final HashMap<MIPTarget, MIPTargetSequenceStatistics> addStatisticsBlock(ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>> ts) {
        HashMap<MIPTarget, MIPTargetSequenceStatistics> t_stat = this.getNewStatisticsBlock();
        ts.add(t_stat);
        return t_stat;
    }
            
    final void setFilteredReadFileNames() throws Exception {
        if (this.mipAnalysisParameters.outputMatched) {
            this.matchedReadsOutputFiles = new BufferedWriter[2];
            this.matchedReadsOutputFiles[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.mipAnalysisParameters.forwardFileMatched))));
            this.matchedReadsOutputFiles[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.mipAnalysisParameters.reverseFileMatched))));
            System.out.printf("NOTE: Writing matched read pairs to %s,%s\n", this.mipAnalysisParameters.forwardFileMatched, this.mipAnalysisParameters.reverseFileMatched);
        } else {
            this.matchedReadsOutputFiles = null;
        }
        if (this.mipAnalysisParameters.outputUnmatched) {
            this.unmatchedReadsOutputFiles = new BufferedWriter[2];
            this.unmatchedReadsOutputFiles[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.mipAnalysisParameters.forwardFileUnmatched))));
            this.unmatchedReadsOutputFiles[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.mipAnalysisParameters.reverseFileUnmatched))));
            System.out.printf("NOTE: Writing matched read pairs to %s,%s\n", this.mipAnalysisParameters.forwardFileUnmatched, this.mipAnalysisParameters.reverseFileUnmatched);
        } else {
            this.unmatchedReadsOutputFiles = null;
        }
    }
    
    int countTotalUMIs( HashMap<MIPTarget, MIPTargetSequenceStatistics> tStat) {
        int numUMIs = 0;
        for (MIPTargetSequenceStatistics stat : tStat.values()) {
            if (stat != null) {
                numUMIs += stat.umis.size();
            }
        }
        return numUMIs;
    }
    
    int countTotalUMIWithNs( HashMap<MIPTarget, MIPTargetSequenceStatistics> tStat) {
        int numUMIs = 0;
        for (MIPTargetSequenceStatistics stat : tStat.values()) {
            if (stat != null) {
                numUMIs += stat.numUMIsWithNs;
            }
        }
        return numUMIs;
    }
    
    public void analyse(boolean suppressFastqOutput) throws Exception {
        // init
        int blockCount = 0, numBlocks = 0;
        HashMap<MIPTarget, MIPTargetSequenceStatistics> tStat = this.addStatisticsBlock(this.targetStatistics);
        
        boolean done = false;
        int numRead = 0, numMatch = 0, numDimer = 0;
        ProcessFastqPairs.FastqRecord[] readPair = new ProcessFastqPairs.FastqRecord[2];
        
        int numKmerRemoved = 0;
        
        while (!(done || numRead < -1000) ) {            
            int res = pfq.getReadFromFastqFiles(readPair);
            if (res == 0) break; // EOF
            else if (res<0) {
                System.err.printf("There was an error: %d\n", res);
            } else if (res == 1) {
                String ligationRead = readPair[0].getSeq();
                String extensionRead = readPair[1].getSeq();
                
                MIPExperimentProbeSequenceHash.LookupResult mipRes = this.mipExperiment.mipExperimentProbeSequenceHash.getMIPTarget(ligationRead, extensionRead, this.mipExperiment.mipExperimentProperties.umiLengthLigation, this.mipExperiment.mipExperimentProperties.umiLengthExtension, this.mipAnalysisParameters);
                MIPTargetSequenceStatistics stat = null, kmerStat = null;
                
                if (mipRes.getResult() >= 0) {
                    numMatch++;
                    stat = tStat.get(mipRes.mip);
                    stat.incNumCount();
                    
                    kmerStat = this.targetStatisticsForKmerCounting.get(mipRes.mip);
                    kmerStat.incNumCount();
                    
                    // detect dimers
                    boolean isDimer = false;
                    if (mipRes.dimerDistanceScore[0] < mipAnalysisParameters.thresholdExtLigDimerDistance || mipRes.dimerDistanceScore[1] < mipAnalysisParameters.thresholdExtLigDimerDistance) {
                        isDimer = true;
                        numDimer++;
                    }
                    
                    String umi = mipRes.umiExtension + mipRes.umiLigation;
                    if (!umi.contains("Nn")) {
                        stat.addUMI(umi,1, isDimer);  
                        kmerStat.addUMI(umi, 1, isDimer);
                        
                        // we are happy with this read.
                        if (this.mipAnalysisParameters.createMIPKmerHash) {
                            kmerStat.updateKMERHash(ligationRead, extensionRead);                                 
                        } else if (this.mipAnalysisParameters.filterReadsUsingKMERs) {
                            boolean remove = this.mipAnalysisParameters.statsWithKMERFreqsForFiltering.get(mipRes.mip).mipKmerCounts.removeReadPair(ligationRead, extensionRead, 0);
                            if (remove) numKmerRemoved++;
                        }
                    } else {
                        stat.incNumUMIsWithNs();
                        kmerStat.incNumUMIsWithNs();
                    }
                    
                    if (!suppressFastqOutput && this.mipAnalysisParameters.outputMatched) {
                        if (this.mipAnalysisParameters.removeUMIs) {
                            UMIRemover.trimUMI(readPair, this.mipExperiment.mipExperimentProperties.umiLengthLigation, this.mipExperiment.mipExperimentProperties.umiLengthExtension);
                        }
                        readPair[0].outputToFile(this.matchedReadsOutputFiles[0]);
                        readPair[1].outputToFile(this.matchedReadsOutputFiles[1]);
                    }
                    
                } else {
                    // write unmatched readpairs to fastq file.
                    if (!suppressFastqOutput && this.mipAnalysisParameters.outputUnmatched) {
                            if (this.mipAnalysisParameters.removeUMIs) {
                                UMIRemover.trimUMI(readPair, this.mipExperiment.mipExperimentProperties.umiLengthLigation, this.mipExperiment.mipExperimentProperties.umiLengthExtension);
                            }

                            readPair[0].outputToFile(this.unmatchedReadsOutputFiles[0]);
                            readPair[1].outputToFile(this.unmatchedReadsOutputFiles[1]);
                    }                    
                }
                if (this.mipAnalysisParameters.estimateFragmentLength) {
                    int fragmentLength = this.mipExperiment.estimateFragmentLength(ligationRead, extensionRead); 
                    if (stat != null) {
                        stat.addFragLengthObservation(fragmentLength,1);
                    }
                }
            } else {
                throw new Exception("HUH?");
            }
            numRead++;
            if (this.mipAnalysisParameters.blockCountUMIs) {
                blockCount = this.countTotalUMIs(tStat);
                // System.out.print("blockCount: "+blockCount+ " numReads: "+numRead + "\n");
            } else {
                blockCount++;
            }
            if (blockCount == this.mipAnalysisParameters.blockSize) {
                int tn = this.countTotalUMIWithNs(tStat);
                tStat = this.addStatisticsBlock(this.targetStatistics);
                blockCount = 0;
                numBlocks++;
                //System.out.printf("Read block %d\n", numBlocks);
                System.out.printf("BLOCK %d Read %d fastq records and matched %d reads to MIP target (%1.2f percent). %d Lig-Ext dimers (%1.2f percent of matched read pairs at distance %d). %d reads had Ns in UMI. %d blocks.\n", numBlocks, numRead, numMatch, 100.*( (double) numMatch)/( (double) numRead ), numDimer, 100.*( (double) numDimer)/( (double) numMatch ),mipAnalysisParameters.thresholdExtLigDimerDistance, tn, numBlocks);
                if (this.mipAnalysisParameters.read_num_blocks != -1 && this.mipAnalysisParameters.read_num_blocks == numBlocks) {
                    System.out.printf("Stopping early because user option readNumBlocks=%d was specified\n", this.mipAnalysisParameters.read_num_blocks);
                    break;
                }
            }
            
            
        }
        if (!suppressFastqOutput) {
            if (this.mipAnalysisParameters.outputMatched) {
                this.matchedReadsOutputFiles[0].close();
                this.matchedReadsOutputFiles[1].close();
            }
            if (this.mipAnalysisParameters.outputUnmatched) {
                this.unmatchedReadsOutputFiles[0].close();
                this.unmatchedReadsOutputFiles[1].close();
            }          
        }
        System.out.printf("ALL Read %d fastq records and matched %d reads to MIP target (%1.2f percent). %d Lig-Ext dimers (%1.2f percent of matched read pairs at distance %d). Read %d blocks.\nLast blockCount: %d\n", numRead, numMatch, 100.*( (double) numMatch)/( (double) numRead ), numDimer, 100.*( (double) numDimer)/( (double) numMatch ), mipAnalysisParameters.thresholdExtLigDimerDistance, numBlocks, blockCount);
        if (this.mipAnalysisParameters.filterReadsUsingKMERs) {
            System.out.printf("Removed %d read pairs using KMER-based filtering\n", numKmerRemoved);
        }
        if (this.mipAnalysisParameters.createMIPKmerHash) {
            this.printKmerHashSummary();
        }
    }
    public void outputToTables(String prefix) throws Exception {
        
        // output individual blocks
        if (this.mipAnalysisParameters.blockCountUMIs) {
            this.outputToTable(prefix+".blocks.statistics.txt", this.targetStatistics);
        }
        
        // output 10x merged blocks
        ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>> tenBlockStatistics = new ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>>();
        int numTen = this.targetStatistics.size()/10;
        
        // Don't add the last bit, we only want full blocks of 10.
        // if (this.targetStatistics.size() % 10 !=0 ) numTen++;
        int i = 0;
        for (int x = 0; x < numTen; x++) {
            HashMap<MIPTarget, MIPTargetSequenceStatistics> hm = this.addStatisticsBlock(tenBlockStatistics);
            for (int k = 0; k < 10; k++) {
                for (MIPTarget mip : hm.keySet()) {
                    hm.get(mip).addStatistics(this.targetStatistics.get(i).get(mip));
                }
                i++;
                if (i==this.targetStatistics.size()) {
                    break;
                }
            }
        }
        
        if (this.mipAnalysisParameters.blockCountUMIs) {
            this.outputToTable(prefix+".ten_blocks.statistics.txt", tenBlockStatistics);
        }
        
        // merge into single block
        ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>> mergedStatistics = new ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>>();
        HashMap<MIPTarget, MIPTargetSequenceStatistics> merged = this.addStatisticsBlock(mergedStatistics);
        for (int k = 0; k < this.targetStatistics.size(); k++) {
            for (MIPTarget mip : merged.keySet()) {
                    merged.get(mip).addStatistics(this.targetStatistics.get(k).get(mip));
                }
        }
        this.outputToTable(prefix+".merged.statistics.txt", mergedStatistics);
        
    }
    
    public void printKmerHashSummary() throws Exception {
        if (!this.mipAnalysisParameters.kmerHistogramOutputFile.isEmpty()) {
            System.out.println("Writing kmer histogram file.");
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter((new FileOutputStream(this.mipAnalysisParameters.kmerHistogramOutputFile))));

            for (MIPTarget mip : this.targetStatisticsForKmerCounting.keySet()) {
                MIPTargetSequenceStatistics ms = this.targetStatisticsForKmerCounting.get(mip);            
                ms.mipKmerCounts.printSummary(mip.uniqueID, bw);
            }
            bw.close();
        } else {
            System.out.println("NOT Writing kmer histogram file because kmer output file was not specified.");
        }
        
        
    }

    public void outputToTable(String fileName, ArrayList< HashMap<MIPTarget, MIPTargetSequenceStatistics>> blockTargetStatistics) throws Exception {
        if (blockTargetStatistics == null) {
            throw new Exception("Need to call analyse first.");
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
        TableData table = this.mipExperiment.getTable();
        // write new header
        bw.write(table.getHeaderString("\t")+"\t"+"blockIndex" + "\t" + MIPTargetSequenceStatistics.getHeaderString()+"\n");
        for (int bidx = 0; bidx < blockTargetStatistics.size(); bidx++) {
            for (int i = 0; i < table.numRows();++i) {
                
                String old = table.getRowAsString(i, "\t");
                String uniqueID = table.getData(i, this.mipExperiment.columnNames.get("uniqueID")).toString();
                MIPTarget mip = this.mipExperiment.getMIP(uniqueID);
                MIPTargetSequenceStatistics stat = blockTargetStatistics.get(bidx).get(mip);
                // System.out.print("i,bidx" + i  + " " + bidx + "mip:" + mip + " id: " + uniqueID + "\n");
                String statString = stat.getRow(this.mipAnalysisParameters);
                bw.write(old + "\t" + bidx + "\t" + statString + "\n");
                
// System.out.printf("UMIS: %s\n",stat.umis.toString());
            }
        }
        bw.close();
    }
}



/*
                    try {
                    if ( (mipRes.dimerDistanceScore[0]<5 || mipRes.dimerDistanceScore[1]<5)) {
                        System.out.printf("\n");
                        System.out.printf("ligationReadRC: %s\n ", SeqUtil.reverseComplement(ligationRead));
                        System.out.printf("extensionRead: %s\n", SeqUtil.reverseComplement(extensionRead));
                        System.out.printf("Ligation/extension ps: %s %s\n", mipRes.getMip().ligationProbeSequence, mipRes.getMip().extensionProbeSequence);
                        System.out.printf("LigationRC/extensionRC ps: %s %s\n", SeqUtil.reverseComplement(mipRes.getMip().ligationProbeSequence), SeqUtil.reverseComplement(mipRes.getMip().extensionProbeSequence));
                        System.out.printf("Distance score: %d %d\n", mipRes.dimerDistanceScore[0],mipRes.dimerDistanceScore[1] );      
                        System.exit(0);
                    }    
                    } catch (Exception e) {
                        System.out.print("Error:" + e.getMessage());
                    }
                    */