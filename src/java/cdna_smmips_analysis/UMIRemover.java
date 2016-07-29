/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Z009157
 */
public class UMIRemover {
    ProcessFastqPairs fq;
    BufferedWriter[] matchedReadsOutputFiles;
    int forwardUMILength, reverseUMILength;

    public UMIRemover(String fq1, 
                      String fq2, 
                      String outputFq1, 
                      String outputFq2, 
                      int forwardUMILength, 
                      int reverseUMILength ) throws Exception {
        fq = new ProcessFastqPairs(fq1, fq2, false);
        this.setFilteredReadFileNames(outputFq1, outputFq2);
        this.forwardUMILength = forwardUMILength;
        this.reverseUMILength = reverseUMILength;
    }
    
    final void setFilteredReadFileNames(String forwardFile, String reverseFile) throws Exception {
        this.matchedReadsOutputFiles = new BufferedWriter[2];
        this.matchedReadsOutputFiles[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(forwardFile))));
        this.matchedReadsOutputFiles[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(reverseFile))));
        System.out.printf("NOTE: Writing matched read pairs to %s,%s\n", forwardFile, reverseFile);
    }
    /***
     * Modifies ProcessFastqPairs.FastqRecord[] readPair instance so that UMIs are append to read ID
     * @param readPair
     * @param forwardUMILength
     * @param reverseUMILength
     * @throws Exception 
     */
    static void trimUMI(ProcessFastqPairs.FastqRecord[] readPair, int forwardUMILength, int reverseUMILength) throws Exception {
        String forwardRead = readPair[0].getSeq();
        String reverseRead = readPair[1].getSeq();

        readPair[0].setSeq(forwardRead.substring(forwardUMILength, forwardRead.length()));
        readPair[1].setSeq(reverseRead.substring(reverseUMILength, reverseRead.length()));
        readPair[0].setQual(readPair[0].getQual().substring(forwardUMILength, readPair[0].getQual().length()));
        readPair[1].setQual(readPair[1].getQual().substring(reverseUMILength, readPair[1].getQual().length()));

        String umi = forwardRead.substring(0,forwardUMILength) + ":" + reverseRead.substring(0,reverseUMILength);
        String idForward = readPair[0].getId() + ":"+umi;                
        String idReverse = readPair[1].getId() + ":"+umi;
        if (idForward.length()>254 || idReverse.length() > 254) {
            System.err.print("Read ID: " + idForward + "\n");
            throw new Exception("Read ID combined with UMI is too long.");
        }
        readPair[0].setId(idForward);
        readPair[1].setId(idReverse);
    }
    
    public void run() throws Exception {
        boolean done = false;
        int numRead = 0;
        ProcessFastqPairs.FastqRecord[] readPair = new ProcessFastqPairs.FastqRecord[2];
        
        while (!(done || numRead < -1000) ) {            
            int res = this.fq.getReadFromFastqFiles(readPair);
            if (res == 0) break; // EOF
            else if (res<0) {
                System.err.printf("There was an error: %d\n", res);
            } else if (res == 1) {
                UMIRemover.trimUMI(readPair, this.forwardUMILength, this.reverseUMILength);                
                readPair[0].outputToFile(this.matchedReadsOutputFiles[0]);
                readPair[1].outputToFile(this.matchedReadsOutputFiles[1]);
                numRead++;
            }
        }
        System.out.printf("%d fastq record read and written\n" , numRead);
        this.matchedReadsOutputFiles[0].close();
        this.matchedReadsOutputFiles[1].close();
    }
}
