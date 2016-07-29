/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
/**
/**
 *
 * @author Z009157
 */
public class ProcessFastqPairs {
    String[] forwardFileNames, reverseFileNames;
    int forwardFileIndex, reverseFileIndex;
    BufferedReader forwardR, reverseR;
    boolean isGzipped;
    public class FastqRecord {
        public String line1, line2, line3, line4;
        public String id;
        public String seq;
        public String qual;
        public int status;
        String[] id_array;
        public FastqRecord(String line1, String line2, String line3, String line4) {
            this.line1 = line1;
            this.line2 = line2;
            this.line3 = line3;
            this.line4 = line4;
            status = init();
        }
        
        final int init() {
            if ((line2 == null || line3 == null || line4 == null) && line1 != null) return -1;
            // check for EOF
            if (line1 == null) return 0; // EOF
            // check if read id starts with @
            if (line1.charAt(0) != '@') return -2; 
            if (line3.charAt(0) != '+') return -3; 
            // we are ok
            this.id_array = line1.split(" ");
            if (id_array.length<1) return -4;
            id = id_array[0];
            seq = line2;            
            qual = line4;
            return 1; // 1 record read from fastq file
        }

        void outputToFile(BufferedWriter bw) throws Exception {
            if (qual.length() != seq.length()) {
                throw new Exception("Unequal length of sequence and quality.");
            }
            if (this.status!=1) {
                throw new Exception("Cannot write Fastq record with status != 1.");
            }
            bw.write(id_array[0]); for (int i = 1; i < id_array.length; i++) bw.write(" " + id_array[i]); bw.write("\n");
            bw.write(seq+"\n");
            bw.write(line3+"\n");
            bw.write(qual+"\n");
        }
        
        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
            this.id_array[0] = id; // note that there may be additional fields in the id line of the fastq file.
        }

        public String getSeq() {
            return seq;
        }

        public void setSeq(String seq) {
            this.seq = seq;
        }

        public int getStatus() {
            return status;
        }

        public void setStatus(int status) {
            this.status = status;
        }
        
        public String getQual() {
            return this.qual;
        }
        
        public void setQual(String newQual) {
            qual = newQual;
        }
        
    }
    
    public ProcessFastqPairs(String forwardFileNames, String reverseFileNames, boolean isGzipped) throws Exception {
        this.forwardFileNames = forwardFileNames.split(",");
        this.reverseFileNames = reverseFileNames.split(",");
        if (this.forwardFileNames.length != this.reverseFileNames.length) {
            throw new Exception("Number of forward and reverse files must be the same.");
        }        
        this.isGzipped = isGzipped;
        for (String fn : this.forwardFileNames) {
            if (fn.endsWith(".gz")) {
                this.isGzipped = true;
                System.out.println("Detected gzipped input files.");
                break;
            }
        }
        for (String fn : this.reverseFileNames) {
            if (fn.endsWith(".gz")) {
                this.isGzipped = true;
                System.out.println("Detected gzipped input files.");
                break;
            }
        }
        initReader();
    }
    
    final void initReader() throws Exception {
        this.forwardFileIndex = 0;
        this.reverseFileIndex = 0;
        forwardR = getReader(forwardFileNames[this.forwardFileIndex]);
        reverseR = getReader(reverseFileNames[this.reverseFileIndex]);
    }
    
    String getCurrentForwardFileName() {
        return forwardFileNames[this.forwardFileIndex];
    }
    
    String getCurrentReverseFileName() {
        return reverseFileNames[this.reverseFileIndex];
    }
    
    BufferedReader getReader(String fileName) throws Exception {
        if (this.isGzipped) {
            FileInputStream fileStream = new FileInputStream(fileName);
            GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
            BufferedReader buffered = new BufferedReader(new InputStreamReader(gzipStream));
            return buffered;
        } else {
            FileInputStream fileStream = new FileInputStream(fileName);            
            BufferedReader buffered = new BufferedReader(new InputStreamReader(fileStream));
            return buffered;
        }
    }
    
    private FastqRecord getFastqRecord(BufferedReader br) throws Exception {
        String line1 = br.readLine();
        String line2 = br.readLine();
        String line3 = br.readLine();
        String line4 = br.readLine();        
        FastqRecord fr = new FastqRecord(line1,line2, line3, line4);
        return fr;
    }
    
    public int getReadFromFastqFiles(FastqRecord[] recordPair) throws Exception {
        if (recordPair.length != 2) { throw new Exception("readPair should be initialized."); }
        
        FastqRecord forward = getFastqRecord(forwardR);
        FastqRecord reverse = getFastqRecord(reverseR);
                
        int fStatus = forward.status;
        int rStatus = reverse.status;
        if (fStatus != rStatus) {
            throw new Exception(String.format("Error: the two fastq %s and %s files were not in sync: %d %d\n", this.getCurrentForwardFileName(), this.getCurrentReverseFileName(), fStatus, rStatus));
        }
        
        if (fStatus == 1 && !forward.id.equals(reverse.id)) {
            throw new Exception(String.format("Error: the IDs between the two fastq files %s and %s are not in sync: %s %s\n", this.getCurrentForwardFileName(), this.getCurrentReverseFileName(), forward.id, reverse.id));
        }
        if (fStatus == 1) {
            recordPair[0] = forward;
            recordPair[1] = reverse;
        }
        if (fStatus == 0) {
            // this is EOF. Check if there are any files left to process.
            this.forwardFileIndex++;
            this.reverseFileIndex++;
            if (this.forwardFileIndex == this.forwardFileNames.length) {
                return 0; // we are done, return EOF
            } else {
                forwardR = getReader(forwardFileNames[this.forwardFileIndex]);
                reverseR = getReader(reverseFileNames[this.reverseFileIndex]);
                System.out.println("Starting to process next pair of fastq files: "+this.forwardFileNames[this.forwardFileIndex]);
                return getReadFromFastqFiles(recordPair);
            }
        }
        
        return fStatus; // should be the same as rStatus
    }
    
   
}
