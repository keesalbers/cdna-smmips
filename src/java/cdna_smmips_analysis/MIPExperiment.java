/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.util.HashMap;
import java.util.HashSet;
import java.util.ArrayList;

/**
 *
 * @author Z009157
 */
public class MIPExperiment {

    
    protected TableData table;
    ArrayList<MIPTarget> mipTargets;
    MIPExperimentProperties mipExperimentProperties;
    MIPExperimentProbeSequenceHash mipExperimentProbeSequenceHash;
    
    protected HashMap<String, MIPTarget> uniqueIDs;
    
    protected HashMap <String, String> columnNames;
    protected HashMap <String, Integer> optionalColumnNames;
    protected int verbose;
    
    public MIPExperiment(String fileName, MIPExperimentProperties _mipExperimentProperties) throws Exception {        
        
        this.verbose = 1;        
        this.mipExperimentProperties = _mipExperimentProperties;
        this.initFromTable(fileName);
    }
    
    final void checkRow(Row row) throws Exception {
         // check for empty required columns
        for (String key : columnNames.keySet()) {
            
            boolean isOptional = columnNames.containsKey(key);
            if (!isOptional && row.get(columnNames.get(key)).toString().length() == 0 ) {
                throw new Exception("Empty field for " + key + " for " + columnNames.get(key) + " in row " + row.asString());
            }
        }
    }
    
    final protected void initFromTable(String fileName) throws Exception {
        mipTargets = new ArrayList<MIPTarget>();        
        this.uniqueIDs = new HashMap<String, MIPTarget>();
        
        this.table = new TableData(fileName);
        this.setColumnLabels();
        this.checkColumnLabels();
        // extract targets from table and create MIPTargets
        
        for (int r = 0; r < table.numRows(); r++) {
            Row row = table.getRow(r);
            checkRow(row);
            MIPTarget mip = new MIPTarget(row, columnNames);
            if (mip.isOK()) {
                this.addMIPTarget(mip);
                if (verbose>1) mip.print();
            } else {
                System.out.print("Error reading MIP:");
                System.out.println(row.asString());
            }
        }
        if (verbose>0) System.out.printf("Number of rows in design file: %d. Number of MIP targets added: %d\n", table.numRows(), mipTargets.size());
        this.mipExperimentProbeSequenceHash = new MIPExperimentProbeSequenceHash(mipTargets, mipExperimentProperties.seedSequenceHashLength);
    }
    
    void addMIPTarget(MIPTarget mip) throws Exception {
        if (this.uniqueIDs.containsKey(mip.uniqueID)) throw new Exception("The MIP IDs are not unique!");
        mip.indexInMIPTargetArray = mipTargets.size();
        mipTargets.add(mip);
        this.uniqueIDs.put(mip.uniqueID, mip);
    }
    
    
    
    
    
    
    final protected void checkColumnLabels() throws Exception {
        if (columnNames == null) {
            throw new Exception("columnNames is null");
        }
        
        for (String key : columnNames.keySet()) {
            boolean hasKey = table.hasColumn(columnNames.get(key));
            boolean isOptional = columnNames.containsKey(key);
            if (!hasKey && !isOptional) {
                throw new Exception("Design file does not have label " + key + " for " + columnNames.get(key));
            }
        }
    }
    public MIPTarget getMIP(String uniqueID) {
        return this.uniqueIDs.get(uniqueID);
    }
    
    public int estimateFragmentLength(String ligationRead, String extensionRead) throws Exception {
        String ligReadStripped = ligationRead.substring(this.mipExperimentProperties.getUmiLengthLigation(), ligationRead.length());
        String extReadStrippedRC = null;
        try {
             extReadStrippedRC = SeqUtil.reverseComplement(extensionRead.substring(this.mipExperimentProperties.getUmiLengthExtension(), extensionRead.length()));
        } catch (Exception e) {
            return -1;
        }
        MakeFragment frag = new MakeFragment(ligReadStripped, extReadStrippedRC);
        if (frag.confidence>0.25) {
            return frag.fragmentLength;
        } else {
            return -1;
        }
    }
    /*
    maps variable names to labels in design file. Can be changed
    */
    final protected void setColumnLabels() {
        columnNames = new HashMap<String, String>();
        optionalColumnNames = new HashMap<String, Integer>();
        columnNames.put("uniqueID", "unique_probe_id");
        columnNames.put("geneName", "hgnc_gene");
        columnNames.put("targetStrand", "probe_is_on_target_strand");
        columnNames.put("exon", "exon");
        columnNames.put("designStrand", "probe_strand");
        columnNames.put("spansEJB", "spans_ejb");
        columnNames.put("isIn3UTR", "in_3utr");
        columnNames.put("extensionProbeSequence", "ext_probe_sequence");
        columnNames.put("ligationProbeSequence", "lig_probe_sequence");
        columnNames.put("targetSequence", "scan_target_sequence");
        columnNames.put("mipSequence", "mip_sequence");
        columnNames.put("chrom", "chr");
        columnNames.put("extensionProbeStart", "ext_probe_start");
        columnNames.put("extensionProbeEnd", "ext_probe_stop");
        columnNames.put("ligationProbeStart", "lig_probe_start");
        columnNames.put("ligationProbeEnd", "lig_probe_stop");        
        
        this.optionalColumnNames.put("exon",1);
        this.optionalColumnNames.put("spansEJB",1);
        this.optionalColumnNames.put("isIn3UTR",1);
        this.optionalColumnNames.put("targetStrand",1);
        this.optionalColumnNames.put("geneName",1);
        
        
    }
    public TableData getTable() {
        return table;
    }

    public ArrayList<MIPTarget> getMipTargets() {
        return mipTargets;
    }

    public MIPExperimentProperties getMipExperimentProperties() {
        return mipExperimentProperties;
    }

    public MIPExperimentProbeSequenceHash getMipExperimentProbeSequenceHash() {
        return mipExperimentProbeSequenceHash;
    }

    public HashMap<String, String> getColumnNames() {
        return columnNames;
    }
              
}
