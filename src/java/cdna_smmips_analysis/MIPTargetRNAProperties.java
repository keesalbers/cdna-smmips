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
import java.lang.NumberFormatException;
public class MIPTargetRNAProperties {
    Boolean spansEJB, isIn3UTR;
    Integer exon; 
    String geneName;
    
    public MIPTargetRNAProperties(Row row, HashMap <String, String> columnNames) throws Exception {
        this.init(row, columnNames);
    }
    
    final void init(Row row, HashMap <String, String> columnNames) throws Exception {
        this.spansEJB = row.get(columnNames.get("spansEJB")).toString().toUpperCase().equals("TRUE");
        this.spansEJB = row.get(columnNames.get("isIn3UTR")).toString().toUpperCase().equals("isIn3UTR");
        try {
        this.exon = Integer.decode(row.get(columnNames.get("exon")).toString());
        } catch (NumberFormatException e) {
            this.exon = 0;
        }
        this.geneName = row.get(columnNames.get("geneName")).toString();
    }
}
