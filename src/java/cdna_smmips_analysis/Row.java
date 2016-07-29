/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

/**
 *
 * @author caa
 */
import java.util.HashMap;
import java.util.Set;

public class Row {
        Object data[];
        HashMap <String, Integer> labelToColumn;
        public Row(int ncols, HashMap <String, Integer> _labelToColumn) {
            data = new Object[ncols];
            labelToColumn = _labelToColumn;
        }
        public Object get(int col) {
            return data[col];
        }

        public boolean hasColumn(String label) {
            return labelToColumn.containsKey(label);
        }
        public Set<String> getColumnLabels() {
            return labelToColumn.keySet();
        }

        public void set(int col, Object o) {
            data[col]=o;
        }
        public Object get(String key) {
            return data[labelToColumn.get(key)];
        }
        
        public String asString() {
            String x = "";
            for (int  i =0; i < data.length; i++) {
                if (i>0) {
                    x += "\t";
                }
                x += data[i].toString();
            }
            return x;
        }
    }
