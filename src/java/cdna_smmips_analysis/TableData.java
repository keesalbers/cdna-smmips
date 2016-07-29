/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.io.*;
/**
 *
 * @author caa
 */


public class TableData {
    String fileName = "";
    String name = "";

    private int numColumns = 0, numRows = 0;

    ArrayList< Row > data;
    String[] headerLabels;
    HashMap <String, Integer> labelToColumn;

 
    public TableData(String _fileName) throws Exception {
        fileName = _fileName;
        name = fileName;
        data = new ArrayList< Row >();

        loadFromFile();

    }

    public TableData( String[] _headerLabels, String _name) {
        headerLabels = new String[_headerLabels.length];
        numColumns = headerLabels.length;

        System.arraycopy(_headerLabels, 0, headerLabels, 0, _headerLabels.length);
        name = _name;
        labelToColumn = new HashMap <String, Integer>();
        for (int i=0;i<headerLabels.length;i++) labelToColumn.put(headerLabels[i],i);

        data = new ArrayList< Row >();

    }
    public String[] getHeaderLabels() {
        return headerLabels;
    }

    public int numColumns() {
        return numColumns;
    }

    public int numRows() {
        return data.size();
    }

    public boolean hasColumn(String label) {
        return labelToColumn.containsKey(label);
    }

    public Object getData(int row, String label) {
        return data.get(row).get(label);
    }

    public Object getData(int row, int col) {
        return data.get(row).get(col);
    }

    public Row getRow(int row) {
        return data.get(row);
    }

    public String getRowAsString(int row, String sep) {
        String s = "";
        int n=0;
        for (String lab : headerLabels) {
            if (n>0) s += sep;
            n++;
            s += getData(row, lab).toString();
        }
        return s;
     }
    public String getHeaderString(String sep) {
        String s = "";
        int n=0;
        for (String lab : headerLabels) {
            if (n>0) s += sep;
            n++;
            s += lab;
        }
        return s;
     }
    public String getName() {
        return name;
    }

    public void addRow(HashMap<String, Object> row_data) {
        Row row = new Row(numColumns, labelToColumn);
        for (String k : row_data.keySet()) {
            if (labelToColumn.containsKey(k)) {

                int col = labelToColumn.get(k);
                //System.out.println("Data " + k + " " + String.valueOf(col));
                row.set(col, row_data.get(k));
            }
        }

        data.add(row);
    }

    public void addRow(Object[] row_data) {
        if (row_data.length == numColumns) {
            Row row = new Row(numColumns, labelToColumn);

            for (int k=0;k<row_data.length; k++) {
                row.set(k, row_data[k]);
            }
            data.add(row);
        }
    }

    // loading function
    private void loadFromFile() throws Exception {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            int numRows = 0;
            try {
                String line = null; //not declared within while loop
                /*
                * readLine is a bit quirky :
                * it returns the content of a line MINUS the newline.
                * it returns null only for the END of the stream.
                * it returns an empty String if two newlines appear in a row.
                */

                // First line is header

                boolean foundHeader = false;
                while (!foundHeader) {
                    line = br.readLine();
                    if ((line.charAt(0)=='#' && line.charAt(1) != '#') || (line.charAt(0)!='#') ) {
                        foundHeader = true;
                    }
                }
                if (line != null) { 
                    headerLabels = line.split("\t");
                    int i=0;
                    System.out.print("Labels: ");
                    labelToColumn = new HashMap<String, Integer> ();
                    for (String lab : headerLabels) {
                        String flab = lab.replace("#", "");
                        labelToColumn.put(flab, new Integer(i));
                        System.out.print(flab+",");
                        i++;
                    }
                    System.out.println("\n");

                    numColumns = headerLabels.length;
                    System.out.println("Number of columns in header: " + numColumns);

                    int shown = 0;
                    while (( line = br.readLine()) != null){
                        String dat[] = line.split("\t",-1);
                        if (shown<2) {
                            for (Object d : dat) {
                                System.out.print(d.toString()+",");
                            }
                            System.out.print("\n");
                            shown++;
                        }
                        //System.out.println("dat.length: " + dat.length + " LINE: " + line);
                        if (dat.length == numColumns) {
                            addRow(dat);
                        } else {
                            System.err.println("Skipping row. Num columns: " + dat.length  + " line: " + line);
                        }
                    }
                }
              }
              finally {
                br.close();
              }
       



    }
}
