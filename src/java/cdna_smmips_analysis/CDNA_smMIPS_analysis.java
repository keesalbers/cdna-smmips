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
public class CDNA_smMIPS_analysis {
    /**
     * @param args the command line arguments
     */
    
    public static void mainCountMIPS(String[] args) throws Exception {
        HashMap<String, String> parsedArgs = new HashMap<String, String>();
        String[] required = {"Fastq1","Fastq2","inputDesignFile","outputFile","ligationUmiLength","extensionUmiLength","blockSize","seedSequenceHashLength"};
        
        String[] optional = {"BLOCKCOUNTUMIS",
            "thresholdExtLigDimerDistance",
            "readNumBlocks",
            "OUTPUTMATCHED",
            "OUTPUTUNMATCHED",
            "removeUMIs",
            "kmerHistogramOutputFile",
            "filterReadsUsingKmers",
            "histogramKmerLength",
            "fractionReadsExplained" };
        
        for (String s : args) {
            String[] pars = s.split("=");
            if (pars.length==2) {
                parsedArgs.put(pars[0].toUpperCase(), pars[1]);
            } else if (pars.length == 1) {
                parsedArgs.put(pars[0].toUpperCase(),"PRESENT");
            } else {
                System.err.print("Did not understand argument "+s);
                System.err.print("Required arguments: "+required.toString());
                System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz");
                System.err.print("Options: " + optional.toString());
                System.exit(1);
            }
        }
        System.out.printf("Read arguments: %s\n", parsedArgs.keySet().toString());
        for (String s : required) {
            if (!parsedArgs.containsKey(s.toUpperCase())) {
                System.err.println("Required arguments: ");
                for (String k : required) {
                    System.err.println("\t"+k);
                }
                System.err.println("Optional arguments: ");
                for (String k : optional) {
                    System.err.println("\t"+k);
                }
                System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz\n");
                System.err.print("Use GZIP option to force reading gzipped input.\n");
                System.exit(1);
            }
        }


        boolean blockCountUMIs = false;
        boolean isGzipped = false;
        
        String fq1 = parsedArgs.get("FASTQ1");
        String fq2 = parsedArgs.get("FASTQ2");
        String inputFile = parsedArgs.get("INPUTDESIGNFILE");
        String outputFile = parsedArgs.get("OUTPUTFILE");
        int ligationUMILength = Integer.decode(parsedArgs.get("ligationUMILength".toUpperCase()));
        int extensionUMILength = Integer.decode(parsedArgs.get("extensionUMILength".toUpperCase()));
        int blockSize = Integer.decode(parsedArgs.get("blockSize".toUpperCase()));
        int seedSequenceHashLength = Integer.decode(parsedArgs.get("seedSequenceHashLength".toUpperCase()));
        if (seedSequenceHashLength > 16) {
            seedSequenceHashLength = 16;
            System.err.print("seedSequenceHashLength out of bounds. Setting to "+seedSequenceHashLength +"\n");
        }
        if (seedSequenceHashLength < 0) {
            seedSequenceHashLength = -1;
            System.out.print("Doing exhaustive matching of probe sequences\n");
        }

        System.out.printf("Read arguments:\n");
        for (String s : parsedArgs.keySet()) {
            System.out.printf("\t%s : %s\n",s, parsedArgs.get(s));
        }
     
        if (parsedArgs.containsKey("BLOCKCOUNTUMIS") && Integer.decode(parsedArgs.get("BLOCKCOUNTUMIS"))==1) {
            blockCountUMIs = true;
            System.out.print("Counting total number of UMIs to determine blocks.\n");
        } else {
            System.out.print("Counting total number of sequenced reads to determine blocks.\n");
        }

        
        if (parsedArgs.containsKey("GZIP")) isGzipped = true;
        
        


        MIPExperiment mipExp = new MIPExperiment(inputFile, new MIPExperimentProperties(ligationUMILength,extensionUMILength, seedSequenceHashLength));
        ProcessFastqPairs fq = new ProcessFastqPairs(fq1, fq2, isGzipped);
        MIPAnalysisParameters mipAnalysisParameters = new MIPAnalysisParameters();
        
        if (parsedArgs.containsKey("thresholdExtLigDimerDistance".toUpperCase())) {
            mipAnalysisParameters.thresholdExtLigDimerDistance = Integer.decode(parsedArgs.get("thresholdExtLigDimerDistance"));
        }
        
        if (parsedArgs.containsKey("readNumBlocks".toUpperCase())) {
            mipAnalysisParameters.read_num_blocks = Integer.decode(parsedArgs.get("readNumBlocks".toUpperCase()));
        }
        
        if (parsedArgs.containsKey("fractionReadsExplained".toUpperCase())) {
            mipAnalysisParameters.fractionReadsExplained = Double.valueOf(parsedArgs.get("fractionReadsExplained".toUpperCase()));
            if (mipAnalysisParameters.fractionReadsExplained < 0.90) mipAnalysisParameters.fractionReadsExplained = 0.90;
            if (mipAnalysisParameters.fractionReadsExplained > 1.00) mipAnalysisParameters.fractionReadsExplained = 1.00;           
        }
        System.out.printf("fractionReadsExplained: %1.2f\n", mipAnalysisParameters.fractionReadsExplained);       
         
        mipAnalysisParameters.blockSize = blockSize;
        mipAnalysisParameters.blockCountUMIs = blockCountUMIs;        
        mipAnalysisParameters.outputMatched = false;
        mipAnalysisParameters.outputUnmatched = false;
        
        if (parsedArgs.containsKey("OUTPUTMATCHED")) {
            String[] outputFilteredFiles = parsedArgs.get("OUTPUTMATCHED").split(",");
            mipAnalysisParameters.forwardFileMatched = outputFilteredFiles[0];
            mipAnalysisParameters.reverseFileMatched = outputFilteredFiles[1];       
            mipAnalysisParameters.outputMatched = true;
        } else if (parsedArgs.containsKey("OUTPUTUNMATCHED")) {
            String[] outputFilteredFiles = parsedArgs.get("OUTPUTUNMATCHED").split(",");            
            mipAnalysisParameters.outputUnmatched = true;
            mipAnalysisParameters.forwardFileUnmatched = outputFilteredFiles[0];
            mipAnalysisParameters.reverseFileUnmatched = outputFilteredFiles[1];
        }
        
        if (parsedArgs.containsKey("kmerHistogramOutputFile".toUpperCase())) {
            mipAnalysisParameters.kmerHistogramOutputFile = parsedArgs.get("kmerHistogramOutputFile".toUpperCase());
            mipAnalysisParameters.createMIPKmerHash = true;
            System.out.printf("kmerHistogramOutputFile specified: %s\n",  mipAnalysisParameters.kmerHistogramOutputFile );
        }
                
        if (parsedArgs.containsKey("histogramKmerLength".toUpperCase())) {
            mipAnalysisParameters.mipTargetKMERSize = Integer.decode(parsedArgs.get("histogramKmerLength".toUpperCase()));
            if (!parsedArgs.containsKey("kmerHistogramOutputFile".toUpperCase()) && !parsedArgs.containsKey("filterReadsUsingKmers".toUpperCase())) {
                System.err.println("You probably want to specify histogramKmerLength or filterReadsUsingKmers option.");
                System.exit(1);
            }
        }       
        
        if (parsedArgs.containsKey(("removeUMIs".toUpperCase()))) {
            mipAnalysisParameters.removeUMIs = true;
            System.out.println("Removing UMIs in output to Fastq file");
        }
        
        if (parsedArgs.containsKey("filterReadsUsingKmers".toUpperCase())) {     
            mipAnalysisParameters.filterReadsUsingKMERs = true;
            mipAnalysisParameters.createMIPKmerHash = true;
        }
        
        MIPProbeCoverageAnalysisBlock analysis = new MIPProbeCoverageAnalysisBlock(mipExp, fq,  mipAnalysisParameters);
        
        boolean suppressFastqOutput = false;
        if (mipAnalysisParameters.filterReadsUsingKMERs) {
            suppressFastqOutput = true;
            // the following options need to be set so that a kmer has is created using all the reads in the fastq files,
            // and not just using the first block.
        }   
        analysis.analyse(suppressFastqOutput);
        //analysis.printKmerHashSummary();
        
        if (mipAnalysisParameters.filterReadsUsingKMERs) {
            mipAnalysisParameters.blockCountUMIs = blockCountUMIs;
            if (parsedArgs.containsKey("readNumBlocks".toUpperCase())) {
                mipAnalysisParameters.read_num_blocks = Integer.decode(parsedArgs.get("readNumBlocks".toUpperCase()));
            }
            System.out.printf("Rerunning and applying kmer-based filtering\n");

            // set parameters for kmer-based filtering
            mipAnalysisParameters.createMIPKmerHash = false;
            mipAnalysisParameters.filterReadsUsingKMERs = true;
            mipAnalysisParameters.statsWithKMERFreqsForFiltering = analysis.targetStatisticsForKmerCounting;
            suppressFastqOutput = false;
            // create new analysis instance
            fq = new ProcessFastqPairs(fq1, fq2, isGzipped);
            MIPProbeCoverageAnalysisBlock kmerAnalysis = new MIPProbeCoverageAnalysisBlock(mipExp, fq,  mipAnalysisParameters);
            kmerAnalysis.analyse(suppressFastqOutput);
            kmerAnalysis.outputToTables(outputFile);
        } else {
            analysis.outputToTables(outputFile);
        }
        
    }
    
    public static void mainRemoveUMIs(String[] args) throws Exception {
        HashMap<String, String> parsedArgs = new HashMap<String, String>();
        String[] required = {"Fastq1","Fastq2","outputFastq1","outputFastq2","ligationUmiLength","extensionUmiLength"};
        for (String s : args) {
            String[] pars = s.split("=");
            if (pars.length==2) {
                parsedArgs.put(pars[0].toUpperCase(), pars[1]);
            } else if (pars.length == 1) {
                parsedArgs.put(pars[0],"PRESENT");
            } else {
                System.err.print("Did not understand argument "+s);
                System.err.print("Required arguments: "+required.toString());
                System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz");
                System.exit(1);
            }
        }
        System.out.printf("Read arguments: %s\n", parsedArgs.keySet().toString());
        for (String s : required) {
            if (!parsedArgs.containsKey(s.toUpperCase())) {
                System.err.println("Required arguments: ");
                for (String k : required) {
                    System.err.println("\t"+k);
                }
                System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz\n");
                System.err.print("Use GZIP option to force reading gzipped input.\n");
                System.exit(1);
            }
        }


        
        String fq1 = parsedArgs.get("FASTQ1");
        String fq2 = parsedArgs.get("FASTQ2");
        String outputFq1 = parsedArgs.get("OUTPUTFASTQ1");
        String outputFq2 = parsedArgs.get("OUTPUTFASTQ2");
        int ligationUMILength = Integer.decode(parsedArgs.get("ligationUMILength".toUpperCase()));
        int extensionUMILength = Integer.decode(parsedArgs.get("extensionUMILength".toUpperCase()));
       
        System.out.printf("Read arguments:\n");
        for (String s : parsedArgs.keySet()) {
            System.out.printf("\t%s : %s\n",s, parsedArgs.get(s));
        }

        ProcessFastqPairs fq = new ProcessFastqPairs(fq1, fq2, false);
       
        UMIRemover umiRemover = new UMIRemover(fq1, fq2, outputFq1, outputFq2, ligationUMILength, extensionUMILength );
        umiRemover.run();
        
    }
    public static void main(String[] args) {
        // TODO code application logic here
        
        
        try {
            HashMap<String, String> parsedArgs = new HashMap<String, String>();
            String[] required = {"PROGRAM"};
            for (String s : args) {
                String[] pars = s.split("=");
                if (pars.length==2) {
                    parsedArgs.put(pars[0].toUpperCase(), pars[1]);
                } else if (pars.length == 1) {
                    parsedArgs.put(pars[0],"PRESENT");
                } else {
                    System.err.print("Did not understand argument "+s);
                    System.err.print("Required arguments: "+required.toString());
                    System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz");
                    System.exit(1);
                }
            }
            System.out.printf("Read arguments: %s\n", parsedArgs.keySet().toString());
            for (String s : required) {
                if (!parsedArgs.containsKey(s.toUpperCase())) {
                    System.err.println("Required arguments: ");
                    for (String k : required) {
                        System.err.println("\t"+k);
                    }
                    System.err.print("Specify each option as follows:  Fastq1=/data/file.R1.fastq.gz\n");
                    System.err.print("Use GZIP option to force reading gzipped input.\n");
                    System.exit(1);
                }
            }
            
         
            String program = parsedArgs.get("PROGRAM");
            
            if (program.equals("COUNTUMIS")) {
                mainCountMIPS(args);
            } else if (program.equals("REMOVEUMIS")) {
                mainRemoveUMIs(args);
            } else {
                System.err.println("Unrecognized program. Options are COUNTUMIS");
                System.exit(1);
            }
            
        }
        catch (Exception e) {
            System.err.print("Exception:\n");
            e.printStackTrace();
        }
    }
    
}
