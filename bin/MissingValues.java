/**
 * This version was in sync with SVN troilkatt repository (svn://gen-svn.princeton.edu/hefalmp) as of 2011.05.27
   package edu.princeton.function.troilkatt.tools;
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author larsab
 *
 */
public class MissingValues {
	// genes with # missing values > (# expression values / GENE_CUTOFF_FACTOR) are discarded
	public float geneCutoffFactor;
	// datasets with #missing values > (# expression values / DATASET_CUTOFF_FACTOR) are discarded
	public float datasetCutoffFactor;
	// datasets with less than SAMPLE_CUTOFF_FACTOR are discarded
	public int sampleCutoffFactor;
	
	/*
	 * INFO file values
	 */
	// Are zero values treated as missing?
	public boolean zerosAsMissingValues;
	// Values less than the cutoff are set as missing
	public float missingValueCutoff;
	
	/**
	 * Constructor. 
	 * 
	 * @param inputFilename PCL file to convert
	 * @param infoFilename meta data file with information about the dataset to convert
	 * @param geneCutoff maximum allowed missing values per gene (in percentage of total values). Genes
	 * with more missing values are discarded.
	 * @param sampleCutoff minimum samples requires for a dataset. Datasets with fewer samples are
	 * discarded.
	 * @param datasetCutoff maximum allowed missing values per dataset (in percentage of total values). Genes
	 * with more missing values are discarded.
	 * @throws IOException if .info file cannot be read.
	 */
	public MissingValues(String inputFilename, String infoFilename,
			float geneCutoff, int sampleCutoff, float datasetCutoff) throws IOException {
		String dsetID = getId(inputFilename);
		if (! parseInfoFile(infoFilename, dsetID)) {
			System.err.println("Dataset ID not found in info file: " + dsetID);
			System.exit(-1);
		}
		
		geneCutoffFactor = geneCutoff;
		sampleCutoffFactor = sampleCutoff;
		datasetCutoffFactor = datasetCutoff;
		
		System.out.println("Discard genes with more than " + geneCutoff + "% missing values.");
		System.out.println("Discard datasets with more than " + datasetCutoff + "% missing values.");
		System.out.println("Discard datasets with less than " + sampleCutoff + " samples.");
	}

	/**
	 * Get dataset ID by parsing filename. It is assumed that the filename without extensions
	 * is the dataset ID.
	 * 
	 * @param inputFilename: input filename.
	 * @param return dataset ID
	 */
	public String getId(String inputFilename) {
		String basename = new File(inputFilename).getName();
		return basename.split("\\.")[0];
	}
	
	
	/**
	 * Parse the info file to initialize global variables: numberOfChannels, zeroAsmissingValues,
	 * and missingValueCutoff
	 * 
	 * @param infoFilename: info filename
	 * @param id: dataset id
	 * @return true if global variables where successfully initialized, false otherwise
	 * @throws IOException 
	 */
	public boolean parseInfoFile(String infoFilename, String id) throws IOException {
		// Important configuration columns
		int iDCol = 1;		
		int zeroMVCol = 22;
		int mvCutoffCol = 23;
				
		BufferedReader ib = new BufferedReader(new FileReader(infoFilename));
		
		// Read the configuration file and locate the needed info
		String line;
		boolean entryFound = false;
		while ((line = ib.readLine()) != null) {
			String[] cols = line.trim().split("\t");
			if (cols[iDCol].equals(id)) {				
				entryFound = true;				
				if (Integer.valueOf(cols[zeroMVCol]) == 0) {
					zerosAsMissingValues = false;
				}
				else {
					zerosAsMissingValues = true;
				}
				try {
					missingValueCutoff = Float.valueOf(cols[mvCutoffCol]);
				}
				catch (NumberFormatException e) {
					missingValueCutoff = Float.NaN;
				}
				
				System.out.println("Zeros as missing values: " + zerosAsMissingValues);
				System.out.println("Missing value cutoff: " + missingValueCutoff);				
				break;
			}
		}
		
		ib.close();
		return entryFound;
	}
	
	/**
	 * Insert missing values and filter out rows with too many missing values.
	 * @throws IOException 
	 * 
	 */
	public void process(String inputFilename, String outputFilename) throws IOException {
		String tmpFilename = inputFilename + ".tmp";
		
		int maxAllowedMissingGenes = -1; // Maximum number of values that can be missing in a row		
		int nCols = -1;                  // number of columns per row (set after parsing first line)
		
		BufferedReader ib = new BufferedReader(new FileReader(inputFilename));
		// Write data first to temporary file
		PrintWriter os = new PrintWriter(new FileWriter(tmpFilename));
		
		String line;		
		int lineCnt = 0;
		int missingTotal = 0;   // number of missing values in included rows
		int totalValues = 0;    // nubmer of expression values in included rows
		int discardedRows = 0;  // number of discarded rows
		
		while ((line = ib.readLine()) != null) {
			lineCnt++;
			if (lineCnt == 1) {
				nCols = line.split("\t").length;	
				maxAllowedMissingGenes = (int) Math.ceil( (nCols - 3) * (geneCutoffFactor / 100));
				// Make sure there are enough samples
				if (nCols < sampleCutoffFactor) {
					System.err.printf("Too few samples in dataset: %d (%d required)\n", nCols, sampleCutoffFactor);
				}
			}			
			
			// Skip the two header lines
			if (lineCnt < 3) { 				
				os.write(line + "\n");				
				continue;
			}
			
			String[] cols = line.trim().split("\t");
			// The three first columns are headers
			if (cols.length < 3) {
				System.err.println("Warning: invalid PCL table row: " + line);
				continue;
			}
			
			int missingInRow = 0;
			String outputLine = cols[0] + "\t" + cols[1] + "\t" + cols[2];
			for (int i = 3; i < cols.length; i++) {
				try {
					Float val = Float.valueOf(cols[i]);
					if (zerosAsMissingValues && (val == 0)) {
						missingInRow++;							
						outputLine = outputLine + "\t";
					}										
					else if ((Float.isNaN(missingValueCutoff)) || (val >= missingValueCutoff)) {
						outputLine = outputLine + "\t" + val;							
					}
					else {
						missingInRow++;							
						outputLine = outputLine + "\t";
					}
				} catch (NumberFormatException e) {
					missingInRow++;						
					outputLine = outputLine + "\t";
				}
			}
			
			// Set non-existing columns to missing
			for (int i = cols.length; i < nCols; i++) {
				missingInRow++;							
				outputLine = outputLine + "\t";
			}
				
			// Discard lines with too many missing values
			if (missingInRow >= maxAllowedMissingGenes) { 
				// More than 50% of values in row are missing: discard line
				discardedRows++;
			}
			else {
				// Do not include discarded lines in the missing count
				missingTotal += missingInRow;
				totalValues += nCols - 3;
				os.write(outputLine + "\n");
			}	
		} // while more lines
		
		os.close();
		ib.close();
		
		int maxAllowedMissingInDataset = (int) Math.ceil(totalValues * (datasetCutoffFactor / 100));		
		if (missingTotal >= maxAllowedMissingInDataset) {
			System.err.printf("Discarding dataset: too many missing values: %d of %d\n",
					missingTotal, totalValues);
		}
		else {
			System.err.printf("Discarded %d of %d genes with too many missing values\n",
					discardedRows, lineCnt);
			System.err.printf("Of remaining, %d of %d values are missing\n", missingTotal, totalValues);
			File tmpFile = new File(tmpFilename);
			tmpFile.renameTo(new File(outputFilename));
		}

	}

	public static void usage() {	  			
		System.out.println("Required arguments:");
		System.out.println("0 - pcl file to insert missing values into");
		System.out.println("1 - configuration info file");
		System.out.println("2 - output filename");
		System.out.println("3 - max allowed missing values per gene (percentage of values for gene)");
		System.out.println("4 - minimum required samples (integer)");
		System.out.println("5 - max allowed missing values in dataset(percentage of total values)");		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 6) {
			usage();
			System.exit(-1);
		}
		String inputFilename = args[0];
		String infoFilename = args[1];
		String outputFilename = args[2];
		float geneCutoff = Float.valueOf(args[3]);
		int sampleCutoff = Integer.valueOf(args[4]);
		float datasetCutoff = Float.valueOf(args[5]);
				
		MissingValues mv = new MissingValues(inputFilename, infoFilename, 
				geneCutoff, sampleCutoff, datasetCutoff);
		mv.process(inputFilename, outputFilename);
	}
}
