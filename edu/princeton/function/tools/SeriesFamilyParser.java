// Based on Troilkatt version that was current as of Feb 4, 2011
// with mods by peak@princeton.edu (spelling corrections, and changes at //peak)

//peak:
package edu.princeton.function.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Convert a soft file to a pcl file (per-file)
 */
public class SeriesFamilyParser {	
	// Keywords used to identify columns in the series file
	protected HashMap<String, String[]> keywords = new HashMap<String, String[]>();	

	// During zeros-as-missing-values estimation, zeros with fewer than this number decimals
	// are considered to be true values (and hence not missing values)
	protected int valueDecimals = 5;
	
	// Column indexes
	private int platformProbeIDCol = -1;
	private int geneNamesCol = -1;
	private int sampleProbeIDCol = -1;
	// Expression values column in sample table
	private int expressionValueCol = -1;
	// Column header for expression values column
	private String expressionValueHeader = null;
	
	// Whether to treat zeros as missing values
	private boolean zerosAsMissingVals;
	
	// Filenames
	private String geneMapFilename;
	private String softFilename;
	private String pclFilename;
	private String infoFilename; 

	/**
	 * Constructor
	 *       
	 * @param softFile: series file to convert
	 * @param pclFile: output filename
	 * @param infoFile: output file for meta data
	 * @param geneMap: file with gene names (gene name<tab>gene ID). If null, gene name column
	 * is not predicted based on gene names.
	 */
	public SeriesFamilyParser(String softFile, String pclFile, String infoFile, String geneMap) {	      	     	
		// These keywords are used to identify the various columns of interest in the .soft file.
		// The columns of interest are the column numbers used as arguments in the script that 
		// converts the .soft file to a .pcl file
		keywords.put("Platform probe ID", new String[] {"ID"}); //# List of keywords used to identify the platform probe column number	                          
		// List of keywords used to detect the gene names column number (in prioritized order)                
		//keywords.put("ngn",  {"ORF", "GB_ACC", "GENOME_ACC", "RANGE_GB", "GB_LIST", "Gene_ID", "CLONE_ID", "Genename", "GENE_NAME", "Common Name", "Gene Symbol", "GENE_SYMBOL"});
		// Conservative list
		keywords.put("Gene names", new String[] {"GENE_SYMBOL", "GENE SYMBOL", "GENE_NAME", "GENE NAME", "GENE ID", "ORF", "DDB"});		
		//keywords.put("Gene names", new String[] {"GENE_SYMBOL", "GENE SYMBOL", "ORF", "GB_ACC", "GENOME_ACC", "RANGE_GB", "GB_LIST", "DDB"});
		// List of keywords used to detect the sample probe ID    
		keywords.put("Sample probe ID",  new String[] {"ID_REF"});    
		// List of keywords used to detect the expression values    
		keywords.put("Expression values",  new String[] {"VALUE"});
	
		softFilename = softFile;
		geneMapFilename = geneMap;
		pclFilename = pclFile;
		infoFilename = infoFile;
	}
	
	/**
	 * Constructor
	 * 
	 * See above for arguments
	 */
	public SeriesFamilyParser(String softFile, String pclFile, String infoFile) {
		this(softFile, pclFile, infoFile, null);
	}

	/**
	  * Find columns with gene names, expression values and so on
	  *  
	  * @return: none, but global variables are initialized
	  *
	  * @throws IOException 
	  */    
	protected void findColumns() throws IOException {
		//self.logger.info('Parse SOFT file: %s' % (filename))
		
		BufferedReader fin = new BufferedReader(new FileReader(softFilename));	
	         
		// Number of zero expression values with at least N decimals, and no decimals. 
		// These are used to determine whether to treat zero's as missing values
		int zeroAsValue = 0;
		int zeroAsMissing = 0;
		// Set to true when parsing sample table rows
		boolean countZeros = false;       
		
		// Information about type of zero expression values found in file
		boolean missingExpressionValueCol = false;
		boolean nullExpressionValue = false;
		boolean intExpressionValue = false;
		boolean emptyExpressionValue = false;
		
		// Number of columns found
		int nPlatfomProbeIDColFound = 0;
		int nGeneNameColFound = 0;
		int nSampleProbeIDFound = 0;
		int nExpressionValueColFound = 0;

		/*
		 *   Search for columns with platform probe number and gene names
		 */
		while (true) {
			String line = fin.readLine();			
			if ((line == "") || (line == null)) {
				break;
			}
			line.trim();	        	                        	        
			
	        if (line.contains("!platform_table_begin")) {
	        	// The next line contains the platform header columns
	        	line = fin.readLine();
	        	System.out.print(line);
	        	String[] parts = line.split("\t");

	        	for (String s: keywords.get("Platform probe ID")) {
	        		for (int i = 0 ; i < parts.length; i++) {
	        			String header = parts[i].trim().toUpperCase();
	        			if (header.equals(s)) {
	        				if (platformProbeIDCol == -1) {	               
	        					System.out.println("Platform probe ID col: " + header + " (col: " + i + ")");
	        					platformProbeIDCol = i;	                				
	        				}
	        				else {
	        					System.err.println("Multiple probe IDs col: " + header + " (col: " + i + ")");
	        				}

	        				nPlatfomProbeIDColFound++;
	        			}	        			
	        		}
	        	}
	        	for (String s: keywords.get("Gene names")) {
	        		for (int i = 0 ; i < parts.length; i++) {
	        			String header = parts[i].trim().toUpperCase();
	        			if (header.equals(s)) {
	        				if (geneNamesCol == -1) {	               
	        					System.out.println("Gene names col: " + header + " (col: " + i + ")");
	        					geneNamesCol = i;	                				
	        				}
	        				else {
	        					System.err.println("Multiple gene name columns: " + header + " (col: " + i + ", old: " + geneNamesCol + ")" );
	        				}	        				
	        				nGeneNameColFound ++;
	        			}
	        		}	                	
	        	}
	        }
	        else if (line.contains("!sample_table_begin")) {                
	        	countZeros = true;

	        	// The next line contains the platform header columns
	        	line = fin.readLine();
	        	String[] parts = line.split("\t");	        	
	        	
	        	for (String s: keywords.get("Sample probe ID")) {
	        		for (int i = 0; i < parts.length; i++) {
	        			String header = parts[i].trim().toUpperCase();

	        			if (header.equals(s)) {
	        				if (sampleProbeIDCol == -1) {	               
	        					System.out.println("Sample probe ID col: " + header + " (col: " + i + ")");
	        					sampleProbeIDCol = i;
	        					nSampleProbeIDFound++;
	        				}
	        				else if (sampleProbeIDCol != i) {
	        					System.err.println("Inconsistent sample probe IDs cols: new: " + i + " previous: " + sampleProbeIDCol + ")");
	        					nSampleProbeIDFound++;
	        				}	        											
	        			}
	        		}
	        	}

	        	for (String s: keywords.get("Expression values")) {
	        		for (int i = 0; i < parts.length; i++) {
	        			String header = parts[i].trim().toUpperCase();

	        			if (header.equals(s)) {
	        				if (expressionValueCol == -1) {	               
	        					System.out.println("Expression value col: " + header + " (col: " + i + ")");
	        					expressionValueCol = i;
	        					expressionValueHeader = header;
	        					nExpressionValueColFound++;
	        				}
	        				else if (expressionValueCol != i) {
	        					System.err.println("Inconsistent expression value cols: new: " + i + " previous: " + expressionValueCol + ")");
	        					nExpressionValueColFound++;
	        				}	        											
	        			}
	        		}
	        	}	        	
	        }
	        else if (line.contains("!sample_table_end")) {
	        	countZeros = false;	        	
	        }
	        else if (countZeros) {
	        	// Column header not set
	        	String[] parts = line.split("\t");
	        	if (expressionValueCol == -1) { // First line after table_begin                    
	        		for (int i = 0; i < parts.length; i++) {
	        			if (parts[i].equals(expressionValueHeader)) {
	        				expressionValueCol = i;
	        				nExpressionValueColFound++;
	        				break;
	        			}
	        		}
	        		if (expressionValueCol == -1) {                        
	        			missingExpressionValueCol = true;                    
	        			countZeros = false; // Ignore sample table
	        		}
	        	}
	        	else { // Expression value row
	        		if (parts.length <= expressionValueCol) {
	        			// Ignore
	        		}
	        		else {         	        				        			
	        			String ev = parts[expressionValueCol].trim();
	        			if (ev.equals("")) {
	        				// Is a missing value
	        				emptyExpressionValue = true;
	        			}
	        			else if (ev.equals("null")) {                            
	        				nullExpressionValue = true;
	        			}
	        			else if (ev.equals("0")) {
	        				zeroAsMissing += 1;
	        				intExpressionValue = true;
	        			}
	        			else if (Float.valueOf(ev) == 0.0) {
	        				if ((ev.charAt(1) == '.') && (ev.length() >= (valueDecimals - 1))) {
	        					zeroAsValue += 1;
	        				}
	        				else {
	        					zeroAsMissing += 1;
	        				}
	        			}
	        		}
	        	}
	        } // else if count zero
		} // while true
		
		fin.close();
	        
		if ((geneNamesCol == -1) && (geneMapFilename != null)) {
			HashSet<String> geneSet = initGeneSet(geneMapFilename);
			geneNamesCol = predictGeneNameCol(geneSet);
		}
		                    
		// Create zero-as-missing-value messages
		if (missingExpressionValueCol) {
			System.err.println("Missing expression value column");
		}
		if (nullExpressionValue) {
			System.err.println("'null' expression values found");
		}
		if (intExpressionValue) {
			System.err.println("'0' (no decimals) expression values found");
		}
		if (emptyExpressionValue) {
			System.err.println("Empty expression values found");
		}
		if ((zeroAsValue == 0) && (zeroAsMissing == 0) && (! emptyExpressionValue) && (! nullExpressionValue)) {
			System.out.println("No zero expression values found");
		}

		// Determine if zero values should be treated as missing
		if ((zeroAsValue == 0) && (zeroAsMissing == 0)) {
			zerosAsMissingVals = false;
		}
		else if (zeroAsValue >= zeroAsMissing) {
			zerosAsMissingVals = false;
		}
		else {
			zerosAsMissingVals = true;
		}
		
        if ((nPlatfomProbeIDColFound == 1) && (nGeneNameColFound == 1) &&
        		(nSampleProbeIDFound == 1) && (nExpressionValueColFound == 1)) {            
            System.out.println("Found exactly one probe ID, gene name, sample ID, and expression value column");
        }
        else if ((nPlatfomProbeIDColFound == 0) || (nGeneNameColFound == 0) || 
        		(nSampleProbeIDFound == 0) || (nExpressionValueColFound == 0)) {            
            System.err.printf("One or more columns are missing (pp: %d, gn: %d, sp: %d, ev: %d\n",
            		nPlatfomProbeIDColFound, nGeneNameColFound, nSampleProbeIDFound, nExpressionValueColFound);
            // TODO
            throw new IOException("Could not parse SOFT file");
        }
        else if ((nSampleProbeIDFound > 1) || (nExpressionValueColFound > 1)) {
        	System.err.printf("Inconsistent column names: (sp: %d, ev: %d)\n",
        			nSampleProbeIDFound, nExpressionValueColFound);
        	throw new IOException("Could not parse SOFT file");
        }
        
        System.out.println("Platform probe ID column: " + platformProbeIDCol);
        System.out.println("Gene name column: " + geneNamesCol);
        System.out.println("Sample probe ID column: " + sampleProbeIDCol);
        System.out.println("Expression value column: " + expressionValueCol);
        System.out.println("Treat zeros as missing values: " + zerosAsMissingVals);        
	}        
	 
	/**
	 * Do the soft to pcl file conversion.
	 * 
	 * Note! The parseSOFT must have been called before calling this function.
	 * @throws IOException 
	 */
	public void convertSoft2PCL() throws IOException {
		boolean inPlatTable = false;
		boolean inSampTable = false;
		boolean inSample = false;

		HashMap<String, Integer> platIdIdxs = new HashMap<String, Integer>(); //Per-platform, keep which column has unique IDs
		HashMap<String, Integer> platNameIdxs = new HashMap<String, Integer>(); //Per-platform, keep which column has gene names

		HashMap<String, Integer>sampIdIdxs = new  HashMap<String, Integer>();		//Per-platform, which column of the sample table has unique IDs
		HashMap<String, Integer>sampValueIdxs = new HashMap<String, Integer>();	//Per-platform, which column of the sample table has desired value

		HashMap<String, String> sampPlatform = new HashMap<String, String>();	//per-sample, which platform to use
		HashMap<String, String> sampTitle = new HashMap<String, String>();		//per-sample, title string for sample

		String platID = null;
		String sampID = null;

		int missingValueCount = 0;
		int zeroValueCount = 0;		     
		
		findColumns();
		
		BufferedReader fin = new BufferedReader(new FileReader(softFilename));			
		
		//Do an initial pass through to determine information about platforms, samples, etc.
		while (true) {
			String line = fin.readLine();
			if ((line == "") || (line == null)) {
				break;
			}
			
			line.trim();
			if (line.contains("^PLATFORM")) {
				platID = line.substring(12, line.length());
				//platLines.clear();
			}
	
			if (line.contains("!platform_table_begin")) {
				inPlatTable = true;
			}
			else if (inPlatTable == true) {
				if (line.contains("!platform_table_end")) {
					inPlatTable = false;
				}
				if (platIdIdxs.get(platID) == null) {
					platIdIdxs.put(platID, platformProbeIDCol);
					platNameIdxs.put(platID, geneNamesCol);               	
				}
			}

			if (line.contains("^SAMPLE")) {
				inSample = true;
				sampID = line.substring(10,line.length());
			}
			else if (inSample == true) {
				if (line.contains("!Sample_title")) {
					sampTitle.put(sampID, line.substring(16, line.length()));
				}
				else if (line.contains("!Sample_platform_id")) {
					sampPlatform.put(sampID, line.substring(22,line.length()));
				}
				else if (line.contains("!sample_table_begin")) {
					inSampTable = true;
				}
				else if (line.contains("!sample_table_end")) {
					inSampTable = false;
					inSample = false;
				}
				else if (inSampTable) {															
					if (sampIdIdxs.get(sampPlatform.get(sampID)) == null) {
						sampIdIdxs.put( sampPlatform.get(sampID), sampleProbeIDCol);
						sampValueIdxs.put( sampPlatform.get(sampID), expressionValueCol);
					}
					if (sampIdIdxs.get(sampPlatform.get(sampID)) != null) {
						String val = "";
						String[] parts = line.split("\t");
						int valIndex = sampValueIdxs.get(sampPlatform.get(sampID));

						//peak: a terminating tab indicates a missing value, but java ignores it, so:

						if (valIndex < parts.length) {
						    val = parts[valIndex];
						} else if (valIndex == parts.length && line.endsWith("\t") ) {
						    val = "";
						} else {
						    System.err.println("Invalid line: " + line);
						    continue;
						}

						if (val.equals("")) {
							missingValueCount += 1;
						}						
						else if (! isNumeric(val)) {
							missingValueCount += 1;
						}
						else if (Float.valueOf(val) == 0) {
							zeroValueCount += 1;
						}
					}
				}
			}
		}

		/*
		 * Convert file into PCL format
		 */
		HashMap<String, HashMap<String, String>> samples = new HashMap<String, HashMap<String, String>>(); //Hash from sampID to a Hash from geneName to value
		HashMap<String, HashMap<String, String>> platMaps = new HashMap<String, HashMap<String, String>>();	// Hash from platID to a Hash from uniqueID to geneName
		HashMap<String, String> idToName = null;
		String plat = null;
		HashSet<String >genes = new HashSet<String>();
		boolean first = false;

		
		// No seek in this Java close, so the file is closed and reopened
		fin.close();
		fin = new BufferedReader(new FileReader(softFilename));
		
		//	On a second pass through, actually parse out everything
		while (true) {
			String line = fin.readLine();
			if ((line == "") || (line == null)) {
				break;
			}			
			line.trim();

			if (line.contains("^PLATFORM")) {
				platID = line.substring(12, line.length());
				platMaps.put(platID, new HashMap<String, String>());
			}
			if (line.contains("!platform_table_begin")) {		
				inPlatTable = true;
			}
			else if (inPlatTable == true) {				
				if (line.contains("!platform_table_end")) {
					inPlatTable = false;
				}
				else {
					String[] parts = line.split("\t");
					HashMap<String, String> pm = platMaps.get(platID);
					String pn = "N/A";
					if (platNameIdxs.get(platID) < parts.length) {						
						pn = parts[platNameIdxs.get(platID)];
					}
					else {
						//System.err.println("Warning: platform name not found in line: " + line);
					}
					pm.put(parts[platIdIdxs.get(platID)], pn);					
				}
			}

			if (line.contains("^SAMPLE")) {
				inSample = true;
				sampID = line.substring(10,line.length());
				samples.put(sampID, new HashMap<String, String>());
				plat = sampPlatform.get(sampID);
				idToName = platMaps.get(plat);
			}
			else if (inSample) {
				if (line.contains("!sample_table_begin")) {
					inSampTable = true;
					first = true;
				}
				else if (inSampTable) {
					if (first == true) {
						first = false;
					}					
					else if (line.contains("!sample_table_end")) {
						inSampTable = false;
						inSample = false;
						plat = null;
						idToName = null;
					}
					else {
						String[] parts = line.split("\t");
						String gname = idToName.get(parts[sampIdIdxs.get(plat)]);
						int valueIndex = sampValueIdxs.get(plat);
						if ((gname != null) && (gname != "") && (valueIndex < parts.length)) {
							HashMap<String, String> sm = samples.get(sampID);
							sm.put(gname, parts[valueIndex]);
							genes.add(gname);
						}
						else {
							//$stderr.puts "WARNING: For sample " + sampID + ", " + parts[sampIdIdxs[plat]] + " was not present in the " + sampPlatform[sampID] + " platform table"
						}
					}
				}
			}
		}

		/*
		 * Output the pcl file
		 */
		PrintStream fout = null;
		try {
			fout = new PrintStream(new FileOutputStream(pclFilename));
		} catch (FileNotFoundException e) {
			System.err.println("Could not open PCL file: " + pclFilename + ": " + e.getMessage());
		}

		String[] sampArr = samples.keySet().toArray(new String[samples.size()]);
		Arrays.sort(sampArr);

		fout.print("ID_REF\tIDENTIFIER\tGWEIGHT");
		for (String s: sampArr) {
			fout.print("\t" + sampTitle.get(s));
		}
		fout.println();
		fout.print("EWEIGHT\t\t");
		for (int i = 0; i < sampArr.length; i++) {
			fout.print("\t1");
		}
		fout.println();

		String[] geneArr = genes.toArray(new String[genes.size()]);
		Arrays.sort(geneArr);
		for (String gene: geneArr) {
		        //peak per larsab: do not set gene=gene.toUpperCase() just yet
			String geneParts[] = gene.split(",");
			// TODO: how to handle this situation?
			if (geneParts.length > 1) {				
				gene = geneParts[0];
				for (String p: geneParts) {	
					p = p.toUpperCase();
					if (! p.contains("NO OVERLAP") && ! p.contains("NO HIT") && ! p.isEmpty()) {
						gene = p;
						break;
					}
				}				
			}
			
			if (gene.isEmpty() || gene.equals("")) {
			    continue; //peak (was: gene = "NOT SPECIFIED";)
			}
			//peak per larsab: toUpperCase only on output:
			fout.print(gene.toUpperCase() + "\t" + gene.toUpperCase() + "\t1");
			for (int i = 0; i < sampArr.length; i++) {
				String val = samples.get(sampArr[i]).get(gene);
				//peak: handle "null"
				if (val == null || "null".equals(val) ) {
					fout.print("\t");
				}
				else {
					if (zerosAsMissingVals && (Float.valueOf(val) == 0)) {
						fout.print("\t");
					}
					else {
						fout.print("\t" + val);
					}
				}
			}
			fout.println();

		}
		fout.close();
	}
	
	/**
	 * Create info file
	 * @throws IOException 
	 */
	public void createInfoFile() throws IOException {
		String tagOrg = "!Platform_organism = ";
		String tagPlat = "!Series_platform_id = ";
		String tagDset = "!Series_geo_accession = ";
		String tagTitle = "!Series_title = ";
		String tagDesc = "!Series_summary = "; // may be multiple
		String tagPmid = "!Series_pubmed_id = ";
		String tagFeatCount = "!Platform_data_row_count = ";
		String tagChanCount = "!Sample_channel_count = "; // may be multiple	
		//String tagValType = "!dataset_value_type = "; NOT IN GSE FILES
		String tagModDate = "!Series_last_update_date = ";		
		String tagPlatformTitle = "!Platform_title = ";
		
		//HashMap<String, String>labels = new HashMap<String, String>();
		//boolean inDataSection = false;
		//boolean atLabels = false;

		/*
		 * Parse soft file to get headers
		 */
		BufferedReader fin = new BufferedReader(new FileReader(softFilename));
		
		String org = "";
		String platform = "";
		String dset = "";
		String title = "";
		String desc = "";
		String pmid = "";
		String featCount = "";
		String chanCount = "";
		int sampCount = 0;
		String valType = "";
		String date = "";
		
		while (true) {
			String line = fin.readLine();			
			if (line == null) {
				break;
			}			
			line.trim();
			
			// Parse header
			if (line.contains(tagOrg)) {					
				org = line.substring(tagOrg.length()).replace("\t"," ");
			}
			else if (line.contains(tagPlat)) {
				platform = line.substring(tagPlat.length()).replace("\t", " ");
			}
			else if (line.contains(tagDset)) {
				dset = line.substring(tagDset.length()).replace("\t", " ");
			}
			else if (line.contains(tagTitle)) {
				title = line.substring(tagTitle.length()).replace("\t", " ");
			}
			else if (line.contains(tagDesc)) {
				desc = line.substring(tagDesc.length()).replace("\t", " ");
			}
			else if (line.contains(tagPmid)) {
				pmid = line.substring(tagPmid.length()).replace("\t", " ");
			}
			else if (line.contains(tagFeatCount)) {
				featCount = line.substring(tagFeatCount.length()).replace("\t", " ");
			}
			else if (line.contains(tagChanCount)) {
				chanCount = line.substring(tagChanCount.length()).replace("\t", " ");
			}
			else if (line.contains("^SAMPLE")) {
				sampCount++;
			}
			//else if (line.contains(tagValType)) {
			//	valType = line.substring(tagValType.length()).replace("\t", " ");
			//}
			else if (line.contains(tagModDate)) {
				date = line.substring(tagModDate.length()).replace("\t", " ");
			}
			else if (line.contains(tagPlatformTitle)) {
				System.err.println("Platform title: " + line.substring(tagPlatformTitle.length()));
			}
		}
		fin.close();
		
		/*
		 * Parse pcl file to get statistics
		 */
		BufferedReader fin2 = new BufferedReader(new FileReader(pclFilename));
		
		double min = 1e10;
		double max = -1e10;
		double mean = 0.0;
		double numPos = 0.0;
		double numNeg = 0.0;
		double numZero = 0.0;
		double numMissing = 0.0;
		double numTotal = 0.0;
		
		// skip header line
		fin2.readLine();
		// skip GWEIGHT line
		fin2.readLine();
		while (true) {
			String line = fin2.readLine();			
			if (line == null) {
				break;
			}			
			line.trim();			
			
			String[] parts = line.split("\t");

			// 3 first columns are gene name, gene name, and EWEIGHT
			for (int i = 2; i < parts.length; i++) {
				if (! isNumeric(parts[i])) {
					numMissing += 1;							
				}
				else {
					numTotal += 1;							
					Float val = Float.valueOf(parts[i]);
					mean += val;
					if (val > max) {
						max = val;
					}
					if (val < min) {
						min = val;
					}
					if (val > 0) {
						numPos += 1;
					}
					else if (val < 0) {
						numNeg += 1;
					}
					else {
						numZero += 1;
					}
				}
			}
		}		

		int tested_numChans = 2;
		int tested_logXformed = 1;
		int tested_zerosMVs = 0;
		String tested_MVcutoff = "NA";

		if (numTotal != 0) {
		  mean /= numTotal;
		}
		else {
		  mean = 0;
		}

		if ((numNeg == 0) || ((numPos / numNeg) > 7.5)) {
			tested_numChans = 1;
		}
		if (max > 100) {
			tested_logXformed = 0;
		}
		if (numZero > (5 * numMissing)) {
			tested_zerosMVs = 1;
		}
		if (tested_numChans == 1) {
			if (tested_logXformed == 1) {
				tested_MVcutoff = "0";
			}
			else if (min > -500) {
				tested_MVcutoff = "2";
			}
		}

		PrintStream finfo = new PrintStream(new FileOutputStream(infoFilename));
		finfo.println("File\tDatasetID\tOrganism\tPlatform\tValueType\t#channels\tTitle\tDescription\tPubMedID\t" +
				"#features\t#samples\tdate\tMin\tMax\tMean\t#Neg\t#Pos\t#Zero\t#MV\t#Total\t#Channels\t" +
				"logged\tzerosAreMVs\tMVcutoff");

                //peak 2011.06.28: Support the GSEnnn_setA_family.soft file-naming convention.
		// That is, if the filename begins with the value of "Series_geo_accession" (e.g. GSE123)
		// and contains the string "_family.", then set DatasetID to the substring of the
		// filename prior to "_family.".
		// For standard GEO "soft" files, this prefix is the "Series_geo_accession" anyway.
		// Example: GSE1986_setA_family.soft => GSE1986_setA
		// Rationale:
		// The "DatasetID" is used for matching, and some values in an .info file 
		// are not determined simply by the value of dset
                String softname = new File(softFilename).getName();
		int idx = softname.indexOf("_family.");
		if ( (idx > -1) && softname.startsWith(dset + "_")) {
		    dset = softname.substring(0, idx);
		}
		finfo.print(softname);
		finfo.print("\t" + dset + "\t" + org + "\t" + platform + "\t" + valType + "\t" + chanCount + "\t" + title);
		finfo.print("\t" + desc + "\t" + pmid + "\t" + featCount + "\t" + sampCount + "\t" + date);
		finfo.print("\t" + min + "\t" + max + "\t" + mean + "\t" + numNeg + "\t" + numPos + "\t");
		finfo.print(numZero + "\t" + numMissing + "\t" + numTotal + "\t" + tested_numChans + "\t");
		finfo.println(tested_logXformed + "\t" + tested_zerosMVs + "\t" + tested_MVcutoff);
		finfo.close();
	}
	
	/**
	 * Test if the given string is a valid number
	 * 
	 * @param x: string to test
	 * @return true if a valid floating point number, false otherwise
	 */
	public static boolean isNumeric(String x) {
		try {
			Float.valueOf(x);
			return true;
		} catch (NumberFormatException e) {
			return false;
		}		
	}

	/*
	 * Initialize hashtable with gene names.
	 *
	 * @param genesFile: file with gene names.
	 *
	 * @return: frozenset with gene names. The frozen set provides O(1) search operation.
	 */
	protected HashSet<String> initGeneSet(String mapFile) throws IOException {
		//self.logger.info('Read set of gene names from: %s' % (genesFile))
		HashSet<String> geneSet = new HashSet<String>();
		
		BufferedReader fin = new BufferedReader(new FileReader(mapFile));
	
		while (true) {
			String line = fin.readLine();			
			if (line == null) {
				break;
			}
	
			String[] parts = line.split("\t");
			String gene = parts[0].trim();
			if (! geneSet.contains(gene)) {
				geneSet.add(gene);
			}
		}                  
		fin.close();
		
		return geneSet;
	}

	/**
	 * Predict the column number of the 'gene names' column by comparing the columns in the platform 
	 * table against a set of known gene names. The column with most gene names is assumed to be the
	 * gene names column, 
	 *
	 * @return: column number for predicted gene name column
	 * @throws IOException 
	 */     	    
	protected int predictGeneNameCol(HashSet<String> geneSet) throws IOException {     
		//self.logger.info('Predict gene names column in file: %s' % (inputFile))	        
	
		/*
		 * Find the header for the column in platform_table that has gene names
		*/
		BufferedReader fin = new BufferedReader(new FileReader(softFilename));				
		String headerLine = null;
	
		// Goto platform_table_begin line
		while (true) {
			String line = fin.readLine();
			if ((line == null) || (line == "")) { 
				throw new IOException("Error: no platform table in file");				
			}
	
			if (line.contains("!platform_table_begin")) {            
				headerLine = fin.readLine().trim();                 
				break;
			}
		}
	
		// Parse header line to find number of columns and column headers
		String[] colHeaders = headerLine.split("\t");    
		int nColHeaders = colHeaders.length;
	
		// List with number of genes per column
		int[] geneCount = new int[nColHeaders];
		for (int i = 0; i < nColHeaders; i++) {
			geneCount[i] = 0;
		}
	
		// Count the number of genes per column
		int nRows = 0;
		while (true) {
			String line = fin.readLine();
			if ((line == null) || (line == "")) {
				// TODO: throw specialized exception
				throw new IOException("Error: unexpected end of platform table in file");				
			}
	
			if (line.contains("!platform_table_end")) {
				break;
			}
	
			String[] cols = line.trim().split("\t");        
			for (int i = 0; i < cols.length; i++) {
				if (geneSet.contains(cols[i])) {
					geneCount[i] += 1;
				}
			}
	
			nRows += 1;
		}	        
	
		// Find the header for the column with most gene names
		int maxVal = 0;
		int maxIndex = 0;
		for (int i = 0; i < nColHeaders; i++) {
			if (geneCount[i] > maxVal) {
				maxVal = geneCount[i];
				maxIndex = i;	    		
			}
		}	        
		String maxHeader = colHeaders[maxIndex];
	
		if (maxVal == 0) {
			throw new IOException("Error: no gene names found in platform table");			
		}
		
		System.out.printf("Column %s has max gene count: %d of %d lines (%2.2f%%)\n", maxHeader, maxVal, nRows, (100.0 * maxVal) / nRows);
		return maxIndex;
	}

        //peak
	public static void usage() {
	    System.err.println("Usage: input.soft output.pcl output.info [gene.map]");
	}

        //peak
	public static void help() {
	    System.out.println("Usage: input.soft output.pcl output.info [gene.map]");
	    System.out.println("If the input filename begins with the string \"{ACCESSION}_{NAME}_family.\" where");
	    System.out.println("{ACCESSION} is the value of 'Series_geo_accession' and {NAME} is a character string,");
            System.out.println("then the DatasetID is set to {ACCESSION}_{NAME}");
	    System.out.println("<gene.map> is optional:");
	    System.out.println("  if given, it should be the filename of gene names (gene name<tab>gene ID);");
	    System.out.println("  otherwise, the gene name column is not predicted based on gene names.");
	}

	/**
	 * Debug
	 * 
	 * @param argv: soft-filename pcl-filename info-filename <gene-map>
	 * @throws IOException 
	 */
	public static void main(String[] argv) throws IOException {
		SeriesFamilyParser dset = null;
		if (argv.length > 0
                    && ( argv[0].equals("-h") || argv[0].equals("--help"))) {
		    help();
		    System.exit(0);
		}
		if (argv.length < 3) {
		    usage();
		    System.exit(-1);
		}
		else if (argv.length == 3) {
			dset = new SeriesFamilyParser(argv[0], argv[1], argv[2], null);
		}
		else {
			dset = new SeriesFamilyParser(argv[0], argv[1], argv[2], argv[3]);
		}
		
		dset.convertSoft2PCL();
		dset.createInfoFile();
	}
}
	        
