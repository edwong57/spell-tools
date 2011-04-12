import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

/*
 * Copyright (C) 2005-7 Matt Hibbs
 * License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
 * See http://creativecommons.org/licenses/by-nc/3.0/
 * Attribution shall include the copyright notice above.
 *
 * Created on Jun 9, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/**
 * @author mhibbs
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class DivLogNorm {

	public class Gene  {
		public String 	orf;
		public String 	name;
	    public double[] exprLvls;
	    public int		id;
	    public ArrayList	utilArray;

	    public Gene (String _orf, String _name) {
	    	orf = _orf;
	    	name = _name;
	    	utilArray = new ArrayList();
	    }   
	}
	
	
	public ArrayList genes;
	public ArrayList exprNames;
	int numExprs;

	public DivLogNorm() {
		genes = new ArrayList();
		exprNames = new ArrayList();
		numExprs = 0;
	}

	/** Clears out all data, and reads in new data in a pcl file from the
	 *  buffered reader provided
	 */
	public boolean readDataFromPCL(BufferedReader fin) throws IOException {
		// Clear out all arrays and re-initialize
		genes.clear();
		numExprs = 0;
		int geneID = 0;

		// Try to read in the data from the buffered reader
		try {
			int row = 0;
			int gweightCol = -1;
			int colOffset = 2;

			//While there's more data to read in...
			String line = fin.readLine();
			while ((line != null) && (!line.equals(""))) {
				boolean processRow = true;

				//If this is the first row, this determines number of experiments
				if (row == 0) {
					String[] names = line.split("\t");
					for (int i = 2; i < names.length; i++) {
						//If there is a GWEIGHT column, make a note
						if (names[i].equalsIgnoreCase("GWEIGHT")) {
							gweightCol = i;
							colOffset = 3;
						}
						//Otherwise, record the experiment name
						else {
							exprNames.add(names[i]);
							numExprs++;
						}
					}
					if (names.length < 4) {
						return false;
					}
					processRow = false;
				}

				//If this is the second row, it may contain weights, or it may
				//contain data
				else if (row == 1) {
					String[] lineParts = line.split("\t");

					//If the first column is EWEIGHT, then just ignore this row
					if (lineParts[0].equalsIgnoreCase("EWEIGHT")) {
						processRow = false;
					}
				}

				//If the line is up for processing...
				if (processRow) {
					String[] lineParts = line.split("\t");

					if (lineParts.length > 1) {
						//Columns 0 and 1 are the SpotIDs and Names
						Gene gene = new Gene(lineParts[0].toUpperCase(), lineParts[1].toUpperCase());
						gene.id = geneID++;
						gene.exprLvls = new double[numExprs];
						genes.add(gene);
					
						//Columns from 2 to the end are data points (except the
						// gweightCol)
						for (int i = 2; i < lineParts.length; i++) {
							if (i != gweightCol) {
								//Try to convert the string to a number, if you
								// fail, just fill in with NaNs;
								try {
									double val = Double.parseDouble(lineParts[i]);
									if (!Double.isNaN(val)) {
										gene.exprLvls[i - colOffset] = val;
									}
									else {
										gene.exprLvls[i - colOffset] = Double.NaN;
									}

								}
								catch (NumberFormatException except) {
									gene.exprLvls[i - colOffset] = Float.NaN;
								}
							}
						}
						//If there are remaining holes, they should be NaNs
						for (int i=(lineParts.length - colOffset); i<numExprs; i++) {
							gene.exprLvls[i] = Float.NaN;
						}
					}
				}

				//Increment row number and repeat
				row++;
				line = fin.readLine();
			}

		}
		catch (IOException except) {
			throw (except);
		}
		
		//Finally, return success
		return true;
	}

	public void logTransformAll() {
		for (int i = 0; i < genes.size(); i++) {
			Gene g = (Gene) genes.get(i);

			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j])) {
					if (g.exprLvls[j] < 0)
						g.exprLvls[j] = -Math.log(-g.exprLvls[j]) / Math.log(2);
					else if (g.exprLvls[j] > 0)
						g.exprLvls[j] = Math.log(g.exprLvls[j]) / Math.log(2);
					else { //exactly zero entries
						g.exprLvls[j] = Math.log(0.001) / Math.log(2);
					}
				}
			}
		}
	}

	public void normalizeAllGenes() {
		for (int i = 0; i < genes.size(); i++) {
			Gene g = (Gene) genes.get(i);
			//Find the mean and standard deviation
			double mean = 0;
			double stdDev = 0;
			int count = 0;
			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j])) {
					mean += g.exprLvls[j];
					count++;
				}
			}
			if (count > 0) mean /= (float)count;
			
			count = 0;
			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j])) {
					stdDev += (g.exprLvls[j] - mean) * (g.exprLvls[j] - mean);
					count++;
				}
			}
			if (count > 0) stdDev = Math.sqrt(stdDev / (float)count);

			//Mean shift to 0, make deviaiton 1
			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j])) {
					g.exprLvls[j] -= mean;
					if (stdDev > 0) {
						g.exprLvls[j] /= stdDev;
					}
				}
			}
		}
	}

	public void divideAllGenesByIndividualMedians() {
		for (int i=0; i<genes.size(); i++) {
			Gene g = (Gene) genes.get(i);
			ArrayList exprs = new ArrayList();
			for (int j=0; j<g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j]))
					exprs.add(g.exprLvls[j]);
			}
			Double[] expArr = new Double [exprs.size()];
			expArr = (Double[])exprs.toArray(expArr);
			Arrays.sort(expArr);
			
			double median = expArr[(int)(expArr.length/2.0)];
			
			for (int j=0; j<g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j]))
					g.exprLvls[j] /= median;
			}
		}
	}

	public void writeNewPclToStdOut() {
		System.out.print("YORF\tNAME\tGWEIGHT");
		for (int i = 0; i < exprNames.size(); i++) {
			System.out.print("\t" + (String) exprNames.get(i));
		}
		System.out.println();
		
		System.out.print("EWEIGHT\t\t");
		for (int i = 0; i < numExprs; i++) {
			System.out.print("\t1");
		}
		System.out.println();

		for (int i = 0; i < genes.size(); i++) {
			Gene g = (Gene) genes.get(i);
			System.out.print(g.orf + "\t" + g.name + "\t1");
			DecimalFormat threePlaces = new DecimalFormat("0.000");
			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Double.isNaN(g.exprLvls[j])) {
					System.out.print("\t"+threePlaces.format(g.exprLvls[j]));
					//System.out.print("\t" + g.exprLvls[j]);
				}
				else {
					System.out.print("\t");
				}
			}
			System.out.println();
		}
	}

	public static void usage() {

		System.out.println("Usage: java -jar DivLogNorm.jar <dataset> <avg> <div> <log> <norm>");
		System.out.println("\tdataset - .pcl formatted file");
		System.out.println("\tdiv - 0 (don't divide by gene median, 1 (divide by per gene median)");
		System.out.println("\tlog  - 0 (don't log2 transform), 1 (do log2 xform)");
		System.out.println("\tnorm - 0 (don't normalize), 1 (normalize)");

	}

	public static void main(String[] args) {
		if (args.length < 4) {
			usage();
			System.exit(0);
		}

		DivLogNorm an = new DivLogNorm();

		//Load in the .pcl file
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(args[0]));
			an.readDataFromPCL(br);
		}
		catch (FileNotFoundException e) {
			System.err.println("FileNotFound: " + args[0] + "\nStopping execution.");
			System.exit(0);
		}
		catch (IOException e) {
			System.err.println("Problem reading file: " + args[0]
					+ "\nStopping execution.");
			System.exit(0);
		}

		if (args[1].equalsIgnoreCase("1")) {
			an.divideAllGenesByIndividualMedians();
		}
		if (args[2].equalsIgnoreCase("1")) {
			an.logTransformAll();
		}
		if (args[3].equalsIgnoreCase("1")) {
			an.normalizeAllGenes();
		}

		an.writeNewPclToStdOut();

	}
}
