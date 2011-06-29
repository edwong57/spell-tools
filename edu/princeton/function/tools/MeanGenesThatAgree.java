package edu.princeton.function.tools;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * @author mhibbs
 */
public class MeanGenesThatAgree {
	// this is used for probes as well as combined genes
	public class Gene {
		public String probe;
		public String name;
		public Integer numProbes;
		public float[] exprLvls;

		public Gene(String _probe, String _name) {
			probe = _probe;
			name = _name;
			numProbes = 1;
		};

		public String toString() {
			return probe + "__" + name;
		}
	}

	public HashMap<String, HashSet<Gene>> geneSets;
	public HashMap<String, Gene> genes;
	public HashMap<String, Gene> probe2gene;
	public ArrayList<String> exprNames;
	int numExprs;
	public Float fCutoff;

	public MeanGenesThatAgree() {
		geneSets = new HashMap<String, HashSet<Gene>>();
		genes = new HashMap<String, Gene>();
		probe2gene = new HashMap<String, Gene>();
		exprNames = new ArrayList<String>();
		numExprs = 0;
		fCutoff = Float.NaN;
	}

	/**
	 * Clears out all data, and reads in new data in a pcl file from the
	 * buffered reader provided
	 */
	public boolean readDataFromPCL(BufferedReader fin) throws IOException {
		// Clear out all arrays and re-initialize
		geneSets.clear();
		numExprs = 0;

		// Try to read in the data from the buffered reader
		try {
			int row = 0;
			int gweightCol = -1;
			int colOffset = 2;

			// While there's more data to read in...
			String line = fin.readLine();
			while ((line != null) && (!line.equals(""))) {
				boolean processRow = true;

				// If this is the first row, this determines number of
				// experiments
				if (row == 0) {
					String[] names = line.split("\t");
					for (int i = 2; i < names.length; i++) {
						// If there is a GWEIGHT column, make a note
						if (names[i].equalsIgnoreCase("GWEIGHT")) {
							gweightCol = i;
							colOffset = 3;
						}
						// Otherwise, record the experiment name
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

				// If this is the second row, it may contain weights, or it may
				// contain data
				else if (row == 1) {
					String[] lineParts = line.split("\t");

					// If the first column is EWEIGHT, then just ignore this row
					if (lineParts[0].equalsIgnoreCase("EWEIGHT")) {
						processRow = false;
					}
				}

				// If the line is up for processing...
				if (processRow) {
					String[] lineParts = line.split("\t");

					if (lineParts.length > 1) {
						// Columns 0 and 1 are the SpotIDs and Names
						HashSet<Gene> geneSet = geneSets.get(lineParts[1]);
						if (geneSet == null) {							
							geneSet = new HashSet<Gene>();
							geneSets.put(lineParts[1], geneSet);
						}
						Gene gene = new Gene(lineParts[0], lineParts[1]);
						gene.exprLvls = new float[numExprs];
						geneSet.add(gene);

						if (probe2gene.containsKey(lineParts[0]) == false) {						
							probe2gene.put(lineParts[0], gene);
						}

						// Columns from 2 to the end are data points (except the
						// gweightCol)
						for (int i = 2; i < lineParts.length; i++) {
							if (i != gweightCol) {
								// Try to convert the string to a number, if you
								// fail, just fill in with NaNs
								try {
									float val = Float.parseFloat(lineParts[i]);
									if (!Float.isNaN(val)) {
										gene.exprLvls[i - colOffset] = val;
									}

								} catch (NumberFormatException except) {
									gene.exprLvls[i - colOffset] = Float.NaN;
								}
							}
						}
						// If there are remaining holes, they should be NaNs
						for (int i = (lineParts.length - colOffset); i < numExprs; i++) {
							gene.exprLvls[i] = Float.NaN;
						}

					}
				}

				// Increment row number and repeat
				row++;
				line = fin.readLine();
			}
		} catch (IOException except) {
			throw (except);
		}

		// Finally, return success
		return true;
	}

	public float normalizedDistance(float[] g1, float[] g2) {
		// fisher transform correlation

		float dMX, dMY, dRet, dDX, dDY, dX, dY;
		int i;
		int iN = g1.length;
		dMX = dMY = 0;
		for (i = 0; i < iN; ++i) {
			dMX += g1[i];
			dMY += g2[i];
		}
		dMX /= iN;
		dMY /= iN;

		dRet = dDX = dDY = 0;
		for (i = 0; i < iN; ++i) {
			dX = g1[i] - dMX;
			dY = g2[i] - dMY;
			dRet += dX * dY;
			dDX += dX * dX;
			dDY += dY * dY;
		}
		if (dDX != 0)
			dRet /= Math.sqrt(dDX);
		if (dDY != 0)
			dRet /= Math.sqrt(dDY);
		dRet = (float) Math.log((1 + dRet) / (1 - dRet));
		return dRet;

		// return 0;
		/*
		 * float result = 0; int count = 0; for(int i=0; i<g1.length; i++) { if
		 * (!Float.isNaN(g1[i]) && !Float.isNaN(g2[i])) { result += (g1[i] -
		 * g2[i]) * (g1[i] - g2[i]); count++; } } return (float)Math.sqrt(result
		 * / (float)count);
		 */
	}

	public void makeMeanGeneIfAgree(String key) {
		HashSet<Gene> gSet = geneSets.get(key);
		Gene[] gArray = new Gene[gSet.size()];
		gArray = (Gene[]) gSet.toArray(gArray);
		if (gArray.length == 1) {
			genes.put(key, gArray[0]);
		} else {
			int[] numAgree = new int[gArray.length];
			for (int i = 0; i < gArray.length; i++) {
				for (int j = i + 1; j < gArray.length; j++) {
					if (genesMLagree(gArray[i], gArray[j])) {
						numAgree[i] += 1;
						numAgree[j] += 1;
						// System.err.println("AGREED "+gArray[i].toString()+"\t"+gArray[j].toString());
					} else {
						// System.err.println("Did not agree "+gArray[i].toString()+"\t"+gArray[j].toString());
					}
				}
			}
			Gene meanGene = new Gene(key, key);
			meanGene.exprLvls = new float[numExprs];
			int[] counts = new int[numExprs];
			for (int c = 0; c < numExprs; c++)
				counts[c] = 0;

			int numGood = 0;
			// Include a gene if it agrees with at least half of the others
			for (int i = 0; i < gArray.length; i++) {
				if (numAgree[i] >= ((numAgree.length - 1) / 2.0)) {
					numGood += 1;
					for (int c = 0; c < numExprs; c++) {
						if (!Float.isNaN(gArray[i].exprLvls[c])) {
							meanGene.exprLvls[c] += gArray[i].exprLvls[c];
							counts[c]++;
						}
					}
				}
			}
			meanGene.numProbes = numGood;
			// If at least half of the genes are included, mean them
			if (numGood > (numAgree.length / 2.0)) {
				for (int c = 0; c < numExprs; c++) {
					if (counts[c] > 0)
						meanGene.exprLvls[c] /= (float) counts[c];
					else
						meanGene.exprLvls[c] = Float.NaN;
				}
				genes.put(key, meanGene);
				System.err.println(key + " expression AVERAGED from " + numGood
						+ " probes (" + (gArray.length - numGood)
						+ " probes were excluded)");
			}
			// Otherwise, toss out this gene
			else {
				System.err.println(key + " did not agree, and was discarded");
			}
		}
	}

	public void meanAllThatAgree() {
		Set<String> keySet = geneSets.keySet();
		String[] keys = new String[keySet.size()];
		keys = (String[]) keySet.toArray(keys);
		for (int i = 0; i < keys.length; i++) {
			makeMeanGeneIfAgree(keys[i]);
		}
	}

	public float[] sampleDists(float rate) {
		Set<String> keys = probe2gene.keySet();
		String[] probes = new String[keys.size()];
		probes = (String[]) keys.toArray(probes);

		ArrayList<Float> dists = new ArrayList<Float>();
		for (int i = 0; i < probes.length; i++) {
			for (int j = i + 1; j < probes.length; j++) {
				float r = (float) Math.random();
				if (r < rate) {
					Gene g1 = (Gene) probe2gene.get(probes[i]);
					Gene g2 = (Gene) probe2gene.get(probes[j]);
					if (!g1.name.equals(g2.name)) {
						dists.add(new Float(normalizedDistance(g1.exprLvls,
								g2.exprLvls)));
					}
				}
			}
		}

		float[] result = new float[dists.size()];
		for (int i = 0; i < dists.size(); i++) {
			result[i] = ((Float) dists.get(i)).floatValue();
		}
		return result;
	}

	public float[] mappedDists() {
		Set<String> keys = geneSets.keySet();
		String[] names = new String[keys.size()];
		names = (String[]) keys.toArray(names);

		ArrayList<Float> dists = new ArrayList<Float>();
		for (int i = 0; i < names.length; i++) {
			HashSet<Gene> gSet = geneSets.get(names[i]);
			if (gSet.size() > 1) {
				Gene[] gArray = new Gene[gSet.size()];
				gArray = (Gene[]) gSet.toArray(gArray);
				for (int x = 0; x < gArray.length; x++) {
					for (int y = x + 1; y < gArray.length; y++) {
						dists.add(new Float(normalizedDistance(
								gArray[x].exprLvls, gArray[y].exprLvls)));
					}
				}
			}
		}

		float[] result = new float[dists.size()];
		for (int i = 0; i < dists.size(); i++) {
			result[i] = ((Float) dists.get(i)).floatValue();
		}
		return result;
	}

	public float[] bg_pdf;
	public float[] map_pdf;
	public float bg_min;
	public float bg_max;
	public float bin_size;
	public int num_bins;

	public static float Mean(float[] fArray) {
		float fMean = 0;
		for (int i = 0; i < fArray.length; i++)
			fMean += fArray[i];
		fMean /= fArray.length;
		return fMean;
	}

	/**
	 *
	 * @param distributionFilename: filename where distributions are written. Can be null.
	 *  If null, no output is written.
	 */
	public void createPDFs(String distributionFilename) {
		float[] bg = sampleDists(0.02f);
		float[] map = mappedDists();
		float bgMean = Mean(bg);
		float mapMean = Mean(map);
		fCutoff = bgMean + (float) 0.5 * (mapMean - bgMean);
		System.err.println("Mapped correlation mean = " + mapMean);
		System.err.println("Background correlation mean = " + bgMean);
		System.err.println("Cutoff = " + fCutoff);
		
		if (distributionFilename != null) {
			FileOutputStream outS;
			PrintStream pS;
			try {
				outS = new FileOutputStream(distributionFilename);
				pS = new PrintStream(outS);
				int i;
				pS.print("Background");
				for (i = 0; i < bg.length; i++)
					pS.print("\t" + bg[i]);
				pS.print('\n');

				pS.print("Mapped");
				for (i = 0; i < map.length; i++)
					pS.print("\t" + map[i]);
				pS.print('\n');

			} catch (Exception e) {
				System.err.println("Error writing to file");

			}
		}
		
		fCutoff = bgMean + (float) 0.5 * (mapMean - bgMean);
		// old code for binning
		/*
		 * bg_min = Float.POSITIVE_INFINITY; bg_max = Float.NEGATIVE_INFINITY;
		 * for (int i = 0; i < bg.length; i++) { if (bg[i] < bg_min) bg_min =
		 * bg[i]; if (bg[i] > bg_max) bg_max = bg[i]; } num_bins = Math.max(30,
		 * (int) ((float) bg.length / 300)); bin_size = (float) (bg_max -
		 * bg_min) / (float) num_bins;
		 * 
		 * bg_pdf = new float[num_bins]; for (int i = 0; i < bg.length; i++) {
		 * bg_pdf[bin(bg[i])] += 1; }
		 * 
		 * 
		 * map_pdf = new float[num_bins]; for (int i = 0; i < map.length; i++) {
		 * map_pdf[bin(map[i])] += 1; }
		 * 
		 * for (int i = 0; i < num_bins; i++) { bg_pdf[i] /= bg.length;
		 * map_pdf[i] /= map.length; }
		 */
	}

	public int bin(float value) {
		int result = (int) ((float) (value - bg_min) / bin_size);
		if (result < 0)
			return 0;
		if (result >= num_bins)
			return num_bins - 1;
		return result;
	}

	public boolean genesMLagree(Gene g1, Gene g2) {
		float dist = normalizedDistance(g1.exprLvls, g2.exprLvls);
		// if (map_pdf[bin(dist)] >= bg_pdf[bin(dist)])
		if (dist > fCutoff)
			return true;
		else
			return false;
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

		Set<String> keySet = genes.keySet();
		String[] keys = new String[keySet.size()];
		keys = (String[]) keySet.toArray(keys);

		for (int i = 0; i < keys.length; i++) {
			Gene g = (Gene) genes.get(keys[i]);
			System.out.print(g.probe + "\t" + g.name + "\t" + g.numProbes);
			DecimalFormat threePlaces = new DecimalFormat("0.000000");
			for (int j = 0; j < g.exprLvls.length; j++) {
				if (!Float.isNaN(g.exprLvls[j]))
					System.out.print("\t" + threePlaces.format(g.exprLvls[j]));
				else
					System.out.print("\t");
			}
			System.out.println();
		}
	}

	public void writeNewPclToFile(String OutFile) {
		FileOutputStream outS;
		PrintStream pS;
		try {
			outS = new FileOutputStream(OutFile);
			pS = new PrintStream(outS);

			pS.print("YORF\tNAME\tGWEIGHT");
			for (int i = 0; i < exprNames.size(); i++) {
				pS.print("\t" + (String) exprNames.get(i));
			}
			pS.println();

			pS.print("EWEIGHT\t\t");
			for (int i = 0; i < numExprs; i++) {
				pS.print("\t1");
			}
			pS.println();

			Set<String> keySet = genes.keySet();
			String[] keys = new String[keySet.size()];
			keys = (String[]) keySet.toArray(keys);

			for (int i = 0; i < keys.length; i++) {
				Gene g = (Gene) genes.get(keys[i]);
				pS.print(g.probe + "\t" + g.name + "\t" + g.numProbes);
				DecimalFormat threePlaces = new DecimalFormat("0.000000");
				for (int j = 0; j < g.exprLvls.length; j++) {
					if (!Float.isNaN(g.exprLvls[j]))
						pS.print("\t" + threePlaces.format(g.exprLvls[j]));
					else
						pS.print("\t");
				}
				pS.println();
			}
		} catch (Exception e) {
			System.err.println("Error writing to file");
		}

	}

	public static void usage() {
		System.out.println("Usage: java -jar MeanGenesThatAgree.jar <input-file> <output-file> <optional: dist file>");
		System.out.println("\tinput-file  - pcl formatted file");
		System.out.println("\toutput-file - pcl formatted averaged file");
		System.out.println("\tdist-file   - optional file where distances are written");
	}

	/**
	 * @param args:
	 *   [0]: input file
	 *   [1]: Output file
	 *   [2]: Distance file (optional)
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			usage();
			System.exit(0);
		}

		MeanGenesThatAgree mm = new MeanGenesThatAgree();

		// Load in the .pcl file
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(args[0]));
			mm.readDataFromPCL(br);
		} catch (FileNotFoundException e) {
			System.err.println("FileNotFound: " + args[0]
					+ "\nStopping execution.");
			System.exit(0);
		} catch (IOException e) {
			System.err.println("Problem reading file: " + args[0]
					+ "\nStopping execution.");
			System.exit(0);
		}
		
		String tmpOutFile = null;
		if (args.length > 2) {
			tmpOutFile = args[2];			
		}
		mm.createPDFs(tmpOutFile);
				
		mm.meanAllThatAgree();		
		String OutFile = args[1];
		mm.writeNewPclToFile(OutFile);
		// mm.writeNewPclToStdOut();
	}
}
