import proteinGrouper.DTASelectProteinGrouper;
import psmPostProcess.MissingPSMRetriever;

import java.io.*;
import java.util.*;

//DTA file selector
//David L. Tabb
//August 10, 2000
//modified April 4, 2017, Titus Jung titusjung@gmail.com
/*
 * Summary of program flow:
 * DTASelect is run in the directory containing .dat and sequest.params
 * files.  It anticipates that .dta files and .out files will be found
 * in subdirectories of the current directory.  The first time
 * DTASelect is run in any particular directory, it visits each .out
 * file in each subdirectory, cataloging the key measures from each
 * .out file in a tab-delimited flat file database entitled
 * DTASelect.txt.  Information from the FASTA database indicated in
 * sequest.params is integrated into the file as well, providing
 * descriptive information about each locus.

 * Once all .out files have been cataloged, the selection process can
 * begin.  The software first applies filters to individual .dta
 * files, removing those that fail to meet criteria of XCorr, DeltCN,
 * or other measures.  Once this process has been undertaken, the
 * loci are filtered; duplicate sequences for each locus can be
 * eliminated, and loci with too few representative peptides are
 * removed.

 * Once the filtering is complete, an HTML document summarizing the
 * loci and supporting peptides is produced.  DTASelect.html links to
 * a CGI for spectrum display, a CGI for .out file display, and a
 * BLAST search CGI.
 */

public class DTASelect {
	DataSet Data = new DataSet();
	IniFile Configuration = new IniFile();
	String CommandLineOptions;

	public static void usage() {
		SelectCriteria Cutoffs = new SelectCriteria();
		System.out.println(Protein.Version()
				+ ", a program for selecting specific SEQUEST results");
		System.out
				.println("Copyright 2001, Yates Lab, The Scripps Research Institute");
		System.out.println("Run in directory containing .sqt or .out files.");
		System.out.println("Creates the file DTASelect.html");
		System.out.println();
		System.out.println("INDIVIDUAL SPECTRUM FILTERS:");
		System.out.println("--pr #\tSet lowest peptide probability");
		System.out.println("--fp #\tSet peptide false discovery rate");
		// added by Howard Choi
		System.out.println("--sfp #\tSet peptide false discovety rate");
		System.out.println("--pfp #\tSet protein false discovery rate");

		System.out.println("-m 0\tRequire peptides to be modified");
		System.out.println("-m 1\tInclude peptides regardless of modification");
		System.out.println("-m 2\tExclude modified peptides");
		System.out.println("--nm #\tSet minimum number of modification sites");
		System.out.println("--NM #\tSet maximum number of modification sites");
		System.out
				.println("-y 0\tInclude peptides regardless of cleavage status");
		System.out
				.println("-y 1\tInclude only half- or fully cleaved peptides");
		System.out.println("-y 2\tInclude only fully cleaved peptides");
		System.out
				.println("--mc #\tSet minimum number of missed tryptic sites");
		System.out
				.println("--MC #\tSet maximum number of missed tryptic sites");
		System.out
				.println("-a true | false\tAllow peptides with ambiguous IDs");
		System.out.println("-v -1\tKeep N peptides, discard all others");
		System.out.println("-v 0\tIgnore manual validation info");
		System.out.println("-v 1\tKeep Y peptides, discard N peptides");
		System.out.println("-v 2\tKeep Y and M peptides, discard N peptides");
		System.out
				.println("-v 3\tKeep Y and M peptides, discard N and U peptides");
		System.out.println();
		System.out.println("SPECIALIZED FILTERS (off by default):");
		System.out.println("-1 #\tSet lowest +1 XCorr");
		System.out.println("-2 #\tSet lowest +2 XCorr");
		System.out.println("-3 #\tSet lowest +3 XCorr");
		System.out.println("-4 #\tSet lowest +4 XCorr");
		System.out.println("-X1 #\tSet highest +1 XCorr");
		System.out.println("-X2 #\tSet highest +2 XCorr");
		System.out.println("-X3 #\tSet highest +3 XCorr");
		System.out.println("-X4 #\tSet highest +4 XCorr");
		System.out.println("-d #\tSet lowest DeltCN");
		System.out.println("-c #\tSet lowest charge state");
		System.out.println("-C #\tSet highest charge state");
		System.out.println("--mz #\tSet lowest peptide m/z");
		System.out.println("--MZ #\tSet highest peptide m/z");
		System.out.println("-dm #\tSet lowest PPM delta mass");
		System.out.println("-DM #\tSet highest PPM delta mass");
		System.out
				.println("-i #\tSet lowest proportion of fragment ions observed");
		System.out.println("-s #\tSet maximum Sp ranking");
		System.out.println("-S #\tSet minimum Sp score");
		System.out
				.println("-Sic $\tSequences must contain all of these residues");
		System.out.println("-Sip $\tSequences must contain this pattern"); // maybe
																			// change
																			// to
																			// something
																			// more
																			// user
																			// friendly
																			// "any of these residues"
		System.out
				.println("-Sec $\tSequences must not contain any of these residues");
		System.out
				.println("-Si2 $\tSequences must contain exactly two of these residues");
		System.out
				.println("-Si1 $\tSequences must contain exactly one of these residues");
		System.out
				.println("-SiN $\tSequences must contain at least one of these residues");
		System.out.println("-Stn $\tPreceding residue must be one of these");
		System.out.println("-Stc $\tC terminal residue must be one of these");
		System.out
				.println("-Smn #\tSequence must be at least this long (default = 6)");
		System.out
				.println("-Smx #\tSequence must be no longer than this (default = 100)");
		System.out
				.println("-Ser #\tSequence must have this many complete ends");
		System.out
				.println("-Ngly 0\tFilter out N-linked glycosylated peptides");
		System.out.println("-Ngly 1\tKeep only N-linked glycosylated peptides");
		System.out
				.println("-file $\tOnly show spectra from MS2 files matching this string");
		System.out.println();
		System.out.println("LOCUS FILTERS:");
		System.out.println("-t 0\tShow all spectra for each sequence");
		System.out
				.println("-t 1\tShow only one spectrum per salt step for each sequence");
		System.out.println("-t 2\tShow only one spectrum for each sequence");
		System.out.println("--CL #\tSet lowest protein confidence level");
		System.out.println("--FP #\tSet protein false discovery rate");
		System.out.println("-V -1\tKeep N proteins, discard all others");
		System.out.println("-V 0\tIgnore manual validation info");
		System.out.println("-V 1\tKeep Y proteins, discard N proteins");
		System.out.println("-V 2\tKeep Y and M proteins, discard N proteins");
		System.out
				.println("-V 3\tKeep Y and M proteins, discard N and U proteins");
		System.out
				.println("-u\tInclude only loci with uniquely matching peptides");
		System.out.println("-in\tInclude proteins that are subsets of others");
		System.out
				.println("-e $\tRemove proteins with IDs matching this string");
		System.out.println("--Seq #\tSet lowest protein sequence coverage");
		System.out.println("--mw #\tSet lowest protein molecular weight");
		System.out.println("--MW #\tSet highest protein molecular weight");
		System.out.println("-U \t Print only unique peptides");
		System.out
				.println("-StN $\tInclude only proteins that start with one of these residues");
		System.out
				.println("-StC $\tInclude only proteins that end with one of these residues");
		System.out
				.println("-E $\tInclude only proteins with IDs matching this string");
		System.out
				.println("-l $\tRemove proteins with descriptions including this word");
		System.out
				.println("-L $\tInclude only proteins with descriptions including this word");
		System.out
				.println("-M #\tSet minimum modified peptides per locus criterion");
		System.out
				.println("-tp #\tSet minimum fully cleaved peptides per locus criterion");
		System.out
				.println("-ty #\tSet minimum cleavage status per locus criterion");
		System.out.println("-tS #\tSet minimum Sp Score per locus criterion");
		System.out.println("-tX #\tSet minimum XCorr per locus criterion");
		System.out.println("-tfp #\tSet maximum FP per locus criterion");
		System.out
				.println("-tpr #\tSet minimum probability per locus criterion");
		System.out
				.println("-tDM #\tSet maximum PPM delta mass per locus criterion");
		System.out.println("-p #\tSet minimum peptides per locus criterion");
		System.out
				.println("-N #\tSet minimum spectral count per locus criterion");
		System.out.println();
		System.out.println("OPTIONS:");
		System.out
				.println("--short\t\tUse actual (short form) accession numbers");
		System.out
				.println("--dcn 0\t\tDefine DeltCN with respect to the second-best peptide");
		System.out
				.println("--dcn 1\t\tDefine DeltCN with respect to the next best DIFFERENT peptide");
		System.out
				.println("--iso\t\tDefine delta mass with respect to nearest isotope");
		System.out
				.println("--noiso\t\tDefine delta mass with respect to monoisotopic/average mass");
		System.out.println("--enzyme 1\tTryptic digest (after K or R)");
		System.out.println("--enzyme 2\tLysC digest (after K)");
		System.out.println("--enzyme 3\tAspN digest (before D)");
		System.out.println("--enzyme 99\tIntact proteins (no digest)");
		System.out
				.println("--before $\tEnzyme cuts before one of these residues");
		System.out
				.println("--after $\tEnzyme cuts after one of these residues");
		System.out.println();
		System.out.println("STATISTICS:");
		System.out
				.println("--decoy $\tDatabase protein decoy label (default: \"Reverse_\")");
		System.out.println("--hidedecoy\tHide decoy database hits");
		System.out
				.println("--trueprot $\tTrue proteins have IDs matching this string");
		System.out.println("--nostats\tNo statistics");
		// System.out.println("--noemp\t\tNo empirical coefficients used");
		System.out
				.println("--cmax #\tMaximum charge state for statistics options");
		System.out
				.println("--raw\t\tUse raw primary score (XCorr) values, no normalization");
		System.out
				.println("--nochgstat\tDon't use models for each charge state");
		System.out
				.println("--trypstat\tUse cleavage status when calculating probabilities");
		System.out
				.println("--modstat\tUse separate models for modified peptides");
		System.out.println("--noxc\t\tXCorr not used for statistics");
		System.out.println("--nodcn\t\tDeltCN not used for statistics");
		System.out.println("--sp\t\tSp used for statistics");
		System.out.println("--mass\t\tDelta mass used for statistics");
		System.out.println("--logmass\tLog delta mass used for statistics");
		System.out.println("--logspr\tLog SpRank used for statistics");
		System.out.println("--ionstat\tIon proportion used for statistics");
		// System.out.println("--fit 0\t\tHistogram-based statistics");
		// System.out.println("--fit 1\t\tMaximum likelihood-based statistics");
		// System.out.println("--fit 2\t\tKernel function-based statistics");
		System.out
				.println("--fptol #\tInitial FP rate used in determining linear coefficients");
		System.out
				.println("--plot 0\tDo not plot score distributions (default)");
		System.out.println("--plot 1\tPlot 2D peptide score distributions");
		System.out.println("--plot 2\tAlso plot ROC curves for each parameter");
		System.out
				.println("--plot 3\tAlso plot normal and decoy score histograms");
		System.out.println();
		System.out.println("UTILITIES:");
		System.out.println("--omssa\t\tUse OMSSA scoring parameters");
		System.out.println("--dm\t\tDisplay delta mass information");
		System.out.println("--pI\t\tDisplay peptide pI information");
		System.out.println("--KD\t\tDisplay Kyte Doolittle hydropathy score");
		System.out.println("--BB\t\tDisplay Bull Breese hydrophobicity index");
		System.out
				.println("--HPLC\t\tDisplay HPLC index adjusted for pH3.4 conditions");
		System.out
				.println("--seqpos\tDisplay sequence position in the DTASelect-filter.txt file");
		System.out
				.println("--extra\t\tOutput extra columns in DTASelect-filter.txt and DB-Peptides.txt files");
		System.out.println("--class\t\tClassify proteins by supplied list");
		System.out.println("--aux\t\tIncorporate protein info from aux file");
		System.out.println("--nofilter, -n\tDo not apply any criteria");
		System.out.println("--CGI\t\tUse new accessory programs");
		System.out
				.println("--GUI\t\tReport through GUI instead of output files");
		System.out.println("--BE\t\tDisplay Bird's Eye View of proteins");
		System.out
				.println("--compress\tCreate .IDX and .SPM files from spectra");
		System.out.println("--copy\t\tPrint commands to copy selected IDs");
		System.out.println("--XML\t\tSave XML report of filtered results");
		System.out.println("--DB\t\tSave in format for database import");
		System.out.println("--chroma\tSave chromatography report");
		System.out.println("--similar\tSave protein similarity table");
		System.out.println("--align\t\tSave sequence alignment report");
		System.out.println("--mods\t\tSave modification report");
		System.out.println("--help, -h\tPrint this help");
		System.out
				.println("--Mascot\tExpect Mascot output rather than SEQUEST");
		System.out.println("--NoBackground\tDo not use a background graphic");
		System.out.println("--here, -.\tInclude only IDs in current directory");
		System.out.println("--quiet\t\tReduced terminal output");
		System.out.println("--brief\t\tReduced complexity of html file output");
		System.out
				.println("--path $\tSpecify location of the DTASelect folder");
		System.out.println("-ug\t\tRuns Protein Grouper:  this take non unique peptides from existing proteins" +
				" and creates protein groups consisting of only unique peptides");
		System.out.println("\t\tPlease check http://manual.integratedproteomics.com/1/en/topic/16-2-unique-group-filter" +
				" for more information ");
		System.out.println("-addPSM\t\tinclude all spectra from unfiltered data for identified peptides");
		System.out.println("Defaults:");
		System.out.println(Cutoffs.PrintCriteria("", "\n", "\t"));
		System.exit(0);
	}

	// Parse command line arguments and create the SpecBrowser object.
	// If an extension is specified, create a FilenameFilter for it.
	// If no directory is specified, use the current directory.
	public static void main(String args[]) throws IOException {
		int i;
		String ParamsFileOptions[];
		boolean OkayToProceed = true;
		DTASelect AppObject = new DTASelect();
		// added by Howard Choi
		boolean run = true;
		boolean isPfp = false;
		boolean isSfp = false;
		boolean iterationStarts = false;
		boolean poApplyed = false;
		float previousOutputFPValue = -1.0f;
		float targetFPValue = -1.0f;
		float holdingFPValue = 0.0f;
		float iterationCount = 0.0f;
		// added by Howard Choi
		// take DTAParam file
		if (AppObject.Data.Cutoffs.ProcessParamsFile() == null) { //need to change directory chooser for debug
		} else {
			ParamsFileOptions = AppObject.Data.Cutoffs.ProcessParamsFile();
			String[] tempArgs = concatArray(ParamsFileOptions, args);
			args = tempArgs;
		}
		// read args[] and define boolean
		// change args[] --fp
		String additionalFP = "";
		for (int u = 0; u < args.length; u++) {
			if (args[u].equals("--fp")) {
				targetFPValue = Float.valueOf(args[u + 1]);
			} else if (args[u].equals("--pfp")) {
				targetFPValue = Float.valueOf(args[u + 1]);
				isPfp = true;
				args[u] = "--fp";
				additionalFP = "--pfp " + args[u + 1];
			} else if (args[u].equals("--sfp")) {
				targetFPValue = Float.valueOf(args[u + 1]);
				isSfp = true;
				args[u] = "--fp";
				additionalFP = "--sfp " + args[u + 1];
			}
		}

		// added by Howard Choi
		// ---------while run is true
		
		int no_improvement_iterations=0;
		
		while (run) {
			AppObject.Configuration.Initialize(); //change

			// for debug //TODO
			AppObject.CommandLineOptions = "";
			System.out.println(Protein.Version());
			// discarded by Howard Choi
			/*
			 * // Set criteria according to DTASelect.params ParamsFileOptions =
			 * AppObject.Data.Cutoffs.ProcessParamsFile(); if (ParamsFileOptions
			 * != null) { for(i = 0; i < ParamsFileOptions.length; i++) {
			 * AppObject.CommandLineOptions += " " + ParamsFileOptions[i]; } }
			 */
			// Read in command-line arguments
			for (i = 0; i < args.length; i++) {
				AppObject.CommandLineOptions += " " + args[i];
				if ((args[i].equals("-h")) || (args[i].equals("/?"))
						|| (args[i].equals("--help"))) {
					usage();
				}
			}
			// Set criteria according to command-line arguments
			AppObject.Data.Cutoffs.SetCriteria(args); //TODO to change for debug
			// Report criteria to screen
			// System.out.println(AppObject.Data.Cutoffs.PrintCriteria("", "\n",
			// "\t"));
			// Proceed through the full run only if we can write to this
			// directory or we know we're using the GUI.
			// ------if
			if (AppObject.Data.Cutoffs.UseGUIInstead || Protein.PrintTest()) {
				// Read .out or .sqt files into memory
				AppObject.ReadAndPreprocessProteins(); //TODO to change for debug

				// Apply the user-selected criteria
				System.out.println("Applying criteria to spectra and loci...");
				if(AppObject.Data.Cutoffs.dms)
				{
					AppObject.Data.ApplyCriteriaDMS();
				}
				else
				{
					AppObject.Data.ApplyCriteria();
				}
				// Compress spectra to IDX and SPM files
				if (AppObject.Data.Cutoffs.CompressDTAs)
					AppObject.CompressSpectra();
				AppObject.Data.SetColumnHeadings();
				// Use GUI?
				if (AppObject.Data.Cutoffs.UseGUIInstead) {
					AppObject.Data.LocusList.SortByClass();
					AppObject.Data.LocusList.DisplaySortDTALists();
					AppObject.PresentGUI();
				} else {
					System.out.println("Creating selected reports...");
					// Report results
					AppObject.Data.PrintReports(AppObject.CommandLineOptions,
							AppObject.Configuration, "DTASelect", additionalFP);
				}
				System.out.println("DTASelect is completed.");

				// added by Howard Choi
				// no --fp paramater
				// ---Howard's
				if (!(hasFPParam(args))) {
					run = false;
				}
				// has --fp parameter
				else {
					float currentInputFPValue = 0.0f;
					float currentOutputFPValue = 0.0f;
					float newInputFPValue = 0.0f;
					// --pfp
					if (isPfp) {
						currentInputFPValue = Float.valueOf(AppObject.Data
								.getProtFPInput());
						currentOutputFPValue = Float.valueOf(AppObject.Data
								.getProtFPOutput());
						run = !(isTargeted(currentOutputFPValue, targetFPValue,
								0.001f));
						// --sfp
					} else if (isSfp) {
						currentInputFPValue = Float.valueOf(AppObject.Data
								.getPepFPInput());
						currentOutputFPValue = Float.valueOf(AppObject.Data
								.getPepFPOutput());
						run = !(isTargeted(currentOutputFPValue, targetFPValue,
								0.001f)); // 0.0005f
						// --fp
					} else {
						// so it
						run = false;
						// currentOutputFPValue = -1.0f;
					}
					System.out.println("previous     " + previousOutputFPValue);
					System.out.println("current    " + currentOutputFPValue);
					
					if (previousOutputFPValue == currentOutputFPValue) {
						no_improvement_iterations++;
					} else {
						no_improvement_iterations=0;
					}
					
					if (no_improvement_iterations>10) {
						run = false;
					}
					
					// run = !(isTargeted(previousOutputFPValue,
					// currentOutputFPValue, targetFPValue, 0.002f,
					// currentOutputFPValue<targetFPValue, iterationStarts));
					previousOutputFPValue = currentOutputFPValue;
					if (currentOutputFPValue <= targetFPValue) {
						holdingFPValue = currentInputFPValue;
					}
					// output is not close to target
					if (run) {
						// need up
						if (currentOutputFPValue < targetFPValue) {
							// no down yet
							if (iterationCount == 0.0f) {
								newInputFPValue = currentInputFPValue + 0.1f;
								poApplyed = true;
							} else {
								iterationCount += 1.0f;
								newInputFPValue = binarySearch(
										currentInputFPValue, targetFPValue,
										true, poApplyed, iterationCount);
							}
						}
						// need down
						else {
							iterationCount += 1.0f;
							iterationStarts = true;
							newInputFPValue = binarySearch(currentInputFPValue,
									targetFPValue, false, poApplyed,
									iterationCount);
						}
						System.out.println("new Input Value    "
								+ newInputFPValue);
					}
					// terminate if currentFPOutputValue is NaN
					boolean nanCheck = Float.isNaN(currentOutputFPValue);
					if (nanCheck) {
						run = false;
					}
					// terminate at --fp 1.0
					if (newInputFPValue >= 1.0f || iterationCount > 10) {
						// terminate
						if (currentOutputFPValue <= targetFPValue) {
							run = false;
						}
						// make new arguments with holding fp value
						else {
							for (int iu = 0; iu < args.length; iu++) {
								if (args[iu].equals("--fp")) {
									args[iu + 1] = String
											.valueOf(holdingFPValue);
								}
							}
						}
					}
					// make new argument
					else {
						for (int y = 0; y < args.length; y++) {
							if (args[y].equals("--fp")) {
								args[y + 1] = String.valueOf(newInputFPValue);
							}
						}
					}
				}// ---Howard's end
			}// ------if end
			else {
				System.out
						.println("I can't write DTASelect.html in this directory.");
				run = false;
			}
		}// ---------while end'
		if(AppObject.Data.Cutoffs.UseProteinGrouper)
		{
			DTASelectProteinGrouper.process( "DTASelect-filter.txt", "DTASelect-filter.txt.clustered");
			File oldDTASelect = new File("DTASelect-filter.txt");
			oldDTASelect.delete();
			File clusteredDTASelect = new File("DTASelect-filter.txt.clustered");
			clusteredDTASelect.renameTo(new File("DTASelect-filter.txt"));
		}
		if(AppObject.Data.Cutoffs.addPSM)
		{
			MissingPSMRetriever retriever = new MissingPSMRetriever();
			retriever.setIsotopes(AppObject.Data.Cutoffs.Isotopes);
			retriever.collectMissedPSM("DTASelect-filter.txt",System.getProperty("user.dir"),
					"DTASelect-filter.txt.newPSM",!AppObject.Data.Cutoffs.IncludeOnlyUniques);
			File oldDTASelect = new File("DTASelect-filter.txt");
			oldDTASelect.delete();
			File clusteredDTASelect = new File("DTASelect-filter.txt.newPSM");
			clusteredDTASelect.renameTo(new File("DTASelect-filter.txt"));
		}
	}

	private void CompressSpectra() {
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		EndsWithFilter FilenameFilter = new EndsWithFilter(".dta");
		String[] DirectoryListing;
		int FileIterator;
		File CurrentSubDir;
		DirectoryFilter SubDirFilter = new DirectoryFilter();
		String[] SubDirList = CurrentDirectory.list(SubDirFilter);
		int DirIterator;
		IDXFile Index = new IDXFile();
		SPMFile Hog = new SPMFile();
		RandomAccessFile RAHandle;
		StringTokenizer Parser;
		long FileOffset;
		short LastLowScan = 0;
		short LastHighScan = 0;
		short CurrentLowScan;
		short CurrentHighScan;
		BubbleSort(SubDirList);
		try {
			// Remove any existing SPM and IDX files
			CurrentSubDir = new File(CurrentDirectory, "DTASelect.SPM");
			if (CurrentSubDir.exists()) {
				CurrentSubDir.delete();
			}
			CurrentSubDir = new File(CurrentDirectory, "DTASelect.IDX");
			if (CurrentSubDir.exists()) {
				CurrentSubDir.delete();
			}
			RAHandle = new RandomAccessFile(new File(CurrentDirectory,
					"DTASelect.SPM"), "rw");
			for (DirIterator = 0; DirIterator < SubDirList.length; DirIterator++) {
				CurrentSubDir = new File(CurrentDirectory,
						SubDirList[DirIterator]);
				DirectoryListing = CurrentSubDir.list(FilenameFilter);
				BubbleSort(DirectoryListing);
				Index.AddSubdir(SubDirList[DirIterator]);
				// Let users know that we're starting to read files in this
				// subdirectory
				System.out.println("Now reading spectra in "
						+ SubDirList[DirIterator]);
				for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
					Parser = new StringTokenizer(
							DirectoryListing[FileIterator], ".");
					Parser.nextToken();
					CurrentLowScan = new Short(Parser.nextToken()).shortValue();
					CurrentHighScan = new Short(Parser.nextToken())
							.shortValue();
					if ((CurrentLowScan != LastLowScan)
							|| (CurrentHighScan != LastHighScan)) {
						FileOffset = Hog.AddSpectrumToQueue(CurrentSubDir,
								DirectoryListing[FileIterator]);
						Index.AddDTAToIndex(CurrentLowScan, CurrentHighScan,
								Hog.LastAddedMOverZ, FileOffset);
						LastLowScan = CurrentLowScan;
						LastHighScan = CurrentHighScan;
					}
				}
				System.out.println("Compressing spectra for "
						+ SubDirList[DirIterator]);
				Hog.AddSpectraToFile(RAHandle);
			}
			System.out.println("Writing spectrum index to disk...");
			Index.WriteToDisk();
		} catch (IOException failure) {
			System.out.println("Error while compressing .dta files.");
			System.out.println(failure.getMessage());
			System.out.println(failure.getLocalizedMessage());
			System.out.println(failure.toString());
			System.exit(0);
		}
	}

	// added by Howard Choi
	// has fp parameter?
	public static boolean hasFPParam(String[] arguments) {
		boolean has = false;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals("--fp") || arguments[i].equals("--pfp")
					|| arguments[i].equals("--sfp")) {
				has = true;
			}
		}
		return has;
	}

	// added by Howard Choi
	public static boolean hasPSFPParam(String[] arguments) {
		boolean hasFP = false;
		for (int i = 0; i < arguments.length; i++) {
			if (arguments[i].equals("--pfp") || arguments[i].equals("--sfp")) {
				hasFP = true;
			}
		}
		return hasFP;
	}

	/*
	 * //added by Howard Choi public static boolean isTargeted(float
	 * currentOutputPFP, float targetPFP, float gap){ boolean targeted = false;
	 * if(Math.abs(currentOutputPFP-targetPFP) <= gap){ targeted = true; }
	 * return targeted; }
	 */

	// added by Howard Choi
	public static boolean isTargeted(float currentOutputPFP, float targetPFP,
			float gap) {
		boolean targeted = false;
		if ((targetPFP - currentOutputPFP <= gap)
				&& (targetPFP >= currentOutputPFP)) {
			targeted = true;
		}
		return targeted;
	}

	// added by Howard Choi
	public static float binarySearch(float currentInputVal, float targetVal,
			boolean isUp, boolean pointOneApplyed, float iterationNumber) {
		float searchedVal = 0.0f;
		float movingVal = 0.0f;
		if (!pointOneApplyed) {
			movingVal = targetVal * (float) (Math.pow(0.5f, iterationNumber));
		} else {
			movingVal = 0.1f * (float) (Math.pow(0.5f, iterationNumber));
		}
		if (isUp) {
			searchedVal = currentInputVal + movingVal;
		} else {
			searchedVal = currentInputVal - movingVal;
		}
		return searchedVal;
	}

	// added by Howard Choi
	public static float RoundTo(float Value, int Places) {
		// Converts a value to a rounded value
		if (Places == 0) {
			return new Integer(Math.round(Value)).intValue();
		} else {
			double Multiplier = Math.pow(10, Places);
			return new Double(Math.rint(Value * Multiplier) / Multiplier)
					.floatValue();
		}
	}

	// added by Howard Choi
	public static <String> String[] concatArray(String[] A, String[] B) {
		String[] C = Arrays.copyOf(A, A.length + B.length);
		System.arraycopy(B, 0, C, A.length, B.length);
		return C;
	}

	public static void BubbleSort(String[] list) {
		int pass;
		int i;
		int ListLength = list.length;
		int Limit;
		String hold;
		for (pass = 1; pass < ListLength; pass++) {
			Limit = ListLength - pass;
			for (i = 0; i < Limit; i++) {
				if (list[i].compareTo(list[i + 1]) > 0) {
					hold = list[i];
					list[i] = list[i + 1];
					list[i + 1] = hold;
				}
			}
		}
	}

	private void PresentGUI() {
		DTASelectGUI GUI = new DTASelectGUI(Data);
	}

	private void ReadAndPreprocessProteins() {
		String StringBuffer = System.getProperty("user.dir");
		
		//TODO
		//File CurrentDirectory = new File("/home/rpark/projects2014_11_20_00_103_yaoyang/"); // for local debug
	//	File CurrentDirectory = new File("/home/diego/Desktop/jolene/"); // for local debug
		File CurrentDirectory = new File(StringBuffer); // for deployment
		// //TODO
		EndsWithFilter FilenameFilter = new EndsWithFilter(".out");
		EndsWithFilter SQTFilter = new EndsWithFilter("qt");
		EndsWithFilter DATFilter = new EndsWithFilter(".dat");
		String[] DirectoryListing = null;
		int FileIterator;
		Protein TempList;
		boolean IgnoreProteinFiles = true;
		File CurrentSubDir;
		File DTASelectTXT;
		DirectoryFilter SubDirFilter = new DirectoryFilter();
		String[] SubDirList = CurrentDirectory.list(SubDirFilter);
		int DirIterator;
		OUTFile OutFileList = new OUTFile();
		OUTFile Buffer;
		OUTFile PlaceHolder;
		LDAConfidence PeptideConfidence, ProteinConfidence;
		/*
		 * STAGE 1: Read files from disk whether from .out files or from .txt
		 * file. If DTASelect.txt exists, read it. Otherwise resort to parsing
		 * the .out files and match their sequences to loci listed in FASTA
		 * database.
		 */
		// Check if txt file has already been written
		DTASelectTXT = new File(CurrentDirectory, "DTASelect.txt");
		if (DTASelectTXT.exists()) {
			// if so, read the txt file
			System.out.println("Reading DTASelect.txt...");
			try {
				if(!Data.Cutoffs.noDB)
				{
					try {
						// Read sequest.params file
						Data.SequestParams = ParamsFile
								.ReadFile(CurrentDirectory);
					} catch (IOException failure) {
						System.out
								.println("Couldn't read sequest.params file.");
						System.exit(0);
					}

				}
				Data.ReadFromFile(DTASelectTXT);
				if(Data.Cutoffs.readFasta)
				{
					Data.LocusList.LookUpLoci(Data.SequestParams, Configuration,
							Data.Cutoffs);
				}
			} catch (IOException failure) {
				System.out.println(failure.getMessage());
				System.out
						.println("Couldn't read DTASelect.txt file.  You might try deleting it and re-running.");
				System.exit(0);
			}
		}
		/*
		 * If DTASelect.txt doesn't exist, read information in from .out files,
		 * sequest.params file, and database.
		 */
		else {
			/*
			 * If DTASelect.txt doesn't exist, read data: 0) If --Mascot, read
			 * .dat files 1) If .sqt files exist, read them. 2) If no .sqt,
			 * read: A) .outs in current directory (if --here is specified) B)
			 * .outs in subdirectories
			 */
			if (Data.Cutoffs.UseMascotOutput) {
				CurrentSubDir = CurrentDirectory;
				// Check to see if DTASelect.params lists the .dats to be
				// included
				try {
					File ParamsFile = new File(CurrentDirectory,
							"DTASelect.params");
					// Read a DTASelect.params file, grabbing the DAT list from
					// there.
					if (ParamsFile.exists()) {
						if (ParamsFile.canRead()) {
							FileReader InputFileReader = new FileReader(
									ParamsFile);
							BufferedReader Incoming = new BufferedReader(
									InputFileReader);
							String LineBuffer;
							// Skip the cutoffs listed.
							LineBuffer = Incoming.readLine();
							if (LineBuffer != null) {
								// At least the first line was there!
								LineBuffer = Incoming.readLine();
								// Queue up the DAT list
								while ((LineBuffer != null)
										&& (!LineBuffer
												.startsWith("[DAT List]"))) {
									LineBuffer = Incoming.readLine();
								}
								// Read each directory into an AuxInfo object
								if (LineBuffer != null) {
									AuxInfo AIList = new AuxInfo();
									AuxInfo AIRunner = AIList;
									int DATCount = 0;
									LineBuffer = Incoming.readLine();
									while ((LineBuffer != null)
											&& (LineBuffer.length() > 0)
											&& (!LineBuffer.startsWith("["))) {
										if (!LineBuffer.startsWith("#")) {
											AIRunner.Next = new AuxInfo();
											AIRunner = AIRunner.Next;
											AIRunner.Descriptor = LineBuffer;
											DATCount++;
										}
										LineBuffer = Incoming.readLine();
									}
									if (DATCount < 1)
										DirectoryListing = null;
									else {
										DirectoryListing = new String[DATCount];
										AIRunner = AIList.Next;
										DATCount = 0;
										while (AIRunner != null) {
											DirectoryListing[DATCount] = AIRunner.Descriptor;
											DATCount++;
											AIRunner = AIRunner.Next;
										}
									}
								}
							}
						}
					}
				} catch (IOException failure) {
					// If there's a failure, punt, reading all .dats in current
					// directory
					DirectoryListing = null;
				}
				if (DirectoryListing == null) {
					DirectoryListing = CurrentDirectory.list(DATFilter);
				}
				if (DirectoryListing.length > 0) {
					int SQTCount;
					// Mascot .dat files have been found in this
					// directory OR in the DTASelect.params file
					// Get the search configuration and database from the first
					// file.
					Data.SequestParams = ParamsFile.ReadMascotConfig(
							DirectoryListing[0], Configuration);
					Data.IDFileFormat = "Mascot";
					Data.SQTGenerator = "Mascot";
					System.out
							.println("Now reading Mascot results in DAT files:");
					Buffer = OutFileList;
					for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
						System.out.println("Now reading "
								+ DirectoryListing[FileIterator]);
						// Function returns numbers of SEQUEST .outs represented
						// by this file
						try {
							SQTCount = OUTFile.ReadMascotDAT(
									DirectoryListing[FileIterator], Buffer,
									Configuration);
							System.out
									.println("\t"
											+ SQTCount
											+ " identifications read.                              ");
							Data.UnfilteredOUTCount += SQTCount;
							while (Buffer.Next != null) {
								Buffer = Buffer.Next;
							}
						} catch (IOException failure) {
							System.out
									.println("Error while reading .dat files.");
							System.out.println(failure.getMessage());
							System.out.println(failure.getLocalizedMessage());
							System.out.println(failure.toString());
							System.exit(0);
						}
					}
				} else {
					System.out.println("No Mascot .dat file found.");
					System.exit(0);
				}
			} else {
				// Read sequest.params file
				if(!Data.Cutoffs.noDB)
				{
					try {
						Data.SequestParams = ParamsFile.ReadFile(CurrentDirectory);
					} catch (IOException failure) {
						System.out.println("Could not read sequest.params.");
						//System.out.println("Proceeding without db");

						System.exit(0); // TO REABILITATE!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}

				}

				// See if DB is present
				if(!Data.Cutoffs.noDB)
				{
					File Test;
					Test = new File(Data.SequestParams.DBName);
					if (!Test.canRead()) {
						System.out.println("Couldn't find database "
								+ Data.SequestParams.DBName);
						//System.exit(0); // TO REABILITATE!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				else
				{
					System.out.println("proceeding without db");
				}

				// check for presence of .sqt (unified .out) files
				CurrentSubDir = CurrentDirectory;
				DirectoryListing = CurrentDirectory.list(SQTFilter);
				if (DirectoryListing.length > 0) {
					Data.IDFileFormat = "SQT";
					int SQTCount;
					// SQT Files are present
					System.out.println("Now reading results in SQT files:");
					Buffer = OutFileList;
					try {
						BufferedReader Incoming = new BufferedReader(
								new FileReader(new File(CurrentDirectory,
										DirectoryListing[0])));
						String WholeLine = Incoming.readLine();
						String LineBuffer;
						StringTokenizer Parser;
						while ((WholeLine != null)
								&& (!WholeLine.startsWith("S\t"))) {
							Parser = new StringTokenizer(WholeLine);
							if (Parser.hasMoreTokens()) {
								LineBuffer = Parser.nextToken();
								if (LineBuffer.startsWith("SEQUEST")
										|| LineBuffer.equals("GutenTag")
										|| LineBuffer.equals("ProLuCID")
										|| LineBuffer.equals("STATS_PROB")
										|| LineBuffer.equals("Probe_pep")
                                                                        	|| LineBuffer.equals("Comet")
                                                                        ) 
                                                                {
									Data.SQTGenerator = LineBuffer;
									if (Parser.hasMoreTokens()) {
										Data.SQTGeneratorVersion = Parser
												.nextToken();
									}
								} else if (LineBuffer.equals("H")) {
									if (Parser.hasMoreTokens()) {
										LineBuffer = Parser.nextToken();
										if (LineBuffer.equals("SQTGenerator")
												&& (Parser.hasMoreTokens()))
											Data.SQTGenerator = Parser
													.nextToken();
										else if (LineBuffer
												.equals("SQTGeneratorVersion")
												&& (Parser.hasMoreTokens()))
											Data.SQTGeneratorVersion = Parser
													.nextToken();
									}

								}
							}
							WholeLine = Incoming.readLine();
						}
						Incoming.close();
						System.out.println("Identifications come from "
								+ Data.SQTGenerator + ", version "
								+ Data.SQTGeneratorVersion);
					} catch (IOException failure) {
						System.out
								.println("Error while attempting to identify program used to make SQT.");
					}
					for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
						System.out.println("Now reading "
								+ DirectoryListing[FileIterator]);
						// Function returns numbers of SEQUEST .outs represented
						// by this file
						try {
							SQTCount = OUTFile.ReadSQT(CurrentSubDir,
									DirectoryListing[FileIterator], Buffer,
									Configuration, Data.Cutoffs);
							System.out.println("\t" + SQTCount
									+ " identifications read. ");
							Data.UnfilteredOUTCount += SQTCount;

							while (Buffer.Next != null) {
								Buffer = Buffer.Next;
							}
						} catch (IOException failure) {
							System.out
									.println("Error while reading .sqt files.");
							System.out.println(failure.getMessage());
							System.out.println(failure.getLocalizedMessage());
							System.out.println(failure.toString());
							System.exit(0);
						}
					}
				} else {
					// if DTASelect.txt doesn't exist, read the .out files
					System.out.println("Reading .out files...");
					Data.IDFileFormat = "OUT";
					Data.SQTGenerator = "SEQUEST";
					Data.SQTGeneratorVersion = "";
					try {
						/*
						 * This if is triggered when the user has specified the
						 * -. (current directory only) option.
						 */
						if (Data.Cutoffs.LookInCurrentDirectory) {
							// Examine only this directory for .out files
							System.out
									.print("Now reading in current directory: ");
							DirectoryListing = CurrentDirectory
									.list(FilenameFilter);
							Data.UnfilteredOUTCount += DirectoryListing.length;
							System.out.println(new Integer(
									DirectoryListing.length).toString()
									+ " files");
							// Loop through all the .outs in the current
							// directory
							for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
								OUTFile.GetVersion(CurrentSubDir,
										DirectoryListing[0], Data);
								// Create an OUTFile object for each .out file
								// in this subdirectory
								Buffer = OUTFile.ReadOUT(CurrentSubDir,
										DirectoryListing[FileIterator], true,
										Configuration, Data.Cutoffs);
								if (Buffer != null) {
									// Since some .out files are multiple loci,
									// a linked list
									// might have been returned rather than a
									// single object.
									PlaceHolder = OutFileList.Next;
									OutFileList.Next = Buffer;
									while (Buffer.Next != null)
										Buffer = Buffer.Next;
									Buffer.Next = PlaceHolder;
								}
							}
						}
						/*
						 * if the user hasn't specified --here, the
						 * subdirectories of the current directory are searched
						 * for .out files, and a count for each directory is
						 * printed as the reading proceeds. If the SEQUEST box
						 * failed to process all DTAs, a directory of zero files
						 * should flag the user that something's gone wrong.
						 */
						else {
							boolean DontKnowVersionYet = true;
							// Loop through each subdirectory in search of .out
							// files
							for (DirIterator = 0; DirIterator < SubDirList.length; DirIterator++) {
								CurrentSubDir = new File(CurrentDirectory,
										SubDirList[DirIterator]);
								DirectoryListing = CurrentSubDir
										.list(FilenameFilter);
								// Let users know that we're starting to read
								// files in this subdirectory
								System.out.print("Now reading in "
										+ SubDirList[DirIterator] + ": ");
								Data.UnfilteredOUTCount += DirectoryListing.length;
								System.out.println(new Integer(
										DirectoryListing.length).toString()
										+ " files");
								// Loop through all the .outs in this
								// subdirectory
								for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
									if (DontKnowVersionYet) {
										OUTFile.GetVersion(CurrentSubDir,
												DirectoryListing[FileIterator],
												Data);
										DontKnowVersionYet = false;
									}
									// Create an OUTFile object for each .out
									// file in this subdirectory
									Buffer = OUTFile.ReadOUT(CurrentSubDir,
											DirectoryListing[FileIterator],
											false, Configuration, Data.Cutoffs);
									if (Buffer != null) {
										// Since some .out files are multiple
										// loci, a linked list
										// might have been returned rather than
										// a single object.
										PlaceHolder = OutFileList.Next;
										OutFileList.Next = Buffer;
										while (Buffer.Next != null)
											Buffer = Buffer.Next;
										Buffer.Next = PlaceHolder;
									}
								}
							}
						}
					} catch (IOException failure) {
						System.out.println("Error while reading .out files.");
						System.out.println(failure.getMessage());
						System.out.println(failure.getLocalizedMessage());
						System.out.println(failure.toString());
						System.exit(0);
					}
				}
			}
			// OutFileList.DebugPrint();
			System.out
					.println("Finished reading database search identifications.");

			OutFileList.SetOriginalNext();

			if (Data.Cutoffs.UseStatistics) {
				double[] LinearCoefficients;
				int iCharge, iTryptic, iModified;
				boolean Modified;
				int ChargeState;
				byte Tryptic;
				int NChargeStates = (Data.Cutoffs.UseChargeState ? Data.Cutoffs.MaxStatisticsCharge
						: 1);
				int NModified = (Data.Cutoffs.UseModStats ? 2 : 1);
				int NTryptic = (Data.Cutoffs.UseTrypticInfo ? 3 : 1);
				ChargeState = Data.Cutoffs.MinChargeState; //was 1
				for (iCharge = 0; iCharge < NChargeStates; iCharge++) {
					Tryptic = 0;
					Modified = false;
					System.out
							.println("Now computing linear coefficients for the following subset:");
					System.out.println("\tCharge state:\t\t\t"
							+ (Data.Cutoffs.UseChargeState ? ChargeState
									: "all"));

					PeptideConfidence = new LDAConfidence(OutFileList,
							Data.Cutoffs, ChargeState, Modified, Tryptic, true);

					PeptideConfidence.ComputeConfidence(); // THIS IS WHERE THE
															// COEFFICIENT ARE
															// CREATED
					LinearCoefficients = new double[PeptideConfidence.NDim];
					PeptideConfidence
							.SaveLinearCoefficients(LinearCoefficients);
					PeptideConfidence.CloseStatisticsFlag(OutFileList);
					PeptideConfidence = null;

					Tryptic = 0;
					for (iTryptic = 0; iTryptic < NTryptic; iTryptic++) {
						Modified = false;
						for (iModified = 0; iModified < NModified; iModified++) {
							System.out
									.println("Now computing peptide confidence levels for the following subset:");
							System.out
									.println("\tCharge state:\t\t\t"
											+ (Data.Cutoffs.UseChargeState ? ChargeState
													: "all"));
							System.out.println("\tEnzyme status:\t\t\t"
									+ (Data.Cutoffs.UseTrypticInfo ? Tryptic
											: "all"));
							System.out.println("\tModified status:\t\t"
									+ (Data.Cutoffs.UseModStats ? Modified
											: "all"));
							PeptideConfidence = new LDAConfidence(OutFileList,
									Data.Cutoffs, ChargeState, Modified,
									Tryptic, false);
							PeptideConfidence
									.GetLinearCoefficients(LinearCoefficients);
							PeptideConfidence.ComputeConfidence();
							PeptideConfidence.ReportScores(OutFileList);
							PeptideConfidence = null;
							Modified = (!Modified);
						}
						Tryptic++;
					}
					LinearCoefficients = null;
					ChargeState++;
				}
			}
			// Sort all the .out files by the identified locus
			System.out.println("Now grouping and sorting...");
			OutFileList.SortListByLocus();
			// Transmogrify from a list of .out files to a list of loci
			// with associated .out files
			Data.LocusList = OutFileList.CollateByLocus(Data.Cutoffs);
			// We're done with the OUTFile list. Let's free up memory.
			OutFileList = null;
			System.out.println("Now matching loci to database sequences...");
			/*
			 * Match up a database sequence and gene name with each identified
			 * locus. Determine where individual peptides align to protein
			 * sequence. Calculate MW and pI for each protein.
			 */
			if(!Data.Cutoffs.noDB)
			{
				try {
					Data.LocusList.LookUpLoci(Data.SequestParams, Configuration,
							Data.Cutoffs);
				} catch (IOException failure) {
					System.out.println("Error reading database file.");
				}
			}


			// Resort by possibly modified Sequence and charge state
			Data.LocusList.SortDTALists();
			// Compute protein confidence levels
			/*
			 * if (Data.Cutoffs.UseStatistics) {
			 * System.out.println("Now computing protein confidence levels...");
			 * ProteinConfidence = new LDAConfidence(Data.LocusList,
			 * Data.Cutoffs); ProteinConfidence.OldComputeConfidence();
			 * ProteinConfidence.ReportScores(Data.LocusList); ProteinConfidence
			 * = null; }
			 */
			// Write DTASelect.txt database to current directory
			System.out
					.println("Now writing information to DTASelect.txt.            ");
			Data.SourceDirectory = System.getProperty("user.dir");
			Data.CutoffsString = CommandLineOptions;
			Data.PrintTXT();
			// Data.PrintTXT(Data.Cutoffs, CommandLineOptions, "DTASelect.txt",
			// Data.SequestParams);
			System.out
					.println("Finished writing information to DTASelect.txt.");
		}
	}

	/*
	 * This filter allows one to list subdirectories found in another directory.
	 * It is used to list directories containing .out files.
	 */
	class DirectoryFilter implements FilenameFilter {
		public boolean accept(File dir, String name) {
			return (new File(dir, name).isDirectory());
		}
	}
}
