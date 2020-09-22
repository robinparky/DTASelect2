import java.io.*;
import java.util.*;
import java.net.URLEncoder;

// DTASelect DataSet
// David L. Tabb
// September 7, 2000
//modified April 4, 2017, Titus Jung titusjung@gmail.com

/*
 * Object to store Protein / DTAFile list with associated criteria.
 * Forms a linked list of datasets capable of producing
 * meta-analytical summaries.
 */

/* Each DataSet object stores information about which directory a
 * sample can be found in, what criteria are to be applied to it, and
 * what sequest.params are in effect for it.  Once the appropriate
 * files have been read into memory, the DataSet's Proteins are listed
 * under LocusList.  The criteria are used against this list, and then
 * DataSet.LocusList is used to create DTASelect.html output for this
 * combination.  The DataSets are the source of information for
 * Contrast to create its LocusSummary objects.  */

public class DataSet {
	String FriendlyName;
	String SourceDirectory;
	String CriteriaName;
	Protein LocusList = new Protein();
	Classification Classifieds = new Classification();
	SelectCriteria Cutoffs = new SelectCriteria();
	String CutoffsString;
	ParamsFile SequestParams;
	DataSet Next;
	Protein LocusCursor;
	boolean Removed = false;
	int UnfilteredLocusCount = 0;
	int UnfilteredOUTCount = 0;
	int UnfilteredSpectrumCount = 0;
	AuxInfo AddedInfo = new AuxInfo();
	/* The following fields are used only during a Contrast merge */
	BufferedReader Incoming;
	String CurrentLocus;
	String WholeLine;
	// These are the headers used in the various reports
	String hXCorr = "XCorr";
	String hCalcPreMass = "CalcM+H+";
	String hDeltCN = "DeltCN";
	String hSpScore = "SpScore";
	String hSpR = "SpR";
	String SQTGenerator = "ProLuCID";
	String SQTGeneratorVersion = "1.4";
	// The following variable tells us which format the IDs are stored in.
	// Acceptable values include:
	// OUT The original Sequest .out file format (with or without ID field)
	// Mascot Mascot-format MIME .dat files
	// SQT Unified SEQUEST result files
	String IDFileFormat = "OUT";

	// INDIVIDUAL FUNCTIONS

	public void SetColumnHeadings() {
		if (SQTGenerator.equals("GutenTag")) {
			hXCorr = "m*DotP";
			hSpScore = "%TIC";
		} else if (SQTGenerator.equals("ProLuCID")
				|| SQTGenerator.equals("Probe_pep")
				|| SQTGenerator.equals("STATS_PROB")) {
			hSpScore = "Prob Score";
		} else if (SQTGenerator.equals("Mascot")) {
			hXCorr = "IonsScore";
			hDeltCN = "Signif";
		} else if (SQTGenerator.equals("DBDigger")) {
			hXCorr = "Score";
			hSpScore = "Big8";
			hSpR = "Iso";
		} else if (SQTGenerator.equals("OMSSA")
				|| SQTGenerator.equals("OmssaXMLReader")) {
			hXCorr = "-Log(E)";
			hSpScore = "-Log(P)";
		} else if (SQTGenerator.equals("XTANDEM")
				|| SQTGenerator.equals("XTandemXMLReader")) {
			hXCorr = "-Log(E)";
		} else if (SQTGenerator.equals("SEQUEST")
				&& SQTGeneratorVersion.equals("3.0")) {
			hSpScore = "ZScore";
		} 
	}

	/*
	 * Caution! This function is used by both Contrast and DTASelect. Print
	 * reports on the remaining proteins in this dataset
	 */
	// added by Howard Choi
	public void PrintReports(String CommandLineOptions, IniFile Configuration,
			String RootName, String addFP) {
		String DBPrefix = RootName;
		if (RootName.equals("DTASelect"))
			DBPrefix = "DB";
		// Check to see if all proteins were eliminated. If so, skip reports.
		if (LocusList.CountProteins() == 0) {
			System.out.println("\tNo proteins pass criteria!");
		}
		/*
		 * If the new CGIs are used, we need to rerun the sequence coverage
		 * generator because we used the old consensus URL to detect identicals.
		 */
		if (IDFileFormat.equals("SQT"))
			Cutoffs.DisplayType = 1;
		else if (IDFileFormat.equals("Mascot"))
			Cutoffs.DisplayType = 2;
		if (Configuration.UseDepthCGI) {
			LocusList.CalculateCoverageForList(true);
		}
		// If some proteins were retained, print all selected reports.
		if (Cutoffs.PrintSimilarities) {
			System.out.println("\tCreating protein similarity table...");
			this.PrintProteinSimilarity(RootName);
		}
		if (Cutoffs.PrintXML) {
			System.out.println("\tCreating XML output...");
			LocusList.PrintXML(Cutoffs, CommandLineOptions, SequestParams,
					RootName + "-filter.xml", this.SourceDirectory);
		}
		/*
		 * Users may indicate their desire to have a script or batch file
		 * created for them that will allow them to copy the selected .dta and
		 * .out files to a separate directory. If they have used the --copy
		 * parameter, this script will be sent out to the console.
		 */
		if (Cutoffs.PrintCopyList && this.IDFileFormat.equals("OUT")) {
			System.out.println("\tCreating list of files to copy...");
			LocusList.PrintCPList(Configuration);
		}
		// Sort the peptides by sequence position
		LocusList.DisplaySortDTALists();
		// Sort the proteins appropriately
		if (Cutoffs.UseClassifications) {
			System.out.println("\tReading classification list...");
			this.StructureByClassification();
		}
		if (Cutoffs.UseAuxInfo) {
			System.out.println("\tReading auxiliary info...");
			this.IncorporateAuxInfo();
		}
		// Sort by class and then by other measures
		LocusList.SortByClass();
		if (Cutoffs.PrintMods) {
			System.out.println("\tCreating modification report...");
			this.PrintModReport(Configuration, RootName, this.SourceDirectory);
		}
		if (Cutoffs.PrintAlignment) {
			System.out.println("\tCreating alignment report...");
			this.PrintAlignReport(Configuration, RootName);
		}
		if (Cutoffs.PrintDB) {
			System.out.println("\tCreating Proteins database...");
			LocusList.PrintDatabase(Cutoffs, SequestParams, DBPrefix);
		}
		if (Cutoffs.DisplayBirdsEye) {
			BirdsEye GUI = new BirdsEye(this);
		}
		System.out
				.println("\tCreating DTASelect.html and DTASelect-filter.txt...");
		// added by Howard Choi
		this.PrintOutput(CommandLineOptions, RootName, Configuration, DBPrefix,
				addFP);
		if (Cutoffs.UseClassifications && Cutoffs.PrintDB)
			Classifieds.PrintDatabase(DBPrefix);
	}

	// added by Howard Choi
	public String fpGreper(String cmdLine) {
		String[] tmpArr = cmdLine.split(" ");
		String fpVal = "";
		for (int i = 0; i < tmpArr.length; i++) {
			if (tmpArr[i].matches("--fp")) {
				fpVal = tmpArr[i + 1];
			}
		}
		return fpVal;
	}

	// added by Howard Choi
	public String protFPInput = "";
	public String protFPOutput = "";

	// added by Howard Choi
	public void setProtFPInput(String protFPInput) {
		this.protFPInput = protFPInput;
	}

	// added by Howard Choi
	public void setProtFPOutput(String protFPOutput) {
		this.protFPOutput = protFPOutput;
	}

	// added by Howard Choi
	public String getProtFPInput() {
		return this.protFPInput;
	}

	// added by Howard Choi
	public String getProtFPOutput() {
		return this.protFPOutput;
	}

	// added by Howard Choi
	public String pepFPInput = "";
	public String pepFPOutput = "";

	// added by Howard Choi
	public void setPepFPInput(String pepFPInput) {
		this.pepFPInput = pepFPInput;
	}

	// added by Howard Choi
	public void setPepFPOutput(String pepFPOutput) {
		this.pepFPOutput = pepFPOutput;
	}

	// added by Howard Choi
	public String getPepFPInput() {
		return this.pepFPInput;
	}

	// added by Howard Choi
	public String getPepFPOutput() {
		return this.pepFPOutput;
	}

	// added by Howard Choi
	public String getTmpCmdLineOptions(String oriCmdOptions, String addFP) {
		String[] tmpArr = oriCmdOptions.split(" ");
		String tmpCmdOptions = "";
		int i = 0;
		// boolean shouldAdd = false;

		while (i < tmpArr.length) {
			if (tmpArr[i].startsWith("--fp") && addFP.startsWith("--")) {
				tmpCmdOptions += addFP + " ";
				i += 1;
				// shouldAdd = true;
			} else {
				tmpCmdOptions += tmpArr[i] + " ";
			}

			i += 1;
		}

		// if(shouldAdd){
		// tmpCmdOptions += addFP;
		// }
		return tmpCmdOptions;
	}

	/*
	 * Print an HTML web page and TXT file from the loci and DTAFiles remaining
	 * in memory after the filtering process. NOTE: after this funcion is
	 * called, the integrity of the protein list is destroyed.
	 */
	// added by Howard Choi
	public void PrintOutput(String CommandLineOptions, String FileName,
			IniFile Config, String DBPrefix, String additionalFP) {
		String Citation = "Reference: Cociorva D, Tabb DL, Yates JR 3rd, \"Validation of tandem mass spectrometry database search results using DTASelect\", Current Protocols in Bioinformatics, Chapter 13, Unit 13.4 (2007).";
		String Descrip;
		String tab = "\t";
		String newline = "\n";
		String CellSep = "</td><td>";
		String DBName = Cutoffs.noDB ? "null":SequestParams.DBName;
		String Consensus = "<a TARGET=\"Win1\" href=\"" + Config.BasicSeqCov
				+ "?Db=" + DBName + "&Ref=";
		String SeqCov = "<a TARGET=\"Win1\" href=\"" + Config.DepthSeqCov + "?"
				+DBName+ "&";
		String DisplayIons = "</td><td><A TARGET=\"Win1\" HREF=\""
				+ Config.DTASpecDisplay + "?Dta=";
		String ShowOut = "</a></td><td><a TARGET=\"Win1\" href=\""
				+ Config.DTAOUTDisplay + "?OutFile=";
		String NewShow = "</td><td><a target=\"Win1\" HREF=\""
				+ Config.SQTDisplay
				+ "?Dr="
				+ (this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
						: this.SourceDirectory) + "&Da=";
		String BLAST1 = "</td><td><a target=\"Win1\" href=\"" + Config.BLAST;
		String BLAST2 = Config.BLASTArgs + "\">";
		String ProtEval = "<a target=\"Win1\" href=\""
				+ Config.ProtValidation
				+ "?"
				+ URLEncoder
						.encode(this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
								: this.SourceDirectory) + "&";
		// StringBuffer CGIBuild;
		boolean printHTML = this.Cutoffs.printHTML;
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File SpectrumDirectory;
		String SpectrumFile;
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing =null;
		File OutputTXTFile;
		FileWriter OutputTXTFileWriter;
		BufferedWriter OutgoingTXT;
		Protein Runner;
		Protein MatchRunner;
		Protein Mega;
		DTAFile DTARunner;
		int Counter = 0;
		int ForwardCounter = 0;
		int DecoyCounter = 0;
		int fcount  =0;
		int rcount =0;
		float ProteinFPRate = 0f;
		int RProteinCounter = 0;
		float RProteinFPRate = 0f;
		int RForwardCounter=0;
		int RDecoyCounter = 0;
		int NRPeptideCounter = 0;
		int NRForwardPeptideCounter = 0;
		int NRDecoyPeptideCounter = 0;
		float NRPeptideFPRate = 0f;
		int NRSpectrumCounter = 0;
		int NRForwardSpectrumCounter = 0;
		int NRDecoySpectrumCounter = 0;
		float NRSpectrumFPRate = 0f;
		int UnclassifiedNRCount = 0;
		int UnclassifiedRCount = 0;
		boolean MSBug = false;
		StringBuffer StaticModCGI = new StringBuffer();
		String StaticModString;
		String Temp1 = "DTASelect";
		String Temp2;
		StringTokenizer Parser;
		StringBuffer Similarities;
		int MatchingPeptides;
		int PepCount;
		int ScoreBuffer;
		Classification LastClassification = null;
		Classification CurrentClassification = null;
		AuxInfo AIRunner;
		Point PoRunner;
		Runner = LocusList.Next;
		String cgiString =  Cutoffs.noDB ? "null":SequestParams.CGIString();
		StaticModCGI.append(cgiString);
		StaticModString = StaticModCGI.toString() + "&Pep=";
		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			if(printHTML)
			{
				OutputFile = new File(CurrentDirectory, FileName + ".html");
				OutputFileWriter = new FileWriter(OutputFile);
				Outgoing = new BufferedWriter(OutputFileWriter);
			}
			OutputTXTFile = new File(CurrentDirectory, FileName + "-filter.txt");
			OutputTXTFileWriter = new FileWriter(OutputTXTFile);
			OutgoingTXT = new BufferedWriter(OutputTXTFileWriter);
			// added by Howard Choi
			setProtFPInput(fpGreper(CommandLineOptions));
			setPepFPInput(fpGreper(CommandLineOptions));
			// Write HTML header
			if(printHTML) Outgoing.write("<HTML><HEAD><TITLE>\nD ");
			if (FileName.equals("DTASelect")) {
				Parser = new StringTokenizer(System.getProperty("user.dir"),
						System.getProperty("file.separator"));
				while (Parser.hasMoreTokens())
					Temp1 = Parser.nextToken();
				if(printHTML) Outgoing.write(Temp1);
			} else {
				if(printHTML) Outgoing.write(FileName);
			}
			if(printHTML) Outgoing.write("</title>\n<style type=\"text/css\">\n");
			if(printHTML) Outgoing.write("td {\n");
			if(printHTML) Outgoing.write("\tfont-family: Courier, monospace;\n");
			if(printHTML) Outgoing.write("\tfont-weight: normal;\n");
			if(printHTML) Outgoing.write("\tfont-size: 10pt;\n");
			if(printHTML) Outgoing.write("\tpadding: 0px 10px\n");
			if(printHTML) Outgoing.write("\t}\n");
			if(printHTML) Outgoing.write("th {\n");
			if(printHTML) Outgoing.write("\tfont-family: \"MS Arial\", sans-serif;\n");
			if(printHTML) Outgoing.write("\tfont-weight: bold;\n");
			if(printHTML) Outgoing.write("\tfont-size: 10pt\n\t}\n");
			if(printHTML) Outgoing.write("em {\n");
			if(printHTML) Outgoing.write("\tcolor: red;\n");
			if(printHTML) Outgoing.write("\tfont-family: \"MS Arial\", sans-serif;\n");
			if(printHTML) Outgoing.write("\tfont-style: normal;\n");
			if(printHTML) Outgoing.write("\tfont-weight: bold;\n");
			if(printHTML) Outgoing.write("\tfont-size: 12pt\n");
			if(printHTML) Outgoing.write("\t}\n");
			if(printHTML) Outgoing.write("dfn {\n");
			if(printHTML) Outgoing.write("\tfont-family: \"MS Arial\", sans-serif;\n");
			if(printHTML) Outgoing.write("\tfont-style: normal;\n");
			if(printHTML) Outgoing.write("\tfont-weight: normal;\n");
			if(printHTML) Outgoing.write("\tfont-size: 10pt\n");
			if(printHTML) Outgoing.write("\t}\n");
			if(printHTML) Outgoing.write("</style>\n");
			if(printHTML) Outgoing.write("</head>");
			if(printHTML) Outgoing.write("<body>\n");
			if(printHTML) Outgoing.write(Protein.Version() + "<BR>\n");
			if(printHTML) Outgoing.write(Citation + "<BR>\n");
			if(printHTML) Outgoing.write((this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
					: this.SourceDirectory)
					+ "<BR>\n");
			String dbName = Cutoffs.noDB ? "null" : SequestParams.DBName;
			if(printHTML) Outgoing.write(dbName + "<BR>\n");
			if(printHTML) Outgoing.write(this.SQTGenerator + " " + this.SQTGeneratorVersion
					+ " in " + this.IDFileFormat + " format.<BR>\n");
			// added by Howard Choi
			String tmpCommandLineOptions = getTmpCmdLineOptions(
					CommandLineOptions, additionalFP);
			if(printHTML) Outgoing.write(tmpCommandLineOptions + "<BR>\n");
			// If other HTML reports were requested, link to them
			if (Cutoffs.PrintAlignment) {
				if(printHTML) Outgoing.write("<a TARGET=\"mods\" href=\"" + FileName
						+ "-align.html\">View</a> alignment report.<BR>\n");
			}
			if (Cutoffs.PrintMods) {
				if(printHTML) Outgoing.write("<a TARGET=\"mods\" href=\"" + FileName
						+ "-mods.html\">View</a> modification report.<BR>\n");
			}
			if(printHTML) Outgoing.write("<a href=\"#Summary\"> Jump </a> to the summary table.<br><BR>\n");
			String htmlMods = Cutoffs.noDB ? "null":SequestParams.GetHTMLMods();
			if(printHTML) Outgoing.write(htmlMods + newline);
			// Write the criteria in effect for this run
			if(printHTML) Outgoing.write("<table border>");
			if(printHTML) Outgoing.write(Cutoffs.PrintCriteria("<TR><TD>", "</TD></TR>\n",
					CellSep));
			if(printHTML) Outgoing.write("</table>\n");
			if(printHTML) Outgoing.write("<h4>Locus Key:</h4>\n<table border><TR><td><a href=\"\">Validation Status</a></td><TD><Font COLOR=\"red\">"
					+ "Locus</td><td>"
					+ "Sequence Count</td><td>Spectrum Count</td><td>"
					+ "<a href=\"\">Sequence Coverage</a></td><TD>Length</td>"
					+ "<TD>MolWt</td><TD>pI</td>");
			// Create TXT file header, too.
			OutgoingTXT.write(Protein.Version() + newline);
			OutgoingTXT
					.write((this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
							: this.SourceDirectory)
							+ newline);
			String db = Cutoffs.noDB ? "null":SequestParams.DBName;
			OutgoingTXT.write(db+ newline);
			OutgoingTXT.write(this.SQTGenerator + " "
					+ this.SQTGeneratorVersion + " in " + this.IDFileFormat
					+ " format.\n");
			// added by Howard Choi
			OutgoingTXT.write(tmpCommandLineOptions + newline);
			// Write the criteria in effect for this run
			OutgoingTXT
					.write(Cutoffs.PrintCriteria("", newline, tab) + newline);
			OutgoingTXT
					.write("Locus\tSequence Count\tSpectrum Count\tSequence Coverage\tLength\tMolWt\tpI\tValidation Status\t");
			AIRunner = AddedInfo.Next;
			while (AIRunner != null) {
				if(printHTML) Outgoing.write("<td>" + AIRunner.Descriptor + "</td>");
				OutgoingTXT.write(AIRunner.Descriptor + tab);
				AIRunner = AIRunner.Next;
			}
			if(printHTML) Outgoing.write("<TD>NSAF</td><TD>EMPAI</td><TD>Descriptive Name</td></tr></table>\n");
			if(printHTML) Outgoing.write("<h4>Similarity Key:</h4>\n<table border><tr><td>Locus</td><td># of identical peptides</td>"
					+ "<td># of differing peptides</td></tr></table>\n");
			OutgoingTXT
					.write("NSAF\tEMPAI\tDescriptive Name\nUnique\tFileName\t"
							+ this.hXCorr
							+ tab
							+ this.hDeltCN
							+ tab
							+ "Conf%\tM+H+\t"
							+ this.hCalcPreMass
							+ (Cutoffs.DisplayDeltaMass && Cutoffs.ExtraColumns ? "\tPPM"
									: "")
							+ "\tTotalIntensity\t"
							+ this.hSpR
							+ tab
							+ this.hSpScore
							+ (Cutoffs.DisplayPI && Cutoffs.ExtraColumns ? "\tpI"
									: "")
							+ (Cutoffs.DisplayKD && Cutoffs.ExtraColumns ? "\tKD"
									: "")
							+ (Cutoffs.DisplayBB && Cutoffs.ExtraColumns ? "\tBB"
									: "")
							+ (Cutoffs.DisplayHPLC && Cutoffs.ExtraColumns ? "\tHPLC"
									: "") + "\tIonProportion\t"
							+ "Redundancy\tSequence"
							+ (Cutoffs.DisplaySeqPos ? "\tSeqPosition" : ""));
			// Initialize our counter of loci
			Counter = 0;
			DecoyCounter = 0;
			fcount =0;
			rcount =0;
			ForwardCounter = 0;

			float NSAF_sum = 0;

			// calculate NSAF
			while (Runner != null) {
				MatchRunner = Runner;
				
				/*if (MatchRunner.Locus.equals("tr|Q5T6W5|Q5T6W5_HUMAN")) {
					System.out.println();
				}*/
				
				NSAF_sum += MatchRunner.getNSAF_normalizer();
				Runner = Runner.Next;
			}

			if (NSAF_sum == 0) {
				NSAF_sum = 1;
			}

			Runner = LocusList.Next;

			// Print off all fields, continuing down list
			while (Runner != null) {
				Counter++;
				MatchRunner = Runner;
				/*if (MatchRunner.Locus.startsWith(Cutoffs.DecoyLabel)) {
					DecoyCounter++;
					//MatchRunner.IsDecoy = true;
				} else {
					ForwardCounter++;
					//MatchRunner.IsDecoy = false;
				}*/
				//MatchRunner.DTASetDecoy();
				if(printHTML) Outgoing.write("<HR><table>\n");
				// If we're using classifications supplied by the user...
				if (Cutoffs.UseClassifications) {
					// If we haven't a class yet
					if (LastClassification == null) {
						CurrentClassification = Classifieds
								.GetClass(MatchRunner.Classification);
					} else {
						// If we were already on a class, but this one doesn't
						// match
						if (LastClassification.Identifier != MatchRunner.Classification) {
							CurrentClassification = Classifieds
									.GetClass(MatchRunner.Classification);
						}
					}
					// Only print a new label if this is a change of
					// classification
					if (CurrentClassification != LastClassification) {
						LastClassification = CurrentClassification;
						// If this number corresponded to a class
						if (CurrentClassification != null) {
							if(printHTML) Outgoing.write("<h2><a name=\"Class"
									+ LastClassification.Identifier + "\">"
									+ LastClassification.Descriptor
									+ "</h2></table>\n<table>");
							OutgoingTXT.write("\nClass\t"
									+ LastClassification.Descriptor);
						} else {
							if(printHTML) Outgoing.write("<h2><a name=\"UnClass\">Uncategorized Loci</h2></table>\n<table>");
							OutgoingTXT.write("\nClass\tUnclassified");
						}
					}
					// Add to this classification's count
					if (CurrentClassification != null) {
						CurrentClassification.NRCountObserved++;
					} else {
						UnclassifiedNRCount++;
					}
				}
				boolean isForward = false;
				while (MatchRunner != null) {
					if (Cutoffs.UseClassifications) {
						// Add to this classification's count
						if (CurrentClassification != null) {
							CurrentClassification.RCountObserved++;
						} else {
							UnclassifiedRCount++;
						}
					}
					RProteinCounter++;
					OutgoingTXT.write(newline);
					// Write the Locus-specific information
					Temp1 = new Integer(MatchRunner.NPeptides).toString();
					Temp2 = new Integer(MatchRunner.NSpectra).toString();
					if(printHTML) Outgoing.write("<TR><TD>" + ProtEval
							+ URLEncoder.encode(MatchRunner.Locus) + "\">"
							+ new Character(MatchRunner.Validated).toString()
							+ "</a>" + CellSep + "<em><a name=\""
							+ MatchRunner.Locus + "\">" + MatchRunner.Locus
							+ "</a></td><td>" + Temp1 + CellSep + Temp2
							+ CellSep);
					OutgoingTXT.write(MatchRunner.Locus + tab + Temp1 + tab
							+ Temp2 + tab);
					
					if (MatchRunner.Locus.startsWith("Reverse")) {
						RDecoyCounter++;
					} else {
						isForward = true;
						RForwardCounter++;
					}
					
					if (MatchRunner.SequenceLength == 0) {
						if(printHTML) Outgoing.write("<a href=\"\">0%</a></td><td><dfn>Not found in database</td></td></tr>\n");
						OutgoingTXT.write("0%\t0\t0\t0\t0\t0\t0\tNot found in database");
					} else {
						Descrip = MatchRunner.Gene;
						// Remove inappropriate characters for HTML from the
						// description...
						Descrip.replace('&', '*');
						Descrip.replace('<', '*');
						Descrip.replace('>', '*');
						Descrip.replace('\"', '*');
						Descrip.replace('\'', '*');
						if (Config.UseDepthCGI) {
							if(printHTML) Outgoing.write(SeqCov
									+ URLEncoder.encode(MatchRunner.Locus));
						} else {
							if(printHTML) Outgoing.write(Consensus + MatchRunner.Locus
									+ "&Pep=");
						}
						if(printHTML) Outgoing.write(MatchRunner.Coverage
								+ "\">"
								+ MatchRunner.SequenceCoverage
								+ "%"
								+ "</a></td><td>"
								+ new Integer(MatchRunner.SequenceLength)
										.toString()
								+ /*
								 * CellSep +
								 * Protein.RoundTo(100f*MatchRunner.ProtConf,1)
								 * + "%" + CellSep +
								 * Protein.RoundTo(100f*MatchRunner.ProtFP,1) +
								 * "%" +
								 */
								CellSep
								+ new Integer(Math.round(MatchRunner.MolWt))
										.toString()
								+ CellSep
								+ new Float(Protein.RoundTo(MatchRunner.pI, 1))
										.toString() + CellSep);
						OutgoingTXT.write(MatchRunner.SequenceCoverage
								+ "%\t"
								+ new Integer(MatchRunner.SequenceLength)
										.toString()
								/*
								 * + tab +
								 * Protein.RoundTo(100f*MatchRunner.ProtConf,1)
								 * + tab +
								 * Protein.RoundTo(100f*MatchRunner.ProtFP,1)
								 */
								+ tab
								+ new Integer(Math.round(MatchRunner.MolWt))
										.toString()
								+ tab
								+ new Float(Protein.RoundTo(MatchRunner.pI, 1))
										.toString()
								+ tab
								+ new Character(MatchRunner.Validated)
										.toString() + tab);
						if (MatchRunner.AuxInfo != null) {
							AIRunner = AddedInfo.Next;
							PoRunner = MatchRunner.AuxInfo.Next;
							while ((AIRunner != null) && (PoRunner != null)) {
								if (AIRunner.IsFloat) {
									if(printHTML) Outgoing.write(new Float(PoRunner.Intensity)
											.toString());
									OutgoingTXT.write(new Float(
											PoRunner.Intensity).toString());
								} else {
									if(printHTML) Outgoing.write(new Integer(Math
											.round(PoRunner.Intensity))
											.toString());
									OutgoingTXT.write(new Integer(Math
											.round(PoRunner.Intensity))
											.toString());
								}
								if(printHTML) Outgoing.write(CellSep);
								OutgoingTXT.write("\t");
								AIRunner = AIRunner.Next;
								PoRunner = PoRunner.Next;
							}
						}

						// Write the NSAF and EMPAI values of this locus to HTML
						// and TXT output files

						MatchRunner.calculateNSAFandEMPAI(NSAF_sum);

						if(printHTML) Outgoing.write(MatchRunner.NSAF + CellSep
								+ MatchRunner.EMPAI);
						OutgoingTXT.write(MatchRunner.NSAF + "\t"
								+ MatchRunner.EMPAI);

						if(printHTML) Outgoing.write(CellSep);
						OutgoingTXT.write("\t");

						// Write the descriptive name of this locus to HTML and
						// TXT output files
						if(printHTML) Outgoing.write("<dfn>" + Descrip + "</td></tr>\n");
						OutgoingTXT.write(MatchRunner.Gene);

					}
					MatchRunner = MatchRunner.IdenticalLoci; //here is the problsm with similar stuff that has 0
					if (MatchRunner != null) {
						MatchRunner.getNSAF_normalizer();
					}
					
				}
				MatchRunner = Runner;
				if(isForward)
				{
					MatchRunner.IsDecoy = false;
					fcount++;
				}
				else {
					MatchRunner.IsDecoy = true;
					rcount++;
				}
				MatchRunner.DTASetDecoy();
				MatchRunner = null;
				if(printHTML) Outgoing.write("</table>\n");

				// Print DTAFile data for this protein
				if(printHTML) Outgoing.write("<table>\n");
				if(printHTML) Outgoing.write("<tr><TH> <th>Filename<th>" + hXCorr + "<TH>"
						+ hDeltCN + "<th>Conf%<th>ObsM+H+<th>" + hCalcPreMass
						+ (Cutoffs.DisplayDeltaMass ? "<th>PPM" : "") + "<th>"
						+ hSpR + "<th>" + hSpScore
						+ (Cutoffs.DisplayPI ? "<th>pI" : "")
						+ (Cutoffs.DisplayKD ? "<th>KD" : "")
						+ (Cutoffs.DisplayBB ? "<th>BB" : "")
						+ (Cutoffs.DisplayHPLC ? "<th>HPLC" : "")
						+ "<th>Ion%<th>#<th>Sequence<TH></tr>\n");
				DTARunner = Runner.DTAs.Next;
				//if(printHTML) Outgoing.write("<tr><td>");
				//OutgoingTXT.write("\n");
				while (DTARunner != null) {

					/*if(Cutoffs.PrintOnlyUnique && !DTARunner.UniqueToLocus)
					{
						DTARunner = DTARunner.Next;
						continue;
					}*/
					if(printHTML) Outgoing.write("<tr><td>");
					OutgoingTXT.write("\n");
					if (DTARunner.UniqueToLocus) {
						if(printHTML) Outgoing.write("*");
						OutgoingTXT.write("*");
					}
					/*
					else if(Cutoffs.PrintOnlyUnique){
						DTARunner = DTARunner.Next;
						continue;
					}*/
					if (DTARunner.EquivSeq != 1) {
						if(printHTML) Outgoing.write(new Byte(DTARunner.EquivSeq).toString());
						OutgoingTXT.write(new Byte(DTARunner.EquivSeq)
								.toString());
					}
					if (DTARunner.Subdirectory.equals("pta")
							|| DTARunner.Subdirectory.equals("phs")) {
						if(printHTML) Outgoing.write("P");
						OutgoingTXT.write("P");
					}
					switch (DTARunner.Validated) {
					case 'Y':
						if(printHTML) Outgoing.write("<FONT COLOR=\"GREEN\">Y</FONT>");
						OutgoingTXT.write("Y");
						break;
					case 'M':
						if(printHTML) Outgoing.write("<FONT COLOR=\"ORANGE\">M</FONT>");
						OutgoingTXT.write("M");
						break;
					case 'N':
						if(printHTML) Outgoing.write("<FONT COLOR=\"RED\">N</FONT>");
						OutgoingTXT.write("N");
						break;
					}
					double dm = Cutoffs.showCorrectedDmValue ? DTARunner.Shifted_PPM_Offset : DTARunner.Adjusted_PPM_Offset;
					OutgoingTXT
							.write(tab
									+ DTARunner.FileName
									+ tab
									+ DTARunner.XCorr
									+ tab
									+ DTARunner.DeltCN
									+ tab
									+ Protein.RoundTo(100f * (1.0f - new Float(
											DTARunner.PepFP).floatValue()), 1)
									+ tab
									+ DTARunner.PrecursorMass
									+ tab
									+ DTARunner.CalcPreMass
									+ tab
									+ (Cutoffs.DisplayDeltaMass
											&& Cutoffs.ExtraColumns ? Protein
											.RoundTo(
													new Float(
															dm)
															.floatValue(), 1)
											+ tab
											: "")
									+ DTARunner.TotalIntensity
									+ tab
									+ DTARunner.Sp
									+ tab
									+ DTARunner.SpScore
									+ tab
									+ (Cutoffs.DisplayPI
											&& Cutoffs.ExtraColumns ? Protein
											.RoundTo(DTARunner.CalculatepI(), 2)
											+ tab
											: "")
									+ (Cutoffs.DisplayKD
											&& Cutoffs.ExtraColumns ? Protein
											.RoundTo(DTARunner
													.CalculateKyteDoolittle(),
													1)
											+ tab : "")
									+ (Cutoffs.DisplayBB
											&& Cutoffs.ExtraColumns ? Protein
											.RoundTo(DTARunner
													.CalculateBullBreese(), 1)
											+ tab : "")
									+ (Cutoffs.DisplayHPLC
											&& Cutoffs.ExtraColumns ? Protein
											.RoundTo(DTARunner
													.CalculateHPLCpH34(), 1)
											+ tab : "")
									+ Protein.RoundTo(
											100f * DTARunner.IonProportion, 1)
									+ tab
									+ DTARunner.Redundancy
									+ tab
									+ DTARunner.Sequence
									+ (Cutoffs.DisplaySeqPos ? tab
											+ DTARunner.SequencePosition : ""));
					/* Print details for each spectrum */
					SpectrumDirectory = new File(this.SourceDirectory,
							DTARunner.Subdirectory);
					if (Cutoffs.DisplayType > 0) {
						// Note that we'll need to add specific
						// handling for Mascot CGIs when feasible.
						// DisplayType == 2 for Mascot.
						 dm = Cutoffs.showCorrectedDmValue? DTARunner.Shifted_PPM_Offset : DTARunner.Adjusted_PPM_Offset;
						if(printHTML) Outgoing.write(NewShow
								+ DTARunner.ShowString()
								+ CellSep
								+ new Float(DTARunner.XCorr).toString()
								+ CellSep
								+ new Float(DTARunner.DeltCN).toString()
								+ CellSep
								+ new Float(Protein.RoundTo(
										100f - 100f * DTARunner.PepFP, 1))
										.toString()
								+ "%"
								+ CellSep
								+ new Float(DTARunner.PrecursorMass).toString()
								+ CellSep
								+ new Float(DTARunner.CalcPreMass).toString()
								+ (Cutoffs.DisplayDeltaMass ? CellSep
										+ new Float(Protein.RoundTo(new Float(
												dm)
												.floatValue(), 1)).toString()
										: "")
								+ CellSep
								+ new Integer(DTARunner.Sp).toString()
								+ CellSep
								+ new Float(DTARunner.SpScore).toString()
								+ (Cutoffs.DisplayPI ? CellSep
										+ new Float(Protein.RoundTo(
												DTARunner.CalculatepI(), 2))
												.toString() : "")
								+ (Cutoffs.DisplayKD ? CellSep
										+ new Float(Protein.RoundTo(DTARunner
												.CalculateKyteDoolittle(), 1))
												.toString() : "")
								+ (Cutoffs.DisplayBB ? CellSep
										+ new Float(
												Protein.RoundTo(DTARunner
														.CalculateBullBreese(),
														1)).toString() : "")
								+ (Cutoffs.DisplayHPLC ? CellSep
										+ new Float(Protein.RoundTo(
												DTARunner.CalculateHPLCpH34(),
												1)).toString() : "")
								+ CellSep
								+ new Float(Protein.RoundTo(
										100f * DTARunner.IonProportion, 1))
										.toString() + "%</td><td>"
								+ new Integer(DTARunner.Redundancy).toString()
								+ BLAST1 + DTARunner.TrimmedSequence() + BLAST2
								+ DTARunner.Sequence + "</a></td><td>");
					} else {
						try {
							SpectrumFile = new File(SpectrumDirectory,
									DTARunner.FileName).getCanonicalPath();
						} catch (IOException failure) {
							MSBug = true;
							SpectrumFile = "";
						}
						 dm = Cutoffs.showCorrectedDmValue? DTARunner.Shifted_PPM_Offset : DTARunner.Adjusted_PPM_Offset;

						if(printHTML) Outgoing.write(DisplayIons
								+ SpectrumFile
								+ ".dta"
								+ (DTARunner.Modified ? DTARunner
										.DisplayIonsString(SequestParams) : "")
								+ StaticModString
								+ DTARunner.TrimmedSequence()
								+ "\">"
								+ DTARunner.FileName
								+ ShowOut
								+ SpectrumFile
								+ ".out\">"
								+ new Float(DTARunner.XCorr).toString()
								+ "</a></td><td>"
								+ new Float(DTARunner.DeltCN).toString()
								+ /*
								 * CellSep + new
								 * Float(Protein.RoundTo(100f*DTARunner
								 * .PepConf,1)).toString() + "%" +
								 */
								CellSep
								+ new Float(Protein.RoundTo(
										100f - 100f * DTARunner.PepFP, 1))
										.toString()
								+ "%"
								+ CellSep
								+ new Float(DTARunner.PrecursorMass).toString()
								+ CellSep
								+ new Float(DTARunner.CalcPreMass).toString()
								+ (Cutoffs.DisplayDeltaMass ? CellSep
										+ new Float(Protein.RoundTo(new Float(dm).floatValue(), 1)).toString()
										: "")
								+ CellSep
								+ new Integer(DTARunner.Sp).toString()
								+ CellSep
								+ new Float(DTARunner.SpScore).toString()
								+ (Cutoffs.DisplayPI ? CellSep
										+ new Float(Protein.RoundTo(
												DTARunner.CalculatepI(), 2))
												.toString() : "")
								+ (Cutoffs.DisplayKD ? CellSep
										+ new Float(Protein.RoundTo(DTARunner
												.CalculateKyteDoolittle(), 1))
												.toString() : "")
								+ (Cutoffs.DisplayBB ? CellSep
										+ new Float(
												Protein.RoundTo(DTARunner
														.CalculateBullBreese(),
														1)).toString() : "")
								+ (Cutoffs.DisplayHPLC ? CellSep
										+ new Float(Protein.RoundTo(
												DTARunner.CalculateHPLCpH34(),
												1)).toString() : "")
								+ CellSep
								+ new Float(Protein.RoundTo(
										100f * DTARunner.IonProportion, 1))
										.toString() + "%</td><td>"
								+ new Integer(DTARunner.Redundancy).toString()
								+ BLAST1 + DTARunner.TrimmedSequence() + BLAST2
								+ DTARunner.Sequence + "</a></td><td>");
					}
					if (!Cutoffs.Brief) {
						if(printHTML)
						{
							Descrip = LocusList
									.PrintLinksForPeptide(DTARunner.FileName);
							Outgoing.write(Descrip);

						}

					}
					DTARunner = DTARunner.Next;
					if(printHTML) Outgoing.write("</tr>\n");

				}
				if(printHTML) Outgoing.write("</table>\n");
				// Set up Similarities field for DTASelect.html creation
				if (!Cutoffs.Brief) {
					if(printHTML)
					{
						Similarities = new StringBuffer();
						MatchRunner = LocusList.Next;
						PepCount = Runner.DTACount();
						while (MatchRunner != null) {
							ScoreBuffer = Runner.Similarity(MatchRunner);
							MatchingPeptides = (PepCount + MatchRunner.DTACount() + ScoreBuffer) / 4;
							// if ( (ScoreBuffer > -1) && (MatchRunner != Runner) )
							// {
							if ((MatchingPeptides > 0) && (MatchRunner != Runner)) {
								Similarities.append("<a href=\"#");
								Similarities.append(MatchRunner.Locus);
								Similarities.append("\">");
								Similarities.append(MatchRunner.Locus);
								Similarities.append("</a>(");
								Similarities.append(MatchingPeptides);
								Similarities.append(":");
								Similarities.append(PepCount - MatchingPeptides);
								Similarities.append(")<BR>");
							}
							MatchRunner = MatchRunner.Next;
						}
						if (Similarities.length() > 0) {
							if(printHTML) Outgoing.write("Similarities:\n");
							if(printHTML) Outgoing.write(Similarities.toString());
						}
					}

				}
				Runner = Runner.Next;
			}

			// PRINT TAIL OF DTASelect.html
			// If the user has selected database output,
			if (Cutoffs.PrintDB) {
				System.out.println("\tCreating Prot2Pep database...");
				// Create the file
				File PeptideDB = new File(CurrentDirectory, DBPrefix
						+ "-Prot2Pep.txt");
				FileWriter PDBFW = new FileWriter(PeptideDB);
				BufferedWriter PDBBW = new BufferedWriter(PDBFW);
				PDBBW.write("LocusID\tMatchID\tSequence Position\tTryptic\n");
				Runner = LocusList.Next;
				while (Runner != null) {
					MatchRunner = Runner;
					while (MatchRunner != null) {
						DTARunner = MatchRunner.DTAs.Next;
						while (DTARunner != null) {
							PDBBW.write(MatchRunner.Locus + tab
									+ DTARunner.FileName + tab
									+ (DTARunner.SequencePosition + 1) + tab
									+ DTARunner.Tryptic + newline);
							DTARunner = DTARunner.Next;
						}
						MatchRunner = MatchRunner.IdenticalLoci;
					}
					Runner = Runner.Next;
				}
				PDBBW.flush();
				PDBBW.close();
			} /*
			 * Manipulate to get nonredundant peptide counts. NOTE: THIS
			 * DESTROYS THE INTEGRITY OF THE PROTEIN LIST!!! This operation
			 * results in a single protein with all the peptides attached to it.
			 */
			Mega = LocusList.CreateMegaProtein();
			Mega.DTAs.FileSortList();
			Mega.RemoveRedundantDTAs();
			if (Cutoffs.PrintChroma) {
				System.out.println("\tCreating chroma report...");
				Mega.PrintChromaFile(FileName + "-chroma.txt");
			}
			DTARunner = Mega.DTAs.Next;
			while (DTARunner != null) {
				NRPeptideCounter++;
				NRSpectrumCounter += DTARunner.Redundancy;
				if (DTARunner.IsDecoy) {
					NRDecoyPeptideCounter++;
					NRDecoySpectrumCounter += DTARunner.Redundancy;
				} else {
					NRForwardPeptideCounter++;
					NRForwardSpectrumCounter += DTARunner.Redundancy;
				}
				DTARunner = DTARunner.Next;
			}
			DecoyCounter = rcount;
			ForwardCounter = fcount;
			if (DecoyCounter == 0) {
				ProteinFPRate = 0.0f;
			} else if (DecoyCounter >= ForwardCounter) {
				ProteinFPRate = 1.0f;
			} else {
				//ProteinFPRate = 1.0f * ((float)rcount)/((float)fcount);
				ProteinFPRate = 1.0f * DecoyCounter / ForwardCounter;
			}
			if (NRDecoyPeptideCounter == 0) {
				NRPeptideFPRate = 0.0f;
			} else if (NRDecoyPeptideCounter >= NRForwardPeptideCounter) {
				NRPeptideFPRate = 1.0f;
			} else {
				NRPeptideFPRate = 1.0f * NRDecoyPeptideCounter
						/ NRForwardPeptideCounter;
			}
			if (NRDecoySpectrumCounter == 0) {
				NRSpectrumFPRate = 0.0f;
			} else if (NRDecoySpectrumCounter >= NRForwardSpectrumCounter) {
				NRSpectrumFPRate = 1.0f;
			} else {
				NRSpectrumFPRate = 1.0f * NRDecoySpectrumCounter
						/ NRForwardSpectrumCounter;
			}
			// Print counts of loci and DTAFiles
			if(printHTML) Outgoing.write("\n<a name=\"Summary\"><table border>");
			if(printHTML) Outgoing.write("<tr><td></td><td>Proteins</td><td>Peptide IDs</td><td>Spectra</td></tr>\n");
			if(printHTML) Outgoing.write("<tr><td>Unfiltered</td><td>"
					+ new Integer(this.UnfilteredLocusCount).toString()
					+ CellSep + new Integer(this.UnfilteredOUTCount).toString()
					+ CellSep
					+ new Integer(this.UnfilteredSpectrumCount).toString()
					+ "</td></tr>\n");
			/*
			 * if(printHTML) Outgoing.write("<tr><td>Redundant</td><td>" + new
			 * Integer(RProteinCounter).toString() + CellSep + new
			 * Integer(PeptideCounter).toString() + CellSep + new
			 * Integer(SpectrumCounter).toString() + "</td></tr>\n");
			 */
			if(printHTML) Outgoing.write("<tr><td>Filtered</td><td>"
					+ new Integer(Counter).toString() + CellSep
					+ new Integer(NRPeptideCounter).toString() + CellSep
					+ new Integer(NRSpectrumCounter).toString()
					+ "</td></tr>\n");
			if (!Cutoffs.HideDecoy) {
				// added by Howard Choi
				System.out.println("this is decoy counter  "
						+ new Integer(DecoyCounter).toString());
				System.out.println("this is Forward matches  "
						+ new Integer(ForwardCounter).toString());
				float decoyVal = DecoyCounter;
				float forwardVal = ForwardCounter;
				setProtFPOutput(Float.toString(Rounding(
						((decoyVal / forwardVal)), 4)));

				// added by Howard Choi
				System.out.println("this is Spec decoy counter  "
						+ new Integer(NRDecoyPeptideCounter).toString());
				System.out.println("this is Spec Forward matches  "
						+ new Integer(NRForwardPeptideCounter).toString());
				float nrDecoyVal = NRDecoyPeptideCounter;
				float nrForwardVal = NRForwardPeptideCounter;
				setPepFPOutput(Float.toString(Rounding(
						((nrDecoyVal / nrForwardVal)), 4)));

				if(printHTML) Outgoing.write("<tr><td>Forward matches</td><td>"
						+ new Integer(ForwardCounter).toString() + CellSep
						+ new Integer(NRForwardPeptideCounter).toString()
						+ CellSep
						+ new Integer(NRForwardSpectrumCounter).toString()
						+ "</td></tr>\n");
				if(printHTML) Outgoing.write("<tr><td>Decoy matches</td><td>"
						+ new Integer(DecoyCounter).toString() + CellSep
						+ new Integer(NRDecoyPeptideCounter).toString()
						+ CellSep
						+ new Integer(NRDecoySpectrumCounter).toString()
						+ "</td></tr>\n");
				if(printHTML) Outgoing.write("<tr><td>Forward FDR</td><td>"
						+ new Float(Protein.RoundTo(100f * ProteinFPRate, 2))
								.toString()
						+ "%"
						+ CellSep
						+ new Float(Protein.RoundTo(100f * NRPeptideFPRate, 2))
								.toString()
						+ "%"
						+ CellSep
						+ new Float(Protein.RoundTo(100f * NRSpectrumFPRate, 2))
								.toString() + "%" + "</td></tr>\n");
			}
			if(printHTML) Outgoing.write("</table></a><BR>\n");
			if (Cutoffs.UseClassifications) {
				if(printHTML) Outgoing.write("\n<table border>");
				if(printHTML) Outgoing.write("<tr><td>Classification</td><td>Nonredundant Proteins</td><td>Redundant Proteins</td></tr>\n");
				CurrentClassification = Classifieds.Next;
				while (CurrentClassification != null) {
					if(printHTML) Outgoing.write("<tr><td><a href=\"#Class"
							+ CurrentClassification.Identifier + "\">"
							+ CurrentClassification.Descriptor
							+ "</a></td><td>"
							+ CurrentClassification.NRCountObserved
							+ "</td><td>"
							+ CurrentClassification.RCountObserved
							+ "</td></tr>");
					CurrentClassification = CurrentClassification.Next;
				}
				if(printHTML) Outgoing.write("<tr><td><a href=\"#UnClass\">Unclassified</a></td><td>"
						+ UnclassifiedNRCount
						+ "</td><td>"
						+ UnclassifiedRCount + "</td></tr>");
				if(printHTML) Outgoing.write("</table>\n");
			}
			if(printHTML) Outgoing.write((this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
					: this.SourceDirectory)
					+ newline);
			if(printHTML) Outgoing.write("</body></html>");
			if(printHTML )Outgoing.flush();
			if(printHTML) Outgoing.close();
			// Print counts of loci and DTAFiles
			OutgoingTXT.write("\n\tProteins\tPeptide IDs\tSpectra\n");
			OutgoingTXT.write("Unfiltered\t"
					+ new Integer(this.UnfilteredLocusCount).toString() + tab
					+ new Integer(this.UnfilteredOUTCount).toString() + tab
					+ new Integer(this.UnfilteredSpectrumCount).toString()
					+ newline);
			/*
			 * OutgoingTXT.write("Redundant\t" + new
			 * Integer(RProteinCounter).toString() + tab + new
			 * Integer(PeptideCounter).toString() + tab + new
			 * Integer(SpectrumCounter).toString() + newline);
			 */
			OutgoingTXT.write("Filtered\t" + new Integer(Counter).toString()
					+ tab + new Integer(NRPeptideCounter).toString() + tab
					+ new Integer(NRSpectrumCounter).toString() + "\n");
			if (!Cutoffs.HideDecoy) {
				OutgoingTXT.write("Forward matches\t"
						+ new Integer(ForwardCounter).toString() + tab
						+ new Integer(NRForwardPeptideCounter).toString() + tab
						+ new Integer(NRForwardSpectrumCounter).toString()
						+ "\n");
				
				OutgoingTXT.write("Redundant Forward matches\t"
						+ new Integer(RForwardCounter).toString() + tab
						+ new Integer(NRForwardPeptideCounter).toString() + tab
						+ new Integer(NRForwardSpectrumCounter).toString()
						+ "\n");
				
				OutgoingTXT
						.write("Decoy matches\t"
								+ new Integer(DecoyCounter).toString()
								+ tab
								+ new Integer(NRDecoyPeptideCounter).toString()
								+ tab
								+ new Integer(NRDecoySpectrumCounter)
										.toString() + "\n");
				
				OutgoingTXT
				.write("Redundant Decoy matches\t"
						+ new Integer(RDecoyCounter).toString()
						+ tab
						+ new Integer(NRDecoyPeptideCounter).toString()
						+ tab
						+ new Integer(NRDecoySpectrumCounter)
								.toString() + "\n");
				
				OutgoingTXT.write("Forward FDR\t"
						+ Protein.RoundTo(100f * ProteinFPRate, 2) + tab
						+ Protein.RoundTo(100f * NRPeptideFPRate, 2) + tab
						+ Protein.RoundTo(100f * NRSpectrumFPRate, 2) + "\n");
			}
			OutgoingTXT
					.write("\nClassification\tNonredundant Proteins\tRedundant Proteins\n");
			CurrentClassification = Classifieds.Next;
			while (CurrentClassification != null) {
				OutgoingTXT.write(CurrentClassification.Descriptor + tab
						+ CurrentClassification.NRCountObserved + tab
						+ CurrentClassification.RCountObserved + newline);
				CurrentClassification = CurrentClassification.Next;
			}
			OutgoingTXT.write("Unclassified\t" + UnclassifiedNRCount + tab
					+ UnclassifiedRCount + newline);
			OutgoingTXT.flush();
			OutgoingTXT.close();
			// If the user has selected database output,
			if (Cutoffs.PrintDB) {
				Mega.PrintDBPeptides(CurrentDirectory, DBPrefix, hXCorr,
						hCalcPreMass, hDeltCN, hSpScore, Cutoffs);
			}
			/*
			 * Handle --copy option for unified files
			 */
			if (Cutoffs.PrintCopyList && this.IDFileFormat.equals("SQT")) {
				Mega.MakeSQTSubsets(CurrentDirectory);
				Mega.MakeMS2Subsets(CurrentDirectory);
			}
			if (MSBug) {
				System.out
						.println("*** An error occured because of one of these possibilities:");
				System.out
						.println("*** 1) your SEQUEST results have been moved to a new directory");
				System.out.println("*** 2) you need to use the --CGI option.");
				System.out
						.println("*** Some links in the DTASelect.html file are likely to be broken.");
			}
		} catch (IOException failure) {
			System.out
					.println("\tSomething went wrong while writing output files");
			System.out.println(failure);
		}
	}

	/*
	 * Print an HTML web page which reports all sequence alignments within each
	 * protein.
	 */
	public void PrintAlignReport(IniFile Config, String RootName) {
		Protein PRunner = this.LocusList.Next;
		Protein QRunner;
		DTAFile DTARunner;
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		File TOutputFile;
		FileWriter TOutputFileWriter;
		BufferedWriter TOutgoing;
		int LeftMostResidue;
		int RightMostResidue;
		int LastClass = 127;
		Classification CurrentClassification = null;
		String SeqString = "oops.";
		String ModString;
		String NewLine = "\n";
		String Tab = "\t";
		StringTokenizer Parser;
		int SeqLength;
		int Looper;
		int ToAdd;
		int PeptideLengthSum;
		int MaxDepth;
		int CurrentDepth;
		// Because we want to write the endpoint of each contig, we
		// store sequence alignment until each row is done.
		StringBuffer HTMLAlignment;
		boolean NonFirstRow;
		try {
			// Open the output files
			OutputFile = new File(CurrentDirectory, RootName + "-align.html");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			TOutputFile = new File(CurrentDirectory, RootName + "-align.txt");
			TOutputFileWriter = new FileWriter(TOutputFile);
			TOutgoing = new BufferedWriter(TOutputFileWriter);
			// Write html header
			Outgoing.write("<html><head><title>A ");
			if (RootName.equals("DTASelect")) {
				Parser = new StringTokenizer(System.getProperty("user.dir"),
						System.getProperty("file.separator"));
				while (Parser.hasMoreTokens())
					SeqString = Parser.nextToken();
				Outgoing.write(SeqString);
			} else {
				Outgoing.write(RootName);
			}
			Outgoing.write("</title></head>");
			Outgoing.write("<body>\n");
			Outgoing.write("Columns include:<ul>\n");
			Outgoing.write("<LI>Start and stop positions of sequence contig\n");
			Outgoing.write("<LI>Sequence alignment\n");
			Outgoing.write("</ul>\n");
			TOutgoing.write("L\tLocusID\tDescription\n");
			TOutgoing.write("P\tPeptideID\tSequence Assembly\tModified\n");
			TOutgoing.write("C\tLocusID\tStart Position\tEnd Position\n");
			TOutgoing.write("S\tLocusID\tMaximum Depth\tResidues Observed\t"
					+ "Sequence Coverage\tSequence Length\n");
			// Loop through all proteins in list
			while (PRunner != null) {
				if (PRunner.Classification != LastClass) {
					LastClass = PRunner.Classification;
					CurrentClassification = Classifieds
							.GetClass(PRunner.Classification);
					if (CurrentClassification != null) {
						Outgoing.write("<h2>"
								+ CurrentClassification.Descriptor + "</h2>\n");
						TOutgoing.write("\nB\t"
								+ CurrentClassification.Descriptor + "\n");
					} else {
						Outgoing.write("<h2>Uncategorized Loci</h2>\n");
						TOutgoing.write("\nB\tUncategorized Loci\n");
					}
				}
				PeptideLengthSum = 0;
				MaxDepth = 0;
				CurrentDepth = 0;
				// Create a table for this protein, with the locus
				// name at the top.
				Outgoing.write("<table border><tr><td colspan=1>\n");
				Outgoing.write("<a target=\"" + RootName + "\" href=\""
						+ RootName + ".html#" + PRunner.Locus + "\">"
						+ PRunner.Locus + "</a>");
				Outgoing.write("</td><td colspan=4>" + PRunner.Gene);
				TOutgoing.write("\nL\t" + PRunner.Locus + Tab + PRunner.Gene
						+ NewLine);
				QRunner = PRunner.IdenticalLoci;
				while (QRunner != null) {
					Outgoing.write("<tr><td><a target=\"" + RootName
							+ "\" href=\"" + RootName + ".html#"
							+ QRunner.Locus + "\">" + QRunner.Locus
							+ "</a></td></tr>");
					TOutgoing.write("\nL\t" + QRunner.Locus + NewLine);
					QRunner = QRunner.Next;
				}
				NonFirstRow = false;
				// Loop through all peptides remaining for this protein
				DTARunner = PRunner.DTAs.Next;
				RightMostResidue = -2;
				LeftMostResidue = -2;
				HTMLAlignment = new StringBuffer("\n");
				while (DTARunner != null) {
					SeqString = DTARunner.TrimmedSequence();
					SeqLength = SeqString.length();
					PeptideLengthSum += SeqLength;
					ModString = DTARunner.SymbolModString();
					// Determine if the current peptide fits the current
					// alignment
					if (DTARunner.SequencePosition > RightMostResidue) {
						// We must start a new table row; this peptide isn't in
						// the last contig.
						if (NonFirstRow) {
							HTMLAlignment.append("</pre>");
							Outgoing.write(new Integer(RightMostResidue)
									.toString());
							TOutgoing.write("C\t"
									+ PRunner.Locus
									+ Tab
									+ new Integer(LeftMostResidue + 1)
											.toString() + Tab
									+ new Integer(RightMostResidue).toString()
									+ NewLine);
						} else
							NonFirstRow = true;
						Outgoing.write(HTMLAlignment.toString());
						Outgoing.write("</td></tr><tr><td>"
								+ new Integer(DTARunner.SequencePosition + 1)
										.toString() + "<BR>");
						if (CurrentDepth > MaxDepth)
							MaxDepth = CurrentDepth;
						CurrentDepth = 0;
						HTMLAlignment = new StringBuffer("</td><td><pre>\n");
						LeftMostResidue = DTARunner.SequencePosition;
						RightMostResidue = LeftMostResidue + SeqLength;
					} else {
						// See if this peptide extends further toward
						// the C-terminus than any peptide so far.
						if (DTARunner.SequencePosition + SeqLength > RightMostResidue)
							RightMostResidue = DTARunner.SequencePosition
									+ SeqLength;
					}
					// Write the row title (leftmost residue position) to the
					// text file
					TOutgoing.write("P\t" + DTARunner.FileName + Tab);
					// Write the appropriate number of spaces to each of the two
					// files
					ToAdd = DTARunner.SequencePosition - LeftMostResidue;
					for (Looper = 0; Looper < ToAdd; Looper++) {
						HTMLAlignment.append(' ');
						TOutgoing.write(' ');
					}
					CurrentDepth++;
					// TOutgoing.write(SeqString);
					// Modified peptides have to be written
					// differently; modified residues are colored,
					// hopefully in the same colors as used by
					// Display_ions.
					if (DTARunner.Modified) {
						for (Looper = 0; Looper < SeqLength; Looper++) {
							switch (ModString.charAt(Looper)) {
							case ' ':
								HTMLAlignment.append(SeqString.charAt(Looper));
								TOutgoing.write(SeqString.charAt(Looper));
								break;
							default:
								HTMLAlignment.append("<font color=\"red\">");
								HTMLAlignment.append(SeqString.charAt(Looper));
								HTMLAlignment.append("</font>");
								TOutgoing.write(new Character(SeqString
										.charAt(Looper)).toString()
										.toLowerCase());
								break;
							}
						}
						HTMLAlignment.append(NewLine);
						TOutgoing.write("\tY");
					} else {
						TOutgoing.write(SeqString + "\tN");
						HTMLAlignment.append(SeqString + NewLine);
					}
					TOutgoing.write(NewLine);
					DTARunner = DTARunner.Next;
				}
				// Now we've finished the protein. Write the end
				// of the table.
				if (CurrentDepth > MaxDepth)
					MaxDepth = CurrentDepth;
				TOutgoing.write("C\t" + PRunner.Locus + Tab
						+ new Integer(LeftMostResidue + 1).toString() + Tab
						+ new Integer(RightMostResidue).toString() + NewLine);
				TOutgoing
						.write("S\t"
								+ PRunner.Locus
								+ Tab
								+ new Integer(MaxDepth).toString()
								+ Tab
								+ new Integer(PeptideLengthSum).toString()
								+ Tab
								+ PRunner.SequenceCoverage
								+ "%\t"
								+ new Integer(PRunner.SequenceLength)
										.toString() + "\n");
				HTMLAlignment
						.append("</pre></td></tr></table>\n<table border><tr><td>\n"
								+ new Integer(MaxDepth).toString()
								+ "</td><td>Maximum Depth\n</td></tr><tr><td>"
								+ new Integer(PeptideLengthSum).toString()
								+ "</td><td>Peptide Residues Observed\n</td></tr><tr><td>"
								+ PRunner.SequenceCoverage
								+ "%</td><td>Sequence Coverage\n</td></tr><tr><td>"
								+ new Integer(PRunner.SequenceLength)
										.toString()
								+ "</td><td>Sequence Length\n</td></tr></table><P>\n");
				Outgoing.write(new Integer(RightMostResidue).toString());
				Outgoing.write(HTMLAlignment.toString());
				PRunner = PRunner.Next;
			}
			// Close off the HTML file. We're done.
			Outgoing.write("</body></html>");
			Outgoing.flush();
			Outgoing.close();
			TOutgoing.flush();
			TOutgoing.close();
		} catch (IOException failure) {
			System.out.println("Error while printing alignment report:");
			System.out.println(failure);
		}
	}

	/*
	 * Print an HTML web page which reports all modifications found for all
	 * modified proteins.
	 */
	public void PrintModReport(IniFile Config, String RootName, String SourceDir) {
		Protein PRunner = this.LocusList.Next;
		Protein QRunner;
		Modification ModList;
		Modification ModRunner;
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		File TOutputFile;
		FileWriter TOutputFileWriter;
		BufferedWriter TOutgoing;
		char LastModType;
		char LastResidue;
		String MassShift = "Error";
		StringTokenizer Parser;
		int LastPosition;
		int LeftMostResidue;
		int Looper;
		int ToAdd;
		int LastClass = 127;
		Classification CurrentClassification = null;
		StringBuffer SequenceStack;
		StringBuffer ModLine;
		StringBuffer XCorrStack;
		StringBuffer FileNameStack;
		File CurrentSubdir;
		String CanonicalName = "DTASelect";
		StringBuffer ModCountStack;
		StringBuffer CopyCountStack;
		StringBuffer StaticModCGI = new StringBuffer();
		String StaticModString;
		String DisplayIons = "<A TARGET=\"Win1\" HREF=\""
				+ Config.DTASpecDisplay + "?Dta=";
		String NewShow = "<a target=\"Win1\" HREF=\""
				+ Config.SQTDisplay
				+ "?Dr="
				+ (this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
						: this.SourceDirectory) + "&Da=";
		String ShowOut = "<a TARGET=\"Win1\" href=\"" + Config.DTAOUTDisplay
				+ "?OutFile=";
		StaticModCGI.append(SequestParams.CGIString());
		StaticModString = StaticModCGI.toString() + "&Pep=";
		try {
			OutputFile = new File(CurrentDirectory, RootName + "-mods.html");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			TOutputFile = new File(CurrentDirectory, RootName + "-mods.txt");
			TOutputFileWriter = new FileWriter(TOutputFile);
			TOutgoing = new BufferedWriter(TOutputFileWriter);
			// Write html header
			Outgoing.write("<html><head><title>M ");
			if (RootName.equals("DTASelect")) {
				Parser = new StringTokenizer(System.getProperty("user.dir"),
						System.getProperty("file.separator"));
				while (Parser.hasMoreTokens())
					CanonicalName = Parser.nextToken();
				Outgoing.write(CanonicalName);
			} else {
				Outgoing.write(RootName);
			}
			Outgoing.write("</title></head>");
			Outgoing.write("<body>\n");
			Outgoing.write("Columns include:<ul>\n");
			Outgoing.write("<LI>Number of modified residue and modification mass difference\n");
			Outgoing.write("<LI>Sequence alignment of peptides, with symbol indicating modified residue\n");
			Outgoing.write("<LI>Spectrum link from spectrum filename\n");
			Outgoing.write("<LI>SEQUEST output file link from " + this.hXCorr
					+ "\n");
			Outgoing.write("<LI>Number of modifications found in this identification\n");
			Outgoing.write("<LI>Number of copies observed for this identification\n</ul>\n");
			TOutgoing.write("L\tLocus\tDescription\n");
			TOutgoing
					.write("M\tResidue Modified\tMass Difference\tSequence Alignment\tFilename\tXCorr\tModCount\tCopy Count\n");
			TOutgoing
					.write("N\tResidue Modified\tMass Difference\tModification Position\tResidue Modified\n\n");
			while (PRunner != null) {
				if (PRunner.Classification != LastClass) {
					LastClass = PRunner.Classification;
					CurrentClassification = Classifieds
							.GetClass(PRunner.Classification);
					if (CurrentClassification != null) {
						Outgoing.write("<h2>"
								+ CurrentClassification.Descriptor + "</h2>\n");
						TOutgoing.write("C\t"
								+ CurrentClassification.Descriptor + "\n");
					} else {
						Outgoing.write("<h2>Uncategorized Loci</h2>\n");
						TOutgoing.write("C\tUncategorized Loci\n");
					}
				}
				ModList = this.GiveModList(PRunner);
				if (ModList != null) {
					// Create a table for this protein, with the locus
					// name at the top.
					Outgoing.write("<table border>\n");
					Outgoing.write("<tr><td colspan=2>\n<a target=\""
							+ RootName + "\" href=\"" + RootName + ".html#"
							+ PRunner.Locus + "\">" + PRunner.Locus
							+ "</a></td><td colspan=4>" + PRunner.Gene);
					TOutgoing.write("L\t" + PRunner.Locus + "\t" + PRunner.Gene
							+ "\n");
					QRunner = PRunner.IdenticalLoci;
					while (QRunner != null) {
						Outgoing.write("<tr><td colspan=2>\n<a target=\""
								+ RootName + "\" href=\"" + RootName + ".html#"
								+ QRunner.Locus + "\">" + QRunner.Locus
								+ "</a>");
						TOutgoing.write("L\t" + QRunner.Locus + "\n");
						QRunner = QRunner.Next;
					}
					// Sort the modifications by 1) residue modified,
					// 2) modification type, 3) peptide position, 4)
					// precursor mass
					ModList.SortList();
					ModRunner = ModList.Next;
					while (ModRunner != null) {
						// Write the residue at which this
						// modification takes place.
						Outgoing.write("</td></tr>\n<tr><td>");
						Outgoing.write(new Integer(
								ModRunner.ModPositionInProt + 1).toString()
								+ "<BR>");
						// Write the mass difference found for this
						// modification type
						MassShift = new Float(
								SequestParams.DiffMods
										.getMassShiftFor(ModRunner.ModType))
								.toString();
						Outgoing.write(MassShift + "</td><td><pre>");
						// Blank out each column; we're going to fill
						// in the values for this protein.
						SequenceStack = new StringBuffer();
						XCorrStack = new StringBuffer();
						FileNameStack = new StringBuffer();
						ModCountStack = new StringBuffer();
						CopyCountStack = new StringBuffer();
						// Set our markers to the current
						// modification's settings.
						LastModType = ModRunner.ModType;
						LastPosition = ModRunner.ModPositionInProt;
						LeftMostResidue = ModRunner.Peptide.SequencePosition;
						LastResidue = ModRunner.Peptide.TrimmedSequence()
								.charAt(LastPosition - LeftMostResidue);
						// While we're on a batch of modifications at
						// the same position and type...
						while (ModRunner != null
								&& (ModRunner.ModPositionInProt == LastPosition)
								&& (ModRunner.ModType == LastModType)) {
							TOutgoing.write("M\t"
									+ new Integer(LastPosition + 1).toString()
									+ "\t" + MassShift + "\t");
							// Store the data needed for each column
							// in our StringBuffer objects.
							// First, align the sequence against the others.
							ToAdd = ModRunner.Peptide.SequencePosition
									- LeftMostResidue;
							for (Looper = 0; Looper < ToAdd; Looper++) {
								SequenceStack.append(' ');
								TOutgoing.write(' ');
							}
							TOutgoing.write(ModRunner.Peptide.TrimmedSequence()
									+ "\t");
							SequenceStack.append(ModRunner.Peptide
									.TrimmedSequence());
							SequenceStack.append("\n");
							// Next, write the XCorr, linking it to
							// the out file viewer
							CurrentSubdir = new File(SourceDir,
									ModRunner.Peptide.Subdirectory);
							CanonicalName = new File(CurrentSubdir,
									ModRunner.Peptide.FileName)
									.getCanonicalPath();
							// We need a separate link for OUTs only if we're
							// using that format
							if (Cutoffs.DisplayType == 0) {
								XCorrStack.append(ShowOut + CanonicalName
										+ ".out\">");
							}
							TOutgoing.write(ModRunner.Peptide.XCorr + "\t");
							XCorrStack.append(ModRunner.Peptide.XCorr);
							if (Cutoffs.DisplayType == 0) {
								XCorrStack.append("</a>\n");
							} else
								XCorrStack.append("\n");
							// Next, write the file name, linking it
							// to the display_ions viewer.
							if (Cutoffs.DisplayType > 0) {
								// Note that we'll need to add code
								// for the Mascot CGIs when this
								// becomes feasible. DisplayType == 2
								// when Mascot results are
								// incorporated.
								FileNameStack
										.append(NewShow
												+ ModRunner.Peptide
														.ShowString() + "\n");
							} else {
								FileNameStack
										.append(DisplayIons
												+ CanonicalName
												+ ".dta"
												+ ModRunner.Peptide
														.DisplayIonsString(SequestParams)
												+ StaticModString
												+ ModRunner.Peptide
														.TrimmedSequence()
												+ "\">");
								FileNameStack
										.append(ModRunner.Peptide.FileName);
								FileNameStack.append("</a>\n");
							}
							TOutgoing.write(ModRunner.Peptide.FileName + "\t");
							// Penultimately, show how many modifications
							// were present in this identification
							TOutgoing.write(ModRunner.ModCount + "\t");
							ModCountStack.append(ModRunner.ModCount);
							ModCountStack.append("\n");
							// Finally, show how many copies of this
							// spectrum were present
							TOutgoing
									.write(ModRunner.Peptide.Redundancy + "\n");
							CopyCountStack.append(ModRunner.Peptide.Redundancy);
							CopyCountStack.append("\n");
							ModRunner = ModRunner.Next;
						}
						// Add a marker showing the location of the
						// modification in these sequences.
						ToAdd = LastPosition - LeftMostResidue;
						TOutgoing.write("N\t"
								+ new Integer(LastPosition + 1).toString()
								+ "\t" + MassShift + "\t");
						ModLine = new StringBuffer();
						for (Looper = 0; Looper < ToAdd; Looper++) {
							SequenceStack.append(' ');
							ModLine.append(' ');
							TOutgoing.write(' ');
						}
						SequenceStack.append(LastModType);
						ModLine.append(LastModType);
						TOutgoing.write(LastModType);
						TOutgoing.write('\t');
						TOutgoing.write(LastResidue);
						TOutgoing.write('\n');
						ModLine.append('\n');
						// Make sure that browsers render the other
						// columns to be as tall as the SequenceStack
						// column. The space adds a last row, oddly
						// enough.
						XCorrStack.append(' ');
						FileNameStack.append(' ');
						ModCountStack.append(' ');
						CopyCountStack.append(' ');
						Outgoing.write(ModLine.toString()
								+ SequenceStack.toString()
								+ "</pre></td><td><pre> \n"
								+ FileNameStack.toString()
								+ "</pre></td><td><pre> \n"
								+ XCorrStack.toString()
								+ "</pre></td><td><pre> \n"
								+ ModCountStack.toString()
								+ "</pre></td><td><pre> \n"
								+ CopyCountStack.toString() + "</pre>");
					}
					// Now we've finished the protein. Write the end
					// of the table.
					Outgoing.write("</td></tr></table><P>\n");
				}
				// Next protein!
				PRunner = PRunner.Next;
			}
			// Close off the HTML file. We're done.
			Outgoing.write("</body></html>");
			Outgoing.flush();
			Outgoing.close();
			TOutgoing.flush();
			TOutgoing.close();
		} catch (IOException failure) {
			System.out.println("Error while printing modification report:");
			System.out.println(failure);
		}
	}

	/*
	 * Create a list of modification objects that correspond to this protein's
	 * modifications. Each nonzero position in each identication's modification
	 * string should result in a modification object.
	 */
	public Modification GiveModList(Protein TargetProtein) {
		DTAFile Runner = TargetProtein.DTAs.Next;
		Modification ModList = new Modification();
		Modification MRunner = ModList;
		String ModString;
		int SpectrumModCount = 0;
		int SeqLength;
		int Looper;
		while (Runner != null) {
			// For each spectrum in this protein's list,
			// Is there a modification present?
			ModString = Runner.SymbolModString();
			SeqLength = ModString.length();
			// Count how many modifications are present here.
			SpectrumModCount = 0;
			for (Looper = 0; Looper < SeqLength; Looper++) {
				if (ModString.charAt(Looper) != ' ')
					SpectrumModCount++;
			}
			if (SpectrumModCount > 0) {
				// Now create modification objects
				for (Looper = 0; Looper < SeqLength; Looper++) {
					if (ModString.charAt(Looper) != ' ') {
						MRunner.Next = new Modification();
						MRunner = MRunner.Next;
						MRunner.ModPositionInProt = Runner.SequencePosition
								+ Looper;
						MRunner.ModType = ModString.charAt(Looper);
						MRunner.Peptide = Runner;
						MRunner.ModCount = new Integer(SpectrumModCount)
								.toString();
					}
				}
			}
			Runner = Runner.Next;
		}
		if (ModList == MRunner)
			return null;
		else
			return ModList;
	}

	/*
	 * Print a table showing similarity between protein groups.
	 */
	public void PrintProteinSimilarity(String RootName) {
		Protein Runner = this.LocusList.Next;
		Protein InsideRunner;
		int Counter = 1;
		int Looper;
		String StringBuffer = System.getProperty("user.dir");
		File CurrentDirectory = new File(StringBuffer);
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		if (Runner != null) {
			Runner = Runner.Next;
		}
		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			OutputFile = new File(CurrentDirectory, RootName + "-pst.txt");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			// Write header
			Outgoing.write(Protein.Version());
			Outgoing.write("\nName\tSelf\t");
			while (Runner != null) {
				Outgoing.write(Runner.Locus + "\t");
				Runner = Runner.Next;
			}
			Outgoing.write("\n");
			// Write similarity information for each protein in list
			Runner = this.LocusList.Next;
			while (Runner != null) {
				// Write locus name
				Outgoing.write(Runner.Locus + "\t");
				// Write self-similarity in one column
				Outgoing.write(new Integer(Runner.Similarity(Runner))
						.toString());
				// Tab out to the proper location
				for (Looper = 0; Looper < Counter; Looper++)
					Outgoing.write("\t");
				InsideRunner = Runner.Next;
				// For each comparison to another protein
				while (InsideRunner != null) {
					Looper = Runner.Similarity(InsideRunner);
					// Write similarity only if it's significant
					if (Looper > -1)
						Outgoing.write(new Integer(Looper).toString());
					Outgoing.write("\t");
					InsideRunner = InsideRunner.Next;
				}
				Outgoing.write("\n");
				Counter++;
				Runner = Runner.Next;
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Error while writing DTASelect-pst.txt");
			System.out.println(failure);
		}
	}

	/*
	 * Create the DTASelect.txt database from all the loci and DTAFiles read
	 * from .out files
	 */
	public void PrintTXT() {
		String StringBuffer = System.getProperty("user.dir");
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		Protein Runner = this.LocusList.Next;
		DTAFile DTARunner;
		String NewLine = "\n";
		String Tab = "\t";
		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			OutputFile = new File(CurrentDirectory, "DTASelect.txt");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			// Write the header to the file.
			Outgoing.write(Protein.Version() + NewLine);
			Outgoing.write((this.Cutoffs.UseCustomPath ? this.Cutoffs.CustomPath
					: CurrentDirectory.getCanonicalPath())
					+ NewLine);
			Outgoing.write(NewLine);
			// Write the fields of the ParamsFile to the next few lines
//			Outgoing.write(this.SequestParams.getDTASelectTxtString());
			Outgoing.write("Source\t" + this.SQTGenerator + Tab
					+ this.SQTGeneratorVersion + NewLine);
			Outgoing.write("Format\t" + this.IDFileFormat + NewLine);
			// Print off all fields, continuing down list
			while (Runner != null) {
				this.UnfilteredLocusCount++;
				if (Runner.Gene.length() == 0)
					Runner.Gene = "no description";
				Outgoing.write(Runner.GetLLine());
				DTARunner = Runner.DTAs.Next;
				while (DTARunner != null) {
					this.UnfilteredSpectrumCount++;
					Outgoing.write(DTARunner.GetDLine());
					DTARunner = DTARunner.Next;
				}
				Runner = Runner.Next;
			}
			Outgoing.write("C\t"
					+ new Integer(this.UnfilteredLocusCount).toString() + Tab
					+ new Integer(this.UnfilteredOUTCount).toString() + Tab
					+ new Integer(this.UnfilteredSpectrumCount).toString()
					+ NewLine);
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing TXT file");
		}
	}

	/*
	 * If the DTASelect.txt file has already been written, don't bother with the
	 * .out files. Just read the big flat-file database that already exists.
	 */
	public void ReadFromFile(File DataFile) throws IOException {
		// Read in the DTASelect.txt file
		Protein ListToReturn = this.LocusList;
		Protein ListRunner = ListToReturn;
		Protein ProtBuffer = new Protein();
		DTAFile DTARunner = ProtBuffer.DTAs;
		DTAFile TempDTA;
		FileReader InputFilereader = new FileReader(DataFile);
		BufferedReader Incoming = new BufferedReader(InputFilereader);
		String LineBuffer;
		String WholeLine;
		String VersionString = Incoming.readLine();
		boolean IsCurrentVersion = (VersionString.equals(Protein.Version()));
		boolean BuildingLegitProtein = false;
		StringTokenizer Parser;
		float ModMass;
		char ModSymbol;
		String ModResidues;
		float ModPreLoss;
		float ModFragLoss;
		int count =0;
		// Read the path from whence these outs came and command line options
		this.SourceDirectory = Incoming.readLine();
		// Read the command-line options at the first run of DTASelect on this
		// set
		this.CutoffsString = Incoming.readLine();
		WholeLine = Incoming.readLine();
		// Until we get to the end of the DTASelect.txt file...
		while (WholeLine != null) {
			Parser = new StringTokenizer(WholeLine, "\t");
			// if this isn't a blank line...
			if (Parser.hasMoreTokens()) {
				LineBuffer = Parser.nextToken();
				// if this is a line describing a locus...
				if (LineBuffer.equals("L")) {
					try {

						this.UnfilteredLocusCount++;
						// add previous protein to the list (if it's got any
						// peptides)

						/*if(ProtBuffer.Locus!=null && ProtBuffer.Locus.equals("sp|Q06830|PRDX1_HUMAN"))
						{
							System.out.println(">>> "+count);
						}*/
						if (BuildingLegitProtein
								&& ProtBuffer.DTAs.Next != null) {
							ListRunner.Next = ProtBuffer;
							ListRunner = ProtBuffer;
							ProtBuffer = new Protein();
							DTARunner = ProtBuffer.DTAs;
							BuildingLegitProtein = false;
						}
						count=0;
						// Use the Protein object's DTASelect.txt line parser.
						ProtBuffer.SetTo(Parser);
					/*	if(ProtBuffer.Locus!=null && ProtBuffer.Locus.equals("sp|Q06830|PRDX1_HUMAN"))
						{
							System.out.println("... "+count);
						}*/
						if (Cutoffs.Permit(ProtBuffer))
							BuildingLegitProtein = true;
					} catch (NoSuchElementException failure) {
						System.out
								.println("This line in the DTASelect.txt is missing some information:");
						System.out.println(WholeLine);
						System.exit(0);
					}
				}
				// if this is a line describing a DTAFile...
				else if (BuildingLegitProtein && LineBuffer.equals("D")) {
					try {
						this.UnfilteredSpectrumCount++;
						// create a temporary object in which it can sit
						TempDTA = new DTAFile();
						count++;
						// For future versions: depending on the version, call
						// different methods. For now, call CurrentSetTo only.
						if (IsCurrentVersion) {
							TempDTA.CurrentSetTo(Parser, Cutoffs);
						} else {
							TempDTA.CurrentSetTo(Parser, Cutoffs);
						}
						// Calculate emergent fields
						TempDTA.DetermineModified();
						String[] ParsedLine = TempDTA.FileName.split("\\.");
						TempDTA.ChargeState = new Byte(ParsedLine[3])
								.byteValue();
						/*
						 * If this one qualifies, store it in the list of the
						 * current locus.
						 */
						if ((!Cutoffs.UseCriteria) || (Cutoffs.Allow(TempDTA, true))) {
						/*	if(ProtBuffer.Locus.equals("sp|Q06830|PRDX1_HUMAN"))
							{
								System.out.println(WholeLine);
							}*/
							DTARunner.Next = TempDTA;
							DTARunner = TempDTA;
						}
					} catch (NoSuchElementException failure) {
						System.out
								.println("This line in the DTASelect.txt is missing some information:");
						System.out.println(WholeLine);
						System.exit(0);
					}
				} else if (LineBuffer.equals("S")) {
					// This line describes information from sequest.params
					this.SequestParams = new ParamsFile(Parser, VersionString);
				} else if (LineBuffer.equals("DM")) {
					// This line describes a differential modification
					ModMass = new Float(Parser.nextToken()).floatValue();
					ModSymbol = Parser.nextToken().charAt(0);
					ModResidues = Parser.nextToken();
					ModPreLoss = new Float(Parser.nextToken()).floatValue();
					ModFragLoss = new Float(Parser.nextToken()).floatValue();
					this.SequestParams.addDiffMod(ModMass, ModSymbol,
							ModResidues, ModPreLoss, ModFragLoss);
				} else if (LineBuffer.equals("SM")) {
					// This line describes a static modification
					ModMass = new Float(Parser.nextToken()).floatValue();
					ModSymbol = Parser.nextToken().charAt(0);
					this.SequestParams.addStaticMod(ModSymbol, ModMass);
				} else if (LineBuffer.equals("Source")) {
					this.SQTGenerator = Parser.nextToken();
					this.SQTGeneratorVersion = Parser.nextToken();
				} else if (LineBuffer.equals("Format")) {
					this.IDFileFormat = Parser.nextToken();
				} else if (LineBuffer.equals("C")) {
					// This line contains counts describing this dataset
					this.UnfilteredLocusCount = new Integer(Parser.nextToken())
							.intValue();
					this.UnfilteredOUTCount = new Integer(Parser.nextToken())
							.intValue();
					this.UnfilteredSpectrumCount = new Integer(
							Parser.nextToken()).intValue();
				}
			}
			// Move to the next line of the DTASelect.txt database
			WholeLine = Incoming.readLine();
		}
		if (ProtBuffer.DTAs.Next != null) {
			ListRunner.Next = ProtBuffer;
			ListRunner = ProtBuffer;
		}
		if(!Cutoffs.noDB)this.SequestParams.ApplyStaticMods();
	}

	/*
	 * When naming the columns in the Contrast.html file, the first part of the
	 * name comes from the sample's alias, and the second part comes from the
	 * criteria set's alias.
	 */
	public String HappyName() {
		return FriendlyName + "-" + CriteriaName;
	}

	public void DebugPrint() {
		System.out.println(HappyName());
		System.out.println(FriendlyName);
		System.out.println(SourceDirectory);
		System.out.println(CriteriaName);
		System.out.println(CutoffsString);
	}

	/*
	 * Read the list of Protein objects from the DTASelect.txt file.
	 */
	public void LoadLociFromFile() {
		File CurrentDTASelectTxt;
		try {
			CurrentDTASelectTxt = new File(this.SourceDirectory,
					"DTASelect.txt");
			this.ReadFromFile(CurrentDTASelectTxt);
			if (!this.Cutoffs.noDB) {
				try {
					// Read sequest.params file
					this.SequestParams = ParamsFile.ReadFile(new File(
							this.SourceDirectory));
				} catch (IOException failure) {
					System.out.println("Couldn't read sequest.params file in "
							+ this.SourceDirectory + ".");
					//System.exit(0);
				}
			}
		} catch (IOException failure) {
			System.out.println("Couldn't read DTASelect.txt in "
					+ this.SourceDirectory + ".");
			//System.exit(0);
		}
	}



	/*
	 * Caution! This function is used by both Contrast and DTASelect. This
	 * function applies spectrum-specific criteria, determines the redundancy of
	 * individual spectra, and then picks out the best spectrum from each
	 * redundant set. Then the program can determine which proteins to include
	 * and which to exclude.
	 */
	public void ApplyCriteria() {
		/*
		 * If the user has specified the -n option, no criteria will be applied
		 * to individual spectra. Otherwise, cut out the spectra that don't pass
		 * the appropriate criteria.
		 */


		if (Cutoffs.UseCriteria)
			LocusList.DumpUnqualifiedDTAs(Cutoffs);


		/*
		 * Determine how many copies of each spectrum remain after the spectrum
		 * filtering.
		 */
		LocusList.CalculateRedundancyForList(Cutoffs);
		if (Cutoffs.UseCriteria) {
			/*
			 * Apply locus-specific criteria, first removing redundant copies of
			 * spectra if specified and then removing proteins with insufficient
			 * redundancy, uniqueness, or peptide variety.
			 */
			// If -t 1 or -t 2 is in place, do that.
			switch (Cutoffs.PurgeDuplicateSequences) {
			case 1:
				LocusList.DitchDuplicateDTAsBySaltStep();
				break;
			case 2:
				LocusList.DitchDuplicateDTAsByXCorr();
				break;
			}
			/* Apply the -p, -r, -M, and -u filters */
			LocusList.DitchProteinsWithoutSufficientDTAs(this.Cutoffs);
			/* Apply the -e, -E, -l, and -L filters */
			/* Apply the --mw, --MW, --pi, and --PI filters */
			/* Apply the -V filter */
			Protein PRunner = LocusList;
			while (PRunner.Next != null) {
				if (Cutoffs.Permit(PRunner.Next)) {
					PRunner = PRunner.Next;
				} else {
					PRunner.Next = PRunner.Next.Next;
				}
			}
		}

		/*
		 * Now that we almost have our final list of DTAs and Loci, calculate
		 * the sequence coverage for each
		 */
		LocusList.CalculateCoverageForList(false);
		if (Cutoffs.UseCriteria) {
			LocusList.DitchProteinsWithLowSequenceCoverage(this.Cutoffs);
		}
		/*
		 * To remove subset proteins, we must first group the identicals. Since
		 * Contrast determines data set similarity after this step, we must
		 * ungroup when we're running that rather than DTASelect.
		 */
		LocusList.GroupIdenticalsOld();
		if (Cutoffs.UseCriteria && Cutoffs.RemoveSubsets) {
			LocusList.RemoveSubsetsQuick();
		}
		//LocusList.CalculateRedundancyForList(Cutoffs.PurgeDuplicateSequences);
	}

	public void ApplyCriteriaDMS() {
		/*
		 * If the user has specified the -n option, no criteria will be applied
		 * to individual spectra. Otherwise, cut out the spectra that don't pass
		 * the appropriate criteria.
		 */


		if (Cutoffs.UseCriteria)
			LocusList.DumpUnqualifiedDTAs(Cutoffs);


		/*
		 * Determine how many copies of each spectrum remain after the spectrum
		 * filtering.
		 */


		/*
		 * Now that we almost have our final list of DTAs and Loci, calculate
		 * the sequence coverage for each
		 */


		if(Cutoffs.UseCriteria)
		{
			/* Apply the -e, -E, -l, and -L filters */
			/* Apply the --mw, --MW, --pi, and --PI filters */
			/* Apply the -V filter */
			Protein PRunner = LocusList;
			while (PRunner.Next != null) {
				if (Cutoffs.Permit(PRunner.Next)) {
					PRunner = PRunner.Next;
				} else {
					PRunner.Next = PRunner.Next.Next;
				}
			}
			LocusList.CalculateFilterMedianAdjustedDeltaMass(Cutoffs);
			LocusList.CalculateCoverageForList(false);
			if (Cutoffs.UseCriteria) {
				LocusList.DitchProteinsWithLowSequenceCoverage(this.Cutoffs);
			}

			LocusList.CalculateRedundancyForList(Cutoffs);

			if (Cutoffs.UseCriteria) {
				/*
				 * Apply locus-specific criteria, first removing redundant copies of
				 * spectra if specified and then removing proteins with insufficient
				 * redundancy, uniqueness, or peptide variety.
				 */
				// If -t 1 or -t 2 is in place, do that.
				switch (Cutoffs.PurgeDuplicateSequences) {
					case 1:
						LocusList.DitchDuplicateDTAsBySaltStep();
						break;
					case 2:
						LocusList.DitchDuplicateDTAsByXCorr();
						break;
				}
				/* Apply the -p, -r, -M, and -u filters */

			}

			LocusList.DitchProteinsWithoutSufficientDTAs(this.Cutoffs);
		}
		else
		{
			LocusList.CalculateFilterMedianAdjustedDeltaMass(Cutoffs);
			LocusList.CalculateRedundancyForList(Cutoffs);
			LocusList.CalculateCoverageForList(false);
			if (Cutoffs.UseCriteria) {
				LocusList.DitchProteinsWithLowSequenceCoverage(this.Cutoffs);
			}

		}
		/*
		 * To remove subset proteins, we must first group the identicals. Since
		 * Contrast determines data set similarity after this step, we must
		 * ungroup when we're running that rather than DTASelect.
		 */
		LocusList.GroupIdenticalsOld();
		if (Cutoffs.UseCriteria && Cutoffs.RemoveSubsets) {
			LocusList.RemoveSubsetsQuick();
		}
		//LocusList.CalculateRedundancyForList(Cutoffs.PurgeDuplicateSequences);

	}


	/*
	 * This method assigns a clssification to each protein listed in the
	 * "Classifications.txt" file. When proteins are sorted as selected by the
	 * user, classification is always the first order applied.
	 */
	public void StructureByClassification() {
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File ParamsFile = new File(CurrentDirectory, "Classifications.txt");
		try {
			FileReader InputFileReader = new FileReader(ParamsFile);
			BufferedReader Incoming = new BufferedReader(InputFileReader);
			String WholeLine = Incoming.readLine();
			String LineBuffer;
			StringTokenizer Parser;
			Classification CRunner = Classifieds;
			StringBuffer Lump;
			Protein PRunner;
			while (WholeLine != null) {
				Parser = new StringTokenizer(WholeLine);
				if (Parser.hasMoreTokens()) {
					LineBuffer = Parser.nextToken();
					if (LineBuffer.equals("class")) {
						// Create a new Classification for this line.
						CRunner.Next = new Classification();
						CRunner = CRunner.Next;
						CRunner.Identifier = new Byte(Parser.nextToken())
								.byteValue();
						Lump = new StringBuffer(Parser.nextToken());
						while (Parser.hasMoreTokens()) {
							Lump.append(" ");
							Lump.append(Parser.nextToken());
						}
						CRunner.Descriptor = Lump.toString();
					} else {
						// Look for this protein descriptor in
						// memory. If it exists, assign it to the
						// appropriate Classification.
						PRunner = LocusList;
						while ((PRunner.Next != null)
								&& (!PRunner.Next.Locus.equals(LineBuffer))) {
							PRunner = PRunner.Next;
						}
						if (PRunner.Next != null) {
							// this is a match.
							try {
								PRunner.Next.Classification = new Byte(
										Parser.nextToken()).byteValue();
							} catch (NoSuchElementException failure) {
								System.out
										.println("\t\t"
												+ LineBuffer
												+ "\tis listed without a class.  Ignoring...");
							}
						}
					}
				}
				WholeLine = Incoming.readLine();
			}
			Incoming.close();
		} catch (IOException failure) {
			System.out
					.println("Error while trying to read Classifications.txt.");
			System.out.println(failure);
		}
	}

	/*
	 * This method assigns auxiliary information to each protein listed in the
	 * "AuxInfo.txt" file.
	 */
	public void IncorporateAuxInfo() {
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File ParamsFile = new File(CurrentDirectory, "AuxInfo.txt");
		try {
			FileReader InputFileReader = new FileReader(ParamsFile);
			BufferedReader Incoming = new BufferedReader(InputFileReader);
			String WholeLine = Incoming.readLine();
			String LineBuffer;
			StringTokenizer Parser;
			AuxInfo CRunner = AddedInfo;
			StringBuffer Lump;
			Protein PRunner;
			Point PoRunner;
			while (WholeLine != null) {
				Parser = new StringTokenizer(WholeLine);
				if (Parser.hasMoreTokens()) {
					LineBuffer = Parser.nextToken();
					if (LineBuffer.equals("float") || LineBuffer.equals("int")) {
						// Create a new Classification for this line.
						CRunner.Next = new AuxInfo();
						CRunner = CRunner.Next;
						CRunner.IsFloat = LineBuffer.equals("float");
						Lump = new StringBuffer(Parser.nextToken());
						while (Parser.hasMoreTokens()) {
							Lump.append(" ");
							Lump.append(Parser.nextToken());
						}
						CRunner.Descriptor = Lump.toString();
					} else {
						// Look for this protein descriptor in
						// memory. If it exists, assign it to the
						// appropriate Classification.
						PRunner = LocusList;
						while ((PRunner.Next != null)
								&& (!PRunner.Next.Locus.equals(LineBuffer))) {
							PRunner = PRunner.Next;
						}
						if (PRunner.Next != null) {
							// this is a match.
							try {
								PRunner.Next.AuxInfo = new Point();
								PoRunner = PRunner.Next.AuxInfo;
								while (Parser.hasMoreTokens()) {
									PoRunner.Next = new Point();
									PoRunner = PoRunner.Next;
									PoRunner.Intensity = new Float(
											Parser.nextToken()).floatValue();
								}
							} catch (NoSuchElementException failure) {
								System.out
										.println("\t\t"
												+ LineBuffer
												+ "\tis listed without aux info.  Ignoring...");
							}
						}
					}
				}
				WholeLine = Incoming.readLine();
			}
			Incoming.close();
		} catch (IOException failure) {
			System.out.println("Error while trying to read AuxInfo.txt.");
			System.out.println(failure);
		}
	}

	/*
	 * If the user has specified that a DataSet is the master, all others are
	 * compared to the master set and any protein not found in the master set is
	 * removed.
	 * 
	 * Assumes that proteins are alphabetically ordered!
	 */
	public void CullLociNotOnList(DataSet Master) {
		Protein ThisPRunner = LocusList;
		Protein MasterPRunner = Master.LocusList.Next;
		while ((ThisPRunner.Next != null) && (MasterPRunner != null)) {
			if (ThisPRunner.Next.Locus.equals(MasterPRunner.Locus)) {
				// The protein in question is on the list and can be
				// skipped over.
				ThisPRunner = ThisPRunner.Next;
				MasterPRunner = MasterPRunner.Next;
			} else if (ThisPRunner.Next.Locus.compareTo(MasterPRunner.Locus) > 0) {
				// "this" is further along in the alphabet than the
				// "master" list's protein.
				MasterPRunner = MasterPRunner.Next;
			} else {
				// "master" is further along in the alphabet than
				// "this" protein. We need to remove this protein
				// because it's not on the list!
				ThisPRunner.Next = ThisPRunner.Next.Next;
			}
		}
		if (MasterPRunner == null) {
			ThisPRunner.Next = null;
		}
	}

	/*
	 * This function is really the heart of Contrast. A LocusSummary is an
	 * object representing a Locus across multiple DataSets. It's like a
	 * meta-Protein; the difference is that it stores where it can be found in
	 * the FoundInPattern field. The FoundInPattern field is not just an
	 * integer; it's a pattern of bits. Each bit represents whether that Locus
	 * can be found in a particular DataSet. Which bit is affected is passed in
	 * as "CurrentPower." A zero means absence, and a one means presence.
	 */
	public void AddToLocusSummary(LocusSummary LSHead, long CurrentPower) {
		Protein Cursor = LocusList.Next;
		LocusSummary LSRunner = LSHead;
		LocusSummary LSBuffer;
		while (Cursor != null) {
			/*
			 * If the next Locus in this DataSet falls after the next in the LS
			 * list, advance the LS list
			 */
			while ((LSRunner.Next != null)
					&& (Cursor.Locus.compareTo(LSRunner.Next.Name) > 0)) {
				LSRunner = LSRunner.Next;
			}
			// if we've reached the end of the LS List, add this to its end
			if (LSRunner.Next == null) {
				LSRunner.Next = new LocusSummary();
				LSRunner = LSRunner.Next;
				LSRunner.Name = Cursor.Locus;
				LSRunner.Description = Cursor.Gene;
				LSRunner.FoundInPattern = CurrentPower;
				LSRunner.FoundCount = 1;
			}
			/*
			 * If the next Locus in this DataSet precedes the next in the LS
			 * list and doesn't yet appear in our LS list, copy it into the LS
			 * list.
			 */
			else if (Cursor.Locus.compareTo(LSRunner.Next.Name) < 0) {
				LSBuffer = new LocusSummary();
				LSBuffer.Name = Cursor.Locus;
				LSBuffer.Description = Cursor.Gene;
				LSBuffer.FoundInPattern = CurrentPower;
				LSBuffer.FoundCount = 1;
				LSBuffer.Next = LSRunner.Next;
				LSRunner.Next = LSBuffer;
			}
			/*
			 * If the next Locus in this DataSet is the same name as the next in
			 * the LS list, update its FoundInPattern
			 */
			else {
				LSRunner.Next.FoundInPattern += CurrentPower;
				LSRunner.Next.FoundCount++;
			}
			Cursor = Cursor.Next;
		}
	}

	// LIST FUNCTIONS

	/*
	 * Read the earliest locus from this collection of DataSets. Sort together
	 * peptides for this locus from each dataset, and return the StringBuffer
	 * Contrast should put into the merged DTASelect.txt
	 */
	public StringBuffer PullTopLocus() {
		boolean LociRemain = false;
		Protein PBuffer = new Protein();
		StringTokenizer Parser;
		String CurrentProtein = "zzzzzzzzzzzzzzz";
		StringBuffer Returned = new StringBuffer();
		DataSet DSRunner = this.Next;
		DTAFile DTARunner;
		// See if any proteins remain in this group of DataSets and
		// set the CurrentProtein string to the correct name. Also,
		// create a protein object for the earliest protein.
		while (DSRunner != null) {
			/*
			 * WholeLine will be the first line in this DTASelect that followed
			 * the last read "D" line. It should be either a "C" [count] line or
			 * a "L" [locus] line. If WholeLine is null, this is an older class
			 * of DTASelect.txt that included no counts.
			 */
			if ((DSRunner.WholeLine != null)
					&& (!DSRunner.WholeLine.startsWith("C\t"))) {
				LociRemain = true;
				if (DSRunner.LocusCursor != null) {
					// is this DSRunner's protein earlier in the alphabet?
					if (CurrentProtein.compareTo(DSRunner.LocusCursor.Locus) > 0) {
						CurrentProtein = DSRunner.LocusCursor.Locus;
						Parser = new StringTokenizer(DSRunner.WholeLine, "\t");
						Parser.nextToken();
						PBuffer = new Protein();
						PBuffer.SetTo(Parser);
					}
				}
			}
			DSRunner = DSRunner.Next;
		}
		if (!LociRemain) {
			return null;
		} else {
			// Build a protein with all the DTAs from each file describing
			// that protein.
			DSRunner = this.Next;
			while (DSRunner != null) {
				if ((DSRunner.LocusCursor != null)
						&& (DSRunner.LocusCursor.Locus.equals(CurrentProtein))) {
					// This DataSet contains this protein. Move past the
					// L line and grab all following D lines, constructing
					// this DataSet's Protein object as we go.
					try {
						DSRunner.WholeLine = DSRunner.Incoming.readLine();
					} catch (IOException failure) {
						System.out
								.println("Could not keep reading from DTASelect.txt in "
										+ DSRunner.SourceDirectory);
						System.exit(0);
					}
					DTARunner = DSRunner.LocusCursor.DTAs;
					while ((DSRunner.WholeLine != null)
							&& (DSRunner.WholeLine.startsWith("D\t"))) {
						DTARunner.Next = new DTAFile();
						DTARunner = DTARunner.Next;
						Parser = new StringTokenizer(DSRunner.WholeLine, "\t");
						Parser.nextToken();
						try {
							DTARunner.CurrentSetTo(Parser, Cutoffs);
						} catch (NoSuchElementException failure) {
							System.out.println("Corrupted DTA line in "
									+ DSRunner.SourceDirectory + "'s locus "
									+ CurrentProtein);
							System.out.println(WholeLine);
							System.exit(0);
						}
						try {
							DSRunner.WholeLine = DSRunner.Incoming.readLine();
						} catch (IOException failure) {
							System.out
									.println("Could not keep reading from DTASelect.txt in "
											+ DSRunner.SourceDirectory);
							System.exit(0);
						}
					}
					// Grab the DTAs from the protein we read from this
					// file, and merge them into the DTAS we've already
					// got for this protein.
					PBuffer.AccumulateDTAsFrom(DSRunner.LocusCursor);
					if (DSRunner.WholeLine.startsWith("L\t")) {
						DSRunner.LocusCursor = new Protein();
						Parser = new StringTokenizer(DSRunner.WholeLine, "\t");
						Parser.nextToken();
						try {
							DSRunner.LocusCursor.SetTo(Parser);
						} catch (NoSuchElementException failure) {
							System.out.println("Corrupted locus line in "
									+ DSRunner.SourceDirectory);
							System.out.println(WholeLine);
							System.exit(0);
						}
					} else {
						DSRunner.LocusCursor = null;
					}
				}
				DSRunner = DSRunner.Next;
			}
			System.out.print("\t" + CurrentProtein + "                   \r");
			Returned.append(PBuffer.GetLLine());
			DTARunner = PBuffer.DTAs.Next;
			while (DTARunner != null) {
				Returned.append(DTARunner.GetDLine());
				DTARunner = DTARunner.Next;
			}
			return Returned;
		}
	}

	/*
	 * How many loci are in the specified DataSet?
	 */
	public int LocusCount(int WhichOne) {
		DataSet Runner = this.Next;
		int Looper;
		int Count = 0;
		Protein Cursor;
		for (Looper = 0; Looper < WhichOne; Looper++) {
			Runner = Runner.Next;
		}
		Cursor = Runner.LocusList.Next;
		while (Cursor != null) {
			Cursor = Cursor.Next;
			Count++;
		}
		return Count;
	}

	/*
	 * Is the specified DataSet hidden?
	 */
	public boolean IsHidden(int WhichOne) {
		DataSet Runner = this.Next;
		for (int Looper = 0; Looper < WhichOne; Looper++) {
			Runner = Runner.Next;
		}
		return Runner.Removed;
	}

	/*
	 * What is the composite (sample and criteria set) name for the specified
	 * DataSet?
	 */
	public String ShortName(int WhichOne) {
		DataSet Runner = this.Next;
		for (int Looper = 0; Looper < WhichOne; Looper++) {
			Runner = Runner.Next;
		}
		return Runner.HappyName();
	}

	/*
	 * Return the HTML string for use in the header of a table for the
	 * Contrast.html file. Create the header for the specified DataSet.
	 */
	public String TableHead(int WhichOne) {
		DataSet Runner = this.Next;
		for (int Looper = 0; Looper < WhichOne; Looper++) {
			Runner = Runner.Next;
		}
		return Runner.FriendlyName + "<br>" + Runner.CriteriaName;
	}

	/*
	 * Pass a Locus name to the list of DataSets. In each DataSet, set its
	 * Cursor to point at the Protein matching that locus name. If the locus
	 * doesn't exist, point to null.
	 */
	public void SetLocusCursors(String Locus) {
		DataSet Runner = this.Next;
		Protein PRunner;
		while (Runner != null) {
			PRunner = Runner.LocusList.Next;
			while ((PRunner != null) && (!PRunner.Locus.equals(Locus))) {
				PRunner = PRunner.Next;
			}
			// Note! If locus isn't found, Cursor is set to null!
			Runner.LocusCursor = PRunner;
			Runner = Runner.Next;
		}
	}

	/* What's the length of the current protein's sequence? */
	public int CurrentSeqLength() {
		DataSet Runner = this.Next;
		int Length = 0;
		while (Runner != null && Length == 0) {
			if (Runner.LocusCursor != null) {
				Length = Runner.LocusCursor.SequenceLength;
			}
			Runner = Runner.Next;
		}
		return Length;
	}

	/*
	 * Return a pointer to the Protein object pointed to by the specified
	 * DataSet's cursor. Send back a null if the cursor is currently unassigned.
	 */
	public Protein GetLocusCursor(int WhichSet) {
		DataSet Runner = this.Next;
		for (int Looper = 0; Looper < WhichSet; Looper++) {
			Runner = Runner.Next;
		}
		return Runner.LocusCursor;
	}

	/*
	 * Set the sequence coverage percentage and consensus fields for the
	 * specified protein in the specified data set.
	 */
	public void SetLocusCoverage(int WhichSet, String Consensus,
			float PercentCoverage) {
		DataSet Runner = this.Next;
		for (int Looper = 0; Looper < WhichSet; Looper++) {
			Runner = Runner.Next;
		}
		if (Runner.LocusCursor == null) {
			System.out
					.println("Attempted to set coverage for a nonexistent locus");
		} else {
			Runner.LocusCursor.SequenceCoverage = PercentCoverage;
			Runner.LocusCursor.Coverage = Consensus;
		}
	}

	/*
	 * This function is called when the protein's sequence is null (of length
	 * 0). Sequence coverage is undefined in this case.
	 */
	public void SetLocusCoverageUndefined() {
		DataSet Runner = this.Next;
		while (Runner != null) {
			if (Runner.LocusCursor != null) {
				Runner.LocusCursor.SequenceCoverage = 0f;
				Runner.LocusCursor.Coverage = "";
			}
			Runner = Runner.Next;
		}
	}

	/*
	 * If Contrast's verbose mode is in use, determine the XCorr for the first
	 * peptide in the specified DataSet's current protein that matches the
	 * passed TargetSequence. This allows a comparison of XCorr scores for the
	 * same peptide across multiple samples.
	 */
	public String GetLocusCursorXCorr(int WhichSet, String TargetSequence,
			byte TargetZ) {
		DataSet Runner = this.Next;
		Protein PRunner;
		boolean WrongOne = true;
		for (int Looper = 0; Looper < WhichSet; Looper++) {
			Runner = Runner.Next;
		}
		if (Runner.LocusCursor == null) {
			return "Error";
		} else {
			DTAFile DRunner = Runner.LocusCursor.DTAs.Next;
			while ((DRunner != null) && (WrongOne)) {
				if ((DRunner.Sequence.equals(TargetSequence))
						&& (DRunner.ChargeState == TargetZ)) {
					WrongOne = false;
				} else {
					DRunner = DRunner.Next;
				}
			}
			if (WrongOne) {
				return "";
			} else {
				return new Float(DRunner.XCorr).toString();
			}
		}
	}

	/*
	 * Determine the current Protein's sequence coverage percentage for the
	 * specified DataSet.
	 */
	public CoverageZone GetLocusCoverageZones(int WhichSet) {
		DataSet Runner = this.Next;
		Protein PRunner;
		CoverageZone CZBuffer;
		for (int Looper = 0; Looper < WhichSet; Looper++) {
			Runner = Runner.Next;
		}
		if (Runner.LocusCursor == null) {
			// Return an error value
			return null;
		} else {
			CZBuffer = Runner.LocusCursor.GenerateZones();
			return CZBuffer;
		}
	}

	// added by Howard Choi
	public static float Rounding(float Value, int Places) {
		// Converts a value to a rounded value
		if (Places == 0) {
			return new Integer(Math.round(Value)).intValue();
		} else {
			double Multiplier = Math.pow(10, Places);
			return new Double(Math.rint(Value * Multiplier) / Multiplier)
					.floatValue();
		}
	}

	/*
	 * How many DataSets are currently in memory?
	 */
	public int CountSets() {
		DataSet Runner = this.Next;
		int Counter = 0;
		while (Runner != null) {
			Counter++;
			Runner = Runner.Next;
		}
		return Counter;
	}

	class Modification {
		int ModPositionInProt = -1;
		char ModType = 0;
		DTAFile Peptide;
		String ModCount;
		Modification Next;

		public void DebugPrint() {
			Modification Runner = this.Next;
			while (Runner != null) {
				System.out.println(new Integer(Runner.ModPositionInProt)
						.toString()
						+ "\t"
						+ new Character(Runner.ModType).toString() + "\t");
				Runner = Runner.Next;
			}
		}

		public void SortList() {
			if (this.Next != null)
				this.Next = this.Next.Sort(null);
		}

		private Modification Sort(Modification Follower) {
			// Recursive quicksorter
			// Returns first item of sorted list (from those starting at this
			// element)
			// Accepts first item to follow
			Modification ListAbove = null;
			Modification ListBelow = null;
			Modification PlaceHolder;
			Modification PlaceHolder2;
			PlaceHolder = this.Next;
			// Partition all remaining points of this linked list
			while (PlaceHolder != null) {
				PlaceHolder2 = PlaceHolder.Next;
				if (this.ModPositionInProt > PlaceHolder.ModPositionInProt) {
					// Move this item to list above this
					PlaceHolder.Next = ListAbove;
					ListAbove = PlaceHolder;
				} else if (this.ModPositionInProt < PlaceHolder.ModPositionInProt) {
					// Move this item to list below this point
					PlaceHolder.Next = ListBelow;
					ListBelow = PlaceHolder;
				} else {
					if (this.ModType > PlaceHolder.ModType) {
						// Move this item to list above this
						PlaceHolder.Next = ListAbove;
						ListAbove = PlaceHolder;
					} else if (this.ModType < PlaceHolder.ModType) {
						// Move this item to list below this point
						PlaceHolder.Next = ListBelow;
						ListBelow = PlaceHolder;
					} else {
						if (this.Peptide.SequencePosition > PlaceHolder.Peptide.SequencePosition) {
							// Move this item to list above this
							PlaceHolder.Next = ListAbove;
							ListAbove = PlaceHolder;
						} else if (this.Peptide.SequencePosition < PlaceHolder.Peptide.SequencePosition) {
							// Move this item to list below this point
							PlaceHolder.Next = ListBelow;
							ListBelow = PlaceHolder;
						} else {
							if (this.Peptide.PrecursorMass > PlaceHolder.Peptide.PrecursorMass) {
								// Move this item to list above this
								PlaceHolder.Next = ListAbove;
								ListAbove = PlaceHolder;
							} else {
								// Move this item to list below this point
								PlaceHolder.Next = ListBelow;
								ListBelow = PlaceHolder;
							}
						}
					}
				}
				// Move to next item to be partitioned
				PlaceHolder = PlaceHolder2;
			}
			if (ListBelow == null)
				this.Next = Follower;
			else
				this.Next = ListBelow.Sort(Follower);
			if (ListAbove == null)
				return this;
			else
				return ListAbove.Sort(this);
		}
	}
}
