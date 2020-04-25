
import java.io.*;
import java.util.*;

//Protein object
//For DTASelect
//Created by Dave Tabb
//Begun 8/11/2000

public class Protein {
	Protein Next;
	String Locus;
	String Gene = "";
	String Coverage = "";
	float MolWt;
	float pI;
	int SequenceLength;
	float SequenceCoverage = 0f;
	float ProtConf = 0.0f;
	float ProtFP = 1.0f;
	DTAFile DTAs = new DTAFile();
	Protein IdenticalLoci = null;
	char Validated = 'U';
	char NChar = 'X';
	char CChar = 'X';
	byte Classification = 127;
	Point AuxInfo;
	boolean IsDecoy = false;
	boolean HasGreatPeptide = true;
	int NPeptides = 0;
	int NSpectra = 0;
	int NModPeptides = 0;
	int NTrypticPeptides = 0;
	// Added by Diego Calzolari, 11/05/2012
	float SAF;
	float NSAF;
	float EMPAI;
	int id;
	boolean toRemove = false;
	private static int ID_COUNTER = 0;
	Set<String> fileNameSet = null;
	/*
	 * INDIVIDUAL FUNCTIONS These are functions to be called on an individual
	 * Protein in a list rather than on the list header itself
	 */

	public Protein()
	{
		id = ID_COUNTER ++ ;
	}

	public Protein Clone() {
		Protein Copy = new Protein();
		DTAFile DTARunner = this.DTAs.Next;
		DTAFile CopyRunner = Copy.DTAs;
		Copy.Locus = this.Locus;
		Copy.Gene = this.Gene;
		Copy.Coverage = this.Coverage;
		Copy.MolWt = this.MolWt;
		Copy.ProtConf = this.ProtConf;
		Copy.ProtFP = this.ProtFP;
		Copy.pI = this.pI;
		Copy.SequenceLength = this.SequenceLength;
		Copy.SequenceCoverage = this.SequenceCoverage;
		Copy.IsDecoy = this.IsDecoy;
		Copy.HasGreatPeptide = this.HasGreatPeptide;
		Copy.NPeptides = this.NPeptides;
		Copy.NSpectra = this.NSpectra;
		Copy.NModPeptides = this.NModPeptides;
		Copy.NTrypticPeptides = this.NTrypticPeptides;
		while (DTARunner != null) {
			CopyRunner.Next = DTARunner.Clone();
			CopyRunner = CopyRunner.Next;
			DTARunner = DTARunner.Next;
		}
		return Copy;
	}

	/*
	 * Set the fields of this protein to the contents of this StringTokenizer
	 * (which delimits by tabs and has the initial "L" truncated), drawn from a
	 * DTASelect.txt file. Called by both Contrast Merger and
	 * DataSet.ReadFromFile
	 */
	public void SetTo(StringTokenizer Parser) throws NoSuchElementException {
		this.Locus = Parser.nextToken();
		// if this file includes information about the sequence...
		if (Parser.hasMoreTokens()) {
			this.SequenceLength = new Integer(Parser.nextToken()).intValue();
			// this.ProtConf = new Float(Parser.nextToken()).floatValue();
			// this.ProtFP = new Float(Parser.nextToken()).floatValue();
			this.MolWt = new Float(Parser.nextToken()).floatValue();
			this.pI = new Float(Parser.nextToken()).floatValue();
			// if we're so lucky as to have a descriptive gene name...
			if (Parser.hasMoreTokens()) {
				this.Gene = Parser.nextToken();
				// Maybe we even have validation info?
				if (Parser.hasMoreTokens()) {
					this.Validated = Parser.nextToken().charAt(0);
				}
			} else
				this.Gene = "";
		} else {
			this.SequenceLength = 0;
			this.MolWt = 0f;
			this.ProtConf = 0.0f;
			this.ProtFP = 1.0f;
			this.pI = 0f;
			this.Gene = "";
		}
	}

	/*
	 * Compare this protein to the passed one and produce a measure of its
	 * similarity.
	 * 
	 * (9/26/02) I noticed that this code relied upon the peptides being ordered
	 * alphabetically by sequence. Now this code presumes nothing about the
	 * order.
	 */
	public int Similarity(Protein Other) {
		DTAFile ThisRunner = this.DTAs.Next;
		DTAFile OtherRunner;
		int SameCount = 0;
		int ThisUnmatched = 0;
		int OtherCount = 0;
		OtherRunner = Other.DTAs.Next;
		while (OtherRunner != null) {
			OtherCount++;
			OtherRunner = OtherRunner.Next;
		}
		while (ThisRunner != null) {
			OtherRunner = Other.DTAs.Next;
			// Try to find a match for This filename in the Other list.
			while ((OtherRunner != null)
					&& (!OtherRunner.FileName.equals(ThisRunner.FileName))) {
				OtherRunner = OtherRunner.Next;
			}
			if (OtherRunner != null)
				SameCount++;
			else
				ThisUnmatched++;
			ThisRunner = ThisRunner.Next;
		}
		return (SameCount * 3 - (ThisUnmatched + OtherCount));
	}

	/*
	 * Return the salient features of this Protein, separating with tabs. Used
	 * for GUI output
	 */

	public double GetRawProbs() {
		DTAFile DTARunner;
		double Negative = 1.0;

		DTARunner = DTAs.Next;
		while (DTARunner != null) {
			Negative = Negative * (1.0 - DTARunner.Probability);
			DTARunner = DTARunner.Next;
		}
		return (1.0 - Negative);
	}

	public String GetLocusTabbed() {
		return (Locus + "\t" + NPeptides + "\t" + NSpectra + "\t"
				+ SequenceCoverage + "%\t"
				+ new Integer(SequenceLength).toString()
				+ "\t"
				+
				// new Float(RoundTo(ProtConf,2)).toString() + "\t" +
				// new Float(RoundTo(ProtFP,2)).toString() + "\t" +
				new Integer(Math.round(MolWt)).toString() + "\t"
				+ new Float(RoundTo(pI, 1)).toString() + "\t" + Gene);
	}

	/*
	 * Return the "L" line for DTASelect.txt corresponding to this protein. Used
	 * both by DataSet.PrintTXT and DataSet.PullTopLocus
	 */
	public String GetLLine() {
		String Tab = "\t";
		return ("L\t"
				+ this.Locus
				+ Tab
				+ new Integer(this.SequenceLength).toString()
				+ Tab
				+
				// new Float(this.ProtConf).toString() + Tab +
				// new Float(this.ProtFP).toString() + Tab +
				new Float(this.MolWt).toString() + Tab
				+ new Float(this.pI).toString() + Tab + this.Gene + Tab
				+ this.Validated + "\n");
	}

	/*
	 * For Contrast verbose mode, be able to create a Protein object that
	 * includes all peptides from Proteins found in each combination of sample
	 * and criteria set. This function is used to accumulate a master Protein by
	 * including all the peptides of Proteins passed to it. It's just a straight
	 * cosequential algorithm like one would use in a merge sort.
	 */
	public void AccumulateDTAsFrom(Protein Other) {
		DTAFile Runner = Other.DTAs.Next;
		DTAFile Buffer;
		DTAFile ThisRunner;
		while (Runner != null) {
			Buffer = Runner.Clone();
			ThisRunner = this.DTAs;
			while ((ThisRunner.Next != null)
					&& (ThisRunner.Next.Sequence.compareTo(Buffer.Sequence) < 0)) {
				ThisRunner = ThisRunner.Next;
			}
			if ((ThisRunner.Next != null)
					&& (ThisRunner.Next.Sequence.equals(Buffer.Sequence))) {
				while (ThisRunner.Next != null
						&& ThisRunner.Next.Sequence.equals(Buffer.Sequence)
						&& ThisRunner.Next.ChargeState < Buffer.ChargeState) {
					ThisRunner = ThisRunner.Next;
				}
			}
			Buffer.Next = ThisRunner.Next;
			ThisRunner.Next = Buffer;
			Runner = Runner.Next;
		}
	}

	/* Implement the --copy option for unified SQT files. */
	public void MakeSQTSubsets(File CurrentDirectory) {
		System.out.println("\tSubsetting SQT files...");
		// Create a Subsets directory containing SQT files with just
		// the peptides remaining after filtering. Implements the
		// --copy feature for unified files
		File TargetDir = new File(CurrentDirectory, "Subsets");
		File CurrentSQT;
		File TargetFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		FileReader InputFilereader;
		BufferedReader Incoming;
		String[] DirectoryListing;
		String[] UniListing;
		StringTokenizer Parser;
		int FileIterator;
		int RealCount; // How many IDs were in the original file? Spectra?
		int SQTCount; // How many are we keeping?
		int SoughtIndex;
		String CurrentBaseName;
		String LineBuffer = "";
		String WholeLine;
		StringBuffer Header;
		StringBuffer CurrentSection;
		StringBuffer Report = new StringBuffer("SQTFile\tIDsKept\tIDsFound\n");
		String newline = "\n";
		TextBlock Contents;
		TextBlock TBRunner;
		DTAFile DTARunner = DTAs.Next;
		EndsWithFilter FilenameFilter = new EndsWithFilter("qt");
		UniListing = CurrentDirectory.list(FilenameFilter);
		if (UniListing.length == 0) {
			System.out.println("\tNo SQT Files in directory");
		} else {
			// Set up the Subsets directory for our new SQT and MS2 files
			if (TargetDir.exists()) {
				// Delete all SQTs currently in Subsets directory
				DirectoryListing = TargetDir.list(FilenameFilter);
				for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
					CurrentSQT = new File(TargetDir,
							DirectoryListing[FileIterator]);
					CurrentSQT.delete();
				}
			} else {
				// Make the subsets directory
				TargetDir.mkdir();
			}
			/*
			 * Visit each SQT file. Read it into memory, sort it, and compare it
			 * to the list of IDs passing filters to determine which are kept.
			 */
			FilenameFilter = new EndsWithFilter("qt");
			UniListing = CurrentDirectory.list(FilenameFilter);
			for (FileIterator = 0; FileIterator < UniListing.length; FileIterator++) {
				CurrentBaseName = UniListing[FileIterator];
				System.out.println("\t\tReading " + CurrentBaseName
						+ " into memory...");
				CurrentSQT = new File(CurrentDirectory, CurrentBaseName);
				RealCount = 0;
				SQTCount = 0;
				TargetFile = new File(TargetDir, CurrentBaseName);
				CurrentBaseName = CurrentBaseName.substring(0,
						CurrentBaseName.length() - 4);
				try {
					InputFilereader = new FileReader(CurrentSQT);
					Incoming = new BufferedReader(InputFilereader);
					Header = new StringBuffer();
					Contents = new TextBlock();
					TBRunner = Contents;
					WholeLine = Incoming.readLine();
					// slap the header stuff into Header StringBuffer
					while ((WholeLine != null)
							&& (!WholeLine.startsWith("S\t"))) {
						Header.append(WholeLine + newline);
						WholeLine = Incoming.readLine();
					}
					// slap the IDs into Contents TextBlocks.
					while (WholeLine != null) {
						TBRunner.Next = new TextBlock();
						RealCount++;
						TBRunner = TBRunner.Next;
						CurrentSection = new StringBuffer(WholeLine + newline);
						Parser = new StringTokenizer(WholeLine);
						try {
							LineBuffer = WholeLine;
							// Skip the S
							Parser.nextToken();
							TBRunner.Index = new Integer(Parser.nextToken())
									.intValue() * 10;
							// Skip the second scan
							Parser.nextToken();
							TBRunner.Index += new Integer(Parser.nextToken())
									.intValue();
							WholeLine = Incoming.readLine();
						} catch (NumberFormatException failure) {
							System.out.println("This line of "
									+ CurrentBaseName + ".sqt was garbled:");
							System.out.println(LineBuffer);
							// Skip to the next identification
							WholeLine = Incoming.readLine();
							while ((WholeLine != null)
									&& (!WholeLine.startsWith("S\t"))) {
								WholeLine = Incoming.readLine();
							}
						} catch (NoSuchElementException oopsie) {
							System.out.println("This line of "
									+ CurrentBaseName + ".sqt was too short:");
							System.out.println(LineBuffer);
							// Skip to the next identification
							WholeLine = Incoming.readLine();
							while ((WholeLine != null)
									&& (!WholeLine.startsWith("S\t"))) {
								WholeLine = Incoming.readLine();
							}
						}
						while ((WholeLine != null)
								&& (!WholeLine.startsWith("S\t"))) {
							CurrentSection.append(WholeLine + newline);
							WholeLine = Incoming.readLine();
						}
						TBRunner.Text = CurrentSection.toString();
					}
					// Sort the IDs stored in memory.
					Contents.Next = MergeSort(Contents.Next);
					// Compare each file identification to the list in memory.
					DTARunner = this.DTAs.Next;
					// Queue up the first DTA with this root name
					while ((DTARunner != null)
							&& (!DTARunner.FileName.startsWith(CurrentBaseName
									+ ".")))
						DTARunner = DTARunner.Next;
					if (DTARunner == null) {
						System.out.println("\t\tNo retained IDs for this file");
					} else {
						// Mark each in the unified file which is retained in
						// memory
						TBRunner = Contents.Next;
						while ((DTARunner != null)
								&& (DTARunner.FileName
										.startsWith(CurrentBaseName + "."))) {
							Parser = new StringTokenizer(DTARunner.FileName,
									".");
							// Skip root name
							Parser.nextToken();
							// multiply ten by loscan number
							SoughtIndex = new Integer(Parser.nextToken())
									.intValue() * 10;
							// skip hiscan
							Parser.nextToken();
							// use charge state
							SoughtIndex += new Integer(Parser.nextToken())
									.intValue();
							TBRunner = Contents.Next;
							while ((TBRunner != null)
									&& (TBRunner.Index != SoughtIndex))
								TBRunner = TBRunner.Next;
							if (TBRunner != null) {
								if (TBRunner.Index == SoughtIndex) {
									TBRunner.Keep = true;
									SQTCount++;
								} else {
									System.out
											.println("SQT: SoughtIndex "
													+ SoughtIndex
													+ " is different from TBRunner.Index "
													+ TBRunner.Index);
								}
							} else {
								System.out
										.println(DTARunner.FileName
												+ " could not be found in the sqt file.");
							}
							DTARunner = DTARunner.Next;
						}
						// Write new copy of unified file
						System.out.println("\t\tWriting " + SQTCount + " / "
								+ RealCount + " IDs to new file...");
						Report.append(CurrentBaseName + "\t" + SQTCount + "\t"
								+ RealCount + newline);
						try {
							OutputFileWriter = new FileWriter(TargetFile);
							Outgoing = new BufferedWriter(OutputFileWriter);
							Outgoing.write(Header.toString());
							TBRunner = Contents.Next;
							while (TBRunner != null) {
								if (TBRunner.Keep)
									Outgoing.write(TBRunner.Text);
								TBRunner = TBRunner.Next;
							}
							Outgoing.flush();
							Outgoing.close();
						} catch (IOException failure) {
							System.out
									.println("An IO Error occured while writing "
											+ TargetFile.toString());
							System.out.println(failure);
						}
					}
				} catch (IOException failure) {
					System.out.println("IO Error while reading "
							+ CurrentBaseName);
					System.out.println(failure);
				}
			}
			TargetFile = new File(TargetDir, "SQTCopyReport.txt");
			try {
				OutputFileWriter = new FileWriter(TargetFile);
				Outgoing = new BufferedWriter(OutputFileWriter);
				Outgoing.write(Report.toString());
				Outgoing.flush();
				Outgoing.close();
			} catch (IOException failure) {
				System.out.println("IO error while writing --copy report");
				System.out.println(failure);
			}
		}
	}

	/* Implement the --copy option for unified MS2 files. */
	public void MakeMS2Subsets(File CurrentDirectory) {
		System.out.println("\tSubsetting MS2 files...");
		// Create a Subsets directory containing MS2 files with just
		// the peptides remaining after filtering. Implements the
		// --copy feature for unified files
		File TargetDir = new File(CurrentDirectory, "Subsets");
		File CurrentMS2;
		File TargetFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		FileReader InputFilereader;
		BufferedReader Incoming;
		String[] DirectoryListing;
		String[] UniListing;
		StringTokenizer Parser;
		int FileIterator;
		int RealCount; // How many IDs were in the original file? Spectra?
		int MS2Count; // How many are we keeping?
		int SoughtIndex;
		String CurrentBaseName;
		String LineBuffer = "";
		String WholeLine;
		StringBuffer Header;
		StringBuffer CurrentSection;
		StringBuffer Report = new StringBuffer(
				"MS2File\tSpectraKept\tSpectraFound\n");
		String newline = "\n";
		TextBlock Contents;
		TextBlock TBRunner;
		DTAFile DTARunner = DTAs.Next;
		EndsWithFilter FilenameFilter = new EndsWithFilter(".ms2");
		UniListing = CurrentDirectory.list(FilenameFilter);
		if (UniListing.length == 0) {
			System.out.println("\tNo MS2 Files in directory");
		} else {
			// Set up the Subsets directory for our new MS2 files
			if (TargetDir.exists()) {
				// Delete all MS2 currently in Subsets directory
				DirectoryListing = TargetDir.list(FilenameFilter);
				for (FileIterator = 0; FileIterator < DirectoryListing.length; FileIterator++) {
					CurrentMS2 = new File(TargetDir,
							DirectoryListing[FileIterator]);
					CurrentMS2.delete();
				}
			} else {
				// Make the subsets directory
				TargetDir.mkdir();
			}
			/*
			 * Visit each MS2 file. Read it into memory, sort it, and compare it
			 * to the list of IDs passing filters to determine which are kept.
			 */
			FilenameFilter = new EndsWithFilter(".ms2");
			UniListing = CurrentDirectory.list(FilenameFilter);
			for (FileIterator = 0; FileIterator < UniListing.length; FileIterator++) {
				CurrentBaseName = UniListing[FileIterator];
				System.out.println("\t\tReading " + CurrentBaseName
						+ " into memory...");
				CurrentMS2 = new File(CurrentDirectory, CurrentBaseName);
				RealCount = 0;
				MS2Count = 0;
				TargetFile = new File(TargetDir, CurrentBaseName);
				CurrentBaseName = CurrentBaseName.substring(0,
						CurrentBaseName.length() - 4);
				try {
					InputFilereader = new FileReader(CurrentMS2);
					Incoming = new BufferedReader(InputFilereader);
					Header = new StringBuffer();
					Contents = new TextBlock();
					TBRunner = Contents;
					WholeLine = Incoming.readLine();
					// slap the header stuff into Header StringBuffer
					while ((WholeLine != null)
							&& (!WholeLine.startsWith("S\t"))) {
						Header.append(WholeLine + newline);
						WholeLine = Incoming.readLine();
					}
					// slap the spectra into Contents TextBlocks.
					while (WholeLine != null) {
						TBRunner.Next = new TextBlock();
						RealCount++;
						TBRunner = TBRunner.Next;
						CurrentSection = new StringBuffer(WholeLine + newline);
						Parser = new StringTokenizer(WholeLine);
						try {
							LineBuffer = WholeLine;
							// Skip the S
							Parser.nextToken();
							TBRunner.Index = new Integer(Parser.nextToken())
									.intValue();
							// Skip the high scan number and charge
							// state!
							WholeLine = Incoming.readLine();
						} catch (NumberFormatException failure) {
							System.out.println("This line of "
									+ CurrentBaseName + ".ms2 was garbled:");
							System.out.println(LineBuffer);
							// Skip to the next identification
							WholeLine = Incoming.readLine();
							while ((WholeLine != null)
									&& (!WholeLine.startsWith("S\t"))) {
								WholeLine = Incoming.readLine();
							}
						} catch (NoSuchElementException oopsie) {
							System.out.println("This line of "
									+ CurrentBaseName + ".ms2 was too short:");
							System.out.println(LineBuffer);
							// Skip to the next identification
							WholeLine = Incoming.readLine();
							while ((WholeLine != null)
									&& (!WholeLine.startsWith("S\t"))) {
								WholeLine = Incoming.readLine();
							}
						}
						while ((WholeLine != null)
								&& (!WholeLine.startsWith("S\t"))) {
							CurrentSection.append(WholeLine + newline);
							WholeLine = Incoming.readLine();
						}
						TBRunner.Text = CurrentSection.toString();
					}
					// Sort the IDs stored in memory.
					Contents.Next = MergeSort(Contents.Next);
					// Compare each file identification to the list in memory.
					DTARunner = this.DTAs.Next;
					// Queue up the first DTA with this root name
					while ((DTARunner != null)
							&& (!DTARunner.FileName.startsWith(CurrentBaseName
									+ ".")))
						DTARunner = DTARunner.Next;
					if (DTARunner == null) {
						System.out
								.println("\t\tNo retained spectra for this file");
					} else {
						// Mark each in the unified file which is retained in
						// memory
						TBRunner = Contents.Next;
						while ((DTARunner != null)
								&& (DTARunner.FileName
										.startsWith(CurrentBaseName + "."))) {
							Parser = new StringTokenizer(DTARunner.FileName,
									".");
							// Skip root name
							Parser.nextToken();
							// Set index to scan number
							SoughtIndex = new Integer(Parser.nextToken())
									.intValue();
							TBRunner = Contents.Next;
							while ((TBRunner != null)
									&& (TBRunner.Index != SoughtIndex))
								TBRunner = TBRunner.Next;
							if (TBRunner != null) {
								if (TBRunner.Index == SoughtIndex) {
									TBRunner.Keep = true;
									MS2Count++;
								} else {
									System.out
											.println("MS2: SoughtIndex "
													+ SoughtIndex
													+ " is different from TBRunner.Index "
													+ TBRunner.Index);
								}
							} else {
								System.out
										.println(DTARunner.FileName
												+ " could not be found in the ms2 file.");
							}
							DTARunner = DTARunner.Next;
						}
						// Write new copy of unified file
						System.out.println("\t\tWriting " + MS2Count + " / "
								+ RealCount + " spectra to new file...");
						Report.append(CurrentBaseName + "\t" + MS2Count + "\t"
								+ RealCount + newline);
						try {
							OutputFileWriter = new FileWriter(TargetFile);
							Outgoing = new BufferedWriter(OutputFileWriter);
							Outgoing.write(Header.toString());
							TBRunner = Contents.Next;
							while (TBRunner != null) {
								if (TBRunner.Keep)
									Outgoing.write(TBRunner.Text);
								TBRunner = TBRunner.Next;
							}
							Outgoing.flush();
							Outgoing.close();
						} catch (IOException failure) {
							System.out
									.println("An IO Error occured while writing "
											+ TargetFile.toString());
							System.out.println(failure);
						}
					}
				} catch (IOException failure) {
					System.out.println("IO Error while reading "
							+ CurrentBaseName);
					System.out.println(failure);
				}
			}
			TargetFile = new File(TargetDir, "MS2CopyReport.txt");
			try {
				OutputFileWriter = new FileWriter(TargetFile);
				Outgoing = new BufferedWriter(OutputFileWriter);
				Outgoing.write(Report.toString());
				Outgoing.flush();
				Outgoing.close();
			} catch (IOException failure) {
				System.out.println("IO error while writing --copy report");
				System.out.println(failure);
			}
		}
	}

	public static TextBlock MergeSort(TextBlock List) {
		TextBlock SecondList;
		if (List == null) {
			return null;
		} else if (List.Next == null) {
			return List;
		} else {
			SecondList = Split(List);
			return Merge(MergeSort(List), MergeSort(SecondList));
		}
	}

	public static TextBlock Merge(TextBlock List1, TextBlock List2) {
		if (List1 == null) {
			return List2;
		} else if (List2 == null) {
			return List1;
		} else if (List1.Index < List2.Index) {
			List1.Next = Merge(List1.Next, List2);
			return List1;
		} else {
			List2.Next = Merge(List1, List2.Next);
			return List2;
		}
	}

	public static TextBlock Split(TextBlock FirstList) {
		// Part of merge sort. Returns the even members of this
		// list as another list. Assumes no header cells!
		TextBlock SecondList;
		if ((FirstList == null) || (FirstList.Next == null)) {
			return null;
		} else {
			SecondList = FirstList.Next;
			FirstList.Next = SecondList.Next;
			SecondList.Next = Split(SecondList.Next);
			return SecondList;
		}
	}

	/*
	 * Create the DB-Peptides.txt file
	 */
	public void PrintDBPeptides(File CurrentDirectory, String DBPrefix,
			String hXCorr, String hCalcPreMass, String hDeltCN,
			String hSpScore, SelectCriteria Cutoffs) {
		System.out.println("\tCreating Peptides database...");
		// Create the file
		DTAFile DTARunner = this.DTAs.Next;
		File PeptideDB = new File(CurrentDirectory, DBPrefix + "-Peptides.txt");
		try {
			FileWriter PDBFW = new FileWriter(PeptideDB);
			BufferedWriter PDBBW = new BufferedWriter(PDBFW);
			PDBBW.write("MatchID\tZ\t"
					+ hXCorr
					+ "\t"
					+ hDeltCN
					+ "\tM+H+\t"
					+ hCalcPreMass
					+ (Cutoffs.DisplayDeltaMass && Cutoffs.ExtraColumns ? "\tPPM"
							: "")
					+ "\tSpR\t"
					+ hSpScore
					+ (Cutoffs.DisplayPI && Cutoffs.ExtraColumns ? "\tpI" : "")
					+ (Cutoffs.DisplayKD && Cutoffs.ExtraColumns ? "\tKD" : "")
					+ (Cutoffs.DisplayBB && Cutoffs.ExtraColumns ? "\tBB" : "")
					+ (Cutoffs.DisplayHPLC && Cutoffs.ExtraColumns ? "\tHPLC"
							: "")
					+ "\tIonProportion\tCount\tSequence\tUnique\n");
			while (DTARunner != null) {
				// If this is a new filename
				PDBBW.write(DTARunner.FileName
						+ "\t"
						+ new Byte(DTARunner.ChargeState).toString()
						+ "\t"
						+ new Float(DTARunner.XCorr).toString()
						+ "\t"
						+ new Float(DTARunner.DeltCN).toString()
						+ "\t"
						+ new Float(DTARunner.PrecursorMass).toString()
						+ "\t"
						+ new Float(DTARunner.CalcPreMass).toString()
						+ (Cutoffs.DisplayDeltaMass && Cutoffs.ExtraColumns ? "\t"
								+ Protein.RoundTo(new Float(
										DTARunner.Adjusted_PPM_Offset)
										.floatValue(), 1)
								: "")
						+ "\t"
						+ new Integer(DTARunner.Sp).toString()
						+ "\t"
						+ new Float(DTARunner.SpScore).toString()
						+ (Cutoffs.DisplayPI && Cutoffs.ExtraColumns ? "\t"
								+ Protein.RoundTo(DTARunner.CalculatepI(), 2)
								: "")
						+ (Cutoffs.DisplayKD && Cutoffs.ExtraColumns ? "\t"
								+ Protein.RoundTo(
										DTARunner.CalculateKyteDoolittle(), 1)
								: "")
						+ (Cutoffs.DisplayBB && Cutoffs.ExtraColumns ? "\t"
								+ Protein.RoundTo(
										DTARunner.CalculateBullBreese(), 1)
								: "")
						+ (Cutoffs.DisplayHPLC && Cutoffs.ExtraColumns ? "\t"
								+ Protein.RoundTo(
										DTARunner.CalculateHPLCpH34(), 1) : "")
						+ "\t"
						+ new Float(100f * DTARunner.IonProportion).toString()
						+ "\t" + new Integer(DTARunner.Redundancy).toString()
						+ "\t" + DTARunner.Sequence + "\t"
						+ (DTARunner.UniqueToLocus ? "*" : " ") + "\n");
				DTARunner = DTARunner.Next;
			}
			PDBBW.flush();
			PDBBW.close();
		} catch (IOException oopsie) {
			System.out.println("Error while writing " + DBPrefix
					+ "-Peptides.txt");
		}
	}

	/*
	 * Return list of peptides for this locus. Used for Contrast's verbose mode
	 * during output.
	 */
	public String ReportPeptideList() {
		StringBuffer List = new StringBuffer();
		DTAFile Runner = this.DTAs.Next;
		while (Runner != null) {
			if (Runner.UniqueToLocus)
				List.append("Y");
			else
				List.append("N");
			List.append("\t");
			List.append(Runner.Sequence);
			List.append("\t");
			List.append(new Integer(Runner.ChargeState).toString());
			List.append("\n");
			Runner = Runner.Next;
		}
		return List.toString();
	}

	/*
	 * Determine whether this Protein has at least one unique peptide associated
	 * with it. This corresponds to the -u criterion.
	 */
	public boolean HasAtLeastOneUnique() {
		DTAFile Runner = this.DTAs.Next;
		while (Runner != null) {
			if (Runner.UniqueToLocus)
				return true;
			Runner = Runner.Next;
		}
		return false;
	}

	/*
	 * Determine how many times each sequence shows up among the peptides for
	 * this locus. Store this value with each peptide. This must be generated
	 * each time DTASelect is run since the criteria may eliminate some
	 * duplicate spectra and this number is generated after the
	 * spectrum-specific criteria are imposed.
	 */
	public void CalculateRedundancy(SelectCriteria criteria) {
		DTAFile Runner = this.DTAs.Next;
		DTAFile SubRunner;
		int Purge = criteria.PurgeDuplicateSequences;
		short Count;

		NPeptides = 0;
		NSpectra = 0;
		NModPeptides = 0;
		NTrypticPeptides = 0;
		// Assume Runner is on the first of identicals
		while (Runner != null) {
			/*if(criteria.PrintOnlyUnique && Runner.UniqueToLocus) NPeptides++;
			else if(!criteria.PrintOnlyUnique)NPeptides++;*/
			NPeptides++;
			if (Runner.Modified) {
				NModPeptides++;
			}
			if (Runner.Tryptic == 2) {
				NTrypticPeptides++;
			}
			SubRunner = Runner.Next;
			Count = 1;
			// changed by Howard Choi
			if (Purge == 1) {
				while ((SubRunner != null)
						&& SubRunner.Sequence.equals(Runner.Sequence)
						&& SubRunner.ChargeState == Runner.ChargeState
						&& SubRunner.FileName.substring(0,
								SubRunner.FileName.indexOf(".")).equals(
								Runner.FileName.substring(0,
										Runner.FileName.indexOf(".")))) {
					Count++;
					SubRunner = SubRunner.Next;
				}
			} else {
				while ((SubRunner != null)
						&& SubRunner.Sequence.equals(Runner.Sequence)
						&& SubRunner.ChargeState == Runner.ChargeState) {
					Count++;
					SubRunner = SubRunner.Next;
				}
			}
			NSpectra = NSpectra + Count;
			while ((Runner != null) && (Runner != SubRunner)) {
				if (Purge == 0) {
					Runner.Redundancy = 1;
				} else if (Purge == 1) {
					Runner.Redundancy = Count;
				} else {
					Runner.Redundancy = Count;
				}
				Runner = Runner.Next;
			}
		}
	}
	public boolean PassScores(SelectCriteria Cutoffs)
	{
		return PassScores(Cutoffs, false);
	}

	public boolean PassScores(SelectCriteria Cutoffs, boolean fromMedianfilter) {
		DTAFile Runner = this.DTAs.Next;

		while (Runner != null) {
			boolean dmFilter;
			if(Cutoffs.dms)
			{
				if(fromMedianfilter)
				{

					double dm =  Runner.Shifted_PPM_Offset;
					dmFilter = Math.abs(dm) <= Cutoffs.MaxProtDM;
				}
				else
				{
					dmFilter = true;
				}
			}
			else
			{
				double dm = Runner.Adjusted_PPM_Offset;
				dmFilter = Math.abs(dm) <= Cutoffs.MaxProtDM;
			}

			if (Runner.Tryptic >= Cutoffs.MinProtTryptic
					&& Runner.SpScore >= Cutoffs.MinProtSpScore
					&& Runner.XCorr >= Cutoffs.MinProtXCorr
					&& Runner.PepFP <= Cutoffs.MaxProtMinFP
					&& Runner.PepConf >= Cutoffs.MinProtPepConf
					&& dmFilter) {
				return true;
			}
			Runner = Runner.Next;
		}
		return false;
	}

	/*
	 * GenerateZones helps determine how many residues of the locus sequence are
	 * represented by residues in the DTAFiles. A CoverageZone object is created
	 * for each peptide remaining for this protein.
	 */
	public CoverageZone GenerateZones() {
		DTAFile Runner = this.DTAs.Next;
		CoverageZone Zones = new CoverageZone();
		CoverageZone ZRunner = Zones;
		/*
		 * Create a list of CoverageZones that corresponds to all DTAFiles that
		 * may map to this locus. Each CZ stores the first and last string
		 * indecies covered by an individual peptide sequence.
		 */
		while (Runner != null) {
			if (Runner.SequencePosition != -1) {
				if(Runner.sequencePositionList.size()>1)
				{
					for(int seqPos: Runner.sequencePositionList)
					{
						ZRunner = new CoverageZone(seqPos,
								Runner.TrimmedSequence());
						ZRunner.Next = Zones.Next;
						Zones.Next = ZRunner;
					}
				}
				else
				{
					ZRunner = new CoverageZone(Runner.SequencePosition,
							Runner.TrimmedSequence());
					ZRunner.Next = Zones.Next;
					Zones.Next = ZRunner;
				}

			}
			Runner = Runner.Next;
		}
		return Zones;
	}

	/*
	 * Determine where each DTAFile sequence maps to the full protein sequence.
	 * Since SEQUEST incorrectly implies that DTAFiles with multiple locus
	 * matches have the same sequence context, fix the sequence contexts for all
	 * peptides (by changing the first and last characters of the Sequence
	 * field). If a subsequence doesn't match the database sequence to which it
	 * is mapped, send a message to the user and prepare to panic.
	 */
	public void FindDTAPositions(String Sequence, String Locus,
			SelectCriteria Cutoffs) {
		DTAFile Runner = this.DTAs.Next;
		if (Sequence.length() > 0) {
			String Subsequence;
			char Leader;
			char Trailer;
			while (Runner != null) {
				Subsequence = Runner.TrimmedSequence();
				if (Subsequence.length() == 0) {
					System.out.println(Runner.FileName
							+ " reports no sequence!");
				} else {
					Runner.SequencePosition = Sequence.indexOf(Subsequence);
					findAllSequenceLocations(Sequence,Subsequence,Runner);
					if (Runner.SequencePosition == -1) {
						System.out.println("Sequence for " + Runner.FileName
								+ " does not match " + Locus + " sequence.");
					} else {
						// Make sequence context match that of locus sequence
						if (Runner.SequencePosition == 0)
							Leader = '-';
						else
							Leader = Sequence
									.charAt(Runner.SequencePosition - 1);
						if (Runner.SequencePosition + Subsequence.length() == Sequence
								.length())
							Trailer = '-';
						else
							Trailer = Sequence.charAt(Runner.SequencePosition
									+ Subsequence.length());
						Runner.Sequence = Leader
								+ "."
								+ Runner.Sequence.substring(2,
										Runner.Sequence.length() - 2) + "."
								+ Trailer;
						byte TempTryptic = Runner.Tryptic;
						Runner.DetermineTryptic(Cutoffs);
						/*
						 * if (Runner.Tryptic != TempTryptic)
						 * System.out.println(
						 * "Warning: tryptic status is actually " +
						 * Runner.Tryptic + " rather than " + TempTryptic +
						 * " for " + Runner.FileName + "!");
						 */
					}
				}
				Runner = Runner.Next;
			}
		}
	}

	public void findAllSequenceLocations(String proteinSeq,String peptideSeq, DTAFile r)
	{
		int index = r.SequencePosition;

		do
		{
			r.sequencePositionList.add(index);
			index = proteinSeq.indexOf(peptideSeq,index+1);
		}
		while(index!=-1);
	}

	/*
	 * Loop through each residue in the database sequence for this locus and
	 * accumulate molecular weight. Use a simple algorithm to estimate pI.
	 */
	public void CalculateMWAndpI(String Sequence) {
		// Get the default amino acid masses...
		ParamsFile SequestParams = new ParamsFile();
		int Length = Sequence.length();
		int Looper;
		int CountLys = 0;
		int CountArg = 0;
		int CountHis = 0;
		int CountAsp = 0;
		int CountGlu = 0;
		int CountCys = 0;
		int CountTyr = 0;
		float CurrentPH = 7.0f;
		float CurrentJump = 3.5f;
		float CurrentCharge;
		float LastCharge = 0;
		float MWAccum = 18.0f;
		char CurrentResidue;
		if (Length > 0) {
			for (Looper = 0; Looper < Length; Looper++) {
				CurrentResidue = Character.toUpperCase(Sequence.charAt(Looper));
				MWAccum += SequestParams.AvgMasses[CurrentResidue - 65];
				switch (CurrentResidue) {
				case 'C':
					CountCys++;
					break;
				case 'D':
					CountAsp++;
					break;
				case 'E':
					CountGlu++;
					break;
				case 'H':
					CountHis++;
					break;
				case 'K':
					CountLys++;
					break;
				case 'R':
					CountArg++;
					break;
				case 'Y':
					CountTyr++;
					break;
				}
			}
			/*
			 * Use a bracketing strategy to identify the isoelectric point.
			 * Calculate charge at pH of 7, and then move up 3.5 or down 3.5 as
			 * necessary. Make each successive move up or down only half as
			 * large. Keep going until two successive charges reported match to
			 * one place past the decimal point.
			 */
			CurrentCharge = Protein.ChargeAtPH(CurrentPH, CountLys, CountArg,
					CountHis, CountAsp, CountGlu, CountCys, CountTyr);
			while (RoundTo(CurrentCharge, 1) != RoundTo(LastCharge, 1)) {
				// System.out.println("pH:\t" + new Float(CurrentPH).toString()
				// + "\tCharge\t" + new Float(CurrentCharge).toString());
				if (CurrentCharge > 0)
					CurrentPH += CurrentJump;
				else
					CurrentPH -= CurrentJump;
				CurrentJump /= 2;
				LastCharge = CurrentCharge;
				CurrentCharge = Protein.ChargeAtPH(CurrentPH, CountLys,
						CountArg, CountHis, CountAsp, CountGlu, CountCys,
						CountTyr);
				if ((CurrentPH > 14) || (CurrentPH < 0)) {
					System.out.println("pI can't be figured for " + Locus);
					System.exit(0);
				}
			}
			pI = CurrentPH;
			MolWt = MWAccum;
		}
	}

	public int DTACount() {
		// Count DTAs belonging to this Protein
		int Count = 0;
		DTAFile Runner = this.DTAs.Next;
		while (Runner != null) {
			Count++;
			Runner = Runner.Next;
		}
		return Count;
	}

	public void DTASetDecoy() {
		// Set IsDecoy value for DTAs belonging to this protein
		DTAFile Runner = this.DTAs.Next;
		while (Runner != null) {
			Runner.IsDecoy = this.IsDecoy;
			Runner = Runner.Next;
		}
	}

	/*
	 * LIST MANAGEMENT FUNCTIONS These functions are those called on the header
	 * node of the Protein list. They always loop through the entire linked list
	 * of Proteins, performing some operation on each.
	 */

	/*
	 * Create one protein with all the peptides from all proteins in this list.
	 * This is useful for getting a nonredundant count of peptides and proteins.
	 * NOTE! This destroys the current list of proteins!!! Do not use if further
	 * printing is planned!
	 */
	public Protein CreateMegaProtein() {
		Protein Mega = new Protein();
		DTAFile DRunner;
		Protein PRunner = this.Next;
		while (PRunner != null) {
			DRunner = PRunner.DTAs;
			while (DRunner.Next != null)
				DRunner = DRunner.Next;
			DRunner.Next = Mega.DTAs.Next;
			Mega.DTAs.Next = PRunner.DTAs.Next;
			PRunner = PRunner.Next;
		}
		return Mega;
	}

	public void RemoveRedundantDTAs() {
		DTAFile DTARunner = this.DTAs.Next;

		while (DTARunner != null) {
			while (DTARunner.Next != null
					&& DTARunner.Next.FileName.equals(DTARunner.FileName)) {
				DTARunner.Next = DTARunner.Next.Next;
			}
			DTARunner = DTARunner.Next;
		}
	}

	/*
	 * Produce chromatography information from a sorted mega protein object.
	 */
	public void PrintChromaFile(String FileName) {
		DTAFile DatStart = this.DTAs.Next;
		DTAFile DatFinish;
		DTAFile InsideRunner;
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		String CurrentPrefix = "";
		String Buffer;
		StringTokenizer Parser;
		int TotalIDs;
		int SumTotalIDs = 0;
		int Totals[] = new int[4];
		int SumTotals[] = new int[4];
		int BinSize = 100;
		int LastBin = 0;
		int LastCount = 0;
		int CurrentBin;
		try {
			OutputFile = new File(CurrentDirectory, FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			Outgoing.write("DAT file\tTotal\t+1s\t+2s\t+3s\t\tHistogram (each is 100 scans)\n");
			SumTotals[1] = 0;
			SumTotals[2] = 0;
			SumTotals[3] = 0;
			while (DatStart != null) {
				/*
				 * Assume DatStart points at the first spectrum ID of a set
				 * created from a single .dat file.
				 */
				Parser = new StringTokenizer(DatStart.FileName, ".");
				CurrentPrefix = Parser.nextToken();
				// Initialize ID counts
				Totals[1] = 0;
				Totals[2] = 0;
				Totals[3] = 0;
				TotalIDs = 0;
				DatFinish = DatStart;
				Buffer = CurrentPrefix;
				// Until we hit the end of the list or the next .dat,
				while ((DatFinish != null) && (Buffer.equals(CurrentPrefix))) {
					if (DatFinish.ChargeState < 4) {
						// Bump the appropriate charge state
						Totals[DatFinish.ChargeState]++;
					}
					// Bump the counter of all Zs
					TotalIDs++;
					DatFinish = DatFinish.Next;
					if (DatFinish != null) {
						Parser = new StringTokenizer(DatFinish.FileName, ".");
						Buffer = Parser.nextToken();
					}
				}
				Outgoing.write(CurrentPrefix + "\t"
						+ new Integer(TotalIDs).toString() + "\t"
						+ new Integer(Totals[1]).toString() + "\t"
						+ new Integer(Totals[2]).toString() + "\t"
						+ new Integer(Totals[3]).toString() + "\t\t");
				SumTotalIDs += TotalIDs;
				SumTotals[1] += Totals[1];
				SumTotals[2] += Totals[2];
				SumTotals[3] += Totals[3];
				/*
				 * InsideRunner = DatStart; LastBin = 0; LastCount = 0; while
				 * (InsideRunner != DatFinish) { Parser = new
				 * StringTokenizer(InsideRunner.FileName,"."); CurrentPrefix =
				 * Parser.nextToken(); CurrentPrefix = Parser.nextToken();
				 * CurrentBin = new Integer(CurrentPrefix).intValue() / BinSize;
				 * while (LastBin != CurrentBin) { Outgoing.write(new
				 * Integer(LastCount).toString() + "\t"); LastBin++; LastCount =
				 * 0; } LastCount++; InsideRunner = InsideRunner.Next; }
				 * Outgoing.write(new Integer(LastCount).toString() + "\n");
				 */
				DatStart = DatFinish;
				Outgoing.write("\n");
			}
			Outgoing.write("Totals\t" + SumTotalIDs + "\t" + SumTotals[1]
					+ "\t" + SumTotals[2] + "\t" + SumTotals[3] + "\n");
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("IO Error while writing chromatography file");
			System.out.println(failure);
		}
	}

	/* Count Proteins in this list */
	public int CountProteins() {
		Protein Runner = this.Next;
		int Count = 0;
		while (Runner != null) {
			Count++;
			Runner = Runner.Next;
		}
		return Count;
	}

	/*
	 * For each protein in the list, examine the other proteins to determine
	 * whether their peptides are supersets of the current protein's peptides.
	 * If so, remove this protein. NOTE: Run after proteins have been grouped
	 * for identical sets of peptides. Note that this function compares peptides
	 * by sequence and charge state only. If peptides have different sequence
	 * context information, they're considered different peptides even if they
	 * have the same observed sequence.
	 */
	public void RemoveSubsetsQuick()
	{
		Protein ThisRunner = this.Next;
		Protein OtherRunner;
		DTAFile ThisDTA;

		// For each protein in the list,
		System.out.println("\tRemoving subset proteins...");
		System.out.println("\tRemoved\tMatched");
		Map<String,List<Protein>> peptideProteinSetMap = new HashMap<>();

		while ((ThisRunner != null)) {

			if(!ThisRunner.toRemove)
			{
			/*	if(ThisRunner.Locus.equals("gi|9910588|ref|NP_064327.1|"))
				{
					System.out.println("DEBUG<<>>< : "+ThisRunner.toRemove);
				}*/

				ThisDTA = ThisRunner.DTAs.Next;
				Map<Integer,Protein> proteinIDMap = new HashMap<>();
				while ((ThisDTA != null)) {
					int peptideId = ThisDTA.id;
					String filename = ThisDTA.FileName;

					List<Protein> proteinList = peptideProteinSetMap.getOrDefault(filename,new ArrayList<>());
					for(Protein protein : proteinList)
					{
						if(!protein.toRemove && protein.id != ThisRunner.id)
						{
							proteinIDMap.put(protein.id, protein);
						}
					}
					proteinList.add(ThisRunner);
					peptideProteinSetMap.put(filename, proteinList);
					ThisDTA = ThisDTA.Next;
				}
				for(Map.Entry<Integer,Protein> entry: proteinIDMap.entrySet())
				{
					OtherRunner = entry.getValue();
					if(OtherRunner.fileNameSet == null)
					{
						OtherRunner.initFileNameSet();
					}
					if(ThisRunner.fileNameSet == null)
					{
						ThisRunner.initFileNameSet();
					}
					Set<String> currentSet = ThisRunner.fileNameSet;
					Set<String> otherSet = OtherRunner.fileNameSet;
					boolean isCurrentSuperSet =  currentSet.containsAll(otherSet);
					boolean isOtherSuperSet  = otherSet.containsAll(currentSet);
					if(isCurrentSuperSet && isOtherSuperSet)
					{
						OtherRunner.toRemove = true;
						OtherRunner.fileNameSet = null;
					}
					else if(isCurrentSuperSet)
					{
						OtherRunner.toRemove = true;
						OtherRunner.fileNameSet = null;
					}
					else if(isOtherSuperSet)
					{
						ThisRunner.toRemove = true;
						ThisRunner.fileNameSet = null;
					}
				}
			}
			ThisRunner = ThisRunner.Next;
		}
	//	System.out.println("DEBUG<<>>> "+count);

		ThisRunner = this;
		boolean IsSuperset = false;
		int RemovedTotal = 0;

		while ((ThisRunner != null) && (ThisRunner.Next != null)) {
			// Start "other" runner at head of list
			IsSuperset = ThisRunner.Next.toRemove;
	/*		if(ThisRunner.Next.Locus.equals("gi|9910588|ref|NP_064327.1|"))
			{
				System.out.println("DEBUG<<>>< Removing: "+ThisRunner.Next.toRemove);
			}*/
			ThisRunner.Next.fileNameSet = null;

			if (IsSuperset) {
				ThisRunner.Next.toRemove = false;
				ThisRunner.Next = ThisRunner.Next.Next;
				IsSuperset = false;
				RemovedTotal++;
			} else {
				ThisRunner = ThisRunner.Next;
			}
		}

		System.out.println("\tRemoved " + new Integer(RemovedTotal).toString()
				+ " subset proteins.");
	}



	public void RemoveSubsets() {
		// Start the runner for "this" protein at head of list
		Protein ThisRunner = this;
		Protein OtherRunner;
		DTAFile ThisDTA;
		DTAFile OtherDTA;
		boolean IsSuperset = false;
		int RemovedTotal = 0;
		// For each protein in the list,
		System.out.println("\tRemoving subset proteins...");
		System.out.println("\tRemoved\tMatched");
		while ((ThisRunner != null) && (ThisRunner.Next != null)) {
			// Start "other" runner at head of list
			OtherRunner = this.Next;
			// For each other protein in the list,
			while ((OtherRunner != null) && (!IsSuperset)) {
				if (OtherRunner != ThisRunner.Next) {
					// Until further notice, this IS a duplicate
					IsSuperset = true;
					// Start with the first peptide of this protein
					ThisDTA = ThisRunner.Next.DTAs.Next;
					/*
					 * Until there are no peptides left to compare or the
					 * proteins are proven incomparable,
					 */
					while ((ThisDTA != null) && (IsSuperset)) {
						// Start with the first peptide of the other protein
						OtherDTA = OtherRunner.DTAs.Next;
						// Find the current peptide in the other protein
						while ((OtherDTA != null)
								&& (!ThisDTA.FileName.equals(OtherDTA.FileName))) {
							OtherDTA = OtherDTA.Next;
						}
						// If we didn't find it,
						if (OtherDTA == null) {
							IsSuperset = false;
						}
						ThisDTA = ThisDTA.Next;
					}
					if (IsSuperset) {
				//		System.out.println("\t" + ThisRunner.Next.Locus + "\t" + ThisRunner.Next.NSpectra + "\t"
				//				+ OtherRunner.Locus + "\t" + OtherRunner.NSpectra);
						RemovedTotal++;
					}
				}
				OtherRunner = OtherRunner.Next;
			}
			if (IsSuperset) {
				ThisRunner.Next = ThisRunner.Next.Next;
				IsSuperset = false;
			} else {
				ThisRunner = ThisRunner.Next;
			}
		}
		System.out.println("\tRemoved " + new Integer(RemovedTotal).toString()
				+ " subset proteins.");
	}

	/*
	 * Grab every protein and peptide from Other list and add it to this list.
	 * Use a basic cosequential algorithm.
	 */
	public void MergeAllFrom(Protein Other) {
		Protein ThisRunner = this;
		Protein OtherRunner = Other.Next;
		Protein Buffer;
		while ((ThisRunner.Next != null) && (OtherRunner != null)) {
			if (ThisRunner.Next.Locus.compareTo(OtherRunner.Locus) > 0) {
				/*
				 * Other's current protein needs to be inserted before
				 * this.Next.
				 */
				Buffer = ThisRunner.Next;
				ThisRunner.Next = OtherRunner.Clone();
				ThisRunner = ThisRunner.Next;
				ThisRunner.Next = Buffer;
				OtherRunner = OtherRunner.Next;
			} else if (ThisRunner.Next.Locus.compareTo(OtherRunner.Locus) < 0) {
				/*
				 * Other's current protein falls after this.Next. Advance to
				 * next in this list and repeat comparison.
				 */
				ThisRunner = ThisRunner.Next;
			} else if (ThisRunner.Next.Locus.compareTo(OtherRunner.Locus) == 0) {
				/*
				 * Other's current protein is the same as this.Next. Append in
				 * the peptides in sorted order.
				 */
				ThisRunner.Next.AccumulateDTAsFrom(OtherRunner);
				OtherRunner = OtherRunner.Next;
			}
		}
		/*
		 * If there are still proteins to be found in the Other list, append
		 * those to the end of this list.
		 */
		while (OtherRunner != null) {
			ThisRunner.Next = OtherRunner.Clone();
			ThisRunner = ThisRunner.Next;
			OtherRunner = OtherRunner.Next;
		}
	}

	/*
	 * Return the longest protein name in this list. Used for DTASelectGUI.
	 */
	public String LongestName() {
		Protein PRunner = this.Next;
		int MaxLength = 0;
		String Longest = "";
		while (PRunner != null) {
			if (PRunner.Locus.length() > MaxLength) {
				MaxLength = PRunner.Locus.length();
				Longest = PRunner.Locus;
			}
			PRunner = PRunner.Next;
		}
		return Longest;
	}

	/*
	 * Return the longest DTA file name in this list. Used for DTASelectGUI.
	 */
	public String LongestFileName() {
		Protein PRunner = this.Next;
		DTAFile DRunner;
		int MaxLength = 0;
		String Longest = "";
		while (PRunner != null) {
			DRunner = PRunner.DTAs.Next;
			while (DRunner != null) {
				if (DRunner.FileName.length() > MaxLength) {
					MaxLength = DRunner.FileName.length();
					Longest = DRunner.FileName;
				}
				DRunner = DRunner.Next;
			}
			PRunner = PRunner.Next;
		}
		return Longest;
	}

	/*
	 * Return the first Protein object in the list matching the supplied locus
	 * name. Used for DTASelectGUI.
	 */
	public Protein FindProtein(String Probe) {
		Protein PRunner = this.Next;
		Protein PRunner2;
		while (PRunner != null) {
			PRunner2 = PRunner;
			while (PRunner2 != null) {
				if (PRunner2.Locus.equals(Probe)) {
					return PRunner2;
				}
				PRunner2 = PRunner2.IdenticalLoci;
			}
			PRunner = PRunner.Next;
		}
		return null;
	}

	/*
	 * Return the first DTA file matching the given sequence. Used for
	 * DTASelectGUI.
	 */
	public DTAFile FindDTA(String Probe) {
		Protein PRunner = this.Next;
		DTAFile DRunner;
		while (PRunner != null) {
			DRunner = PRunner.DTAs.Next;
			while (DRunner != null) {
				if (DRunner.FileName.equals(Probe))
					return DRunner;
				DRunner = DRunner.Next;
			}
			PRunner = PRunner.Next;
		}
		return null;
	}

	/*
	 * Link each protein to all others remaining with identical peptides.
	 */
	public void GroupIdenticalsOld() {
		Protein PRunner = this.Next;
		Protein IdenticalBuffer;
		Protein InsidePRunner;
		DTAFile DRunner;
		DTAFile InsideDRunner;
		boolean FailedToMatch;
		while (PRunner != null) {
			InsidePRunner = PRunner;
			while (InsidePRunner.Next != null) {
				FailedToMatch = false;
				DRunner = PRunner.DTAs.Next;
				InsideDRunner = InsidePRunner.Next.DTAs.Next;
				while (!FailedToMatch && DRunner != null) {
					if (InsideDRunner == null) {
						FailedToMatch = true;
					} else {
						if (!DRunner.FileName.equals(InsideDRunner.FileName)) {
							FailedToMatch = true;
						} else {
							DRunner = DRunner.Next;
							InsideDRunner = InsideDRunner.Next;
						}
					}
				}
				if (!FailedToMatch && InsideDRunner == null) {
					// These two loci have the same DTAs in evidence
					IdenticalBuffer = PRunner.IdenticalLoci;
					PRunner.IdenticalLoci = InsidePRunner.Next;
					InsidePRunner.Next = InsidePRunner.Next.Next;
					PRunner.IdenticalLoci.IdenticalLoci = IdenticalBuffer;
				} else {
					InsidePRunner = InsidePRunner.Next;
				}
			}
			PRunner = PRunner.Next;
		}
	}

	/*
	 * Link each protein to all others remaining with identical sequence
	 * coverage. One must run "CalculateCoverageForList()" before this is run.
	 */
	public void GroupIdenticals() {
		Protein PRunner = this.Next;
		Protein IdenticalBuffer;
		Protein InsidePRunner;
		while (PRunner != null) {
			if (PRunner.Coverage.length() > 0) {
				InsidePRunner = PRunner;
				while (InsidePRunner.Next != null) {
					if (InsidePRunner.Next.Coverage.equals(PRunner.Coverage)) {
						// These two loci have the same DTAs in evidence
						IdenticalBuffer = PRunner.IdenticalLoci;
						PRunner.IdenticalLoci = InsidePRunner.Next;
						InsidePRunner.Next = InsidePRunner.Next.Next;
						PRunner.IdenticalLoci.IdenticalLoci = IdenticalBuffer;
					} else {
						InsidePRunner = InsidePRunner.Next;
					}
				}
			}
			PRunner = PRunner.Next;
		}
	}

	/*
	 * Undo what the above function does. Put all the Proteins attached at
	 * IdenticalLoci pointers back into their proper places in the list.
	 */
	public void UngroupIdenticals() {
		Protein Runner = this.Next;
		Protein SeparateStorage = new Protein();
		Protein SSRunner;
		Protein RRunner;
		Protein TempProt;
		while (Runner != null) {
			if (Runner.IdenticalLoci != null) {
				// Move the proteins attached here to separate storage
				RRunner = Runner.IdenticalLoci;
				while (RRunner != null) {
					SSRunner = SeparateStorage;
					// Advance to proper location for addition in SS list
					while ((SSRunner.Next != null)
							&& (SSRunner.Next.Locus.compareTo(RRunner.Locus) < 0))
						SSRunner = SSRunner.Next;
					// Hook Protein into SS list
					TempProt = SSRunner.Next;
					SSRunner.Next = RRunner;
					RRunner.Next = TempProt;
					// Unlink this Protein's IdenticalLoci chain
					TempProt = RRunner.IdenticalLoci;
					RRunner.IdenticalLoci = null;
					RRunner = TempProt;
				}
				// Unlink the originating Protein's IdenticalLoci chain
				Runner.IdenticalLoci = null;
			}
			Runner = Runner.Next;
		}
		// Merge the SS list with the normal list.
		SSRunner = SeparateStorage.Next;
		Runner = this;
		while (SSRunner != null) {
			// Advance to proper merge location
			while ((Runner.Next != null)
					&& (Runner.Next.Locus.compareTo(SSRunner.Locus) < 0))
				Runner = Runner.Next;
			// Add this Protein into place
			TempProt = Runner.Next;
			Runner.Next = SSRunner;
			Runner = Runner.Next;
			// Advance SS list to next before we ditch the pointer.
			SSRunner = SSRunner.Next;
			Runner.Next = TempProt;
		}
	}

	// Determine how many times each peptide is present for each locus
	public void CalculateRedundancyForList(SelectCriteria criteria) {
		Protein Runner = this.Next;
		while (Runner != null) {
			Runner.CalculateRedundancy(criteria);
			Runner = Runner.Next;
		}
	}

	/*
	 * Print a series of HTML links to other occurrences of this peptide within
	 * the list
	 */
	public String PrintLinksForPeptide(String Probe) {
		Protein PRunner = this.Next;
		DTAFile DRunner;
		StringBuffer Returned = new StringBuffer("");
		String[] ParsedLine = Probe.split("\\.");
		String Trunc = ParsedLine[0] + "." + ParsedLine[1] + "."
				+ ParsedLine[2];
		while (PRunner != null) {
			DRunner = PRunner.DTAs.Next;
			while ((DRunner != null) && (!DRunner.FileName.startsWith(Trunc)))
				DRunner = DRunner.Next;
			if (DRunner != null) {
				String[] ParsedFileName = DRunner.FileName.split("\\.");
				Returned.append("<a href=\"#" + PRunner.Locus + "\">"
						+ ParsedFileName[3] + "</a>");
			}
			PRunner = PRunner.Next;
		}
		return Returned.toString();
	}

	/*
	 * Determine sequence coverage for each Protein in the list. Store both the
	 * percentage and the list of peptides (after merging them together).
	 */
	public void CalculateCoverageForList(boolean UseNewCGI) {
		Protein OuterRunner = this.Next;
		Protein Runner;
		CoverageZone Buffer;
		CoverageZone BRunner;
		float SCAccumulator;
		while (OuterRunner != null) {
			Runner = OuterRunner;
			while (Runner != null) {
				SCAccumulator = 0f;
				Buffer = Runner.GenerateZones();
				if (UseNewCGI) {
					Runner.Coverage = Buffer.GetSeqCovString();
					Buffer.MakeMinimal();
				} else {
					Buffer.MakeMinimal();
					Runner.Coverage = Buffer.GetConsensusList();
				}
				/*
				 * Sum together the remaining CoverageZones' spans and divide by
				 * the full protein length to determine percentage of coverage.
				 */
				BRunner = Buffer.Next;
				while (BRunner != null) {
					SCAccumulator += BRunner.Finish - BRunner.Start + 1;
					BRunner = BRunner.Next;
				}
				Runner.SequenceCoverage = Protein.RoundTo(
						(100.0f * SCAccumulator / Runner.SequenceLength), 1);
				Runner = Runner.IdenticalLoci;
			}
			OuterRunner = OuterRunner.Next;
		}
	}

	/*
	 * Open the FASTA database indicated in sequest.params. Locate each
	 * identified locus in that file, and store the sequence found with each
	 * Protein object. In addition, store the descriptive name for each locus in
	 * the Protein object. FASTA databases are not necessarily sorted by locus!
	 */
	public void LookUpLoci(ParamsFile SEQUESTParams, IniFile Config,
			SelectCriteria Cutoffs) throws IOException {
		// Runner is a pointer into our list of identified proteins
		String DBName = SEQUESTParams.DBName;
		Protein Runner = this.Next;
		File DBFile = new File(DBName);
		FileReader InputFileReader = new FileReader(DBFile);
		BufferedReader Incoming = new BufferedReader(InputFileReader);
		int LongestLocus = Config.LocusLengthCutoff;
		String LineBuffer = Incoming.readLine();
		StringTokenizer Parser;
		StringBuffer SequenceSB;
		String Sequence;
		
		StringBuffer fename;
		
		// Do only if there are any .out results in memory
		if (this.Next == null) {
			System.out.println("No .out files found!");
			System.out
					.println("If .out files are in the current directory, use the --here or -. options");
			System.exit(0);
		} else {
			// Queue up the DB file to the first locus line

			Hashtable<String, FastaEntry> dbtable = new Hashtable<String, FastaEntry>();

			while ((LineBuffer.length() == 0) || (LineBuffer.charAt(0) != '>')) {
				LineBuffer = Incoming.readLine();
			}
			// Until we hit the end of the file
			while ((LineBuffer != null)) {

				Parser = new StringTokenizer(LineBuffer, " \t>");
				if (!Parser.hasMoreTokens()) {
					LineBuffer = Incoming.readLine();
				} else {
					LineBuffer = Parser.nextToken();
				}
				/*
				 * Deal with the fact that SEQUEST .output files list only
				 * first 40 characters of locus names. Binary searches are
				 * even shorter names!
				 */
				if (LineBuffer.length() > LongestLocus) {
					LineBuffer = LineBuffer.substring(0, LongestLocus);
				}
				if (Cutoffs.ShortForm) {
					LineBuffer = NRUtils.getAccession(LineBuffer);
				}
				if (!Cutoffs.Quiet) {
					System.out.print("\t" + LineBuffer + "      \r");
				}
				FastaEntry fe = new FastaEntry(LineBuffer);
				fename = new StringBuffer();
				while (Parser.hasMoreTokens()) {
					fename.append(Parser.nextToken() + " ");
				}
				fe.setName(fename.toString());
				SequenceSB = new StringBuffer();
				LineBuffer = Incoming.readLine();
				while ((LineBuffer != null)
						&& ((LineBuffer.length() == 0) || (LineBuffer
								.charAt(0) != '>'))) {
					SequenceSB.append(LineBuffer);
					LineBuffer = Incoming.readLine();
				}
				// Strip out extraneous characters from sequence
				Sequence = DTAFile.JustLettersFrom(SequenceSB
						.toString());
				fe.setSequence(Sequence);
				dbtable.put(fe.getLocus(), fe);
				
			}

			Runner = this.Next;
			while (Runner != null) {
				
				if (dbtable.containsKey(Runner.Locus)) {
					
					Runner.Gene = dbtable.get(Runner.Locus).getName();
					
					Sequence = dbtable.get(Runner.Locus).getSequence();
					
					Runner.SequenceLength = Sequence.length();
					Runner.NChar = Sequence.charAt(0);
					Runner.CChar = Sequence
							.charAt(Runner.SequenceLength - 1);
					/*
					 * Determine protein's sequence coverage and calculate
					 * molecular weight and pI.
					 */
					Runner.FindDTAPositions(Sequence, Runner.Locus, Cutoffs);
					Runner.CalculateMWAndpI(Sequence);
					
					Runner = Runner.Next;
					
				} else {
					if (!Cutoffs.Quiet) {
						System.out.print("Locus " + Runner.Locus + " not in database");
					}
					Runner = Runner.Next;
				}
				
			}
			
			dbtable.clear();
			dbtable = null;
			
			/*while ((LineBuffer.length() == 0) || (LineBuffer.charAt(0) != '>'))
				LineBuffer = Incoming.readLine();
			// Until we hit the end of the file
			while ((LineBuffer != null)) {
				Parser = new StringTokenizer(LineBuffer, " \t>");
				if (!Parser.hasMoreTokens()) {
					LineBuffer = Incoming.readLine();
				} else {
					LineBuffer = Parser.nextToken();
					
					 * Deal with the fact that SEQUEST .output files list only
					 * first 40 characters of locus names. Binary searches are
					 * even shorter names!
					 
					if (LineBuffer.length() > LongestLocus) {
						LineBuffer = LineBuffer.substring(0, LongestLocus);
					}
					if (Cutoffs.ShortForm) {
						LineBuffer = NRUtils.getAccession(LineBuffer);
					}
					if (!Cutoffs.Quiet) {
						System.out.print("\t" + LineBuffer + "      \r");
					}
					// Start at the head of the list
					Runner = this.Next;
					// Move our identified locus pointer to proper place.
					while ((Runner != null)
							&& (LineBuffer.compareTo(Runner.Locus) > 0))
						Runner = Runner.Next;
					if ((Runner != null)
							&& (LineBuffer.startsWith(Runner.Locus))) {
						
						 * This DB locus corresponds to the current identified
						 * locus. Read gene name and sequence
						 
						while (Parser.hasMoreTokens()) {
							Runner.Gene += Parser.nextToken() + " ";
						}
						SequenceSB = new StringBuffer();
						LineBuffer = Incoming.readLine();
						while ((LineBuffer != null)
								&& ((LineBuffer.length() == 0) || (LineBuffer
										.charAt(0) != '>'))) {
							SequenceSB.append(LineBuffer);
							LineBuffer = Incoming.readLine();
						}
						// Strip out extraneous characters from sequence
						Sequence = DTAFile.JustLettersFrom(SequenceSB
								.toString());
						Runner.SequenceLength = Sequence.length();
						Runner.NChar = Sequence.charAt(0);
						Runner.CChar = Sequence
								.charAt(Runner.SequenceLength - 1);
						
						 * Determine protein's sequence coverage and calculate
						 * molecular weight and pI.
						 
						Runner.FindDTAPositions(Sequence, Runner.Locus, Cutoffs);
						Runner.CalculateMWAndpI(Sequence);
					} else {
						// This DB sequence matches nothing.
						// Advance to next database sequence
						LineBuffer = Incoming.readLine();
						while ((LineBuffer != null)
								&& ((LineBuffer.length() == 0) || (LineBuffer
										.charAt(0) != '>')))
							LineBuffer = Incoming.readLine();
					}*/
			//	}
		//	}
		}
	}

	/*
	 * Run through the entire Protein list. Run through each DTAFile for each
	 * Protein. Eliminate DTAFiles that don't meet the criteria we're using.
	 */
	public void DumpUnqualifiedDTAs(SelectCriteria Cutoffs) {
		Protein Runner = this.Next;
		DTAFile DTARunner;
		while (Runner != null) {
			DTARunner = Runner.DTAs;
			while (DTARunner.Next != null) {
				if (Cutoffs.Allow(DTARunner.Next))
					DTARunner = DTARunner.Next;
				else
					DTARunner.Next = DTARunner.Next.Next;
			}
			if (Cutoffs.UseProteinFilters) {
				Runner.HasGreatPeptide = Runner.PassScores(Cutoffs);
			}
			Runner = Runner.Next;
		}
	}

	public static double getMedian(List<Double> list)
	{
		Collections.sort(list);
		int mid = list.size()/2;
		double median;
		if(list.size()%2 ==0)
		{
			median = (list.get(mid) + list.get(mid-1))/2;
		}
		else
		{
			median = list.get(mid);
		}
		return median;
	}

	public void CalculateProteinFalsePositive()
	{
		Protein Runner = this.Next;
		DTAFile DTARunner;
		List<DTAFile> dtaList = new ArrayList<>();
		boolean isDecoy = true;
		boolean containPeptides = false;
		while (Runner != null) {
			DTARunner = Runner.DTAs.Next;

			if(!Runner.Locus.startsWith("Reverse_"))
			{
				isDecoy =false;
			}
			if(DTARunner !=null)
			{
				containPeptides = true;
			}
			while (DTARunner != null) {
				dtaList.add(DTARunner);
				DTARunner.IsDecoy = isDecoy;
				DTARunner = DTARunner.Next;
			}
			if(containPeptides)
			{
				isDecoy = true;
				containPeptides = false;
			}

			Runner = Runner.Next;
		}
		Collections.sort(dtaList, Comparator.comparing(DTAFile::getPepFP));
		int forwardCount =0;
		int reverseCount =0;
		for(DTAFile dta: dtaList)
		{
			if(dta.IsDecoy)
			{
				reverseCount++;
			}
			else
			{
				forwardCount++;
			}

		}

	}



	public void CalculateFilterMedianAdjustedDeltaMass(SelectCriteria cutoffs)
	{
		Protein Runner = this.Next;
		DTAFile DTARunner;
		List<Double> ppmList = new ArrayList<>();
		while (Runner != null) {
			DTARunner = Runner.DTAs.Next;
			while (DTARunner != null) {
				double ppm = DTARunner.Adjusted_PPM_Offset;
				ppmList.add(ppm);
				DTARunner = DTARunner.Next;
			}
			Runner = Runner.Next;
		}
		double median = getMedian(ppmList);
	//	System.out.println("<<>><>  "+median);
		Runner = this.Next;
		DTARunner =null;
		while (Runner != null) {
			DTARunner = Runner.DTAs;
			while (DTARunner.Next != null) {
				DTARunner.Next.Shifted_PPM_Offset = DTARunner.Next.Adjusted_PPM_Offset - median;
				if(cutoffs.AllowShiftDM(DTARunner.Next))
				{
					DTARunner = DTARunner.Next;
				}
				else
				{
					DTARunner.Next = DTARunner.Next.Next;
				}

			}
			if(cutoffs.UseProteinFilters)
			{
				Runner.HasGreatPeptide = Runner.PassScores(cutoffs, true);
			}
			Runner = Runner.Next;
		}

	}

	/*
	 * Run through the list of Proteins. Count how many DTAFiles are associated
	 * with each locus. Remove Proteins with insufficient numbers of supporting
	 * DTAFiles.
	 */
	public void DitchProteinsWithoutSufficientDTAs(SelectCriteria Cutoffs) {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this;
		boolean Qualifies;
		while (Runner.Next != null) {
			// Allow loci in if they satisfy either allowance rule
			Qualifies = (Runner.Next.NPeptides >= Cutoffs.MinPepsPerLocus)
					&& (Runner.Next.NSpectra >= Cutoffs.MinSpectraPerLocus)
					&& (Runner.Next.NTrypticPeptides >= Cutoffs.MinTrypticPeps)
					&& (Runner.Next.HasGreatPeptide);
			// Apply the minimum modified peptides per locus rule
			if (Cutoffs.MinModPepsPerLocus > 0) {
				Qualifies &= (Runner.Next.NModPeptides >= Cutoffs.MinModPepsPerLocus);
			}
			// Remove loci if they fail required uniqueness
			if (Cutoffs.IncludeOnlyUniques)
				Qualifies &= Runner.Next.HasAtLeastOneUnique();
			if (!Qualifies)
				Runner.Next = Runner.Next.Next;
			else
				Runner = Runner.Next;
		}
	}

	public void DitchProteinsWithLowSequenceCoverage(SelectCriteria Cutoffs) {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this;
		boolean Qualifies;
		while (Runner.Next != null) {
			Qualifies = (Runner.Next.SequenceCoverage >= Cutoffs.MinSequenceCoverage);
			if (!Qualifies)
				Runner.Next = Runner.Next.Next;
			else
				Runner = Runner.Next;
		}
	}

	/*
	 * Run through the list of Proteins. Check DTAFiles for duplicate sequences.
	 * Destroy all but one of the DTAFiles for any given sequence, retaining the
	 * one with the highest XCorr.
	 */
	public void DitchDuplicateDTAsByXCorr() {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this.Next;
		while (Runner != null) {
			Runner.DTAs.KeepOnlyDTAWithHighXCorr();
			Runner = Runner.Next;
		}
	}

	/*
	 * Run through the list of Proteins. Check DTAFiles for duplicate sequences.
	 * Destroy all but one of the DTAFiles for any given sequence, retaining the
	 * one with the highest TIC.
	 */
	public void DitchDuplicateDTAsBySaltStep() {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this.Next;
		while (Runner != null) {
			Runner.DTAs.KeepOnlyDTAInEachSaltStep();
			Runner = Runner.Next;
		}
	}

	public void DebugPrint() {
		Protein PRunner = this.Next;
		while (PRunner != null) {
			System.out.println(PRunner.Locus + "\t" + PRunner.Classification
					+ "\t" + PRunner.MolWt);
			PRunner = PRunner.Next;
		}
	}

	/*
	 * Sort proteins into groups by classification
	 */
	public void SortByClass() {
		if (this.Next != null) {
			this.Next = this.Next.Sort(null);
		}
	}

	/*
	 * Sort by classification, then by sequence coverage
	 */
	private Protein Sort(Protein Follower) {
		// Recursive quicksorter
		// Returns first item of sorted list (from those starting at this
		// element)
		// Accepts first item to follow
		Protein ListAbove = null;
		Protein ListBelow = null;
		Protein PlaceHolder;
		Protein PlaceHolder2;
		PlaceHolder = this.Next;
		// Partition all remaining points of this linked list
		while (PlaceHolder != null) {
			PlaceHolder2 = PlaceHolder.Next;
			if (this.Classification > PlaceHolder.Classification) {
				// Move this item to list above this
				PlaceHolder.Next = ListAbove;
				ListAbove = PlaceHolder;
			} else if (this.Classification < PlaceHolder.Classification) {
				// Move this item to list below this point
				PlaceHolder.Next = ListBelow;
				ListBelow = PlaceHolder;
			} else {
				if (this.SequenceCoverage < PlaceHolder.SequenceCoverage) {
					// Move this item to list above this
					PlaceHolder.Next = ListAbove;
					ListAbove = PlaceHolder;
				} else if (this.SequenceCoverage > PlaceHolder.SequenceCoverage) {
					// Move this item to list below this point
					PlaceHolder.Next = ListBelow;
					ListBelow = PlaceHolder;
				} else {
					if (this.MolWt < PlaceHolder.MolWt) {
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

	/*
	 * Run through the list of Proteins. Sort the DTAFiles for each Protein on
	 * the basis of sequence. Use the nifty quicksorter in DTAFile for this
	 * process.
	 */
	public void SortDTALists() {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this.Next;
		DTAFile DRunner;
		while (Runner != null) {
			// Check each peptide's modification state
			DRunner = Runner.DTAs.Next;
			while (DRunner != null) {
				DRunner.DetermineModified();
				DRunner = DRunner.Next;
			}
			// Now sort the peptides
			Runner.DTAs.SortList();
			Runner = Runner.Next;
		}
	}

	/*
	 * Run through the list of Proteins. Sort the DTAFiles for each Protein on
	 * the basis of sequence. Use the nifty quicksorter in DTAFile for this
	 * process.
	 */
	public void DisplaySortDTALists() {
		// this is just an empty header; the real data starts at the next item
		Protein Runner = this.Next;
		while (Runner != null) {
			Runner.DTAs.DisplaySortList();
			Runner = Runner.Next;
		}
	}

	/*
	 * Print the list of DTAFile filenames as a series of cp commands. We don't
	 * care which directory the files come from. Nor do we bother to write these
	 * to a file but instead splat them out to the screen. We're lazy. That's
	 * why.
	 * 
	 * (9/26/02) Well, clearly I cleaned it up later on. Now it writes a
	 * platform-specific script to do the copying. Progress!
	 */
	public void PrintCPList(IniFile Config) {
		try {
			File OutputFile;
			String CurrentDir = System.getProperty("user.dir");
			if (Config.ServerType.equals("UNIX")) {
				OutputFile = new File(CurrentDir, "copylist.bash");
			} else {
				OutputFile = new File(CurrentDir, "copylist.bat");
			}
			FileWriter OutputFileWriter = new FileWriter(OutputFile);
			BufferedWriter Outgoing = new BufferedWriter(OutputFileWriter);
			Protein Runner = this.Next;
			DTAFile DTAFileRunner;
			while (Runner != null) {
				DTAFileRunner = Runner.DTAs.Next;
				if (Config.ServerType.equals("UNIX")) {
					while (DTAFileRunner != null) {
						Outgoing.write("cp "
								+ DTAFileRunner.CanonicalName(CurrentDir)
								+ " ~/temp/\n");
						Outgoing.write("cp "
								+ DTAFileRunner.DTACanonicalName(CurrentDir)
								+ " ~/temp/\n");
						DTAFileRunner = DTAFileRunner.Next;
					}
				} else {
					while (DTAFileRunner != null) {
						Outgoing.write("copy "
								+ DTAFileRunner.CanonicalName(CurrentDir)
								+ " c:\\temp\\\n");
						Outgoing.write("copy "
								+ DTAFileRunner.DTACanonicalName(CurrentDir)
								+ " c:\\temp\\\n");
						DTAFileRunner = DTAFileRunner.Next;
					}
				}
				Runner = Runner.Next;
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Error writing script to copy spectra");
			System.out.println(failure);
		}
	}

	/* Test to ensure that saving to this directory is possible */
	public static boolean PrintTest() {
		String StringBuffer = System.getProperty("user.dir");
		File CurrentDirectory = new File(StringBuffer);
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			OutputFile = new File(CurrentDirectory, "DTASelect.html");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			return true;
		} catch (IOException failure) {
			return false;
		}
	}

	public void calculateNSAFandEMPAI(float normalizer) {

		float NSAF_sum = normalizer;

		if (NSAF_sum == 0) {
			NSAF_sum = 1;
		}

		this.NSAF = this.SAF / NSAF_sum;

		this.EMPAI = (float) (Math.pow(10, (this.SequenceCoverage / 100))) - 1;

	}

	public void PrintDatabase(SelectCriteria Cutoffs, ParamsFile Params,
			String RootName) {
		Protein Runner = this.Next;
		Protein MatchRunner;
		File CurrentDirectory = new File(System.getProperty("user.dir"));
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		String tab = "\t";

		// Added by Diego Calzolari, 11/05/2012
		float NSAF_sum = 0f;
		// ------------------------------------

		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			OutputFile = new File(CurrentDirectory, RootName + "-Proteins.txt");
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			// modified by Diego Calzolari, 11/05/2012
			// original line
			// :Outgoing.write("LocusID\tSequenceCount\tSpectrumCount\tCoverage\tL\tMW\tpI\tDescription\tDuplicate\tClassID\n");
			Outgoing.write("LocusID\tSequenceCount\tSpectrumCount\tCoverage\tL\tMW\tpI\temPAI\tNSAF\tDescription\tDuplicate\tClassID\n");
			NSAF_sum = getNSAF_normalizer();
			// -----------------------------------------------------------------------------------------------------------------------------------
			while (Runner != null) {
				MatchRunner = Runner;
				while (MatchRunner != null) {
					Outgoing.write(MatchRunner.Locus
							+ tab
							+ MatchRunner.NPeptides
							+ tab
							+ MatchRunner.NSpectra
							+ tab
							+ MatchRunner.SequenceCoverage
							+ tab
							+ new Integer(MatchRunner.SequenceLength)
									.toString()
							+ tab
							+ new Integer(Math.round(MatchRunner.MolWt))
									.toString()
							+ tab
							+ new Float(RoundTo(MatchRunner.pI, 1)).toString()
							+ tab
							+
							// modified by Diego Calzolari, 11/05/2012
							new Float(Math.pow(10,
									(MatchRunner.SequenceCoverage / 100)) - 1)
									.toString() + tab
							+ new Float(MatchRunner.SAF / NSAF_sum).toString()
							+ tab +
							// ---------------------------------------------------------------------------------------
							MatchRunner.Gene);
					if (MatchRunner != Runner)
						Outgoing.write("\tY\t");
					else
						Outgoing.write("\tN\t");
					Outgoing.write(Runner.Classification + "\n");
					MatchRunner = MatchRunner.IdenticalLoci;
				}
				Runner = Runner.Next;
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Failed to write database files.");
			System.out.println(failure);
		}
	}

	// Added by Diego Calzolari, 11/05/2012
	public float getNSAF_normalizer() {

		if (this.SequenceLength != 0) {
			this.SAF = (float) this.NSpectra / (float) this.SequenceLength;
		} else {
			this.SAF = (float) this.NSpectra;
		}

		return this.SAF;
	}

	// ------------------------------------------------------------------------------------------

	public void PrintXML(SelectCriteria Cutoffs, String CommandLineOptions,
			ParamsFile Params, String FileName, String DirName) {
		String StringBuffer = System.getProperty("user.dir");
		String Descrip;
		File CurrentDirectory = new File(StringBuffer);
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		Protein Runner = this.Next;
		Protein MatchRunner;
		DTAFile DTARunner;
		int Counter = 0;
		int PepCounter = 0;
		try {
			CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
			OutputFile = new File(CurrentDirectory, FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			// Write header
			Outgoing.write("<?xml version=\"1.0\" standalone=\"no\"?>\n");
			Outgoing.write("<!DOCTYPE ResultFile SYSTEM \"DTASelect1.7.dtd\">\n");
			Outgoing.write("<DTASelect:ResultFile xmlns:DTASelect=\"http://fields.scripps.edu/DTASelect\">");
			Outgoing.write("\n<DTASelect:Header>\n");
			Outgoing.write("\t<DTASelect:Version>" + Protein.Version()
					+ "</DTASelect:Version>\n");
			Outgoing.write("\t<DTASelect:Directory>" + DirName
					+ "</DTASelect:Directory>\n");
			Outgoing.write("\t<DTASelect:Database>" + Params.DBName
					+ "</DTASelect:Database>\n");
			Outgoing.write("\t<DTASelect:CommandLine>" + CommandLineOptions
					+ "</DTASelect:CommandLine>\n");
			// Write the criteria in effect for this run
			Outgoing.write(Cutoffs.PrintXMLCriteria());
			Outgoing.write("</DTASelect:Header>\n");
			// Initialize our counter of loci
			Counter = 0;
			// Print off all fields, continuing down list
			while (Runner != null) {
				Counter++;
				MatchRunner = Runner;
				Outgoing.write("<DTASelect:Protein>\n");
				while (MatchRunner != null) {
					// Write the Locus-specific information
					Outgoing.write("\t<DTASelect:Locus>\n"
							+ "\t\t<DTASelect:ID>" + MatchRunner.Locus
							+ "</DTASelect:ID>\n" + "\t\t<DTASelect:DTACount>"
							+ new Integer(MatchRunner.NPeptides).toString()
							+ "</DTASelect:DTACount>\n"
							+ "\t\t<DTASelect:SpectrumCount>"
							+ new Integer(MatchRunner.NSpectra).toString()
							+ "</DTASelect:SpectrumCount>\n");
					if (MatchRunner.SequenceLength == 0) {
					} else {
						Descrip = MatchRunner.Gene;
						Descrip.replace('&', '*');
						Descrip.replace('<', '*');
						Descrip.replace('>', '*');
						Descrip.replace('\"', '*');
						Descrip.replace('\'', '*');
						Outgoing.write("\t\t<DTASelect:Coverage>"
								+ MatchRunner.SequenceCoverage
								+ "</DTASelect:Coverage>\n"
								+ "\t\t<DTASelect:SeqLength>"
								+ new Integer(MatchRunner.SequenceLength)
										.toString()
								+ "</DTASelect:SeqLength>\n"
								+ "\t\t<DTASelect:ProtConf>"
								+ new Float(RoundTo(MatchRunner.ProtConf, 2))
										.toString()
								+ "</DTASelect:ProtConf>\n"
								+ "\t\t<DTASelect:ProtFP>"
								+ new Float(RoundTo(MatchRunner.ProtFP, 2))
										.toString()
								+ "</DTASelect:ProtFP>\n"
								+ "\t\t<DTASelect:MW>"
								+ new Integer(Math.round(MatchRunner.MolWt))
										.toString()
								+ "</DTASelect:MW>\n"
								+ "\t\t<DTASelect:pI>"
								+ new Float(RoundTo(MatchRunner.pI, 1))
										.toString() + "</DTASelect:pI>\n"
								+ "\t\t<DTASelect:Description>" + Descrip
								+ "</DTASelect:Description>\n");
					}
					MatchRunner = MatchRunner.IdenticalLoci;
					Outgoing.write("\t</DTASelect:Locus>\n");
				}
				DTARunner = Runner.DTAs.Next;
				while (DTARunner != null) {
					// Update our counter of DTAFiles
					PepCounter++;
					Outgoing.write("\t<DTASelect:DTAFile>\n");
					Outgoing.write(DTARunner.PrintXML(Params, DirName));
					Outgoing.write("\t</DTASelect:DTAFile>\n");
					DTARunner = DTARunner.Next;
				}
				Runner = Runner.Next;
				Outgoing.write("</DTASelect:Protein>\n");
			}
			Outgoing.write("<DTASelect:Counts>\n");
			// Print counts of loci and DTAFiles
			Outgoing.write("\t<DTASelect:ProteinCount>"
					+ new Integer(Counter).toString()
					+ "</DTASelect:ProteinCount>\n");
			Outgoing.write("\t<DTASelect:DTACounts>"
					+ new Integer(PepCounter).toString()
					+ "</DTASelect:DTACounts>\n");
			Outgoing.write("</DTASelect:Counts>\n");
			Outgoing.write("</DTASelect:ResultFile>\n");
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing HTML file");
		}
	}

	// Because I don't like the Math.roundTo function, I wrote this one.
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

	// STATIC FUNCTIONS

	// Return a string containing the version number
	public static String Version() {
		//return "DTASelect v2.0.49";
		return "DTASelect v2.1.8";
	}

	/*
	 * Determine the charge on a theoretical protein at a given pH given its
	 * relevant composition
	 */
	public static float ChargeAtPH(float pH, int CountLys, int CountArg,
			int CountHis, int CountAsp, int CountGlu, int CountCys, int CountTyr) {
		// Start out accumulator with charge of termini
		float Accum = PercentPositive(pH, 8.0f) - PercentNegative(pH, 3.1f);
		Accum += CountLys * PercentPositive(pH, 10.0f);
		Accum += CountArg * PercentPositive(pH, 12.0f);
		Accum += CountHis * PercentPositive(pH, 6.5f);
		Accum -= CountAsp * PercentNegative(pH, 4.4f);
		Accum -= CountGlu * PercentNegative(pH, 4.4f);
		Accum -= CountCys * PercentNegative(pH, 8.5f);
		Accum -= CountTyr * PercentNegative(pH, 10.0f);
		return Accum;
	}

	/*
	 * What percentage of ions of given pK at given pH would be postively
	 * charged?
	 */
	public static float PercentPositive(float pH, float pK) {
		double ConcentrationRatio = Math.pow(10f, pK - pH);
		return new Double(ConcentrationRatio / (ConcentrationRatio + 1))
				.floatValue();
	}

	/*
	 * What percentage of ions at given pK at given pH would be negatively
	 * charged?
	 */
	public static float PercentNegative(float pH, float pK) {
		double ConcentrationRatio = Math.pow(10, pH - pK);
		return new Double(ConcentrationRatio / (ConcentrationRatio + 1))
				.floatValue();
	}

	class TextBlock {
		String Text;
		int Index;
		boolean Keep = false;
		TextBlock Next;
	}

	private class FastaEntry {

		String locus;
		String name;
		String sequence;

		public FastaEntry(String locus) {
			this.locus = locus;
		}

		public void setName(String name) {
			this.name = name;
		}

		public void setSequence(String sequence) {
			this.sequence = sequence;
		}

		public String getLocus() {
			return locus;
		}

		public String getName() {
			return name;
		}

		public String getSequence() {
			return sequence;
		}

	}

	private void initFileNameSet()
	{
		DTAFile dta = this.DTAs.Next;
		fileNameSet = new HashSet<>();
		while ((dta != null)) {
			// Start with the first peptide of the other protein
			fileNameSet.add(dta.FileName);
			dta = dta.Next;
		}
	}

}
