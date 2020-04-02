import javax.swing.plaf.ColorUIResource;
import java.io.*;
import java.util.*;
import java.net.URLDecoder;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//OUTFile object
//For DTASelect
//Created by Dave Tabb
//Begun 8/22/2000

/* An OUTFile represents the information found in an .out file.  Most of
 * that information is stored in a DTAFile object to which the OUTFile
 * points.  The locus identified for this .out, however, is stored
 * here temporarily.  The list of OUTFiles will be amassed and then
 * collated to create the list of Proteins.  The Protein object is
 * better prepared to store more extensive information about each
 * locus.
 */
public class OUTFile {
	String Locus;
	String RawLocus;
	DTAFile DTA = new DTAFile();
	boolean OngoingStatisticsFlag = false;
	OUTFile Next;
	OUTFile OriginalNext;

	// Copy this OUTFile, pointing to a copy of the DTA
	public OUTFile Clone() {
		OUTFile Copy = new OUTFile();
		Copy.Locus = this.Locus;
		Copy.RawLocus = this.RawLocus;
		Copy.DTA = this.DTA.Clone();
		return Copy;
	}

	/*
	 * This is the magical function by which a list of OUTFiles becomes a list
	 * of Proteins. OUTFile lists are one-dimensional; each OUTFile points to
	 * another OUTFile, and each OUTFile points to only one DTAFile. Protein
	 * lists are two-dimensional; each Protein points to another Protein, and
	 * each Protein points to a list of DTAFiles. This function should always be
	 * called on an OUTFile list that has been sorted by locus. For each group
	 * of OUTFiles that have the same locus, the DTAFiles are collated into a
	 * list and that list is attached to a new Protein object that is added to a
	 * growing list of Proteins. To maintain the order in which the OUTFiles
	 * were sorted, the Protein objects are tail added to the growing list.
	 */
	public Protein CollateByLocus(SelectCriteria Cutoffs) {
		Protein ListToReturn = new Protein();
		Protein CurrentProtein = ListToReturn;
		DTAFile CurrentDTAFile = ListToReturn.DTAs;
		OUTFile Runner = this.Next;
		String LastLocus = "";
		while (Runner != null) {
			if (!Runner.Locus.equals(LastLocus)) {
				CurrentProtein.Next = new Protein();
				CurrentProtein = CurrentProtein.Next;
				CurrentProtein.Locus = Runner.Locus;
				CurrentDTAFile = CurrentProtein.DTAs;
				LastLocus = Runner.Locus;
			}
			CurrentDTAFile.Next = Runner.DTA;
			CurrentDTAFile = CurrentDTAFile.Next;
			Runner = Runner.Next;
		}
		return ListToReturn;
	}

	/*
	 * Read the unified SEQUEST result file. Create an OUTFile object for each
	 * peptide for each protein.
	 * 
	 * Format of .sqt file is as follows: For each S line, there should be a
	 * series of M lines. After each M line, there is one or several L lines. S
	 * lines indicate a new spectrum (or a different charge state for a
	 * previously examined spectrum). M lines indicate matching sequences for
	 * the most recent spectrum. L lines indicate database loci containing the
	 * sequence reported in the last M line.
	 * 
	 * This function separates the file into sections starting with S lines.
	 * Each M line with a deltCN of 0.0 is multiplexed with the L lines below to
	 * create OUTFile objects that are appended to the MasterList. If there's an
	 * error in one of these groups, the objects for that spectrum are removed
	 * (by keeping a pointer to the last MasterList until all is clear). Once a
	 * nonzero DeltCN is found, the series of created OUTFile objects for this
	 * spectrum are marked with the DeltCN, and the file is advanced to the next
	 * S line.
	 */
	public static int ReadSQT(File CurrentDirectory, String Filename,
			OUTFile MasterList, IniFile Config, SelectCriteria Cutoffs)
			throws IOException {
		int OUTCount = 0;
		OUTFile MasterListRunner = MasterList;
		OUTFile LastSoFar;
		File OutFile = new File(CurrentDirectory, Filename);
		FileReader InputFilereader = new FileReader(OutFile);
		BufferedReader Incoming = new BufferedReader(InputFilereader);
		String LineBuffer;
		String WholeLine;
		String TempLocus, TempRawLocus;
		StringTokenizer Parser;
		DTAFile CurrentDTA = new DTAFile();
		boolean AbortSpectrum;
		int LocusCount;
		boolean FoundMLine;
		// Storage for S line values
		float bPrecursorMass = 0.0f;
		float bTotalIntensity = 0.0f;
		int SetDeltCN = Cutoffs.SetDeltCN;
		float bDeltCN;
		float SecondDeltCN;
		short bXCorrRank;
		short SecondXCorrRank;
		String bFileName;
		String bSubdirectory = Filename.substring(0, Filename.length() - 4);
		// Storage for M line values
		short bSp;
		float bXCorr;
		float bSpScore;
		float bCalcPreMass;
		float bTertiaryScore=0;
		float bPScore=0;
		float bEScore=0;
		int bIPTop;
		int bIPBottom;
		float bIonProportion;
		String bSequence;
		char bValidated;
		byte bEquivSeq;
		String SecondSequence;
		String LastReadLocus = "";
		// Advance to the first spectrum line
		WholeLine = Incoming.readLine();
		while ((WholeLine != null) && (!WholeLine.startsWith("S\t"))) {
			WholeLine = Incoming.readLine();
		}

		// Proceed until the end of the file
		while (WholeLine != null) {
			// Assume we're on an S line
			LastSoFar = MasterListRunner;
			OUTCount++;
			AbortSpectrum = false;
			bDeltCN = 0.0f;
			SecondDeltCN = 0.01f;
			bXCorrRank = 1;
			SecondXCorrRank = 1;
			bSp = 1;
			FoundMLine = false;
			bSequence = "";
			Parser = new StringTokenizer(WholeLine);
			Parser.nextToken();
			bFileName = bSubdirectory + '.' + Parser.nextToken() + '.'
					+ Parser.nextToken() + '.' + Parser.nextToken();
			bEquivSeq = 0;
			// Skip time to match
			Parser.nextToken();
			// Skip server name
			Parser.nextToken();
			//System.out.println(WholeLine);
			bPrecursorMass = new Float(Parser.nextToken()).floatValue();
			bTotalIntensity = new Float(Parser.nextToken()).floatValue();
			// Ignore lowest Sp and # of matched peptides
			WholeLine = Incoming.readLine();

			try {
				while ((bDeltCN == 0.0f) && (WholeLine != null)
						&& (WholeLine.startsWith("M\t"))) {
					// Assume we're on an M line
					// Determine whether this is a nonzero DeltCN line
					Parser = new StringTokenizer(WholeLine);
					int colCount = Parser.countTokens();
					if(colCount>11)
					{
						Cutoffs.colCount = colCount;
					}
					Parser.nextToken();
					bXCorrRank = new Short(Parser.nextToken()).shortValue();
					bSp = new Short(Parser.nextToken()).shortValue();
					Parser.nextToken();
					bDeltCN = new Float(Parser.nextToken()).floatValue();
					Parser.nextToken();
					Parser.nextToken();
					Parser.nextToken();
					Parser.nextToken();
					SecondSequence = Parser.nextToken().toUpperCase();
					if (bDeltCN == 0.0f) {
						bEquivSeq++;
						SecondXCorrRank = bXCorrRank;
						Parser = new StringTokenizer(WholeLine);
						// Skip the "M" character
						Parser.nextToken();
						// Skip the "XCorr" rank
						Parser.nextToken();
						// Store the information for this M line
						// Store rank by preliminary score
						bSp = new Short(Parser.nextToken()).shortValue();
						// Skip the predicted M/Z
						bCalcPreMass = new Float(Parser.nextToken())
								.floatValue();
						// Skip the DeltCN
						Parser.nextToken();
						// Store XCorr
						bXCorr = new Float(Parser.nextToken()).floatValue();
						// Skip the Sp score
						String scoreString = Parser.nextToken();
						try
						{
							bSpScore = new Float(scoreString).floatValue();
						}
						catch(NumberFormatException nbe)
						{
							bSpScore = 0;
						}
						if(colCount>11)
						{
							bTertiaryScore = new Float(Parser.nextToken()).floatValue();
							if(colCount>12)
							{
								bPScore = new Float(Parser.nextToken()).floatValue();
								bEScore = new Float(Parser.nextToken()).floatValue();
							}
						}
						// Generate the ion proportion
						bIPTop = new Integer(Parser.nextToken()).intValue();
						bIPBottom = new Integer(Parser.nextToken()).intValue();
						bIonProportion = (float) bIPTop / (float) bIPBottom;
						bSequence = Parser.nextToken().toUpperCase();
						if(bSequence.contains("["))
						{
							bSequence = IsoformUtils.getRegularSequence(bSequence);
						}
							
						bValidated = Parser.nextToken().charAt(0);
						WholeLine = Incoming.readLine();
						LocusCount = 0;
						LastReadLocus = "";
						while ((WholeLine != null)
								&& (WholeLine.startsWith("L\t"))) {
							/*
							 * When we get a locus name, we can create an
							 * OUTFile storing the information we've accumulated
							 * in the above S and M lines.
							 */
							Parser = new StringTokenizer(WholeLine);
							// Skip the line label
							Parser.nextToken();
							// Trim this locus to 40 characters
							TempRawLocus = Parser.nextToken();
							if (TempRawLocus.length() > Config.LocusLengthCutoff)
								TempRawLocus = TempRawLocus.substring(0,
										Config.LocusLengthCutoff);
							if (Cutoffs.ShortForm) {
								TempLocus = NRUtils.getAccession(TempRawLocus);
							} else {
								TempLocus = TempRawLocus;
							}
							if (!TempLocus.equals(LastReadLocus)) {
								LastReadLocus = TempLocus;
								LocusCount++;
								// Create a new OUTFile for this locus and this
								// match
								MasterListRunner.Next = new OUTFile();
								MasterListRunner = MasterListRunner.Next;
								FoundMLine = true;
								// Store the information we've gathered in the
								// new object
								MasterListRunner.Locus = TempLocus;
								MasterListRunner.RawLocus = TempRawLocus;
								CurrentDTA = MasterListRunner.DTA;
								CurrentDTA.PrecursorMass = bPrecursorMass;
								CurrentDTA.FileName = bFileName;
								String[] ParsedLine = CurrentDTA.FileName
										.split("\\.");
								CurrentDTA.ChargeState = new Byte(ParsedLine[3])
										.byteValue();

								CurrentDTA.ScanNumber = new Integer(
										ParsedLine[1]).intValue();
								CurrentDTA.RootFileName = ParsedLine[0];
								CurrentDTA.Subdirectory = bSubdirectory;
								CurrentDTA.TotalIntensity = bTotalIntensity;
								CurrentDTA.Sp = bSp;
								CurrentDTA.XCorr = bXCorr;
								CurrentDTA.SpScore = bSpScore;
								CurrentDTA.tertiaryScore = bTertiaryScore;
								CurrentDTA.pvalue = bPScore;
								CurrentDTA.evalue = bEScore;
								CurrentDTA.CalcPreMass = bCalcPreMass;
								CurrentDTA.GetMassOffsets(Cutoffs);
								CurrentDTA.IonProportion = bIonProportion;
								CurrentDTA.Sequence = bSequence;
								CurrentDTA.TrypticSites = CurrentDTA
										.DetermineMissedTryptic(Cutoffs);
								CurrentDTA.DetermineTryptic(Cutoffs);
								CurrentDTA.DetermineModified();
								CurrentDTA.Validated = bValidated;
							}
							WholeLine = Incoming.readLine();
						} // End L line while loop
						if ((CurrentDTA != null) && (LocusCount == 1))
							CurrentDTA.UniqueToLocus = true;
					} // End "if bDeltCN == 0" condition
					else {
						if ((SetDeltCN != 0)
								&& (DTAFile
										.JustLettersFromTrimmed(SecondSequence))
										.equals(DTAFile
												.JustLettersFromTrimmed(bSequence))) {
							SecondDeltCN = bDeltCN;
							bDeltCN = 0.0f;
							SecondXCorrRank = bXCorrRank;
							WholeLine = Incoming.readLine();
							while ((WholeLine != null)
									&& (WholeLine.startsWith("L\t"))) {
								WholeLine = Incoming.readLine();
							}
						} else if ((SetDeltCN != 0) && (bSp == 1)
								&& (bXCorrRank - SecondXCorrRank > 1)) {
							bDeltCN = SecondDeltCN;
							WholeLine = Incoming.readLine();
						} else {
							WholeLine = Incoming.readLine();
						}
					} // End "if bDeltCN != 0" condition
				} // End M line while loop
				if (bDeltCN == 0.0f) {
					if (FoundMLine) {
						bDeltCN = SecondDeltCN;
					}
				}
				if (bDeltCN > 0.0f) {
					// We've hit a nonzero DeltCN
					// Write DeltCN to all newly created OUTFiles
					if (LastSoFar.Next != null) {
						LastSoFar = LastSoFar.Next;
					}
					while (LastSoFar != null) {
						LastSoFar.DTA.DeltCN = bDeltCN;
						LastSoFar.DTA.EquivSeq = bEquivSeq;
						LastSoFar = LastSoFar.Next;
					}
					LastSoFar = MasterListRunner;
				}
			} // End "try" clause
			catch (NumberFormatException failure) {
				System.out.println("\tA number is formatted incorrectly in "
						+ bFileName);
			} catch (NoSuchElementException failure) {
				System.out
						.println("\tA line did not have as many elements as it should in "
								+ bFileName);
			}
			// Queue up the next S line!
			while ((WholeLine != null) && (!WholeLine.startsWith("S\t"))) {
				WholeLine = Incoming.readLine();
			}
		} // End "while (WholeLine != null)"
		return OUTCount;
	}

	public static int ReadMascotDAT(String Filename, OUTFile MasterList,
			IniFile Config) throws IOException {
		int OUTCount = 0;
		OUTFile MasterListRunner = MasterList;
		OUTFile LastSoFar;
		File DatFile = new File(Filename);
		String RootName = Filename.substring(0, Filename.length() - 4);
		FileReader PeptidesFilereader = new FileReader(DatFile);
		BufferedReader PeptidesIncoming = new BufferedReader(PeptidesFilereader);
		String PeptidesWholeLine;
		FileReader QueryFilereader = new FileReader(DatFile);
		BufferedReader QueryIncoming = new BufferedReader(QueryFilereader);
		String QueryWholeLine;
		FileReader SummaryFilereader = new FileReader(DatFile);
		BufferedReader SummaryIncoming = new BufferedReader(SummaryFilereader);
		String SummaryWholeLine;
		String TempLocus;
		String bFileName = "Initial";
		StringTokenizer Parser;
		StringTokenizer Parser2;
		DTAFile CurrentDTA = new DTAFile();
		boolean AbortSpectrum;
		boolean SameAsLast = false;
		int LocusCount;
		int Looper;
		// Mascot-specific storage info
		double Ln10 = Math.log(10);
		int QueryNum;
		int MatchNum;
		byte QueryZ = 0;
		int Line2MatchNum;
		int bQmatch = 0;
		float bObsPrecursorMPlusH;
		int bIonsMatched;
		float bIonProportion;
		String bSequence;
		String bModString;
		StringBuffer SequenceBuilder;
		float bIonsScore;
		float bThreshold;
		float bDeltCN;
		String LocusList;
		// Advance to the first qx_px= line in peptide section
		PeptidesWholeLine = PeptidesIncoming.readLine();
		while ((PeptidesWholeLine != null)
				&& (PeptidesWholeLine.indexOf("name=\"peptides") < 0)) {
			PeptidesWholeLine = PeptidesIncoming.readLine();
		}
		PeptidesWholeLine = PeptidesIncoming.readLine();
		PeptidesWholeLine = PeptidesIncoming.readLine();
		// Advance to the first query section
		QueryWholeLine = QueryIncoming.readLine();
		while ((QueryWholeLine != null)
				&& (QueryWholeLine.indexOf("name=\"query") < 0)) {
			QueryWholeLine = QueryIncoming.readLine();
		}
		// Advance to the summary section
		SummaryWholeLine = SummaryIncoming.readLine();
		while ((SummaryWholeLine != null)
				&& (SummaryWholeLine.indexOf("name=\"summary") < 0)) {
			SummaryWholeLine = SummaryIncoming.readLine();
		}
		if (PeptidesWholeLine == null) {
			System.out.println("This file has no \"peptides\" section!");
			System.exit(0);
		}
		// Go through the peptides section, query by query
		while (PeptidesWholeLine.startsWith("q")) {
			try {
				// Determine which query we're on
				Parser = new StringTokenizer(PeptidesWholeLine, "q_p=,");
				QueryNum = new Integer(Parser.nextToken()).intValue();
				MatchNum = new Integer(Parser.nextToken()).intValue();
				if (Parser.nextToken().equals("-1")) {
					// Skip to next query number
					PeptidesWholeLine = PeptidesIncoming.readLine();
				} else {
					// Get the filename info from "query" section
					TempLocus = "name=\"query"
							+ new Integer(QueryNum).toString() + "\"";
					if (!SameAsLast) {
						while ((QueryWholeLine != null)
								&& (QueryWholeLine.indexOf(TempLocus) < 0)) {
							QueryWholeLine = QueryIncoming.readLine();
						}
						if (QueryWholeLine == null) {
							System.out
									.println("Hit bottom of file while looking up query title "
											+ new Integer(QueryNum).toString());
							QueryFilereader = new FileReader(DatFile);
							QueryIncoming = new BufferedReader(QueryFilereader);
							QueryWholeLine = QueryIncoming.readLine();
							bFileName = Filename + QueryNum;
						} else {
							QueryWholeLine = QueryIncoming.readLine();
							QueryWholeLine = QueryIncoming.readLine();
							Parser2 = new StringTokenizer(QueryWholeLine, "=");
							Parser2.nextToken();
							TempLocus = Parser2.nextToken();
							if (TempLocus.endsWith("dta")) {
								bFileName = URLDecoder.decode(TempLocus);
								bFileName = bFileName.substring(0,
										bFileName.length() - 4);
							} else
								bFileName = Filename + QueryNum;
						}
						// Get the summary info for this query
						TempLocus = "qexp" + new Integer(QueryNum).toString()
								+ "=";
						while ((SummaryWholeLine != null)
								&& (!SummaryWholeLine.startsWith(TempLocus))) {
							SummaryWholeLine = SummaryIncoming.readLine();
						}
						if (SummaryWholeLine == null) {
							System.out
									.println("Hit bottom of file while looking up summary."
											+ TempLocus);
							SummaryFilereader = new FileReader(DatFile);
							SummaryIncoming = new BufferedReader(
									SummaryFilereader);
							SummaryWholeLine = SummaryIncoming.readLine();
							QueryZ = 2;
							bQmatch = 1;
						} else {
							Parser2 = new StringTokenizer(SummaryWholeLine,
									",+-");
							Parser2.nextToken();
							QueryZ = new Byte(Parser2.nextToken()).byteValue();
							TempLocus = "qmatch"
									+ new Integer(QueryNum).toString() + "=";
							SummaryWholeLine = SummaryIncoming.readLine();
							while ((SummaryWholeLine != null)
									&& (!SummaryWholeLine.startsWith(TempLocus))) {
								SummaryWholeLine = SummaryIncoming.readLine();
							}
							if (SummaryWholeLine == null) {
								System.out
										.println("Hit bottom of file while looking up summary."
												+ TempLocus);
								SummaryFilereader = new FileReader(DatFile);
								SummaryIncoming = new BufferedReader(
										SummaryFilereader);
								SummaryWholeLine = SummaryIncoming.readLine();
								bQmatch = 1;
							} else {
								Parser2 = new StringTokenizer(SummaryWholeLine,
										"=");
								Parser2.nextToken();
								bQmatch = new Integer(Parser2.nextToken())
										.intValue();
							}
						}
					}
					/*
					 * System.out.println(new Integer(QueryNum).toString() +
					 * "\t" + new Byte(QueryZ).toString() + "\t" + new
					 * Integer(bQmatch).toString());
					 */
					// Process this query
					bObsPrecursorMPlusH = new Float(Parser.nextToken())
							.floatValue()
							+ new Float(Parser.nextToken()).floatValue() + 1.0f;
					bIonsMatched = new Integer(Parser.nextToken()).intValue();
					bSequence = Parser.nextToken();
					// Skip peaks used from Ions1
					Parser.nextToken();
					// Read in modifications, noting that first and
					// last characters are to be discarded
					bModString = Parser.nextToken();
					bModString = bModString.substring(1,
							bModString.length() - 1);
					bIonsScore = new Float(Parser.nextToken()).floatValue();
					bThreshold = new Float(10.0f * (Math.log(bQmatch)) / (Ln10))
							.floatValue();
					bDeltCN = Protein.RoundTo(
							new Float(Math.pow(10f, bIonsScore / 10.0f)
									/ new Double(bQmatch).doubleValue())
									.floatValue(), 2);
					/*
					 * System.out.println(new Integer(QueryNum).toString() +
					 * "\t" + new Integer(MatchNum).toString() + "\t" +
					 * bFileName + "\t" + bSequence + "\t" + new
					 * Float(bIonsScore).toString() + "\t" + new
					 * Float(bThreshold).toString() + "\t" + new
					 * Integer(bQmatch).toString() + "\t" + new
					 * Float(bDeltCN).toString());
					 */
					/*
					 * System.out.print(new Integer(QueryNum).toString() + "\t"
					 * + new Integer(MatchNum).toString() + "\t" + new
					 * Float(bObsPrecursorMPlusH).toString() + "\t" + new
					 * Integer(bIonsMatched).toString() + "\t" + bSequence +
					 * "\t" + bModString + "\t" + new
					 * Float(bIonsScore).toString() + "\t");
					 */
					// Grab the locus names
					// Set LocusList to include everything after the semicolon
					Parser = new StringTokenizer(PeptidesWholeLine, ";");
					Parser.nextToken();
					LocusList = Parser.nextToken();
					// Now split this string on commas (which separate locus
					// sections)
					Parser = new StringTokenizer(LocusList, ",");
					OUTCount++;
					// Make this sequence like a proper sequence
					SequenceBuilder = new StringBuffer("-.");
					for (Looper = 0; Looper < bSequence.length(); Looper++) {
						SequenceBuilder.append(bSequence.charAt(Looper));
						switch (bModString.charAt(Looper)) {
						case '0':
							break;
						case '1':
							SequenceBuilder.append("*");
							break;
						case '2':
							SequenceBuilder.append("#");
							break;
						case '3':
							SequenceBuilder.append("@");
							break;
						case '4':
							SequenceBuilder.append("$");
							break;
						}
					}
					SequenceBuilder.append(".-");
					Looper = Parser.countTokens();
					while (Parser.hasMoreTokens()) {
						// For each locus, tokenize on doublequotes
						Parser2 = new StringTokenizer(Parser.nextToken(), "\"");
						// Create a new OUTFile for this locus and this match
						MasterListRunner.Next = new OUTFile();
						MasterListRunner = MasterListRunner.Next;
						MasterListRunner.Locus = Parser2.nextToken();
						CurrentDTA = MasterListRunner.DTA;
						CurrentDTA.PrecursorMass = bObsPrecursorMPlusH;
						CurrentDTA.FileName = bFileName;
						CurrentDTA.Subdirectory = RootName;
						// UNFINISHED. Whither TotalIntensity?
						// UNFINISHED. Whither SpScore?
						// UNFINISHED. Whither CalcPreMass?
						CurrentDTA.TotalIntensity = 0.0f;
						CurrentDTA.UniqueToLocus = (Looper == 1);
						CurrentDTA.Sp = 1;
						CurrentDTA.XCorr = bIonsScore;
						CurrentDTA.DeltCN = bDeltCN;
						CurrentDTA.IonProportion = new Float(bIonsMatched)
								.floatValue()
								/ new Float(
										(QueryZ > 2) ? 4 * (bSequence.length() - 1)
												: 2 * (bSequence.length() - 1))
										.floatValue();
						CurrentDTA.Sequence = SequenceBuilder.toString();
					}
					// If the next match's Ions Score is the same, handle it
					// next
					PeptidesWholeLine = PeptidesIncoming.readLine();
					Parser = new StringTokenizer(PeptidesWholeLine, ",");
					// Skip q#_p#=
					Parser.nextToken();
					// Skip calc M
					Parser.nextToken();
					// Skip M offset
					Parser.nextToken();
					// Skip Matched ions
					Parser.nextToken();
					// Skip sequence
					Parser.nextToken();
					// Skip peaks used
					Parser.nextToken();
					// Skip variable mods
					Parser.nextToken();
					// Now we're on Ions Score
					// If the next match has the same score, we process it, too.
					if (new Float(Parser.nextToken()).floatValue() != CurrentDTA.XCorr) {
						// Queue up the next #1
						SameAsLast = false;
						MatchNum = 0;
						while (MatchNum != 1) {
							PeptidesWholeLine = PeptidesIncoming.readLine();
							// Determine which query we're on
							if (PeptidesWholeLine.startsWith("q")) {
								Parser = new StringTokenizer(PeptidesWholeLine,
										"q_p=,");
								QueryNum = new Integer(Parser.nextToken())
										.intValue();
								MatchNum = new Integer(Parser.nextToken())
										.intValue();
							} else {
								MatchNum = 1;
							}
						}
					} else {
						SameAsLast = true;
						OUTCount--;
					}
				}
			} catch (NoSuchElementException oopsie) {
				System.out.println("One of these lines was too short:");
				System.out.println(PeptidesWholeLine);
				System.out.println(QueryWholeLine);
			}
		}
		return OUTCount;
	}

	public void DebugPrint() {
		OUTFile ORunner = this.Next;
		while (ORunner != null) {
			/*System.out.println(ORunner.Locus + "\t" + ORunner.DTA.PrecursorMass
					+ "\t" + ORunner.DTA.FileName + "\t"
					+ ORunner.DTA.Subdirectory + "\t"
					+ ORunner.DTA.TotalIntensity + "\t" + ORunner.DTA.Sp + "\t"
					+ ORunner.DTA.XCorr + "\t" + ORunner.DTA.DeltCN + "\t"
					+ ORunner.DTA.PepConf + "\t" + ORunner.DTA.PepFP + "\t"
					+ ORunner.DTA.IonProportion + "\t" + ORunner.DTA.Sequence);*/
			if (ORunner.DTA.ChargeState==1) {
				System.out.println(ORunner.Locus);
			}
			ORunner = ORunner.Next;
		}
	}

	public static void GetVersion(File CurrentDirectory, String Filename,
			DataSet Data) throws IOException {
		File OutFile = new File(CurrentDirectory, Filename);
		FileReader InputFilereader = new FileReader(OutFile);
		BufferedReader Incoming = new BufferedReader(InputFilereader);
		String LineBuffer = null;
		String WholeLine = Incoming.readLine();
		boolean NotFound = true;
		StringTokenizer Parser;
		while ((WholeLine != null) && (NotFound)) {
			Parser = new StringTokenizer(WholeLine);
			if (Parser.hasMoreTokens()) {
				LineBuffer = Parser.nextToken();
				if (LineBuffer.toUpperCase().startsWith("TURBO")
						|| LineBuffer.toUpperCase().startsWith("SEQUEST")) {
					NotFound = false;
					LineBuffer = WholeLine;
				}
			}
			WholeLine = Incoming.readLine();
		}
		if (NotFound) {
			Data.SQTGenerator = "SEQUEST";
			Data.SQTGeneratorVersion = "?";
		} else {
			Parser = new StringTokenizer(LineBuffer);
			Data.SQTGenerator = Parser.nextToken();
			if (Parser.hasMoreTokens()) {
				Data.SQTGeneratorVersion = Parser.nextToken();
			}
		}
		Incoming.close();
	}

	/*
	 * In an effort to clean up the reading of .out files, I've copied my code
	 * and hope to make this function more flexible in how it reads .out files
	 * as well as more comprehensible. 1) multiple loci may be listed after the
	 * first matching peptide. This is flagged by the presence of a +x before
	 * the first matching sequence. The other matching loci are the first tokens
	 * on the succeeding lines. The number of other matching loci is limited to
	 * 20, although the +x may indicate a larger number of possible loci for
	 * this match.
	 * 
	 * 2) The DeltCN value in line 2 may be 0.0. This could result from an
	 * exchange of L for I or Q for K, but there have been reports of identical
	 * sequences appearing in lines one and 2 as well. Suffice to say that if a
	 * 0.0 DeltCN is found in line 2, the DeltCN on line 3 is sought out.
	 * 
	 * If only one locus maps to this sequence, an isolated OUTFile is returned.
	 * If multiple loci map to this sequence, a linked list of OUTFiles is
	 * returned.
	 */
	public static OUTFile ReadOUT(File CurrentDirectory, String Filename,
			boolean ThisDirOnly, IniFile Config, SelectCriteria Cutoffs)
			throws IOException {
		OUTFile TempOutFile = new OUTFile();
		OUTFile OUTRunner = TempOutFile;
		DTAFile CurrentDTA;
		DTAFile NewDTA;
		File OutFile = new File(CurrentDirectory, Filename);
		FileReader InputFilereader = new FileReader(OutFile);
		BufferedReader Incoming = new BufferedReader(InputFilereader);
		String LineBuffer;
		String WholeLine = Incoming.readLine();
		StringTokenizer Parser;
		boolean NoDeltCNYet = true;
		boolean WorkingOnHeader = true;
		boolean IncludesIDField = false;
		boolean IncludesSfField = false;
		// Config.OutIncludesID;
		float Top;
		short Bottom;
		int OtherMatchCount = 0;
		float bPrecursorMass = 0.0f;
		float bTotalIntensity = 0.0f;
		float bDeltCN = 0.0f;
		byte bEquivSeq = 0;
		String bFileName;
		String bSubdirectory;
		String[] ParsedLine;
		String TempLocus, TempRawLocus;
		String LastReadLocus = "";
		try {
			if (ThisDirOnly)
				bSubdirectory = ".";
			else
				bSubdirectory = CurrentDirectory.getName();
			bFileName = Filename.substring(0, Filename.length() - 4);
			Incoming.readLine();
			WholeLine = Incoming.readLine();
			WholeLine = Incoming.readLine();
			while ((WholeLine != null) && (NoDeltCNYet)) {
				Parser = new StringTokenizer(WholeLine, " \t/=,(");
				// If there is any text in this line...
				if (Parser.hasMoreTokens()) {
					LineBuffer = Parser.nextToken();
					// If this is a sequence line, we're out of the header
					// portion
					if (LineBuffer.equals("#")
							&& !WholeLine.startsWith(" # amino")
							&& !WholeLine.startsWith(" # bases")) {
						WorkingOnHeader = false;
						if (WholeLine.indexOf("Id#") > -1) {
							IncludesIDField = true;
						}
						if (WholeLine.indexOf("Sf") > -1) {
							IncludesSfField = true;
						}
						WholeLine = Incoming.readLine();
						WholeLine = Incoming.readLine();
						if (WholeLine.length() == 0) {
							System.out.println("\t" + Filename
									+ " did not list any matches!");
							return null;
						}
						Parser = new StringTokenizer(WholeLine, " \t/=,(");
						LineBuffer = Parser.nextToken();
					}
					if (WorkingOnHeader) {
						/*
						 * We haven't reached the "1." line yet. Note that if
						 * total intensity and precursor mass are on the same
						 * line, we'll only be able to get the first. In
						 * SEQUEST-SNP files, precursor mass is on the same line
						 * as a mass tolerance, so it's required that we kill
						 * the rest of the line.
						 */
						while (Parser.hasMoreTokens()) {
							/*
							 * We can't be sure which line has the total
							 * intensity and precursor mass information, so
							 * we'll read all of them, moving word by word.
							 */
							if (LineBuffer.equals("total")) {
								/*
								 * The token after next should be total
								 * intensity.
								 */
								LineBuffer = Parser.nextToken();
								bTotalIntensity = new Float(Parser.nextToken())
										.floatValue();
								Parser = new StringTokenizer("");
							} else if (LineBuffer.equals("mass")) {
								/*
								 * The next token should be the observer
								 * precursor mass.
								 */
								bPrecursorMass = new Float(Parser.nextToken())
										.floatValue();
								Parser = new StringTokenizer("");
							} else {
								LineBuffer = Parser.nextToken();
							}
						}
						WholeLine = Incoming.readLine();
					} else {
						/*
						 * We're past the header. This must be a numbered line.
						 */
						// Determine if this is a nonzero DeltCN
						Parser.nextToken();
						Parser.nextToken();
						Parser.nextToken();
						if (IncludesIDField)
							Parser.nextToken();
						bDeltCN = new Float(Parser.nextToken()).floatValue();
						if (bDeltCN == 0.0f) {
							bEquivSeq++;
							// Restore parser
							Parser = new StringTokenizer(WholeLine, " \t/");
							Parser.nextToken();
							// Make a new OUTFile object to store this set of
							// IDs
							OUTRunner.Next = new OUTFile();
							OUTRunner = OUTRunner.Next;
							CurrentDTA = OUTRunner.DTA;
							// Store the invariant fields drawn from the header
							CurrentDTA.PrecursorMass = bPrecursorMass;
							CurrentDTA.FileName = bFileName;
							ParsedLine = CurrentDTA.FileName.split("\\.");
							CurrentDTA.ChargeState = new Byte(ParsedLine[3])
									.byteValue();
							CurrentDTA.ScanNumber = new Integer(ParsedLine[1])
									.intValue();
							CurrentDTA.RootFileName = ParsedLine[0];
							CurrentDTA.Subdirectory = bSubdirectory;
							CurrentDTA.TotalIntensity = bTotalIntensity;
							// Store the characteristics of this ID
							Parser.nextToken();
							CurrentDTA.Sp = new Short(Parser.nextToken())
									.shortValue();
							if (IncludesIDField)
								Parser.nextToken();
							CurrentDTA.CalcPreMass = new Float(
									Parser.nextToken()).floatValue();
							CurrentDTA.GetMassOffsets(Cutoffs);
							Parser.nextToken();
							CurrentDTA.XCorr = new Float(Parser.nextToken())
									.floatValue();
							if (CurrentDTA.XCorr == 0) {
								System.out.println("\t" + Filename
										+ " has XCorr of 0.0!");
								return null;
							}
							CurrentDTA.SpScore = new Float(Parser.nextToken())
									.floatValue();
							if (IncludesSfField)
								Parser.nextToken();
							Top = new Float(Parser.nextToken()).floatValue();
							Bottom = new Short(Parser.nextToken()).shortValue();
							CurrentDTA.IonProportion = Top / Bottom;
							LineBuffer = Parser.nextToken();
							if (LineBuffer.length() > Config.LocusLengthCutoff)
								OUTRunner.RawLocus = LineBuffer.substring(0,
										Config.LocusLengthCutoff);
							else
								OUTRunner.RawLocus = LineBuffer;
							if (Cutoffs.ShortForm) {
								OUTRunner.Locus = NRUtils
										.getAccession(OUTRunner.RawLocus);
							} else {
								OUTRunner.Locus = OUTRunner.RawLocus;
							}
							LastReadLocus = OUTRunner.Locus;
							// Now deal with either +1 or sequence field
							LineBuffer = Parser.nextToken();
							if (LineBuffer.charAt(0) == '+') {
								// There will be more loci containing this
								// sequence
								// Currently, LineBuffer is not on the sequence
								LineBuffer = Parser.nextToken();
							} else {
								// This is a unique ID; LineBuffer contains the
								// sequence.
								CurrentDTA.UniqueToLocus = true;
							}
							// Check for SEQUEST-SNP formatted sequence
							if (LineBuffer.charAt(0) == '(') {
								CurrentDTA.Sequence = new Character(
										LineBuffer.charAt(1)).toString()
										+ "." + LineBuffer.substring(3) + ".-";
							} else {
								// This is a modern SEQUEST formatted sequence
								CurrentDTA.Sequence = LineBuffer.toUpperCase();
							}
							if (CurrentDTA.Sequence.length() == 4) {
								System.out.println("\t" + Filename
										+ " reports no sequence!");
								return null;
							}
							CurrentDTA.DetermineTryptic(Cutoffs);
							CurrentDTA.DetermineModified();
							CurrentDTA.TrypticSites = CurrentDTA
									.DetermineMissedTryptic(Cutoffs);
							WholeLine = Incoming.readLine();
							Parser = new StringTokenizer(WholeLine, " \t");
							// Read in other loci containing this sequence
							while (Parser.hasMoreTokens()) {
								LineBuffer = Parser.nextToken();
								if (LineBuffer.charAt(LineBuffer.length() - 1) == '.') {
									// This is a numbered line! Go back through
									// larger loop.
									Parser = new StringTokenizer("");
								} else {
									if (IncludesIDField)
										LineBuffer = Parser.nextToken();
									/*
									 * Since loci appearing on the line with
									 * XCorr are limited to only 40 characters,
									 * limit this one to that length, too.
									 */
									if (LineBuffer.length() > Config.LocusLengthCutoff)
										TempRawLocus = LineBuffer.substring(0,
												Config.LocusLengthCutoff);
									else
										TempRawLocus = LineBuffer;
									if (Cutoffs.ShortForm) {
										TempLocus = NRUtils
												.getAccession(TempRawLocus);
									} else {
										TempLocus = TempRawLocus;
									}
									if (!TempLocus.equals(LastReadLocus)) {
										// This is an additional locus
										// containing this sequence
										LastReadLocus = TempLocus;
										OUTRunner.Next = new OUTFile();
										OUTRunner = OUTRunner.Next;
										OUTRunner.Locus = TempLocus;
										OUTRunner.RawLocus = TempRawLocus;
										// Write the invariant fields for this
										// .out file
										NewDTA = OUTRunner.DTA;
										NewDTA.PrecursorMass = bPrecursorMass;
										NewDTA.FileName = bFileName;
										ParsedLine = NewDTA.FileName
												.split("\\.");
										NewDTA.ChargeState = new Byte(
												ParsedLine[3]).byteValue();
										NewDTA.ScanNumber = new Integer(
												ParsedLine[1]).intValue();
										NewDTA.RootFileName = ParsedLine[0];
										NewDTA.Subdirectory = bSubdirectory;
										NewDTA.TotalIntensity = bTotalIntensity;
										// Write the sequence-specific
										// information
										NewDTA.Sp = CurrentDTA.Sp;
										NewDTA.XCorr = CurrentDTA.XCorr;
										NewDTA.SpScore = CurrentDTA.SpScore;
										NewDTA.CalcPreMass = CurrentDTA.CalcPreMass;
										NewDTA.GetMassOffsets(Cutoffs);
										NewDTA.IonProportion = CurrentDTA.IonProportion;
										NewDTA.Sequence = CurrentDTA.Sequence;
										NewDTA.TrypticSites = NewDTA
												.DetermineMissedTryptic(Cutoffs);
										NewDTA.Tryptic = CurrentDTA.Tryptic;
										NewDTA.Modified = CurrentDTA.Modified;
									}
									WholeLine = Incoming.readLine();
									Parser = new StringTokenizer(WholeLine,
											" \t");
								}
							}
						} else {
							// We've hit a nonzero DeltCN and are done reading
							// the file.
							NoDeltCNYet = false;
							// Write the DeltCN to every OUTFile object created.
							OUTRunner = TempOutFile.Next;
							while (OUTRunner != null) {
								OUTRunner.DTA.DeltCN = bDeltCN;
								OUTRunner.DTA.EquivSeq = bEquivSeq;
								OUTRunner = OUTRunner.Next;
							}
						}
					}
				} else {
					if (!WorkingOnHeader) {
						NoDeltCNYet = false;
					}
					WholeLine = Incoming.readLine();
				}
			}
			if (TempOutFile.Next == null) {
				System.out.println("\t" + Filename
						+ " did not list any matches!");
			}
			Incoming.close();
			InputFilereader.close();
			return TempOutFile.Next;
		} catch (NumberFormatException failure) {
			System.out.println("\t" + Filename + " is garbled!");
			return null;
		} catch (NoSuchElementException failure) {
			System.out.println("\t" + Filename + " is missing information!");
			return null;
		}
	}

	public void SetOriginalNext() {

		OUTFile OUTRunner = this;
		while (OUTRunner.Next != null) {
			OUTRunner.OriginalNext = OUTRunner.Next;
			OUTRunner = OUTRunner.Next;
		}
	}

	// Function to sort OUTFile lists by locus.
	public void SortListByLocus() {
		// this is just an empty header; the real data starts at the next item
		if (this.Next != null)
			this.Next = this.Next.Sort(null);
	}

	/*
	 * Dave's amazing quicksorter. I leave it as an exercise for the reader to
	 * determine how this works.
	 */
	private OUTFile Sort(OUTFile Follower) {
		// Recursive quicksorter
		// Returns first item of sorted list (from those starting at this
		// element)
		// Accepts first item to follow
		OUTFile ListAbove = null;
		OUTFile ListBelow = null;
		OUTFile PlaceHolder;
		OUTFile PlaceHolder2;
		PlaceHolder = this.Next;
		// Partition all remaining points of this linked list
		while (PlaceHolder != null) {
			PlaceHolder2 = PlaceHolder.Next;
			if (this.Locus.compareTo(PlaceHolder.Locus) > 0) {
				// Move this item to list above this
				PlaceHolder.Next = ListAbove;
				ListAbove = PlaceHolder;
			} else if (this.Locus.compareTo(PlaceHolder.Locus) < 0) {
				// Move this item to list below this point
				PlaceHolder.Next = ListBelow;
				ListBelow = PlaceHolder;
			} else {
				if (this.DTA.PrecursorMass > PlaceHolder.DTA.PrecursorMass) {
					// Move this item to list above this
					PlaceHolder.Next = ListAbove;
					ListAbove = PlaceHolder;
				} else {
					// Move this item to list below this point
					PlaceHolder.Next = ListBelow;
					ListBelow = PlaceHolder;
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
