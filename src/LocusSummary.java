import java.awt.*;
import java.io.*;
import java.util.*;

// DTASelect LocusSummary
// David L. Tabb
// September 7, 2000

/*
 * Object to store meta-analytical results for a particular locus
 * across multiple DataSets including presence information for each
 * DataSet.
 */
public class LocusSummary {
    String              Name = "";
    String              Description = "";
    long                FoundInPattern = 0;
    byte                FoundCount = 0;
    String              CumulativeCoverage;
    String              SequenceCoverage;
    String              DepthSequenceCoverage;
    byte                Classification = 127;
    LocusSummary        Next;
    LocusSummary        IdenticalLS;

    //LIST FUNCTIONS

    /* This method assigns a clssification to each protein listed in
     * the "Classifications.txt" file.  When proteins are sorted as
     * selected by the user, classification is always the first order
     * applied. */
    public void StructureByClassification(Classification Classifieds) {
	File             CurrentDirectory = new File(System.getProperty("user.dir"));
	File             ParamsFile = new File(CurrentDirectory, "Classifications.txt");
	try {
	    FileReader       InputFileReader = new FileReader(ParamsFile);
	    BufferedReader   Incoming = new BufferedReader(InputFileReader);
	    String           WholeLine = Incoming.readLine();
	    String           LineBuffer;
	    StringTokenizer  Parser;
	    Classification   CRunner = Classifieds;
	    StringBuffer     Lump;
	    LocusSummary     PRunner;
	    while (WholeLine != null) {
		Parser = new StringTokenizer(WholeLine);
		if (Parser.hasMoreTokens()) {
		    LineBuffer = Parser.nextToken();
		    if (LineBuffer.equals("class")) {
			// Create a new Classification for this line.
			CRunner.Next = new Classification();
			CRunner = CRunner.Next;
			CRunner.Identifier = new Byte(Parser.nextToken()).byteValue();
			Lump = new StringBuffer(Parser.nextToken());
			while (Parser.hasMoreTokens()) {
			    Lump.append(" ");
			    Lump.append(Parser.nextToken());
			}
			CRunner.Descriptor = Lump.toString();
		    }
		    else {
			// Look for this protein descriptor in
			// memory.  If it exists, assign it to the
			// appropriate Classification.
			PRunner = this;
			while ( (PRunner.Next != null) &&
				(!PRunner.Next.Name.equals(LineBuffer)) ){
			    PRunner = PRunner.Next;
			}
			if (PRunner.Next != null) {
			    // this is a match.
			    try {
				PRunner.Next.Classification = new Byte(Parser.nextToken()).byteValue();
			    }
			    catch (NoSuchElementException failure) {
				System.out.println("\t\t" + LineBuffer +
						   "\tis listed without a class.  Ignoring...");
			    }
			}
		    }
		}
		WholeLine = Incoming.readLine();
	    }
	    Incoming.close();
	}
	catch (IOException failure) {
	    System.out.println("Error while trying to read Classifications.txt.");
	    System.out.println(failure);
	}
    }

    public void DebugPrint() {
	LocusSummary  LSRunner = this.Next;
	while (LSRunner != null) {
	    System.out.println(LSRunner.Name + "\t" + LSRunner.FoundInPattern);
	    LSRunner = LSRunner.Next;
	}
    }

    /* How many proteins match this FoundInPattern? If this is run
     * after the sort, we could make this run a bit faster... */
    public int CountProteins(long ProbePattern) {
	LocusSummary     Runner = this.Next;
	int              Count = 0;
	while (Runner != null) {
	    if ((Runner.FoundInPattern & ProbePattern) == ProbePattern)
		Count++;
	    Runner = Runner.Next;
	}
	return Count;
    }

    /* Before printing HTML or TXT files, run this to generate the
     * cumulative sequence coverage percentages and strings. */
    public void PrepToPrint(DataSet Sets, boolean UseBasicSeqCov) {
	LocusSummary      Cursor = this.Next;
	LocusSummary      LSRunner;
	LocusSummary      LSBuffer;
	CoverageZone      CZBuffer;
	CoverageZone      CZCumulative;
	int               SetCount = Sets.CountSets();
	int               SetLooper;
	int               SequenceLength;
	while (Cursor != null) {
	    /* Each DataSet has a Cursor pointing to a Protein object.
	     * The following line causes these pointers to move to the
	     * protein with a name matching this string.  If the
	     * protein is not present in a DataSet, that DS's pointer
	     * will be null. */
	    Sets.SetLocusCursors(Cursor.Name);
	    SequenceLength = Sets.CurrentSeqLength();
	    if (SequenceLength == 0) {
		Sets.SetLocusCoverageUndefined();
		Cursor.CumulativeCoverage = "?";
		Cursor.SequenceCoverage = "";
	    }
	    else {
		CZCumulative = new CoverageZone();
		for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		    CZBuffer = Sets.GetLocusCoverageZones(SetLooper);
		    /* If CZBuffer is null, that means that this protein
		     * is not found in this data set. */
		    if (CZBuffer != null) {
			// Add all the coverage zones to the growing
			// pile for this locus across data sets.
			CZCumulative.MergeWith(CZBuffer);
			// Find the minimal set of coverage zones to
			// describe sequence coverage for this protein
			// in this particular data set.
			CZBuffer.MakeMinimal();
			// Now set the sequence coverage for finding
			// identical proteins (may be redone later for
			// producing DTASelect.html files for each
			// DataSet).
			Sets.SetLocusCoverage(SetLooper,
					      CZBuffer.GetConsensusList(),
					      Protein.RoundTo(100.0f *
							      CZBuffer.SeqLength() /
							      SequenceLength,1));
		    }
		}
		if (!UseBasicSeqCov)
		    Cursor.DepthSequenceCoverage = CZCumulative.GetSeqCovString();
		CZCumulative.MakeMinimal();
		Cursor.SequenceCoverage = CZCumulative.GetConsensusList();
		Cursor.CumulativeCoverage = new Float(Protein.RoundTo(100.0f *
								      CZCumulative.SeqLength() /
								      SequenceLength,1)
						      ).toString();
	    }
	    Cursor = Cursor.Next;
	}
	/* Now find which LocusSummary objects are identical (here
	 * judged by which have identical SequenceCoverage strings).
	 * These will be grouped in the output. */
	Cursor = this.Next;
	while (Cursor != null) {
	    /* We'll be comparing Cursor against LSRunner.Next.  That
	     * way we can move LSRunner.Next from the main chain to
	     * the Cursor.IdenticalLS chain if it matches.  */
	    LSRunner = Cursor;
	    while ((LSRunner.Next != null) &&
		   (LSRunner.FoundInPattern == Cursor.FoundInPattern)) {
		if (LSRunner.Next.SequenceCoverage.equals(Cursor.SequenceCoverage)) {
		    /* These match; remove LSRunner.Next and add it to
		     * Cursor.IdenticalLS chain. */
		    LSBuffer = Cursor.IdenticalLS;
		    Cursor.IdenticalLS = LSRunner.Next;
		    LSRunner.Next = LSRunner.Next.Next;
		    Cursor.IdenticalLS.IdenticalLS = LSBuffer;
		}
		else {
		    /* No match; continue to next LocusSummary. */
		    LSRunner = LSRunner.Next;
		}
	    }
	    Cursor = Cursor.Next;
	}
    }

    /* Function to print "Contrast.html" and "Contrast.txt" files.
     * Header includes sample and criteria set enumerations.  Proteins
     * are grouped by the data sets in which they appear, and their
     * presence in a data set is indicated by a sequence coverage
     * percentage in the column representing that data set.  A summary
     * table follows.
     * Run PrepToPrint before this function!  */
    public void CreateOutputFiles(DataSet Sets, DataSet Master,
				  DataSet DSList, DataSet CList,
				  boolean Verbose, boolean SpeCount, 
                                  IniFile Config) {
	String          StringBuffer = System.getProperty("user.dir");
	String          Consensus = "<a TARGET=\"Win1\" href=\"" +
	    Config.BasicSeqCov + "?Db=" +
	    Sets.Next.SequestParams.DBName +
	    "&Ref=";
	String                  SeqCov = "<a TARGET=\"Win1\" href=\"" +
	    Config.DepthSeqCov + "?" +
	    Sets.Next.SequestParams.DBName +
	    "&";
	String          DBName;
	File            CurrentDirectory = new File(StringBuffer);
	File            HTMLOutputFile;
	FileWriter      HTMLOutputFileWriter;
	BufferedWriter  HTMLOutgoing;
	File            TXTOutputFile;
	FileWriter      TXTOutputFileWriter;
	BufferedWriter  TXTOutgoing;
	LocusSummary    Cursor = this.Next;
	LocusSummary    IdenticalsRunner;
	int             SetCount = Sets.CountSets();
	int             SetLooper;
	byte            TargetZ;
	long            CurrentPower;
	long            OtherPower;
	long            LastPattern = 0;
	DataSet         DSRunner;
	int             SectionCounter = 0;
	int             RSectionCounter = 0;
	int             TotalCounter = 0;
	int             RTotalCounter = 0;
	Color           BackColor;
	String          HTMLSummary = "";
	String          TXTSummary = "";
	String          Title = "Contrast";
	StringBuffer    PresenceAndAbsence;
	float           PercentBuffer;
	Protein         LocusBuffer;
	StringTokenizer Parser;
	boolean         CurrentGroupHidden = true;
	boolean         PepIsUnique;
	while (Cursor != null) {
	    TotalCounter ++;
	    RTotalCounter ++;
	    Cursor = Cursor.Next;
	}
	Cursor = this.Next;
	try {
	    CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
	    HTMLOutputFile = new File(CurrentDirectory, "Contrast.html");
	    HTMLOutputFileWriter = new FileWriter(HTMLOutputFile);
	    HTMLOutgoing = new BufferedWriter(HTMLOutputFileWriter);
	    TXTOutputFile = new File(CurrentDirectory, "Contrast.txt");
	    TXTOutputFileWriter = new FileWriter(TXTOutputFile);
	    TXTOutgoing = new BufferedWriter(TXTOutputFileWriter);
	    //Write HTML header
	    HTMLOutgoing.write("<HTML><HEAD><TITLE>C ");
	    Parser = new StringTokenizer(System.getProperty("user.dir"), System.getProperty("file.separator"));
	    while (Parser.hasMoreTokens())
		Title = Parser.nextToken();
	    HTMLOutgoing.write(Title);
	    HTMLOutgoing.write("</title></head><body>\n");
	    HTMLOutgoing.write(Protein.Version() + "\n<BR>");
	    HTMLOutgoing.write(CurrentDirectory.getCanonicalPath() +
			       "\n<BR><table border><tr><td>");
	    HTMLOutgoing.write("Included directories:</td><td>Location</td><td>Hidden</td></tr>\n");
	    //Write TXT header
	    TXTOutgoing.write(Protein.Version() + "\n");
	    TXTOutgoing.write(CurrentDirectory.getCanonicalPath() + "\n");
	    TXTOutgoing.write("Included directories:\tLocation\tHidden\n");
	    DSRunner = DSList.Next;
	    //Enumerate samples included in this analysis
	    while (DSRunner != null) {
		HTMLOutgoing.write("<tr><td>" + DSRunner.FriendlyName
				   + "</td><td>" + DSRunner.SourceDirectory
				   + "</td><td>" + DSRunner.Removed
				   + "</td></tr>\n");
		TXTOutgoing.write(DSRunner.FriendlyName + "\t" +
				  DSRunner.SourceDirectory + "\t" +
				  DSRunner.Removed + "\n");
		DSRunner = DSRunner.Next;
	    }
	    //Enumerate criteria sets included in this analysis
	    HTMLOutgoing.write("</table><P>\nCriteria sets:");
	    TXTOutgoing.write("\nCriteria sets:\n");
	    HTMLOutgoing.write("\n<table><tr>\n");
	    DSRunner = CList.Next;
	    while (DSRunner != null) {
		HTMLOutgoing.write("<td>");
		HTMLOutgoing.write("<P>" + DSRunner.CriteriaName);
		if (DSRunner.Removed)
		    HTMLOutgoing.write(" [HIDDEN]");
		HTMLOutgoing.write("<BR>\n" +
				   DSRunner.CutoffsString + "<BR>\n");
		HTMLOutgoing.write("<table border>");
		HTMLOutgoing.write(DSRunner.Cutoffs.PrintCriteria("<TR><TD>", "</TD></TR>\n", "</TD><TD>") + "\n");
		HTMLOutgoing.write("</table>");
		TXTOutgoing.write(DSRunner.CriteriaName + "\t" + DSRunner.CutoffsString + "\n");
		TXTOutgoing.write(DSRunner.Cutoffs.PrintCriteria("", "\n", "\t") + "\n");
		HTMLOutgoing.write("</td>");
		DSRunner = DSRunner.Next;
	    }
	    HTMLOutgoing.write("</tr></table>\n");
	    if (Master != null) {
		HTMLOutgoing.write("<P><b>Master set selected: " +
				   Master.HappyName() + "</b></P>\n");
		TXTOutgoing.write("\nMaster set selected:\t" +
				  Master.HappyName() + "\n");
	    }
	    HTMLOutgoing.write("<P>Numbers in table report the percentage of residues in protein sequence "
			       + "that are represented by at least one peptide passing the criteria.<BR>\n");
            if (SpeCount) {
	        HTMLOutgoing.write("Numbers in parantheses represent (Spectral Count, Spectral Count / Molecular Weight).<BR>\n");
            }
	    HTMLOutgoing.write("<a href=\"#Summary\">Jump</a> to summary table.\n");
	    HTMLOutgoing.write("<br><table border>\n");
	    HTMLOutgoing.write("<tr><td width=\"120\"><a name=\"" +
			       Cursor.FoundInPattern + "\">Locus</a></td>");
	    TXTOutgoing.write("Numbers in table report the percentage of residues in protein sequence "
			      + "that are represented by at least one peptide passing the criteria.\n");
	    TXTOutgoing.write("\nLocus\t");
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		HTMLOutgoing.write("<td>" + Sets.TableHead(SetLooper) + "</td>");
		TXTOutgoing.write(Sets.ShortName(SetLooper) + "\t");
	    }
	    HTMLOutgoing.write("<td>Total</td><td width=\"300\">Description</td></tr>\n");
	    TXTOutgoing.write("Total\tDescription\n");
	    DBName = Sets.Next.SequestParams.DBName;
	    /* While there are still loci in this group to include,
	     * print of their information with appropriate links.
	     * Sequence coverage percentage links to the appropriate
	     * locus in DTASelect.html. */
	    /* Presence or absence is determined by bit analysis of a
	     * locus' FoundInPattern.  Presence in the first data set
	     * adds a 1 to this pattern, presence in the second adds a
	     * 2, and presence in the third adds 4.  In this way, the
	     * presence of loci can be reconstructed by reading the
	     * bits of the FoundInPattern.  Since the LocusSummary
	     * objects are sorted by FoundInPattern, a change in the
	     * FoundInPattern implies that a new table is called for
	     * since a different combination of presence and absence
	     * occured.
	     */
	    while (Cursor != null) {
		if ( (LastPattern > 0) && (LastPattern != Cursor.FoundInPattern) ) {
		    //Add to the summary table string for the previous group
		    //How many of this type?  What percentage does this constitute?
		    HTMLSummary +=
			"<tr><td><a href=\"#" + new Long(LastPattern).toString() + "\">"
			+ new Integer(RSectionCounter).toString() + "</a></td>"
			+ "<td><a href=\"#" + new Long(LastPattern).toString() + "\">"
			+ new Integer(SectionCounter).toString() + "</a></td><td>"
			+ new Float(Protein.RoundTo(
						    100f * new Float(SectionCounter).floatValue() /
						    new	Float(TotalCounter).floatValue(), 1)
				    ).toString() + "%";
		    TXTSummary += new Integer(RSectionCounter).toString() + "\t"
			+ new Integer(SectionCounter).toString() + "\t"
			+ new Float(
				    Protein.RoundTo(
						    100f * new Float(SectionCounter).floatValue() /
						    new	Float(TotalCounter).floatValue(), 1)
				    ).toString() + "%";
		    //Show the X symbols in the proper spots
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			HTMLSummary += "</td><td>" + (((LastPattern & CurrentPower) > 0) ? "X" : " ");
			TXTSummary += "\t" + (((LastPattern & CurrentPower) > 0) ? "X" : " ");
		    }
		    HTMLSummary += "</td></tr>\n";
		    TXTSummary += "\n";
		    //Now actually write the summary for the previous group
		    HTMLOutgoing.write("<tr><td>" +
				       new Integer(SectionCounter).toString());
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper++) {
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			HTMLOutgoing.write("</td><td>" +
					   (((LastPattern & CurrentPower) > 0) ? "X" : ""));
		    }
		    HTMLOutgoing.write("</td><td></td><td></td></tr>");
		    HTMLOutgoing.write("</table><table border>\n<tr><td width=\"120\"><a name=\""
				       + Cursor.FoundInPattern +
				       "\">Locus</a></td>");
		    TXTOutgoing.write(new Integer(SectionCounter).toString());
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper++) {
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			HTMLOutgoing.write("<td>" + Sets.TableHead(SetLooper) +
					   "</td>");
			TXTOutgoing.write("\t" + (((LastPattern & CurrentPower) > 0) ? "X" : " "));
		    }
		    HTMLOutgoing.write("<td>Total</td><td width=\"300\">Description</td></tr>\n");
		    TXTOutgoing.write("\n\n");
		    RSectionCounter = 0;
		    SectionCounter = 0;
		    CurrentGroupHidden = Cursor.IsHidden(Sets);
		}
		if (LastPattern == 0) {
		    CurrentGroupHidden = Cursor.IsHidden(Sets);
		}
		//Move cursors in each DataSet to this locus
		Sets.SetLocusCursors(Cursor.Name);
		SectionCounter++;
		RSectionCounter++;
		if (!CurrentGroupHidden) {
		if (Cursor.IdenticalLS == null) {
		    //Write the name of this locus to the file
		    HTMLOutgoing.write("<tr><td><FONT COLOR=\"red\">"
				       + Cursor.Name + "<FONT COLOR=\"black\">");
		    TXTOutgoing.write(Cursor.Name);
		    //Write sequence coverage percents to the file
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
			HTMLOutgoing.write("</td><td>");
			TXTOutgoing.write("\t");
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			//If this column should appear for this group
			if ( (Cursor.FoundInPattern & CurrentPower) > 0 ) {
			    if (!SpeCount) {
                                HTMLOutgoing.write(
                                               "<a target=\"" + Sets.ShortName(SetLooper) +
                                               "\" href=\"" + Sets.ShortName(SetLooper) +
                                               ".html#" + Cursor.Name +
                                               "\">" +
                                               Sets.GetLocusCursor(SetLooper).SequenceCoverage +
                                               "</a>");
                                TXTOutgoing.write(new Float(Sets.GetLocusCursor(SetLooper).SequenceCoverage).toString());
			    }
			    else {
			        HTMLOutgoing.write(
					       "<table><tr><td><a target=\"" + Sets.ShortName(SetLooper) +
					       "\" href=\"" + Sets.ShortName(SetLooper) +
					       ".html#" + Cursor.Name +
					       "\">" +
					       Sets.GetLocusCursor(SetLooper).SequenceCoverage +
					       "</a></td><td></td><td></td><td>(" + Sets.GetLocusCursor(SetLooper).NSpectra + ",</td><td></td><td>" + Sets.GetLocusCursor(SetLooper).NSpectra/Sets.GetLocusCursor(SetLooper).MolWt + ")</td></tr></table>");
			        TXTOutgoing.write(new Float(Sets.GetLocusCursor(SetLooper).SequenceCoverage).toString());
			    }
			}
		    }
		    //Write cumulative percentage across all DataSets
		    if (Config.UseDepthCGI) {
			HTMLOutgoing.write("</td><td>" +
					   SeqCov + Cursor.Name +
					   Cursor.DepthSequenceCoverage +
					   "\">" + 
					   Cursor.CumulativeCoverage +
					   "</a>");
		    }
		    else {
			HTMLOutgoing.write("</td><td>" +
					   Consensus + Cursor.Name +
					   "&Pep=" +
					   Cursor.SequenceCoverage +
					   "\">" + 
					   Cursor.CumulativeCoverage +
					   "</a>");
		    }
		    TXTOutgoing.write("\t" + Cursor.CumulativeCoverage);
		    //Write the protein's long name to the file
		    HTMLOutgoing.write("</td><td>" + Cursor.Description + "</td></tr>\n");
		    TXTOutgoing.write("\t" + Cursor.Description + "\n");
		}
		/* Alternate version of the above; there are multiple
		 * equivalent loci here.  */
		else {
		    HTMLOutgoing.write("<tr><td><FONT COLOR=\"red\">"
				       + Cursor.Name + "<FONT COLOR=\"black\">");
		    TXTOutgoing.write(Cursor.Name);
		    IdenticalsRunner = Cursor.IdenticalLS;
		    while (IdenticalsRunner != null) {
			RSectionCounter++;
			RTotalCounter++;
			HTMLOutgoing.write("<HR>" + IdenticalsRunner.Name);
			IdenticalsRunner = IdenticalsRunner.IdenticalLS;
		    }
		    PresenceAndAbsence = new StringBuffer();
		    //Write sequence coverage percents to the file
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
			HTMLOutgoing.write("</td><td>");
			TXTOutgoing.write("\t");
			PresenceAndAbsence.append("\t");
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			//If this column should appear for this group
			if ( (Cursor.FoundInPattern & CurrentPower) > 0 ) {
                            if (!SpeCount) {
                                HTMLOutgoing.write(
                                               "<a target=\"" + Sets.ShortName(SetLooper) +
                                               "\" href=\"" + Sets.ShortName(SetLooper) +
                                               ".html#" + Cursor.Name +
                                               "\">" +
                                               Sets.GetLocusCursor(SetLooper).SequenceCoverage +
                                               "</a>");
                                TXTOutgoing.write(new Float(Sets.GetLocusCursor(SetLooper).SequenceCoverage).toString());
                                PresenceAndAbsence.append("X");
                            }
                            else {
			        HTMLOutgoing.write(
					       "<table><tr><td><a target=\"" + Sets.ShortName(SetLooper) +
					       "\" href=\"" + Sets.ShortName(SetLooper) +
					       ".html#" + Cursor.Name +
					       "\">" +
					       Sets.GetLocusCursor(SetLooper).SequenceCoverage +
					       "</a></td><td></td><td></td><td>(" + Sets.GetLocusCursor(SetLooper).NSpectra + ",</td><td></td><td>" + Sets.GetLocusCursor(SetLooper).NSpectra/Sets.GetLocusCursor(SetLooper).MolWt + ")</td></tr></table>");
			        TXTOutgoing.write(new Float(Sets.GetLocusCursor(SetLooper).SequenceCoverage).toString());
			        PresenceAndAbsence.append("X");
                            }
			}
		    }
		    //Write all cumulative percentages across all DataSets
		    IdenticalsRunner = Cursor.IdenticalLS;
		    if (Config.UseDepthCGI) {
			HTMLOutgoing.write("</td><td>" +
					   SeqCov + Cursor.Name +
					   Cursor.DepthSequenceCoverage +
					   "\">" + 
					   Cursor.CumulativeCoverage +
					   "</a>");
		    }
		    else {
			HTMLOutgoing.write("</td><td>" +
					   Consensus + Cursor.Name +
					   "&Pep=" +
					   Cursor.SequenceCoverage +
					   "\">" + 
					   Cursor.CumulativeCoverage +
					   "</a>");
		    }
		    while (IdenticalsRunner != null) {
			if (Config.UseDepthCGI) {
			    HTMLOutgoing.write("<hr>" +
					       SeqCov + IdenticalsRunner.Name +
					       IdenticalsRunner.DepthSequenceCoverage +
					       "\">" + 
					       IdenticalsRunner.CumulativeCoverage +
					       "</a>");
			}
			else {
			    HTMLOutgoing.write("<hr>" +
					       Consensus + IdenticalsRunner.Name +
					       "&Pep=" +
					       IdenticalsRunner.SequenceCoverage +
					       "\">" + 
					       IdenticalsRunner.CumulativeCoverage +
					       "</a>");
			}
			IdenticalsRunner = IdenticalsRunner.IdenticalLS;
		    }
		    TXTOutgoing.write("\t" +
				      Cursor.CumulativeCoverage);
		    //Write all descriptions
		    IdenticalsRunner = Cursor.IdenticalLS;
		    HTMLOutgoing.write("</td><td>" +
				       Cursor.Description);
		    TXTOutgoing.write("\t" + Cursor.Description + "\n");
		    while (IdenticalsRunner != null) {
			HTMLOutgoing.write("<HR>" +
					   IdenticalsRunner.Description);
			TXTOutgoing.write(IdenticalsRunner.Name +
					  PresenceAndAbsence.toString() +
					  "\t" +
					  IdenticalsRunner.CumulativeCoverage +
					  "\t" +
					  IdenticalsRunner.Description
					  + "\n");
			IdenticalsRunner = IdenticalsRunner.IdenticalLS;
		    }
		    HTMLOutgoing.write("</td></tr>\n");
		}
		if (Verbose) {
		    //For verbose output, add peptide information to file
		    LocusBuffer = new Protein();
		    //Accumulate peptide list from each dataset
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper++) {
			CurrentPower = Math.round(Math.pow(2,SetLooper));
			if ( (Cursor.FoundInPattern & CurrentPower) > 0 ) {
			    LocusBuffer.AccumulateDTAsFrom(Sets.GetLocusCursor(SetLooper));
			}
		    }
		    //Ditch duplicate copies of peptides
		    LocusBuffer.DTAs.KeepOnlyDTAWithHighXCorr();
		    //Print list of peptides to file
		    Parser = new StringTokenizer(LocusBuffer.ReportPeptideList());
		    while (Parser.hasMoreTokens()) {
			if (Parser.nextToken().charAt(0) == 'Y')
			    PepIsUnique = true;
			else
			    PepIsUnique = false;
			StringBuffer = Parser.nextToken();
			TargetZ = new Byte(Parser.nextToken()).byteValue();
			HTMLOutgoing.write("<tr><td>");
			if (PepIsUnique) {
			    HTMLOutgoing.write("* ");
			    TXTOutgoing.write("* ");
			}
			HTMLOutgoing.write(StringBuffer +
					    " +" + new Byte(TargetZ).toString());
			TXTOutgoing.write(StringBuffer + " +" + new
			    Byte(TargetZ).toString());
			for (SetLooper = 0; SetLooper < SetCount; SetLooper++) {
			    CurrentPower = Math.round(Math.pow(2,SetLooper));
			    HTMLOutgoing.write("</td><td>");
			    TXTOutgoing.write("\t");
			    if ( (Cursor.FoundInPattern &
				  CurrentPower) > 0 ) {
				HTMLOutgoing.write(Sets.GetLocusCursorXCorr(SetLooper,
									    StringBuffer,
									    TargetZ));
				TXTOutgoing.write(Sets.GetLocusCursorXCorr(SetLooper,
									   StringBuffer,
									   TargetZ));
			    }
			}
			HTMLOutgoing.write("</td></tr>\n");
			TXTOutgoing.write("\n");
		    }
		    HTMLOutgoing.write("</table><table border>\n");
		    HTMLOutgoing.write("<tr><td width=\"375\">Locus</td>");
		    TXTOutgoing.write("\n");
		    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
			HTMLOutgoing.write("<td>" + Sets.TableHead(SetLooper) + "</td>");
		    }
		    HTMLOutgoing.write("<td>Total</td><td>Description</td></tr>\n");
		}
		}
		LastPattern = Cursor.FoundInPattern;
		Cursor = Cursor.Next;
	    }
	    HTMLOutgoing.write("<tr><td>" +
			       new Integer(SectionCounter).toString());
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper++) {
		CurrentPower = Math.round(Math.pow(2,SetLooper));
		HTMLOutgoing.write("</td><td>" +
				   (((LastPattern & CurrentPower) > 0) ? "X" : ""));
	    }
	    HTMLOutgoing.write("</td><td></td><td></td></tr>");
	    HTMLSummary +=
		"<tr><td><a href=\"#" + new Long(LastPattern).toString() + "\">"
		+ new Integer(RSectionCounter).toString() + "</a></td>"
		+ "<td><a href=\"#" + new Long(LastPattern).toString() + "\">"
		+ new Integer(SectionCounter).toString() + "</a></td><td>"
		+ new Float(Protein.RoundTo(
					    100f * new Float(SectionCounter).floatValue() /
					    new	Float(TotalCounter).floatValue(), 1)
			    ).toString() + "%";
	    TXTSummary += new Integer(RSectionCounter).toString() + "\t"
		+ new Integer(SectionCounter).toString() + "\t"
		+ new Float(
			    Protein.RoundTo(
					    100f * new Float(SectionCounter).floatValue() /
					    new	Float(TotalCounter).floatValue(), 1)
			    ).toString() + "%";
	    TXTOutgoing.write(new Integer(SectionCounter).toString());
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		CurrentPower = Math.round(Math.pow(2,SetLooper));
		HTMLSummary += "</td><td>" + (((LastPattern &
						CurrentPower) > 0) ?
					      "X" : " ");
		TXTSummary += "\t" + (((LastPattern & CurrentPower) > 0) ? "X" : " ");
		TXTOutgoing.write("\t" + (((LastPattern & CurrentPower) > 0) ? "X" : " "));
	    }
	    HTMLSummary += "</td></tr>\n";
	    HTMLOutgoing.write("</table><P><table border>\n<a name=\"Summary\"></a><tr><td>Redundant<br>Count</td><td>Nonredundant<BR>Count</td><td>Percent</td>");
	    TXTSummary += "\n";
	    TXTOutgoing.write("\n\nR Count\tNR Count\tPercent\t");
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		HTMLOutgoing.write("<td>" + Sets.TableHead(SetLooper)
				   + "</td>");
		TXTOutgoing.write(Sets.ShortName(SetLooper) + "\t");
	    }
	    HTMLOutgoing.write("</tr>\n");
	    HTMLOutgoing.write(HTMLSummary);
	    HTMLOutgoing.write("<tr><td>" + new
		Integer(RTotalCounter).toString() + "</td><td>"
		+ new Integer(TotalCounter).toString() + "</td><td>");
	    TXTOutgoing.write("\n");
	    TXTOutgoing.write(TXTSummary);
	    TXTOutgoing.write(new Integer(RTotalCounter).toString() + "\t");
	    TXTOutgoing.write(new Integer(TotalCounter).toString() + "\t");
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		CurrentPower = Math.round(Math.pow(2,SetLooper));
		HTMLOutgoing.write("</td><td>" + new Integer(this.CountProteins(CurrentPower)).toString());
		TXTOutgoing.write("\t" + new Integer(this.CountProteins(CurrentPower)).toString());
	    }
	    HTMLOutgoing.write("</td></tr></table>\n");
	    //Show two-way set comparisons in colorful table.
	    HTMLOutgoing.write("<P>The cells of the following table show the percentage of the row sample proteins which are found in the column sample.</P>");
	    HTMLOutgoing.write("<table border><tr><td>Sample</td><td>Proteins</td>");
	    TXTOutgoing.write("\n\nSample\tProteins");
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		HTMLOutgoing.write("<td>" + Sets.ShortName(SetLooper) + "</td>");
		TXTOutgoing.write("\t" + Sets.ShortName(SetLooper));
	    }
	    HTMLOutgoing.write("</tr>\n");
	    TXTOutgoing.write("\n");
	    for (SetLooper = 0; SetLooper < SetCount; SetLooper ++) {
		CurrentPower = Math.round(Math.pow(2,SetLooper));
		HTMLOutgoing.write("<tr><td>");
		// Yes, I realize I'm not using the TotalCounter and
		// RTotalCounter variables correctly below!
		RTotalCounter = new Integer(this.CountProteins(CurrentPower)).intValue();
		HTMLOutgoing.write(Sets.ShortName(SetLooper) + "</td><td>" +
				   new Integer(RTotalCounter).toString() + "</td>");
		TXTOutgoing.write(Sets.ShortName(SetLooper) + "\t" + new Integer(RTotalCounter).toString());
		for (TotalCounter = 0; TotalCounter < SetCount; TotalCounter++) {
		    OtherPower = Math.round(Math.pow(2,TotalCounter));
		    PercentBuffer = new Float(this.CountProteins(OtherPower | CurrentPower)).floatValue() /
			new Float(RTotalCounter).floatValue();
		    BackColor = new Color(Color.HSBtoRGB(0.3f - 0.29f * PercentBuffer, 0.6f, 1.0f));
		    StringBuffer = HTMLColor(BackColor);
		    HTMLOutgoing.write("<td bgcolor=\"#" + StringBuffer +
				       "\">" +
				       new Float(Protein.RoundTo(100f * PercentBuffer, 1)).toString() +
				       "%</td>");
		    TXTOutgoing.write("\t" + new Float(Protein.RoundTo(100f * PercentBuffer, 1)).toString());
		}
		HTMLOutgoing.write("</tr>\n");
		TXTOutgoing.write("\n");
	    }
	    HTMLOutgoing.write("</table></body></html>");
	    HTMLOutgoing.flush();
	    HTMLOutgoing.close();
	    TXTOutgoing.write("\n");
	    TXTOutgoing.flush();
	    TXTOutgoing.close();
	}
	catch (IOException failure) {
	    System.out.println("File error while writing Contrast output files");
	}
    }

    // Get the HTML to describe this color
    public String HTMLColor(Color Target) {
	int            RedComponent = Target.getRed();
	int            GreenComponent = Target.getGreen();
	int            BlueComponent = Target.getBlue();
	StringBuffer   HTML = new StringBuffer();
	String         OneColor;
	OneColor = Integer.toHexString(RedComponent);
	if (OneColor.length() < 2)
	    HTML.append("0");
	HTML.append(OneColor);
	OneColor = Integer.toHexString(GreenComponent);
	if (OneColor.length() < 2)
	    HTML.append("0");
	HTML.append(OneColor);
	OneColor = Integer.toHexString(BlueComponent);
	if (OneColor.length() < 2)
	    HTML.append("0");
	HTML.append(OneColor);
	return HTML.toString();
    }

    /* Add the passed string nonredundantly to this LocusSummary list.
     * We don't care about order, but only add the string if it
     * doesn't already appear in the list. */
    // This is an egregious misuse of this object.  I'm lazy.
    public void AddNonredundantly(String NewItem) {
	LocusSummary    LSRunner = this;
	boolean         NotHereAlready = true;
	while ( (LSRunner.Next != null) && (NotHereAlready) ){
	    if (LSRunner.Next.Name.equals(NewItem)) {
		NotHereAlready = false;
	    }
	    LSRunner = LSRunner.Next;
	}
	if (NotHereAlready) {
	    LSRunner.Next = new LocusSummary();
	    LSRunner.Next.Name = NewItem;
	}
    }

    // Function to sort LocusSummary lists by locus.
    public void	SortListByLocus() {
	//this is just an empty header; the real data starts at the next item
	if (this.Next != null)
	    this.Next = this.Next.Sort(null);
    }

    // INDIVIDUAL FUNCTIONS

    /* Determine whether any of the data sets for this locus should be
     * hidden; if so, this locus should not be shown in Contrast
     * output.  */
    public boolean IsHidden(DataSet Sets) {
	int        Counter;
	long       CurrentPower = 0;
	for (Counter = 0; CurrentPower < this.FoundInPattern;
	     Counter++) {
	    CurrentPower = Math.round(Math.pow(2,Counter));
	    /*
	    System.out.println(new Long(CurrentPower).toString() + "\t" +
			       new Integer(Counter).toString() + "\t" +
			       new Long(FoundInPattern).toString());
	    */
	    if ((this.FoundInPattern & CurrentPower) > 0) {
		if (Sets.IsHidden(Counter)) {
		    return true;
		}
	    }
	}
	return false;
    }

    /* Sort loci so they appear in the proper groups in Contrast.  */
    private LocusSummary Sort(LocusSummary Follower) {
	//Recursive quicksorter
	//Returns first item of sorted list (from those starting at this element)
	//Accepts first item to follow
	LocusSummary		ListAbove = null;
	LocusSummary		ListBelow = null;
	LocusSummary		PlaceHolder;
	LocusSummary		PlaceHolder2;
	PlaceHolder = this.Next;
	//Partition all remaining points of this linked list
	while (PlaceHolder != null) {
	    PlaceHolder2 = PlaceHolder.Next;
	    if (this.FoundCount < PlaceHolder.FoundCount) {
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else if (this.FoundCount > PlaceHolder.FoundCount) {
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    else {
		if (this.FoundInPattern < PlaceHolder.FoundInPattern) {
		    PlaceHolder.Next = ListAbove;
		    ListAbove = PlaceHolder;
		}
		else if (this.FoundInPattern > PlaceHolder.FoundInPattern) {
		    PlaceHolder.Next = ListBelow;
		    ListBelow = PlaceHolder;
		}
		else {
		    if (this.Classification > PlaceHolder.Classification) {
			PlaceHolder.Next = ListAbove;
			ListAbove = PlaceHolder;
		    }
		    else if (this.Classification < PlaceHolder.Classification) {
			PlaceHolder.Next = ListBelow;
			ListBelow = PlaceHolder;
		    }
		    else {
			if (this.Name.compareTo(PlaceHolder.Name) > 0) {
			    PlaceHolder.Next = ListAbove;
			    ListAbove = PlaceHolder;
			}
			else {
			    PlaceHolder.Next = ListBelow;
			    ListBelow = PlaceHolder;
			}
		    }
		}
	    }
	    //Move to next item to be partitioned
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
