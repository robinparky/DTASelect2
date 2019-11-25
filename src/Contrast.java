import java.io.*;
import java.util.*;

//Contrast
//David L. Tabb
//September 11, 2000

/*
 * Contrast program:
 * 1) Read in the Contrast.params file, storing the source directories
 * and criteria sets to be used.
 * 2) Generate a DataSet embodying every combination of source
 * directory and criteria set.
 * 3) Loading each dataset into memory in order, apply the criteria to
 * each set of loci.
 * 4) Generate a LocusSummary summarizing these DataSets.
 * 5) Output the DataSet HTML files and Contrast.html file
 */

public class Contrast {
    IniFile              Configuration = new IniFile();
    LocusSummary         LS = new LocusSummary();
    DataSet              Sets = new DataSet();
    //Store basic info about data sets read from drive
    DataSet              DSList = new DataSet();
    DataSet              MasterSet = null;
    //Store basic info about criteria in use
    DataSet              CList = new DataSet();
    // Store classes of proteins if specified
    Classification       Classifieds = new Classification();
    //Command-line flags
    boolean              VerboseOutput = false;
    boolean              Merge = false;
    boolean              SubsetDatabase = false;
    boolean              Classify = false;
    boolean		 SpeCount = false;

    public static void main(String args[]) throws IOException {
	Contrast             AppObject = new Contrast();
	int                  i;
	AppObject.Configuration.Initialize();
	//Read in command-line arguments
	for(i = 0; i < args.length; i++) {
	    if ( (args[i].equals("-h")) ||
		 (args[i].equals("/?")) ||
		 (args[i].equals("--help")) ) {
		System.out.println(Protein.Version());
		System.out.println("-v\tVerbose output");
		System.out.println("-h\tPrint this help message");
		System.exit(0);
	    }
	    else if (args[i].equals("-v")) {
		//Supply verbose output
		AppObject.VerboseOutput = true;
	    }
	}
	AppObject.Run();
    }

    public void Run() throws IOException {
	DataSet          DSRunner;
	long             CurrentPower = 0;
	String           CurrentName;
	//added by Howard Choi
	String 		 tmpString = "";
	System.out.println("Reading Contrast.params file...");
	try {
	    ReadContrastParams();
	}
	catch (IOException failure) {
	    System.out.println("There was an error reading Contrast.params.");
	    System.exit(0);
	}
	if (Sets.CountSets() > 63) {
	    System.out.println("You have specified too many directories and / or criteria sets.\n");
	    System.out.println("The number of directories multiplied by the number" +
			       " of criteria sets may not exceed 63.\n");
	    System.exit(0);
	}
	DSRunner = Sets.Next;
	while (DSRunner != null) {
	    if (CurrentPower == 0)
		CurrentPower = 1;
	    else CurrentPower *= 2;
	    CurrentName = DSRunner.FriendlyName + "-" + DSRunner.CriteriaName;
	    System.out.println("Reading peptides for " + CurrentName);
	    DSRunner.LoadLociFromFile();
	    System.out.println("\tFiltering at protein level...");
	    DSRunner.ApplyCriteria();
	    DSRunner.LocusList.UngroupIdenticals();
	    DSRunner = DSRunner.Next;
	}
	System.out.println("Comparing protein content...");
	// Find patterns of appearance for each protein
	DSRunner = Sets.Next;
	CurrentPower = 0;
	while (DSRunner != null) {
	    if (CurrentPower == 0)
		CurrentPower = 1;
	    else CurrentPower *= 2;
	    CurrentName = DSRunner.FriendlyName + "-" + DSRunner.CriteriaName;
	    if (MasterSet != null) {
		DSRunner.CullLociNotOnList(MasterSet);
	    }
	    DSRunner.AddToLocusSummary(LS, CurrentPower);
	    DSRunner = DSRunner.Next;
	}
	if (SubsetDatabase) {
	    System.out.println("Subsetting database...");
	    this.SubsetDatabase();
	}
	if (Classify) {
	    System.out.println("Classifying all loci in memory...");
	    LS.StructureByClassification(Classifieds);
	}
	System.out.println("Sorting summary information...");
	LS.SortListByLocus();
	System.out.println("Preparing to print...");
	LS.PrepToPrint(Sets, !Configuration.UseDepthCGI);
	if (LS.Next == null) {
	    System.out.println("No proteins make it past your criteria!\n");
	    System.exit(0);
	}
	System.out.println("Printing Contrast output files...");
	LS.CreateOutputFiles(Sets, MasterSet, DSList,
			     CList, VerboseOutput, SpeCount, Configuration);
	DSRunner = Sets.Next;
	while (DSRunner != null) {
	    if (!DSRunner.Removed) {
		CurrentName = DSRunner.FriendlyName + "-" + DSRunner.CriteriaName;
		System.out.println("Printing " + CurrentName + ".html...");
		DSRunner.SetColumnHeadings();
		DSRunner.LocusList.GroupIdenticalsOld();
		//changed by Howard Choi
		DSRunner.PrintReports(DSRunner.CutoffsString, Configuration, DSRunner.HappyName(), tmpString);
	    }
	    DSRunner = DSRunner.Next;
	}
	System.out.println("Contrast complete.");
    }

    public void MergeDTASelects() {
	SelectCriteria        Cutoffs = new SelectCriteria();
	DataSet               DSRunner = DSList.Next;
	StringBuffer          Header = new StringBuffer();
	StringBuffer          Locus;
	FileReader            InputFilereader;
	File                  OutputFile;
	FileWriter            OutputFileWriter;
	BufferedWriter        Outgoing;
	String                NewLine = "\n";
	String                Tab = "\t";
	long                  LocusCount = 1;
	long                  OUTCount = 0;
	long                  SpecCount = 0;
	StringTokenizer       Parser;
	System.out.println("Merging to create new DTASelect.txt");
	/* Generate a unified DataSet encompassing all proteins and
	 * peptides from this list.  Assume that the list has not yet been
	 * read from disk. */
	// Advance each DataSet's BufferedReader to the first locus.
	// Borrow the header info from the last DataSet.
	System.out.println("\tQueueing up each DTASelect.txt file...");
	try {
	    while (DSRunner != null) {
		Header = new StringBuffer();
		InputFilereader = new FileReader(new File(DSRunner.SourceDirectory,
							  "DTASelect.txt"));
		DSRunner.Incoming = new BufferedReader(InputFilereader);
		DSRunner.WholeLine = DSRunner.Incoming.readLine();
		DSRunner.WholeLine = DSRunner.Incoming.readLine();
		DSRunner.WholeLine = DSRunner.Incoming.readLine();
		while (!DSRunner.WholeLine.startsWith("L\t")) {
		    Header.append(DSRunner.WholeLine + "\n");
		    DSRunner.WholeLine = DSRunner.Incoming.readLine();
		}
		// Create a protein object for this dataset from this L line.
		DSRunner.LocusCursor = new Protein();
		Parser = new StringTokenizer(DSRunner.WholeLine, Tab);
		Parser.nextToken();
		try {
		    DSRunner.LocusCursor.SetTo(Parser);
		}
		catch (NoSuchElementException failure) {
		    System.out.println("The first protein in " +
				       DSRunner.SourceDirectory +
				       "'s DTASelect.txt is garbled.");
		    System.exit(0);
		}
		DSRunner = DSRunner.Next;
	    }
	    System.out.println("\tStarting new DTASelect.txt...");
	    OutputFile = new File("DTASelect.txt");
	    OutputFileWriter = new FileWriter(OutputFile);
	    Outgoing = new BufferedWriter(OutputFileWriter);
	    Outgoing.write(Protein.Version() + NewLine);
	    Outgoing.write(System.getProperty("user.dir") + NewLine);
	    Outgoing.write(Header.toString());
	    // Successively grab the earliest locus from this set of
	    // DataSets and write it to the file
	    Locus = DSList.PullTopLocus();
	    while (Locus != null) {
		Outgoing.write(Locus.toString());
		// System.out.print(Tab + LocusCount + "\r");
		LocusCount++;
		Locus = DSList.PullTopLocus();
	    }
	    // Write count line
	    DSRunner = DSList.Next;
	    while (DSRunner != null) {
		if ( (DSRunner.WholeLine != null) &&
		     (DSRunner.WholeLine.startsWith("C\t")) ) {
		    Parser = new StringTokenizer(DSRunner.WholeLine, "\t");
		    Parser.nextToken();
		    Parser.nextToken();
		    OUTCount += new Long(Parser.nextToken()).longValue();
		    SpecCount += new Long(Parser.nextToken()).longValue();
		}
		DSRunner = DSRunner.Next;
	    }
	    Outgoing.write("C\t" +
			   LocusCount + Tab +
			   OUTCount + Tab +
			   SpecCount + NewLine);
	    // Close each input DTASelect.txt
	    DSRunner = DSList.Next;
	    while (DSRunner != null) {
		DSRunner.Incoming.close();
		DSRunner = DSRunner.Next;
	    }
	    // Close output DTASelect.txt
	    Outgoing.flush();
	    Outgoing.close();
	}
	catch (IOException failure) {
	    System.out.println("IO Error while merging DTASelect.txts:");
	    System.out.println(failure);
	}	
	System.out.println("\tMerging complete: DTASelect.txt created");
    }

    /* Embodies first two stages of Contrast.  This function parses
     * the Contrast.params file, drawing from it the names and
     * particulars of the DTASelect.txt files from which it will draw
     * locus names and the criteria sets which it will apply to each
     * run.  It then combines these to set up the "Sets" field of
     * Contrast.
     */
    public void ReadContrastParams() throws IOException {
        //Read a Contrast.params file, returning the database name
	File             TempDir;
	File             TempFile;
	File             CurrentDirectory = new File(System.getProperty("user.dir"));
	File             ParamsFile = new File(CurrentDirectory, "Contrast.params");
        FileReader       InputFileReader = new FileReader(ParamsFile);
        BufferedReader   Incoming = new BufferedReader(InputFileReader);
        String           LineBuffer = Incoming.readLine();
	String           SecondBuffer;
	String           MasterSpec = null;
	StringBuffer     ExplicitMappings = null;
        StringTokenizer  Parser;
        StringTokenizer  Parser2;
	DataSet          DSBuffer;
	DataSet          SetsRunner = Sets;
	DataSet          DSRunner = DSList;
	DataSet          CRunner = CList;
	int              TokenCounter;
	String           CritParams[];
	StringBuffer     TempBuffer;
	boolean          OkayToAdd;
	// READ INCLUDED DIRECTORIES SECTION OF CONTRAST.PARAMS
	while (!LineBuffer.equals("[Included Directories]")) {
	    LineBuffer = Incoming.readLine();
	}
	System.out.println("Adding directories:");
	LineBuffer = Incoming.readLine();
	Parser = new StringTokenizer(LineBuffer);
	while ( (Parser.hasMoreTokens()) &&
		(!LineBuffer.startsWith("[Criteria Sets]")) ) {
	    LineBuffer = Parser.nextToken();
	    TempBuffer = new StringBuffer();
	    TempBuffer.append(Parser.nextToken());
	    while (Parser.hasMoreTokens()) {
		TempBuffer.append(" ");
		TempBuffer.append(Parser.nextToken());
	    }
	    TempDir = new File(TempBuffer.toString());
	    TempFile = new File(TempDir, "DTASelect.txt");
	    // Check for presence of DTASelect.txt file
	    if (TempFile.canRead()) {
		TempFile = new File(TempDir, "sequest.params");
		DSBuffer = new DataSet();
		DSBuffer.FriendlyName = LineBuffer;
		DSBuffer.SourceDirectory = TempDir.toString();
		DSRunner.Next = DSBuffer;
		DSRunner = DSRunner.Next;
		System.out.println("\tAdded " + TempDir.toString());
	    }
	    else {
		System.out.println(TempFile.toString() + " does not exist.");
		System.exit(0);
	    }
	    LineBuffer = Incoming.readLine();
	    Parser = new StringTokenizer(LineBuffer);
	}
	if (DSRunner == DSList) {
	    System.out.print("Error: no samples included.\nThere cannot be");
	    System.out.print(" a blank line after [Included");
	    System.out.println(" Directories].");
	    System.exit(0);
	}
	// READ CRITERIA SETS SECTION OF CONTRAST.PARAMS
	while (!LineBuffer.startsWith("[Criteria Sets]"))
	    LineBuffer = Incoming.readLine();
	System.out.println("Adding criteria sets:");
	LineBuffer = Incoming.readLine();
	Parser = new StringTokenizer(LineBuffer);
	while ( (LineBuffer != null) && (Parser.hasMoreTokens()) &&
		(!LineBuffer.startsWith("[")) ) {
	    DSBuffer = new DataSet();
	    DSBuffer.CriteriaName = Parser.nextToken();
	    CritParams = new String[Parser.countTokens()];
	    TokenCounter = 0;
	    LineBuffer = "";
	    while (Parser.hasMoreTokens()) {
		CritParams[TokenCounter] = Parser.nextToken();
		LineBuffer += CritParams[TokenCounter] + " ";
		TokenCounter++;
	    }
	    DSBuffer.CutoffsString = LineBuffer;
	    DSBuffer.Cutoffs.SetCriteria(CritParams);
	    CRunner.Next = DSBuffer;
	    CRunner = CRunner.Next;
	    System.out.println("\tAdded " + DSBuffer.CriteriaName);
	    LineBuffer = Incoming.readLine();
	    if (LineBuffer != null)
		Parser = new StringTokenizer(LineBuffer);
	}
	if (CRunner == CList) {
	    System.out.print("Error: no criteria sets included.\nThere cannot be");
	    System.out.print(" a blank line after [Criteria");
	    System.out.println(" Sets].");
	    System.exit(0);
	}
	// READ EXPLICIT MAPPINGS SECTION OF CONTRAST.PARAMS (optional)
	if (LineBuffer != null) {
	    while ( (LineBuffer != null) &&
		    (!LineBuffer.startsWith("[")) ) {
		LineBuffer = Incoming.readLine();
	    }
	    if ((LineBuffer != null) &&
		(LineBuffer.equals("[Explicit Mappings]"))) {
		//Yes, we've found the Explicit Mappings section
		System.out.println("Setting mappings:");
		LineBuffer = Incoming.readLine();
		while ( (LineBuffer != null) &&
			(!LineBuffer.startsWith("[")) ) {
		    Parser = new StringTokenizer(LineBuffer);
		    if (Parser.hasMoreTokens()) {
			SecondBuffer = LineBuffer;
			LineBuffer = Parser.nextToken();
			DSRunner = DSList.Next;
			while ( (DSRunner != null) &&
				(!DSRunner.FriendlyName.equals(LineBuffer)) ) {
			    DSRunner = DSRunner.Next;
			}
			if (DSRunner == null) {
			    System.out.print("Could not find a sample named ");
			    System.out.println(LineBuffer);
			    System.exit(0);
			}
			else {
			    System.out.println("\tAdding mappings for " + LineBuffer);
			    //Okay, this is a legit DS.  Add it to the
			    //list of mappings, and check out the
			    //criteria to be applied to it.
			    if (ExplicitMappings == null) {
				ExplicitMappings = new StringBuffer();
			    }
			    ExplicitMappings.append("\n");
			    ExplicitMappings.append(LineBuffer);
			    while (Parser.hasMoreTokens()) {
				LineBuffer = Parser.nextToken();
				CRunner = CList.Next;
				while ( (CRunner != null) &&
					(!CRunner.CriteriaName.equals(LineBuffer)) ) {
				    CRunner = CRunner.Next;
				}
				if (CRunner == null) {
				    System.out.print("Could not find a criteria set named ");
				    System.out.println(LineBuffer);
				    System.exit(0);
				}
				else {
				    ExplicitMappings.append("\t");
				    ExplicitMappings.append(LineBuffer);
				}
			    }
			}
		    }
		    LineBuffer = Incoming.readLine();
		}
	    }
	}
	// READ OPTIONS SECTION OF CONTRAST.PARAMS (optional)
	if (LineBuffer != null) {
	    while ( (LineBuffer != null) &&
		    (!LineBuffer.startsWith("[Options]")) ) {
		LineBuffer = Incoming.readLine();
	    }
	    if (LineBuffer != null) {
		System.out.println("Setting options:");
		LineBuffer = Incoming.readLine();
		if (LineBuffer != null)
		    Parser = new StringTokenizer(LineBuffer);
	    }
	    while ( (LineBuffer != null) && (Parser.hasMoreTokens() ) ) {
		LineBuffer = Parser.nextToken();
		if (LineBuffer.toUpperCase().equals("VERBOSE")) {
		    System.out.println("\tVerbose mode selected");
		    VerboseOutput = true;
		}
		else if (LineBuffer.toUpperCase().equals("NOBACKGROUND")) {
		    System.out.println("\tNoBackground is an obsolete option.");
		}
		else if (LineBuffer.toUpperCase().equals("MERGE")) {
		    System.out.println("\tDTASelect.txt Merging selected");
		    Merge = true;
		}
		else if (LineBuffer.toUpperCase().equals("DATABASE")) {
		    System.out.println("\tDatabase subsetting selected");
		    SubsetDatabase = true;
		}
		else if (LineBuffer.toUpperCase().equals("CLASS")) {
		    System.out.println("\tContrast locus classification selected");
		    Classify = true;
		}
		else if (LineBuffer.toUpperCase().equals("SPECTRALCOUNT")) {
		    System.out.println("\tSpectral count output selected");
		    SpeCount = true;
		}
		else if (LineBuffer.toUpperCase().equals("MASTER")) {
		    MasterSpec = Parser.nextToken() + "\t" + Parser.nextToken();
		    System.out.println("\tMaster column selected");
		}
		else if (LineBuffer.toUpperCase().equals("HIDE")) {
		    SecondBuffer = Parser.nextToken();
		    System.out.print("\tHiding " + SecondBuffer);
		    TokenCounter = 0;
		    //Mark removed any matching criteria sets
		    CRunner = CList.Next;
		    while (CRunner != null) {
			if (CRunner.CriteriaName.equals(SecondBuffer)) {
			    CRunner.Removed = true;
			    TokenCounter++;
			}
			CRunner = CRunner.Next;
		    }
		    //Mark removed any matching directory
		    DSRunner = DSList.Next;
		    while (DSRunner != null) {
			if (DSRunner.FriendlyName.equals(SecondBuffer)) {
			    DSRunner.Removed = true;
			    TokenCounter++;
			}
			DSRunner = DSRunner.Next;
		    }
		    System.out.println(": " + new Integer(TokenCounter).toString() +
				       " match(es) found.");
		}
		else {
		    System.out.println("\tDid not understand option: " + LineBuffer);
		    System.exit(0);
		}
		LineBuffer = Incoming.readLine();
		if (LineBuffer != null)
		    Parser = new StringTokenizer(LineBuffer);
	    }
	}
	/* Create DataSets for real run.  Create one object for each
           combination of criteria set and directory. */
	CRunner = CList;
	while (CRunner.Next != null) {
	    CRunner = CRunner.Next;
	    DSRunner = DSList.Next;
	    while (DSRunner != null) {
		if (ExplicitMappings != null) {
		    //If mappings have been specified, make sure this
		    //one's okay.
		    OkayToAdd = true;
		    Parser = new StringTokenizer(ExplicitMappings.toString(), "\n");
		    while (Parser.hasMoreTokens()) {
			LineBuffer = Parser.nextToken();
			Parser2 = new StringTokenizer(LineBuffer, "\t");
			LineBuffer = Parser2.nextToken();
			if (LineBuffer.equals(DSRunner.FriendlyName)) {
			    OkayToAdd = false;
			    while (Parser2.hasMoreTokens()) {
				LineBuffer = Parser2.nextToken();
				if (LineBuffer.equals(CRunner.CriteriaName)) {
				    OkayToAdd = true;
				}
			    }
			}
		    }
		}
		else {
		    OkayToAdd = true;
		}
		if (OkayToAdd) {
		    SetsRunner.Next = new DataSet();
		    SetsRunner = SetsRunner.Next;
		    SetsRunner.FriendlyName = DSRunner.FriendlyName;
		    SetsRunner.SourceDirectory = DSRunner.SourceDirectory;
		    SetsRunner.CriteriaName = CRunner.CriteriaName;
		    SetsRunner.Cutoffs = CRunner.Cutoffs;
		    SetsRunner.CutoffsString = CRunner.CutoffsString;
		    SetsRunner.Removed = DSRunner.Removed || CRunner.Removed;
		}
		DSRunner = DSRunner.Next;
	    }
	}
	if (MasterSpec != null) {
	    boolean     FoundIt = false;
	    Parser = new StringTokenizer(MasterSpec);
	    LineBuffer = Parser.nextToken();
	    SecondBuffer = Parser.nextToken();
	    SetsRunner = Sets.Next;
	    while ( (SetsRunner != null) &&
		    (!FoundIt) ) {
		if ( (SetsRunner.FriendlyName.equals(LineBuffer)) &&
		     (SetsRunner.CriteriaName.equals(SecondBuffer)) ) {
		    FoundIt = true;
		    MasterSet = SetsRunner;
		}
		SetsRunner = SetsRunner.Next;
	    }
	    if (!FoundIt) {
		System.out.println("Did not find master column from directory " +
				   LineBuffer + " and criteria " + SecondBuffer);
		System.exit(0);
	    }
	}
	// Create unified DTASelect.txt, if selected
	if (Merge) {
	    this.MergeDTASelects();
	}
    }

    /* Create a sequence database that contains only the proteins
     * observed in the passed LocusSummary. */
    public void SubsetDatabase() {
	try {
	    File                  DBFile;
	    FileReader            InputFileReader;
	    BufferedReader        Incoming;
	    File                  OutputFile = new File("Contrast.fasta");
	    FileWriter            OutputFileWriter = new FileWriter(OutputFile);
	    BufferedWriter        Outgoing = new BufferedWriter(OutputFileWriter);
	    int                   LongestLocus = Configuration.LocusLengthCutoff;
	    String                LineBuffer;
	    StringTokenizer       Parser;
	    String                CurrentLocusName;
	    int                   MemoryLocusCount = 0;
	    int                   OriginalLocusCount = 0;
	    int                   NewLocusCount = 0;
	    LocusSummary          DatabaseList = new LocusSummary();
	    LocusSummary          TempLocList = new LocusSummary();
	    LocusSummary          LSRunner = LS.Next;
	    LocusSummary          DBRunner = TempLocList;
	    char                  NewLine = '\n';
	    DataSet               DSRunner = Sets.Next;
	    // Make a nonredundant list of databases to search
	    while (DSRunner != null) {
		DatabaseList.AddNonredundantly(DSRunner.SequestParams.DBName);
		DSRunner = DSRunner.Next;
	    }
	    // Make a copy of the locus list to our temp copy
	    while (LSRunner != null) {
		DBRunner.Next = new LocusSummary();
		DBRunner = DBRunner.Next;
		DBRunner.Name = LSRunner.Name;
		MemoryLocusCount++;
		LSRunner = LSRunner.Next;
	    }
	    System.out.println("\t" + MemoryLocusCount + " total proteins to be matched.");
	    DBRunner = DatabaseList.Next;
	    // Loop through all listed databases unless we're done matching
	    while ( (DBRunner != null) && (MemoryLocusCount > 0) ) {
		System.out.println("\tNow matching against " + DBRunner.Name);
		OriginalLocusCount = 0;
		NewLocusCount = 0;
		DBFile = new File(DBRunner.Name);
		InputFileReader = new FileReader(DBFile);
		Incoming = new BufferedReader(InputFileReader);
		LineBuffer = Incoming.readLine();
		while (LineBuffer != null) {
		    Parser = new StringTokenizer(LineBuffer, "\t >");
		    if (Parser.hasMoreTokens()) {
			// We're on a > line, and the first token is the locus name
			OriginalLocusCount++;
			CurrentLocusName = Parser.nextToken();
			// Trim this name down to size of loci in memory
			if (CurrentLocusName.length() > LongestLocus) {
			    CurrentLocusName = CurrentLocusName.substring(0, LongestLocus);
			}
			// Search the list in memory to determine if this one matches
			LSRunner = TempLocList;
			while ((LSRunner.Next != null) && (LSRunner.Next.Name.compareTo(CurrentLocusName) < 0)) {
			    LSRunner = LSRunner.Next;
			}
			if ( (LSRunner.Next != null) && (LSRunner.Next.Name.equals(CurrentLocusName)) ) {
			    NewLocusCount++;
			    // Excise this locus from our temporary list
			    LSRunner.Next = LSRunner.Next.Next;
			    // Copy this locus to the new database
			    Outgoing.write(LineBuffer + NewLine);
			    LineBuffer = Incoming.readLine();
			    while ( (LineBuffer != null) && 
				    (LineBuffer.length() > 0) &&
				    (LineBuffer.charAt(0) != '>') ) {
				Outgoing.write(LineBuffer + NewLine);
				LineBuffer = Incoming.readLine();
			    }
			}
			else {
			    LineBuffer = Incoming.readLine();
			    // Skip to the next locus
			    while ( (LineBuffer != null) && 
				    (LineBuffer.length() > 0) &&
				    (LineBuffer.charAt(0) != '>') ) {
				LineBuffer = Incoming.readLine();
			    }
			}
		    }
		    else {
			LineBuffer = Incoming.readLine();
		    }
		}
		Incoming.close();
		System.out.println("\t" + OriginalLocusCount + "\tloci in this database");
		System.out.println("\t" + NewLocusCount + "\tloci matched those left in the list");
		MemoryLocusCount -= NewLocusCount;
		System.out.println("\t" + MemoryLocusCount + "\tmore loci remain to be matched");
		DBRunner = DBRunner.Next;
	    }
	    Outgoing.flush();
	    Outgoing.close();
	}
	catch (IOException failure) {
	    System.out.println("Error while trying to create the subset database.");
	    System.out.println(failure);
	}
    }
    
}
