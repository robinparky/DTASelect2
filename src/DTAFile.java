
import java.io.*;
import java.util.*;

//DTAFile object
//For DTASelect
//Created by Dave Tabb
//Begun 8/11/2000

public class DTAFile {
    DTAFile              Next;
    String               FileName;
    String		 RootFileName;
    String               Subdirectory;
    String               Sequence;
    float                XCorr;
    float                SpScore;
    float				 tertiaryScore;
    float 				 pvalue;
    float 				 evalue;
    float		 PepConf = 0.0f;
    double		 Probability = 0.0;
    float		 PepFP = 1.0f;
    float                CalcPreMass = 0.0f;
    float                PrecursorMass = 0.0f;
    double               Raw_PPM_Offset = 0.0;
    double               Adjusted_PPM_Offset = 0.0;
	double               Shifted_PPM_Offset = 0.0;
    float                TotalIntensity;
    float                IonProportion;
    float                DeltCN = 0.0f;
    boolean              Modified;
    int                  NMods = 0;
    boolean              UniqueToLocus = false;
    char                 Validated = 'U';
    byte                 Tryptic;
    byte                 ChargeState;
    int			 ScanNumber;
    byte                 EquivSeq = 1;
    short                Redundancy = 1;
    short                Sp;
    int                  SequencePosition;
    List<Integer>		sequencePositionList = new ArrayList<>();
    int			 TrypticSites;
    boolean		 IsDecoy = false;
    int id;
	private static int ID_COUNTER = 0;
    /* INDIVIDUAL FUNCTIONS
     * The following functions should be called on individual DTAFiles
     * that may or may not be part of lists
     */

    /* Clone copies all fields of one DTAFile in creating another.
     * Since the Next pointer is not set, the new DTAFile is not part
     * of the chain in which this DTAFile sits.
     */

    public DTAFile()
	{
		id = ID_COUNTER++;
	}
    public DTAFile Clone() {
	DTAFile       Copy = new DTAFile();
	Copy.FileName = this.FileName;
	Copy.RootFileName = this.RootFileName;
	Copy.Subdirectory = this.Subdirectory;
	Copy.XCorr = this.XCorr;
	Copy.PepConf = this.PepConf;
	Copy.PepFP = this.PepFP;
	Copy.SpScore = this.SpScore;
	Copy.CalcPreMass = this.CalcPreMass;
	Copy.PrecursorMass = this.PrecursorMass;
        Copy.Raw_PPM_Offset = this.Raw_PPM_Offset;
        Copy.Adjusted_PPM_Offset = this.Adjusted_PPM_Offset;
	Copy.Sequence = this.Sequence;
	Copy.Modified = this.Modified;
	Copy.NMods = this.NMods;
	Copy.Tryptic = this.Tryptic;
	Copy.UniqueToLocus = this.UniqueToLocus;
	Copy.ChargeState = this.ChargeState;
	Copy.ScanNumber = this.ScanNumber;
	Copy.TotalIntensity = this.TotalIntensity;
	Copy.Sp = this.Sp;
	Copy.IonProportion = this.IonProportion;
	Copy.DeltCN = this.DeltCN;
	Copy.SequencePosition = this.SequencePosition;
	Copy.EquivSeq = this.EquivSeq;
	Copy.Redundancy = this.Redundancy;
	Copy.Validated = this.Validated;
	Copy.IsDecoy = this.IsDecoy;
	Copy.TrypticSites = this.TrypticSites;
	return Copy;
    }

    /* Set the fields of this DTAFile to the contents of this
     * StringTokenizer (which delimits by tabs and has the initial "D"
     * truncated), drawn from a DTASelect.txt file.  Called by both
     * Contrast Merger and DataSet.ReadFromFile */
    public void OldSetTo(StringTokenizer Parser) throws NoSuchElementException {
	this.FileName = Parser.nextToken();
	this.Subdirectory = Parser.nextToken();
	this.XCorr = new Float(Parser.nextToken()).floatValue();
	this.DeltCN = new Float(Parser.nextToken()).floatValue();
	this.PepConf = new Float(Parser.nextToken()).floatValue();
	this.PepFP = new Float(Parser.nextToken()).floatValue();
	this.PrecursorMass = new Float(Parser.nextToken()).floatValue();
	this.TotalIntensity = new Float(Parser.nextToken()).floatValue();
	this.Sp = new Short(Parser.nextToken()).shortValue();
	this.IonProportion = new Float(Parser.nextToken()).floatValue();
	this.Sequence = Parser.nextToken();
	this.SequencePosition = new Integer(Parser.nextToken()).intValue();
	this.Tryptic = new Byte(Parser.nextToken()).byteValue();
	this.UniqueToLocus = new Boolean(Parser.nextToken()).booleanValue();
	if (Parser.hasMoreTokens())
	    this.Validated = Parser.nextToken().charAt(0);
    }

    /* Set the fields of this DTAFile to the contents of this
     * StringTokenizer (which delimits by tabs and has the initial "D"
     * truncated), drawn from a DTASelect.txt file.  Called by both
     * Contrast Merger and DataSet.ReadFromFile */
    public void CurrentSetTo(StringTokenizer Parser, SelectCriteria Cutoffs) throws NoSuchElementException {
	this.FileName = Parser.nextToken();
	String[] ParsedLine = this.FileName.split("\\.");
        this.ChargeState = new Byte(ParsedLine[3]).byteValue();
	this.ScanNumber = new Integer(ParsedLine[1]).intValue();
	this.RootFileName = ParsedLine[0];
	this.Subdirectory = Parser.nextToken();
	this.XCorr = new Float(Parser.nextToken()).floatValue();
	this.DeltCN = new Float(Parser.nextToken()).floatValue();
//	this.PepConf = new Float(Parser.nextToken()).floatValue();
	this.PepFP = 1.0f - new Float(Parser.nextToken()).floatValue();
	this.PrecursorMass = new Float(Parser.nextToken()).floatValue();
	this.CalcPreMass = new Float(Parser.nextToken()).floatValue();
	this.GetMassOffsets(Cutoffs);
	this.TotalIntensity = new Float(Parser.nextToken()).floatValue();
	this.Sp = new Short(Parser.nextToken()).shortValue();
	this.SpScore = new Float(Parser.nextToken()).floatValue();
	this.IonProportion = new Float(Parser.nextToken()).floatValue();
	this.Sequence = Parser.nextToken();
	this.TrypticSites = this.DetermineMissedTryptic(Cutoffs);
	String seqPos = Parser.nextToken();
	if(seqPos.contains(","))
	{
		String [] posStrArr = seqPos.split(",");
		for(String posStr: posStrArr)
		{
			sequencePositionList.add(Integer.parseInt(posStr));
		}
		this.SequencePosition = sequencePositionList.get(0);
	}
	else
	{
		this.SequencePosition = Integer.parseInt(seqPos);
	}
	this.Tryptic = new Byte(Parser.nextToken()).byteValue();
	this.DetermineModified();
	this.UniqueToLocus = new Boolean(Parser.nextToken()).booleanValue();
	this.EquivSeq = new Byte(Parser.nextToken()).byteValue();
	this.Validated = Parser.nextToken().charAt(0);
    }

    /* Since some database sequences include non-letter characters,
     * this function culls out these characters and returns a string
     * containing just the sequence (in upper case characters).
     */
    public static String JustLettersFrom(String Input) {
	StringBuffer     Result = new StringBuffer();
	int              Counter;
	for (Counter = 0; Counter < Input.length(); Counter++) {
	    if (Character.isLetter(Input.charAt(Counter)))
		Result.append(Input.charAt(Counter));
	}
	return (Result.toString()).toUpperCase();
    }

    public static String JustLettersFromTrimmed(String Input) {
	StringBuffer     Result = new StringBuffer();
	int              Counter;
	for (Counter = 2; Counter < Input.length() - 2; Counter++) {
	    if (Character.isLetter(Input.charAt(Counter)))
		Result.append(Input.charAt(Counter));
	}
	return (Result.toString()).toUpperCase();
    }

    //SEQUENCE REPORTING FUNCTIONS

    /* The Sequence field of each DTAFile usually takes on the form
     * "X.XXXX.X," where the character before the first period is the
     * last residue before the N-terminal cleavage site and the final
     * character represents the first residue after the C-terminalF
     * cleavage site.  This function returns only the identified
     * peptide without the context information.  Note that this
     * returned sequence will not include modification characters.  */
    public String TrimmedSequence() {
	StringBuffer  Result = new StringBuffer();
	int           Counter;
	for (Counter = 2; Counter < Sequence.length()-2; Counter ++)
	    if (Character.isUpperCase(Sequence.charAt(Counter)))
		Result.append(Sequence.charAt(Counter));
	return Result.toString();
    }

    public String BetweenDots() {
	return this.Sequence.substring(2, Sequence.length() - 2);
    }

    // get information for this .dta file for use with the new Show CGI
    public String ShowString() {
	StringBuffer     Returned = new StringBuffer();
	StringTokenizer  Parser = new StringTokenizer(this.FileName, ".");
	String           DatFile = Parser.nextToken() + ".";
	String           MyScanNumber;

        while (this.Subdirectory.startsWith(DatFile)) {
            DatFile = DatFile.concat(Parser.nextToken() + ".");
        }
	MyScanNumber = Parser.nextToken();
	DatFile = DatFile.substring(0,DatFile.length() - 1);
	Returned.append(DatFile);
	if (!this.Subdirectory.equals(DatFile)) {
	    Returned.append("&Sd=");
	    Returned.append(this.Subdirectory);
	}
	Returned.append("&Sc=");
	Returned.append(MyScanNumber);
	Returned.append("&Sq=");
	if (this.Modified) {
	    int              Counter;
	    StringBuffer     PassMods = new StringBuffer();
	    int              SeqLength = Sequence.length() - 2;
	    int              CurrentResidue = 0;
	    char             CurrentCharacter;
	    for (Counter = 2; Counter < SeqLength; Counter ++) {
		CurrentCharacter = this.Sequence.charAt(Counter);
		if (Character.isUpperCase(CurrentCharacter)) {
		    // Current character is a amino acid residue!
		    Returned.append(CurrentCharacter);
		    CurrentResidue++;
		}
		else {
		    // Current character is a modification!
		    PassMods.append("&M");
		    PassMods.append(CurrentResidue);
		    PassMods.append('=');
		    PassMods.append(new Integer((int)(CurrentCharacter)).toString());
		}
	    }
	    Returned.append(PassMods.toString());
	}
	else {
	    Returned.append(this.Sequence.substring(2, Sequence.length() - 2));
	}
	Returned.append("&Z=");
	Returned.append(this.ChargeState);
	Returned.append("\">");
	Returned.append(this.FileName);
	Returned.append("</A>");
	return (Returned.toString());
    }

    public String DisplayIonsString(ParamsFile SequestParams) {
	StringBuffer   Returned = new StringBuffer("&DSite=");
	int            MassInUse = -1;
	char           Symbol[];
	float          Mass[];
	String         ModString = this.SymbolModString();
	int            SeqLength = ModString.length();
	int            Looper;
	int            ModLooper;
	int            ThisMod;
	char           ThisSymbol;
	Symbol = new char[3];
	Mass = new float[3];
	for (Looper = 0; Looper < SeqLength; Looper++) {
	    ThisSymbol = ModString.charAt(Looper);
	    if (ThisSymbol != ' ') {
		ThisMod = -1;
		// Check to see if this modification has been used
		for (ModLooper = 0; ModLooper < MassInUse + 1; ModLooper++) {
		    if (Symbol[ModLooper] == ThisSymbol) {
			ThisMod = ModLooper;
			ModLooper = MassInUse;
		    }
		}
		if (ThisMod == -1) {
		    MassInUse++;
		    Symbol[MassInUse] = ThisSymbol;
		    Mass[MassInUse] = SequestParams.DiffMods.getMassShiftFor(ThisSymbol);
		    ThisMod = MassInUse;
		}
		Returned.append(new Integer(MassInUse+1).toString());
	    }
	    else {
		Returned.append('0');
	    }
	}
	for (ModLooper = 0; ModLooper < MassInUse + 1; ModLooper++) {
	    Returned.append("&DMass");
	    Returned.append(ModLooper+1);
	    Returned.append('=');
	    Returned.append(Mass[ModLooper]);
	}
	return Returned.toString();
    }

    /* Prepare a modification string to describe the post-translation
     * modifications found in this identification.  Include a space
     * for each unmodified residue, and include the appropriate symbol
     * for each modified residue (the symbols are enumerated in
     * sequest.params).
     */
    public String SymbolModString() {
	StringBuffer   Result = new StringBuffer();
	int            Counter;
	char           Buffer;
	char           Buffer2;
	for (Counter = 2; Counter < Sequence.length()-2; Counter ++) {
	    Buffer = Sequence.charAt(Counter);
	    if (Character.isLetter(Buffer)) {
		Buffer2 = Sequence.charAt(Counter+1);
		if ( (Character.isLetter(Buffer2)) || (Buffer2 == '.') || (Buffer2 == '-') ) {
		    Result.append(' ');
		}
		else {
		    Result.append(Buffer2);
		}
	    }
	}
	return Result.toString();
    }

    /* If this peptide has a modification character in its sequence,
     * this will set the Modified field to true.  Otherwise it will be
     * set to false.  This is the basis of the "modified only" filter.
     */
    public void DetermineModified() {
	String   ModString = this.SymbolModString();
	int      Looper;
	int      NModsSymbol = 0;
	int      NModsParanthesis = 0;
	int      SeqLength = ModString.length();
	this.Modified = false;
	for (Looper = 0; Looper < SeqLength; Looper++) {
	    if (ModString.charAt(Looper) != ' ') {
		this.Modified = true;
		NModsSymbol++;
	    }
	}
	for (Looper = 2; Looper < Sequence.length()-2; Looper++) {
            if (Sequence.charAt(Looper) == '(') {
                this.Modified = true;
		NModsParanthesis++;
            }
        }
	this.NMods = (NModsParanthesis > NModsSymbol) ? NModsParanthesis : NModsSymbol;
    }

    /* A peptide is said to be tryptic iff both cleavage sites
     * correspond to either termini or cuts after K or R.
     */
    public void DetermineTryptic(SelectCriteria Cutoffs) {
	String   JustSequence = TrimmedSequence();
	if (JustSequence.length() == 0) {
	    System.out.println("Warning: " + this.FileName + " reports no sequence!");
	    return;
	}
	char     LastRealChar = JustSequence.charAt(JustSequence.length()-1);
	char     LastChar = Sequence.charAt(Sequence.length()-1);
	char     FirstChar = Sequence.charAt(0);
	char     FirstRealChar = JustSequence.charAt(0);
	boolean  BeforeNTerminus = false;
	boolean  AfterNTerminus = false;
	boolean  BeforeCTerminus = false;
	boolean  AfterCTerminus = false;
	if (Cutoffs.EnzymeCutsAfter.length() > 0) {
	    if (Cutoffs.EnzymeCutsAfter.indexOf(FirstChar) != -1) {
		BeforeNTerminus = true;
	    }
	    if (Cutoffs.EnzymeCutsAfter.indexOf(LastRealChar) != -1) {
                BeforeCTerminus = true;
            }
	}
	if (Cutoffs.EnzymeCutsBefore.length() > 0) {
            if (Cutoffs.EnzymeCutsBefore.indexOf(FirstRealChar) != -1) {
                AfterNTerminus = true;
            }
            if (Cutoffs.EnzymeCutsBefore.indexOf(LastChar) != -1) {
                AfterCTerminus = true;
            }
        }
	boolean  NTryptic = (BeforeNTerminus ||
			     AfterNTerminus ||
			     (FirstChar == '-'));
	boolean  CTryptic = (BeforeCTerminus ||
			     AfterCTerminus ||
			     (LastChar == '-'));
	Tryptic = 0;
	if (NTryptic && CTryptic)
	    Tryptic = 2;
	else if (NTryptic || CTryptic)
	    Tryptic = 1;
    }

    public int DetermineMissedTryptic(SelectCriteria Cutoffs) {
	String   JustSequence = TrimmedSequence();
	int Looper;
	int NMissedTryptic = 0;

	for (Looper = 0; Looper < JustSequence.length()-1; Looper++) {
	    if ((JustSequence.charAt(Looper) == 'R') ||
	        (JustSequence.charAt(Looper) == 'K'))
		NMissedTryptic++;
	}
	return NMissedTryptic;
    }

    /* Return the full file name of this .out file.  The directory is
     * assumed to be the portion of the filename before its first
     * period.
     */
    public String CanonicalName(String CurrentDir) {
	try {
	    File		CurrentDirectory = new File(CurrentDir, Subdirectory);
	    String 		DTAFile = new File(CurrentDirectory,
						   FileName+ ".out").getCanonicalPath();
	    return (DTAFile);
	}
	catch (IOException failure) {
	    System.out.println("CN: I take exception to that\t" +
			       CurrentDir + "\t" + FileName);
	    return "junk";
	}
    }

    /* Return the full file name of this .dta file.  The directory is
     * assumed to be the portion of the filename before its first
     * period.
     */
    public String DTACanonicalName(String CurrentDir) {
	try {
	    File		CurrentDirectory = new File(CurrentDir, Subdirectory);
	    String 		DTAFile = new File(CurrentDirectory,
						   FileName+ ".dta").getCanonicalPath();
	    return (DTAFile);
	}
	catch (IOException failure) {
	    System.out.println("DCN: I take exception to that" + CurrentDir);
	    return "junk";
	}
    }

    /* Loop through each residue in the database sequence for this
     * locus and accumulate molecular weight.  Use a simple algorithm
     * to estimate pI.
     */
    public float CalculatepI() {
	String         Sequence = this.TrimmedSequence();
	int            Length = Sequence.length();
	int            Looper;
	int            CountLys = 0;
	int            CountArg = 0;
	int            CountHis = 0;
	int            CountAsp = 0;
	int            CountGlu = 0;
	int            CountCys = 0;
	int            CountTyr = 0;
	float          CurrentPH = 7.0f;
	float          CurrentJump = 3.5f;
	float          CurrentCharge;
	float          LastCharge = 0;
	float          MWAccum = 0f;
	char           CurrentResidue;
	if (Length > 0) {
	    for (Looper = 0; Looper < Length; Looper++) {
		CurrentResidue = Sequence.charAt(Looper);
		switch (CurrentResidue) {
		case 'A':
		    MWAccum += 71.0;
		    break;
		case 'C':
		    MWAccum += 103.0;
		    CountCys ++;
		    break;
		case 'D':
		    MWAccum += 115.0;
		    CountAsp ++;
		    break;
		case 'E':
		    MWAccum += 129.0;
		    CountGlu ++;
		    break;
		case 'F':
		    MWAccum += 147.0;
		    break;
		case 'G':
		    MWAccum += 57.0;
		    break;
		case 'H':
		    MWAccum += 137.0;
		    CountHis ++;
		    break;
		case 'I':
		    MWAccum += 113.0;
		    break;
		case 'K':
		    MWAccum += 128.0;
		    CountLys ++;
		    break;
		case 'L':
		    MWAccum += 113.0;
		    break;
		case 'M':
		    MWAccum += 131.0;
		    break;
		case 'N':
		    MWAccum += 114.0;
		    break;
		case 'P':
		    MWAccum += 97.0;
		    break;
		case 'Q':
		    MWAccum += 128.0;
		    break;
		case 'R':
		    MWAccum += 156.0;
		    CountArg ++;
		    break;
		case 'S':
		    MWAccum += 87.0;
		    break;
		case 'T':
		    MWAccum += 101.0;
		    break;
		case 'V':
		    MWAccum += 99.0;
		    break;
		case 'W':
		    MWAccum += 186.0;
		    break;
		case 'Y':
		    MWAccum += 176.0;
		    CountTyr ++;
		    break;
		}
	    }
	    /* Use a bracketing strategy to identify the isoelectric
	     * point.  Calculate charge at pH of 7, and then move up
	     * 3.5 or down 3.5 as necessary.  Make each successive
	     * move up or down only half as large.  Keep going until
	     * two successive charges reported match to one place past
	     * the decimal point.
	     */
	    CurrentCharge = Protein.ChargeAtPH(CurrentPH, CountLys,
					       CountArg, CountHis,
					       CountAsp, CountGlu,
					       CountCys, CountTyr);
	    while (Protein.RoundTo(CurrentCharge,1) != Protein.RoundTo(LastCharge,1)) {
		//		System.out.println("pH:\t" + new Float(CurrentPH).toString() 
		//              + "\tCharge\t" + new Float(CurrentCharge).toString());
		if (CurrentCharge > 0)
		    CurrentPH += CurrentJump;
		else
		    CurrentPH -= CurrentJump;
		CurrentJump /= 2;
		LastCharge = CurrentCharge;
		CurrentCharge = Protein.ChargeAtPH(CurrentPH,
						   CountLys, CountArg,
						   CountHis, CountAsp,
						   CountGlu, CountCys,
						   CountTyr);
		if ( (CurrentPH > 14) || (CurrentPH < 0) ) {
		    System.out.println("pI can't be figured for " + FileName);
		    System.exit(0);
		}
	    }
	}
	return CurrentPH;
    }

    public float CalculateKyteDoolittle() {
        String         Sequence = this.TrimmedSequence();
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                case 'A':
                    MWAccum += 18.0;
                    break;
                case 'C':
                    MWAccum += 25.0;
                    break;
                case 'D':
                    MWAccum -= 35.0;
                    break;
                case 'E':
                    MWAccum -= 35.0;
                    break;
                case 'F':
                    MWAccum += 28.0;
                    break;
                case 'G':
                    MWAccum -= 4.0;
                    break;
                case 'H':
                    MWAccum -= 32.0;
                    break;
                case 'I':
                    MWAccum += 45.0;
                    break;
                case 'K':
                    MWAccum -= 39.0;
                    break;
                case 'L':
                    MWAccum += 38.0;
                    break;
                case 'M':
                    MWAccum += 19.0;
                    break;
                case 'N':
                    MWAccum -= 35.0;
                    break;
                case 'P':
                    MWAccum -= 16.0;
                    break;
                case 'Q':
                    MWAccum -= 35.0;
                    break;
                case 'R':
                    MWAccum -= 45.0;
                    break;
                case 'S':
                    MWAccum -= 8.0;
                    break;
                case 'T':
                    MWAccum -= 7.0;
                    break;
                case 'V':
                    MWAccum += 42.0;
                    break;
                case 'W':
                    MWAccum -= 9.0;
                    break;
                case 'Y':
                    MWAccum -= 13.0;
                    break;
                }
            }
        }
        return MWAccum/10.0f;
    }

    public float CalculateBullBreese() {
        String         Sequence = this.TrimmedSequence();
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0.0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                case 'A':
                    MWAccum += 610;
                    break;
                case 'C':
                    MWAccum += 360;
                    break;
                case 'D':
                    MWAccum += 610;
                    break;
                case 'E':
                    MWAccum += 510;
                    break;
                case 'F':
                    MWAccum += -1520;
                    break;
                case 'G':
                    MWAccum += 810;
                    break;
                case 'H':
                    MWAccum +=690;
                    break;
                case 'I':
                    MWAccum += -1450;
                    break;
                case 'K':
                    MWAccum += 460;
                    break;
                case 'L':
                    MWAccum += -1650;
                    break;
                case 'M':
                    MWAccum += -660;
                    break;
                case 'N':
                    MWAccum += 890;
                    break;
                case 'P':
                    MWAccum += -170;
                    break;
                case 'Q':
                    MWAccum += 970;
                    break;
                case 'R':
                    MWAccum += 690;
                    break;
                case 'S':
                    MWAccum += 420;
                    break;
                case 'T':
                    MWAccum += 290;
                    break;
                case 'V':
                    MWAccum += -750;
                    break;
                case 'W':
                    MWAccum += -1200;
                    break;
                case 'Y':
                    MWAccum += -1430;
                    break;
                }
            }
        }
        return MWAccum;
    }

    public float CalculateHPLCpH34() {
        String         Sequence = this.TrimmedSequence();
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0.0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                case 'A':
                    MWAccum += 42;
                    break;
                case 'C':
                    MWAccum += 84;
                    break;
                case 'D':
                    MWAccum += -51;
                    break;
                case 'E':
                    MWAccum += -37;
                    break;
                case 'F':
                    MWAccum += 174;
                    break;
                case 'G':
                    MWAccum += 0;
                    break;
                case 'H':
                    MWAccum += -228;
                    break;
                case 'I':
                    MWAccum += 181;
                    break;
                case 'K':
                    MWAccum += -203;
                    break;
                case 'L':
                    MWAccum += 180;
                    break;
                case 'M':
                    MWAccum += 118;
                    break;
                case 'N':
                    MWAccum += -103;
                    break;
                case 'P':
                    MWAccum += 86;
                    break;
                case 'Q':
                    MWAccum += -96;
                    break;
                case 'R':
                    MWAccum += -156;
                    break;
                case 'S':
                    MWAccum += -64;
                    break;
                case 'T':
                    MWAccum += -26;
                    break;
                case 'V':
                    MWAccum += 134;
                    break;
                case 'W':
                    MWAccum += 146;
                    break;
                case 'Y':
                    MWAccum += 51;
                    break;
                }
            }
        }
        return MWAccum/100.0f;
    }

    // Print a header for an outdated debug function.
    public static void DebugPrintHeader() {
	System.out.println("FileName\tXCorr\tSequence\tChargeState\tTotalIntensity\tSpRank\tIonProportion\tDeltCN\tConf");
    }

    // Print this DTA's stats for the GUI
    public String GetTabDelimitedFields() {
	return (
		FileName +
		"\t" + new Float(Protein.RoundTo(XCorr,3)).toString() +
		"\t" + new Float(Protein.RoundTo(DeltCN,4)).toString() +
		"\t" + new Float(Protein.RoundTo(PrecursorMass,1)).toString() +
		"\t" + new Float(Protein.RoundTo(CalcPreMass,1)).toString() +
		"\t" + new Integer(Sp).toString() +
		"\t" + new Float(Protein.RoundTo(SpScore,3)).floatValue() +
		"\t" + new Float(Protein.RoundTo(100f*IonProportion,1)).toString() +
		"%\t" + new Integer(Redundancy).toString() +
		"\t" + Sequence +
		"\t" + (UniqueToLocus?"*":" ")
		);
    }

    /* Return the "D" line for DTASelect.txt corresponding to this
     * identification.  Used both by DataSet.PrintTXT and
     * DataSet.PullTopLocus */
    public String GetDLine() {
	String Tab = "\t";
	String sequencePostionStr;
	if(this.sequencePositionList.size()>1)
	{
		StringBuilder sb = new StringBuilder();
		for(int i: sequencePositionList)
		{
			sb.append(i).append(",");
		}
		sequencePostionStr = sb.toString();
	}
	else
	{
		sequencePostionStr = Integer.toString(this.SequencePosition);
	}
	return ("D\t" +
		this.FileName +
		Tab + this.Subdirectory +
		Tab + new Float(this.XCorr).toString() +
		Tab + new Float(this.DeltCN).toString() +
//		Tab + new Float(this.PepConf).toString() +
		Tab + new Float(1.0 - this.PepFP).toString() +
		Tab + new Float(this.PrecursorMass).toString() +
		Tab + new Float(this.CalcPreMass).toString() +
		Tab + new Float(this.TotalIntensity).toString() +
		Tab + new Integer(this.Sp).toString() +
		Tab + new Float(this.SpScore).toString() +
		Tab + new Float(this.IonProportion).toString() +
		Tab + this.Sequence +
		Tab + sequencePostionStr +
		Tab + new Byte(this.Tryptic).toString() +
		Tab + new Boolean(this.UniqueToLocus).toString() +
		Tab + new Byte(this.EquivSeq).toString() +
		Tab + new Character(this.Validated).toString() +
		"\n");
    }

    public String PrintXML(ParamsFile Params, String DirName) {
	return (
		"\t\t<DTASelect:Filename>" +
		FileName +
		"</DTASelect:Filename>\n" +
	        "\t\t<DTASelect:XCorr>" +
		new Float(XCorr).toString() +
		"</DTASelect:XCorr>\n" +
		"\t\t<DTASelect:DeltCN>" +
		new Float(DeltCN).toString() +
		"</DTASelect:DeltCN>\n" +
		"\t\t<DTASelect:PepConf>" +
                new Float(PepConf).toString() +
                "</DTASelect:PepConf>\n" +
		"\t\t<DTASelect:PepFP>" +
                new Float(PepFP).toString() +
                "</DTASelect:PepFP>\n" +
		"\t\t<DTASelect:PrecursorMass>" +
		new Float(PrecursorMass).toString() +
		"</DTASelect:PrecursorMass>\n" +
		"\t\t<DTASelect:TotalIntensity>" +
		new Float(TotalIntensity).toString() +
		"</DTASelect:TotalIntensity>\n" +
		"\t\t<DTASelect:SpRank>" +
		new Integer(Sp).toString() +
		"</DTASelect:SpRank>\n" +
		"\t\t<DTASelect:FragmentIonPercentage>" +
		new Float(Protein.RoundTo(100f*IonProportion,1)).toString() +
		"</DTASelect:FragmentIonPercentage>\n" +
		"\t\t<DTASelect:CopyCount>" +
		new Integer(Redundancy).toString() +
		"</DTASelect:CopyCount>\n" +
	        "\t\t<DTASelect:Sequence>" +
		Sequence +
		"</DTASelect:Sequence>\n" +
		"\t\t<DTASelect:Unique>" +
		UniqueToLocus +
		"</DTASelect:Unique>\n"
		);
    }

    // Yet another crufty debugging print all fields function
    public void DebugPrintThis() {
	System.out.println(Subdirectory + "/" + FileName + "\t" + new Float(XCorr).toString() +
			   "\t" + Sequence + "\t" + new Integer(ChargeState).toString() +
			   "\t" + new Float(TotalIntensity).toString() +
			   "\t" + new Integer(Sp).toString() +
			   "\t" + new Float(IonProportion).toString() +
			   "\t" + new Float(DeltCN).toString() +
			   "\t" + new Float(PepConf).toString() +
			   "\t" + new Float(PepFP).toString()
			       );
    }

    /* LIST FUNCTIONS
     * Each function below is run on the empty header of a linked list
     * of DTA files.  The operation specified is run on each DTAFile
     * in the list.
     */

    // Header call for Quicksort function.  Sorts DTAFiles on basis of Sequence
    public void	SortList() {
	//this is just an empty header; the real data starts at the next item
	if (this.Next != null)
	    this.Next = this.Next.Sort(null);
    }

    // Header call for Quicksort function.  Sorts DTAFiles on basis of sequence alignment
    public void	DisplaySortList() {
	//this is just an empty header; the real data starts at the next item
	if (this.Next != null)
	    this.Next = this.Next.DisplaySort(null);
    }

    // Header call for Quicksort function.  Sorts DTAFiles on basis of filename
    public void	FileSortList() {
	//this is just an empty header; the real data starts at the next item
	if (this.Next != null)
	    this.Next = this.Next.FileSort(null);
    }

    /* Searches for duplicate sequences under this locus.  Removes
     * all but one of duplicate set on basis of XCorr
     */
    public void KeepOnlyDTAWithHighXCorr() {
	//this is just an empty header; the real data starts at the next item.
	//This function removes DTAs that are duplicate sequences on basis of XCorr
	DTAFile         Runner = this;
	while ((Runner.Next != null) && (Runner.Next.Next != null)) {
	    if (Runner.Next.Sequence.equals(Runner.Next.Next.Sequence) &&
		(Runner.Next.ChargeState == Runner.Next.Next.ChargeState)) {
		if (Runner.Next.XCorr > Runner.Next.Next.XCorr) {
		    Runner.Next.Next = Runner.Next.Next.Next;
		}
		else {
		    Runner.Next = Runner.Next.Next;
		}
	    }
	    else {
		Runner = Runner.Next;
	    }
	}
    }

    /* Searches for duplicate sequences under this locus.  Removes
     * all but one of duplicate set on basis of Total Intensity
     */
    public void KeepOnlyDTAInEachSaltStep() {
	//this is just an empty header; the real data starts at the next item
	DTAFile         Runner = this;
	while ((Runner.Next != null) && (Runner.Next.Next != null)) {

	    if (Runner.Next.Sequence.equals(Runner.Next.Next.Sequence) &&
		(Runner.Next.ChargeState == Runner.Next.Next.ChargeState) &&
		(Runner.Next.RootFileName.equals(Runner.Next.Next.RootFileName))) {
		if (Runner.Next.XCorr > Runner.Next.Next.XCorr) {
		    Runner.Next.Next = Runner.Next.Next.Next;
		}
		else {
		    Runner.Next = Runner.Next.Next;
		}
	    }
	    else {
		Runner = Runner.Next;
	    }
	}
    }

    /* Dave's amazingly clever quicksorter.  I leave it as an exercise
     * for the reader to determine how this works.
     */
    private DTAFile Sort(DTAFile Follower) {
	//Recursive quicksorter
	//Returns first item of sorted list (from those starting at this element)
	//Accepts first item to follow
	DTAFile		ListAbove = null;
	DTAFile		ListBelow = null;
	DTAFile		PlaceHolder;
	DTAFile		PlaceHolder2;
	PlaceHolder = this.Next;
	//Partition all remaining points of this linked list
	while (PlaceHolder != null) {
	    PlaceHolder2 = PlaceHolder.Next;
	    if (this.Sequence.compareTo(PlaceHolder.Sequence) > 0) {
				//Move this item to list above this
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else if (this.Sequence.compareTo(PlaceHolder.Sequence) < 0) {
				//Move this item to list below this point
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    else {
		if (this.ChargeState > PlaceHolder.ChargeState) {
				//Move this item to list above this
		    PlaceHolder.Next = ListAbove;
		    ListAbove = PlaceHolder;
		}
		else if (this.ChargeState < PlaceHolder.ChargeState) {
				//Move this item to list below this point
		    PlaceHolder.Next = ListBelow;
		    ListBelow = PlaceHolder;
		}
		else {
		    if (this.FileName.compareTo(PlaceHolder.FileName) > 0) {
				//Move this item to list above this
			PlaceHolder.Next = ListAbove;
			ListAbove = PlaceHolder;
		    }
		    else {
				//Move this item to list below this point
			PlaceHolder.Next = ListBelow;
			ListBelow = PlaceHolder;
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

    /* Sort peptides by alignment position, then by precursor mass
     * Same code as before
     */
    private DTAFile DisplaySort(DTAFile Follower) {
	//Recursive quicksorter
	//Returns first item of sorted list (from those starting at this element)
	//Accepts first item to follow
	DTAFile		ListAbove = null;
	DTAFile		ListBelow = null;
	DTAFile		PlaceHolder;
	DTAFile		PlaceHolder2;
	PlaceHolder = this.Next;
	//Partition all remaining points of this linked list
	while (PlaceHolder != null) {
	    PlaceHolder2 = PlaceHolder.Next;
	    if (this.SequencePosition > PlaceHolder.SequencePosition) {
				//Move this item to list above this
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else if (this.SequencePosition < PlaceHolder.SequencePosition) {
				//Move this item to list below this point
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    else {
		if (this.CalcPreMass > PlaceHolder.CalcPreMass) {
				//Move this item to list above this
		    PlaceHolder.Next = ListAbove;
		    ListAbove = PlaceHolder;
		}
		else if (this.CalcPreMass < PlaceHolder.CalcPreMass) {
				//Move this item to list below this point
		    PlaceHolder.Next = ListBelow;
		    ListBelow = PlaceHolder;
		}
		else {
                    if (this.ChargeState > PlaceHolder.ChargeState) {
                                //Move this item to list above this
                        PlaceHolder.Next = ListAbove;
                        ListAbove = PlaceHolder;
                    }
                    else if (this.ChargeState < PlaceHolder.ChargeState) {
                                //Move this item to list below this point
                        PlaceHolder.Next = ListBelow;
                        ListBelow = PlaceHolder;
                    }
                    else {
                        if (this.FileName.compareTo(PlaceHolder.FileName) > 0) {
                                //Move this item to list above this
                            PlaceHolder.Next = ListAbove;
                            ListAbove = PlaceHolder;
                        }
                        else {
                                //Move this item to list below this point
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
	    this.Next = ListBelow.DisplaySort(Follower);
	if (ListAbove == null)
	    return this;
	else
	    return ListAbove.DisplaySort(this);
    }

    // Sort identifications by 1) FileName, 2) ChargeState
    private DTAFile FileSort(DTAFile Follower) {
	//Recursive quicksorter
	//Returns first item of sorted list (from those starting at this element)
	//Accepts first item to follow
	DTAFile		ListAbove = null;
	DTAFile		ListBelow = null;
	DTAFile		PlaceHolder;
	DTAFile		PlaceHolder2;
	PlaceHolder = this.Next;
	//Partition all remaining points of this linked list
	while (PlaceHolder != null) {
	    PlaceHolder2 = PlaceHolder.Next;
	    if (this.FileName.compareTo(PlaceHolder.FileName) > 0) {
				//Move this item to list above this
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else if (this.FileName.compareTo(PlaceHolder.FileName) < 0) {
				//Move this item to list below this point
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    else {
		if (this.ChargeState > PlaceHolder.ChargeState) {
				//Move this item to list above this
		    PlaceHolder.Next = ListAbove;
		    ListAbove = PlaceHolder;
		}
		else {
				//Move this item to list below this point
		    PlaceHolder.Next = ListBelow;
		    ListBelow = PlaceHolder;
		}
	    }
	    //Move to next item to be partitioned
	    PlaceHolder = PlaceHolder2;
	}
	if (ListBelow == null)
	    this.Next = Follower;
	else
	    this.Next = ListBelow.FileSort(Follower);
	if (ListAbove == null)
	    return this;
	else
	    return ListAbove.FileSort(this);
    }

    public void GetMassOffsets(SelectCriteria Cutoffs) {
        int counter;
        double diffC12C13 = 1.003354826;
        double AbsoluteOffset = PrecursorMass - CalcPreMass;

        Raw_PPM_Offset = AbsoluteOffset;
        if (Cutoffs.Isotopes) {
            for (counter = 1; counter < 2 * ChargeState + 1; counter++) {
                if (Math.abs(Raw_PPM_Offset) > (Math.abs(AbsoluteOffset - counter * diffC12C13)))
                    Raw_PPM_Offset = AbsoluteOffset - counter * diffC12C13;
            }
        }
        Raw_PPM_Offset = Raw_PPM_Offset * 1E6 / PrecursorMass;
        Adjusted_PPM_Offset = Raw_PPM_Offset;
    }

	public float getPepFP() {
		return PepFP;
	}
}
