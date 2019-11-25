import java.io.*;
import java.util.*;

//ParamsFile Object
//For DTASelect
//Created by Dave Tabb
//Begun 8/31/2000

public class ParamsFile {
    String           DBName;
    DiffMod          DiffMods = new DiffMod();
    StaticMod        StaticMods = new StaticMod();
    float            CPepMod = 0;
    float            NPepMod = 0;
    float            CProtMod = 0;
    float            NProtMod = 0;
    float            AvgMasses[];
    float            MonoMasses[];
    boolean          AvgTypeForFragmentIons = false;
    boolean          AvgTypeForParentIon = true;

    public ParamsFile() {
	this.setInitialMasses();
    }

    public ParamsFile(StringTokenizer Parser, String Version) {
	// Is compatible only with 1.8 and later files!
	// If this file is from the current version
	DBName = Parser.nextToken();
	CPepMod = new Float(Parser.nextToken()).floatValue();
	NPepMod = new Float(Parser.nextToken()).floatValue();
	CProtMod = new Float(Parser.nextToken()).floatValue();
	NProtMod = new Float(Parser.nextToken()).floatValue();
	AvgTypeForFragmentIons = new Boolean(Parser.nextToken()).booleanValue();
	AvgTypeForParentIon = new Boolean(Parser.nextToken()).booleanValue();
	this.setInitialMasses();
    }

    public void ApplyStaticMods() {
	// Loop through all the static modifications in this
	// ParamsFile and apply them to the masses of the amino acids.
	StaticMod          SMRunner = this.StaticMods.Next;
	int                Index;
	while (SMRunner != null) {
	    Index = SMRunner.Residue - 'A';
	    AvgMasses[Index] += SMRunner.Mass;
	    MonoMasses[Index] += SMRunner.Mass;
	    SMRunner = SMRunner.Next;
	}
    }

    private void setInitialMasses() {
	AvgMasses = new float[26];
	MonoMasses = new float[26];
	// Set up an array of average residue masses
	AvgMasses['G' - 65] = 57.0519f;
	AvgMasses['A' - 65] = 71.0788f;
	AvgMasses['S' - 65] = 87.0782f;
	AvgMasses['P' - 65] = 97.1167f;
	AvgMasses['V' - 65] = 99.1326f;
	AvgMasses['T' - 65] = 101.1051f;
	AvgMasses['C' - 65] = 103.1388f;
	AvgMasses['L' - 65] = 113.1594f;
	AvgMasses['I' - 65] = 113.1594f;
	AvgMasses['X' - 65] = 113.1594f;
	AvgMasses['N' - 65] = 114.1038f;
	AvgMasses['O' - 65] = 114.1472f;
	AvgMasses['B' - 65] = 114.5962f;
	AvgMasses['D' - 65] = 115.0886f;
	AvgMasses['Q' - 65] = 128.1307f;
	AvgMasses['K' - 65] = 128.1741f;
	AvgMasses['Z' - 65] = 128.6231f;
	AvgMasses['E' - 65] = 129.1155f;
	AvgMasses['M' - 65] = 131.1926f;
	AvgMasses['H' - 65] = 137.1411f;
	AvgMasses['F' - 65] = 147.1766f;
	AvgMasses['R' - 65] = 156.1875f;
	AvgMasses['Y' - 65] = 163.1760f;
	AvgMasses['W' - 65] = 186.2132f;
	// Set up an array of monoisotopic residues masses
	MonoMasses['G' - 65] = 57.02146f;
	MonoMasses['A' - 65] = 71.03711f;
	MonoMasses['S' - 65] = 87.02303f;
	MonoMasses['P' - 65] = 97.05276f;
	MonoMasses['V' - 65] = 99.06841f;
	MonoMasses['T' - 65] = 101.04768f;
	MonoMasses['C' - 65] = 103.00919f;
	MonoMasses['L' - 65] = 113.08406f;
	MonoMasses['I' - 65] = 113.08406f;
	MonoMasses['X' - 65] = 113.08406f;
	MonoMasses['N' - 65] = 114.04293f;
	MonoMasses['O' - 65] = 114.07931f;
	MonoMasses['B' - 65] = 114.53494f;
	MonoMasses['D' - 65] = 115.02694f;
	MonoMasses['Q' - 65] = 128.05858f;
	MonoMasses['K' - 65] = 128.09496f;
	MonoMasses['Z' - 65] = 128.55059f;
	MonoMasses['E' - 65] = 129.04259f;
	MonoMasses['M' - 65] = 131.04049f;
	MonoMasses['H' - 65] = 137.05891f;
	MonoMasses['F' - 65] = 147.06841f;
	MonoMasses['R' - 65] = 156.10111f;
	MonoMasses['Y' - 65] = 163.06333f;
	MonoMasses['W' - 65] = 186.07931f;
    }

    // Add a new static mod to this set of Sequest.params
    public void addStaticMod(char PassedResidue, float PassedMass) {
	StaticMod         SMRunner = this.StaticMods;
	while (SMRunner.Next != null) {
	    SMRunner = SMRunner.Next;
	}
	SMRunner.Next = new StaticMod();
	SMRunner = SMRunner.Next;
	SMRunner.Mass = PassedMass;
	SMRunner.Residue = PassedResidue;
    }

    // Add a new differential mod to this set of Sequest.params
    public void addDiffMod(float PassedMass, char PassedSymbol, String PassedResidues) {
	DiffMod            DMRunner = this.DiffMods;
	while (DMRunner.Next != null)
	    DMRunner = DMRunner.Next;
	DMRunner.Next = new DiffMod();
	DMRunner = DMRunner.Next;
	DMRunner.Mass = PassedMass;
	DMRunner.Symbol = PassedSymbol;
	DMRunner.Residues = PassedResidues;
    }

    // Add a new lossy differential mod to this set of Sequest.params
    public void addDiffMod(float PassedMass, char PassedSymbol, String PassedResidues,
			   float PassedNeuLoss, float PassedFragLoss) {
	DiffMod            DMRunner = this.DiffMods;
	while (DMRunner.Next != null)
	    DMRunner = DMRunner.Next;
	DMRunner.Next = new DiffMod();
	DMRunner = DMRunner.Next;
	DMRunner.Mass = PassedMass;
	DMRunner.Symbol = PassedSymbol;
	DMRunner.Residues = PassedResidues;
	DMRunner.PrecursorNeutralLoss = PassedNeuLoss;
	DMRunner.FragmentNeutralLoss = PassedFragLoss;
    }



    // Returns a string representing all fields of this ParamsFile.
    // Used for writing DTASelect.txt file.
    public String getDTASelectTxtString() {
	int               Index;
	StringBuffer      Returned = new StringBuffer();
	Returned.append("S\t");
	Returned.append(DBName);
	Returned.append("\t");
	Returned.append(CPepMod);
	Returned.append("\t");
	Returned.append(NPepMod);
	Returned.append("\t");
	Returned.append(CProtMod);
	Returned.append("\t");
	Returned.append(NProtMod);
	Returned.append("\t");
	Returned.append(AvgTypeForFragmentIons);
	Returned.append("\t");
	Returned.append(AvgTypeForParentIon);
	Returned.append("\n");
	Returned.append(StaticMods.getDTASelectTxtString());
	Returned.append(DiffMods.getDTASelectTxtString());
	return Returned.toString();
    }

    public String GetHTMLMods() {
	StringBuffer      Returned = new StringBuffer("sequest.params modifications:\n<table border>\n");
	int               Index;
	DiffMod           DMRunner = DiffMods.Next;
	StaticMod         SMRunner = StaticMods.Next;
	while (DMRunner != null) {
	    Returned.append("<tr><td>");
	    Returned.append(DMRunner.Symbol);
	    Returned.append("</td><td>");
	    Returned.append(DMRunner.Residues);
	    Returned.append("</td><td>");
	    Returned.append(DMRunner.Mass);
	    if ( (DMRunner.PrecursorNeutralLoss != 0f) ||
		 (DMRunner.FragmentNeutralLoss != 0f) ) {
		Returned.append("</td><td>");
		Returned.append(DMRunner.PrecursorNeutralLoss);
		Returned.append("</td><td>");
		Returned.append(DMRunner.FragmentNeutralLoss);
	    }
	    Returned.append("</td></tr>\n");
	    DMRunner = DMRunner.Next;
	}
	while (SMRunner != null) {
	    Returned.append("<tr><td>Static</td><td>");
	    Returned.append(SMRunner.Residue);
	    Returned.append("</td><td>");
	    Returned.append(SMRunner.Mass);
	    Returned.append("</td></tr>\n");
	    SMRunner = SMRunner.Next;
	}
	Returned.append("</table>");
	return Returned.toString();
    }

    public String CGIString() {
	StringBuffer     Returned = new StringBuffer();
	StaticMod        SMRunner = StaticMods.Next;
	while (SMRunner != null) {
	    Returned.append("&Mass");
	    Returned.append(SMRunner.Residue);
	    Returned.append("=");
	    if (this.AvgTypeForFragmentIons)
		Returned.append(AvgMasses[SMRunner.Residue - 'A']);
	    else
		Returned.append(MonoMasses[SMRunner.Residue - 'A']);
	    SMRunner = SMRunner.Next;
	}
	return Returned.toString();
    }

    /* Read the configuration info from this Mascot .dat file. */
    public static ParamsFile ReadMascotConfig(String FileName,
					      IniFile Config) {
	try {
	    File	    DatFile = new File(FileName);
	    FileReader	    Filereader = new FileReader(DatFile);
	    BufferedReader  Incoming = new BufferedReader(Filereader);
	    String          WholeLine;
	    String          CurrentToken;
	    String          DBLabel = "";
	    File            DBPath = null;
	    StringTokenizer Parser;
	    ParamsFile      NewOne = new ParamsFile();
	    int             DiffModCount = 0;
	    char            Symbol;
	    String          ModSite;
	    float           ModShift;
	    // Queue to "parameters" tag
	    WholeLine = Incoming.readLine();
	    while ( (WholeLine != null) &&
		    (WholeLine.indexOf("name=\"parameters") < 0) ) {
		WholeLine = Incoming.readLine();
	    }
	    WholeLine = Incoming.readLine();
	    WholeLine = Incoming.readLine();
	    while (!WholeLine.startsWith("--")) {
		Parser = new StringTokenizer(WholeLine,"=");
		CurrentToken = Parser.nextToken();
		if (CurrentToken.equals("MASS")) {
		    if (Parser.nextToken().equals("Monoisotopic")) {
			NewOne.AvgTypeForFragmentIons = false;
		    }
		    else {
			NewOne.AvgTypeForFragmentIons = true;
		    }
		}
		// Grab database label
		else if (CurrentToken.equals("DB")) {
		    DBLabel = Parser.nextToken();
		    // Get path info from mascot.dat
		    try {
			// Find the path to this database in mascot.dat
			File            MascotDir = new File(Config.MascotPath);
			File            CfgDir = new File(MascotDir, "config");
			File	        CfgFile = new File(CfgDir, "mascot.dat");
			FileReader	CfgReader = new FileReader(CfgFile);
			BufferedReader  CfgIncoming = new BufferedReader(CfgReader);
			WholeLine = CfgIncoming.readLine();
			while ( (WholeLine != null) &&
				(!WholeLine.startsWith("Databases")) ) {
			    WholeLine = CfgIncoming.readLine();
			}
			if (WholeLine == null) {
			    System.out.println("No Databases section in mascot.dat.");
			}
			else {
			    CurrentToken = "";
			    while ( (WholeLine != null) &&
				    (!WholeLine.startsWith("end")) &&
				    (!CurrentToken.equals(DBLabel)) ) {
				WholeLine = CfgIncoming.readLine();
				if (WholeLine != null) {
				    Parser = new StringTokenizer(WholeLine,"\t");
				    CurrentToken = Parser.nextToken();
				}
			    }
			    if ( (WholeLine != null) &&
				 (CurrentToken.equals(DBLabel)) ) {
				//Grab the current path
				DBPath = new File(Parser.nextToken());
			    }
			    else {
				System.out.println("Couldn't find matching database path in mascot.dat");
				System.out.println(DBLabel);
				System.out.println(WholeLine);
				DBPath = null;
			    }
			}
		    }
		    catch (IOException failure) {
			System.out.println("Failed to read Mascot configuration file.  Check mascot-path in DTASelect.ini.");
			System.out.println(failure);
			System.exit(0);
		    }
		}
		WholeLine = Incoming.readLine();
	    }
	    // Queue to "header" tag
	    Filereader = new FileReader(DatFile);
	    Incoming = new BufferedReader(Filereader);
	    WholeLine = Incoming.readLine();
	    while ( (WholeLine != null) &&
		    (WholeLine.indexOf("name=\"header") < 0) ) {
		WholeLine = Incoming.readLine();
	    }
	    WholeLine = Incoming.readLine();
	    WholeLine = Incoming.readLine();
	    while (!WholeLine.startsWith("--")) {
		Parser = new StringTokenizer(WholeLine,"=");
		CurrentToken = Parser.nextToken();
		// Grab the filename of the database
		if (CurrentToken.equals("release")) {
		    if (DBPath == null) {
			NewOne.DBName = Parser.nextToken();
		    }
		    else {
			NewOne.DBName = new File(DBPath.getParent(),
						 Parser.nextToken()).toString();
		    }
		}
		WholeLine = Incoming.readLine();
	    }
	    // Queue to "masses" tag
	    Filereader = new FileReader(DatFile);
	    Incoming = new BufferedReader(Filereader);
	    WholeLine = Incoming.readLine();
	    while ( (WholeLine != null) &&
		    (WholeLine.indexOf("name=\"masses") < 0) ) {
		WholeLine = Incoming.readLine();
	    }
	    WholeLine = Incoming.readLine();
	    WholeLine = Incoming.readLine();
	    while (!WholeLine.startsWith("--")) {
		Parser = new StringTokenizer(WholeLine,"=, \t()");
		CurrentToken = Parser.nextToken();
		if (CurrentToken.equals("C_term")) {
		    NewOne.CPepMod = new Float(Parser.nextToken()).floatValue();
		}
		else if (CurrentToken.equals("N_term")) {
		    NewOne.NPepMod = new Float(Parser.nextToken()).floatValue();
		}
		else if (CurrentToken.startsWith("delta")) {
		    DiffModCount++;
		    switch (DiffModCount) {
		    case 1:
			Symbol = '*';
			break;
		    case 2:
			Symbol = '#';
			break;
		    case 3:
			Symbol = '@';
			break;
		    case 4:
			Symbol = '$';
			break;
		    default:
			Symbol = '?';
			break;
		    }
		    ModShift = new Float(Parser.nextToken()).floatValue();
		    ModSite = "";
		    while (Parser.hasMoreTokens()) {
			ModSite = Parser.nextToken();
		    }
		    NewOne.addDiffMod(ModShift,
				      Symbol, 
				      ModSite);
		}
		WholeLine = Incoming.readLine();
	    }
	    return NewOne;
	}
	catch (IOException failure) {
	    System.out.println("An error occured while reading configuration from " + FileName);
	    System.out.println(failure);
	    System.exit(0);
	}
	return null;
    }

    /* Searches current directory for sequest.params file.  Function
     * throws a I/O Exception if that file is not present.  Extracts
     * information about FASTA database file location, differential
     * modification residues, and differential modification mass shifts.
     */
    public static ParamsFile ReadFile(File CurrentDirectory) throws IOException {
	ParamsFile       NewPF = new ParamsFile();
	File             ParamsFile = new File(CurrentDirectory, "sequest.params");
	//Read a sequest.params file, returning the database name
	FileReader       InputFileReader = new FileReader(ParamsFile);
	BufferedReader   Incoming = new BufferedReader(InputFileReader);
	String           LineBuffer = Incoming.readLine();
	String           WholeLine;
	String           TempString;
	StringTokenizer  Parser;
	int              IndexOfResidue;
	float            MassBuffer;
	char             SymbolBuffer;
	//Skip up to the [SEQUEST] line
	while ( (LineBuffer != null) && (!LineBuffer.equals("[SEQUEST]")) ) {
	    LineBuffer = Incoming.readLine();
	}
	if (LineBuffer == null) {
	    System.out.println("sequest.params [SEQUEST] header is missing.");
	    System.exit(0);
	}
	LineBuffer = Incoming.readLine();
	while (LineBuffer != null) {
	    WholeLine = LineBuffer;
	    Parser = new StringTokenizer(LineBuffer, " \t,");
	    if (Parser.hasMoreTokens()) {
		LineBuffer = Parser.nextToken();
		try {
		    if (LineBuffer.equals("database_name") ||
			LineBuffer.equals("first_database_name") ||
			LineBuffer.equals("second_database_name")) {
			//skip equals sign
			LineBuffer = Parser.nextToken();
			if (Parser.hasMoreTokens()) {
			    if (NewPF.DBName != null) {
				System.out.println("Using last listed database rather than " + NewPF.DBName);
			    }
			    NewPF.DBName = Parser.nextToken();
			    if (NewPF.DBName.endsWith(".bin") ||
				NewPF.DBName.endsWith(".hdr")) {
				System.out.println("Cannot use indexed database: " +
						   NewPF.DBName);
				System.out.println("Change sequest.params to indicate FASTA database and re-run.");
				System.exit(0);
			    }
			}
		    }
		    else if (LineBuffer.equals("diff_search_options")) {
			//skip equals sign
			LineBuffer = Parser.nextToken();
			MassBuffer = new Float(Parser.nextToken()).floatValue();
			NewPF.addDiffMod(MassBuffer, '*', Parser.nextToken());
			if (Parser.hasMoreTokens()) {
			    MassBuffer = new Float(Parser.nextToken()).floatValue();
			    NewPF.addDiffMod(MassBuffer, '#', Parser.nextToken());
			}
			if (Parser.hasMoreTokens()) {
			    TempString = Parser.nextToken();
			    if (!TempString.equals("X")) {
				MassBuffer = new Float(TempString).floatValue();
				NewPF.addDiffMod(MassBuffer, '@', Parser.nextToken());
			    }
			}
		    }
		    else if (LineBuffer.equals("dm_standard")) {
			SymbolBuffer = Parser.nextToken().charAt(0);
			NewPF.addDiffMod(new Float(Parser.nextToken()).floatValue(),
					 SymbolBuffer, Parser.nextToken());
		    }
		    else if (LineBuffer.equals("dm_nterm")) {
			SymbolBuffer = Parser.nextToken().charAt(0);
			NewPF.addDiffMod(new Float(Parser.nextToken()).floatValue(),
					 SymbolBuffer, "N-term:" + Parser.nextToken());
		    }
		    else if (LineBuffer.equals("dm_nloss")) {
			SymbolBuffer = Parser.nextToken().charAt(0);
			NewPF.addDiffMod(new Float(Parser.nextToken()).floatValue(),
					 SymbolBuffer, Parser.nextToken(),
					 new Float(Parser.nextToken()).floatValue(),
					 new Float(Parser.nextToken()).floatValue());
		    }
		    else if (LineBuffer.equals("add_C_terminus") ||
			     LineBuffer.equals("add_Cterm_peptide")) {
			Parser.nextToken();
			NewPF.CPepMod = new Float(Parser.nextToken()).floatValue();
		    }
		    else if (LineBuffer.equals("add_Cterm_protein")) {
			Parser.nextToken();
			NewPF.CProtMod = new Float(Parser.nextToken()).floatValue();
		    }
		    else if (LineBuffer.equals("add_N_terminus") ||
			     LineBuffer.equals("add_Nterm_peptide") ) {
			Parser.nextToken();
			NewPF.NPepMod = new Float(Parser.nextToken()).floatValue();
		    }
		    else if (LineBuffer.equals("add_Nterm_protein")) {
			Parser.nextToken();
			NewPF.NProtMod = new Float(Parser.nextToken()).floatValue();
		    }
		    else if (LineBuffer.equals("mass_type_parent")) {
			Parser.nextToken();
			NewPF.AvgTypeForParentIon = (Parser.nextToken().charAt(0)
							== '0') ? true : false;
		    }
		    else if (LineBuffer.equals("mass_type_fragment")) {
			Parser.nextToken();
			NewPF.AvgTypeForFragmentIons = (Parser.nextToken().charAt(0) ==
							'0') ? true : false;
		    }
		    else if (LineBuffer.startsWith("add_")) {
			Parser.nextToken();
			MassBuffer = new Float(Parser.nextToken()).floatValue();
			if (MassBuffer != 0.0f) {
			    NewPF.addStaticMod(LineBuffer.charAt(4), MassBuffer);
			}
		    }
		}
		catch (NumberFormatException failure) {
		    System.out.println("sequest.params has garbled number in this line:");
		    System.out.println(WholeLine);
		    System.exit(0);
		}
		catch (NoSuchElementException failure) {
		    System.out.println("sequest.params is missing information in this line:");
		    System.out.println(WholeLine);
		    System.exit(0);
		}
	    }
	    LineBuffer = Incoming.readLine();
	}
	NewPF.ApplyStaticMods();
	return NewPF;
    }

    public void DebugAcids() {
	int      Looper;
	for (Looper = 0; Looper < 26; Looper++) {
	    System.out.println((char)('A' + Looper) + "\t" +
			       AvgMasses[Looper] + "\t" +
			       MonoMasses[Looper]);
	}
    }

    public void DebugPrint() {
	System.out.println(DBName);
	DiffMods.DebugPrint();
	System.out.println(CPepMod);
	System.out.println(NPepMod);
	System.out.println(CProtMod);
	System.out.println(NProtMod);
	System.out.println(AvgTypeForFragmentIons);
	System.out.println(AvgTypeForParentIon);
    }

    class StaticMod {
	public float          Mass;
	public char           Residue;
	public StaticMod      Next;

	public String getDTASelectTxtString() {
	    StringBuffer      Returned = new StringBuffer();
	    StaticMod         SMRunner = this.Next;
	    while (SMRunner != null) {
		Returned.append("SM\t");
		Returned.append(SMRunner.Mass);
		Returned.append("\t");
		Returned.append(SMRunner.Residue);
		Returned.append("\n");
		SMRunner = SMRunner.Next;
	    }
	    return Returned.toString();
	}
    }

    class DiffMod {
	public float          Mass;
	public char           Symbol;
	public String         Residues;
	public float          PrecursorNeutralLoss = 0f;
	public float          FragmentNeutralLoss = 0f;
	public DiffMod        Next;

	/* Add the passed list of DiffMods to the end of this list.
           If symbols are duplicated, note this to the user.  The
           passed DiffMods are a list with a null header.
	   public void Append(DiffMod NewOnes) {
	   DiffMod    DMRunner = this;
	   DiffMod    DMRunner2;
	   while (DMRunner.Next != null)
	   DMRunner = DMRunner.Next;
	   DMRunner.Next = NewOnes.Next;
	   DMRunner = this.Next;
	   while (DMRunner.Next != null) {
	   DMRunner2 = DMRunner.Next;
	   while (DMRunner2 != null) {
	   if (DMRunner2.Symbol == DMRunner.Symbol) {
	   System.out.println("\tDifferential Modification " +
	   new Character(DMRunner.Symbol).toString() +
	   " appears multiple times.");
	   DMRunner2 = null;
	   }
	   DMRunner2 = DMRunner2.Next;
	   }
	   DMRunner = DMRunner.Next;
	   }
	   }
	*/

	public String getDTASelectTxtString() {
	    StringBuffer      Returned = new StringBuffer();
	    DiffMod           DMRunner = this.Next;
	    while (DMRunner != null) {
		Returned.append("DM\t");
		Returned.append(DMRunner.Mass);
		Returned.append("\t");
		Returned.append(DMRunner.Symbol);
		Returned.append("\t");
		Returned.append(DMRunner.Residues);
		Returned.append("\t");
		Returned.append(DMRunner.PrecursorNeutralLoss);
		Returned.append("\t");
		Returned.append(DMRunner.FragmentNeutralLoss);
		Returned.append("\n");
		DMRunner = DMRunner.Next;
	    }
	    return Returned.toString();
	}

	// Debugging function
	public void DebugPrint() {
	    DiffMod         DMRunner = this.Next;
	    while (DMRunner != null) {
		System.out.println(new Float(DMRunner.Mass).toString() + "\t" +
				   new Character(DMRunner.Symbol).toString() + "\t" +
				   Residues);
		DMRunner = DMRunner.Next;
	    }
	}

	// Return the mass shift for the passed symbol
	public float getMassShiftFor(char PassedSymbol) {
	    DiffMod         DMRunner = this.Next;
	    while (DMRunner != null) {
		if (DMRunner.Symbol == PassedSymbol) {
		    return DMRunner.Mass;
		}
		DMRunner = DMRunner.Next;
	    }
	    return 0.0f;
	}
    }
}
