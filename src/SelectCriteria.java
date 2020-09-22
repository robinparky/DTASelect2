import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;

//SelectCriteria Object
//Begun Aug 11, 2000
//By David Tabb
//modified April 4, 2017, Titus Jung titusjung@gmail.com
public class SelectCriteria {
	boolean UseCriteria = true;
	// Per DTA criteria
	int MinChargeState = 1;
	int MaxChargeState = 50;
	int MaximumSp = 10000;
	float MinSpScore = -1f;
	boolean UseXCorrCriteria = true;
	boolean DisplayXCorrCriteria = false;
	boolean filterReadXorr = false;
	double minReadXorr = 1.0;
	float MinSingleChargeXCorr = 1.0f;
	float MinDoubleChargeXCorr = 1.0f;
	float MinTripleChargeXCorr = 1.0f;
	float MinQuadrupleChargeXCorr = 1.0f;
	float MinDeltCN = 0.00f;
	boolean DisplayDeltaMass = false;
	boolean UseMinDeltaMass = false;
	boolean UseMaxDeltaMass = false;
	float MinDeltaMassPPM = 0.0f;
	float MaxDeltaMassPPM = 1000000.0f;
	float MinPepConf = 0.0f;
	float MaxPepFP = 0.05f;
	float MinIonProportion = -1.0f;
	byte HandleModified = 1;
	byte RequireTryptic = 0;
	byte MinimumMissedTryptic = 0;
	byte MaximumMissedTryptic = 100;
	byte MinimumModifications = 0;
	byte MaximumModifications = 100;
	boolean PermitAmbiguous = false;
	byte PeptideValidation = 0;
	// Statistics options
	String DecoyLabel = "Reverse_";
	String TrueProteinLabel = "";
	boolean HideDecoy = false;
	boolean UseStatistics = true;
	boolean UseTrueProtein = false;
	boolean UseEmpirics = true;
	boolean UseRawXCorr = false;
	boolean UseChargeState = true;
	boolean UseModStats = false;
	boolean UseTrypticInfo = false;
	boolean UseXCorr = true;
	boolean UseDeltCN = true;
	boolean UseSp = false;
	boolean UseMassDifference = false;
	boolean UseLogMass = false;
	boolean UseLogSpRank = false;
	boolean UseIonStat = false;
	boolean UseProteinGrouper = false;
	int FitLevel = 0;
	byte MaxStatisticsCharge = 3;
	int PlotLevel = 0;
	double FPTol = 0.2;
	// Sequence criteria
	boolean UseFilenameFilter = false;
	String FilenameFilter = "";
	boolean UseSequenceCriteria = false;
	String SeqMustIncludeResidues = "";
	String SeqMustIncludePattern = "";
	String SeqMustExcludeBeforeC = "";
	String SeqMustEndWithResidues = "";
	String SeqMustComeAfterResidues = "";
	String SeqMustIncludeTwo = "";
	String SeqMustIncludeOne = "";
	String SeqMustIncludeNone = "";
	String SeqMustIncludeEither = "";
	int SeqMinLength = 6;
	int SeqMaxLength = 10000;
	byte SequenceCompleteness = 0;
	byte NGlycosylation = 2;
	// M/Z criteria
	boolean UseMinMZ = false;
	boolean UseMaxMZ = false;
	float MinMZ = 0f;
	float MaxMZ = 100000f;
	// Maximum XCorr criteria
	boolean UseMaxXCorrs = false;
	float MaxSingleChargeXCorr = 1000.0f;
	float MaxDoubleChargeXCorr = 1000.0f;
	float MaxTripleChargeXCorr = 1000.0f;
	float MaxQuadrupleChargeXCorr = 1000.0f;
	// Per Locus criteria
	float MinProtConf = 0.0f;
	float MaxProtFP = 1.00f;
	boolean UseProteinFilters = false;
	byte MinProtTryptic = 0;
	float MinProtXCorr = 0.00f;
	float MinProtSpScore = 0.00f;
	float MaxProtMinFP = 1.00f;
	float MinProtPepConf = 0.00f;
	float MaxProtDM = 100000.0f;
	int MinPepsPerLocus = 2;
	int MinSpectraPerLocus = 1;
	int MinTrypticPeps = 0;
	int MinModPepsPerLocus = 0;
	boolean IncludeOnlyUniques = false;
	boolean UseMultipleScore = false;
	int PurgeDuplicateSequences = 2;
	boolean RemoveSubsets = true;
	String NameExcludePattern = "";
	String NameIncludePattern = "";
	String DescripExcludePattern = "";
	String DescripIncludePattern = "";
	String ProteinMustStartWithResidues = "";
	String ProteinMustEndWithResidues = "";
	byte LocusValidation = 0;
	boolean UseMinMW = false;
	boolean UseMaxMW = false;
	boolean UseMinPI = false;
	boolean UseMaxPI = false;
	boolean printHTML = false;
	float MinMW = 0f;
	float MaxMW = 10000000.0f;
	float MinSequenceCoverage = 0f;
	float MinPI = 0f;
	float MaxPI = 14f;
	// Options
	String EnzymeCutsAfter = "KR";
	String EnzymeCutsBefore = "";
	int SetDeltCN = 1;
	boolean DisplayPI = false;
	boolean DisplayKD = false;
	boolean DisplayBB = false;
	boolean DisplayHPLC = false;
	boolean DisplaySeqPos = false;
	boolean ExtraColumns = false;
	boolean Isotopes = true;
	boolean FullProteinName = false;
	// Utility settings
	boolean UseMascotOutput = false;
	boolean PrintCopyList = false;
	boolean PrintChroma = false;
	boolean PrintMods = false;
	boolean PrintAlignment = false;
	boolean PrintSimilarities = false;
	boolean PrintXML = false;
	boolean PrintDB = false;
	boolean Quiet = false;
	boolean Brief = false;
	boolean ShortForm = false;
	boolean UseCustomPath = false;
	int colCount = 11;
	//boolean multipleScores = false;
	String CustomPath = "";
	/*
	 * DisplayType meaning: 0: DTAs and OUTs 1: SQT 2: Mascot (not implemented)
	 * Is set on basis of ID file type
	 */
	byte DisplayType = 0;
	boolean UseClassifications = false;
	boolean UseAuxInfo = false;
	boolean LookInCurrentDirectory = false;
	boolean DisplayBirdsEye = false;
	boolean UseGUIInstead = false;
	boolean CompressDTAs = false;

	boolean PrintOnlyUnique = false;
	boolean noDB = false;
	boolean readFasta = false;

	boolean addPSM  = false;
	boolean dms  = false;
	boolean showCorrectedDmValue = false;

	boolean filterByModification = false;
	String modFilterStr = "";
	/*
	 * If the user has commonly used options for DTASelect, he or she can write
	 * them into a file called DTASelect.params and copy them around rather than
	 * continually retype the same command-line options each time.
	 */
	public String[] ProcessParamsFile() {
		
		//TODO
		
		String CurrentDirectory = System.getProperty("user.dir"); // if deployment
	//	File CurrentDirectory = new File("/home/diego/Desktop/jolene/"); // if debug
		
		//TODO
		
		String LineBuffer = null;
		try {
			File ParamsFile = new File(CurrentDirectory, "DTASelect.params");
			// Read a DTASelect.params file, pulling default values from there
			if (ParamsFile.exists()) {
				if (ParamsFile.canRead()) {
					FileReader InputFileReader = new FileReader(ParamsFile);
					BufferedReader Incoming = new BufferedReader(
							InputFileReader);
					LineBuffer = Incoming.readLine();
				}
			}
		} catch (IOException failure) {
			return null;
		}
		if (LineBuffer == null)
			return null;
		else {
			StringTokenizer Parser = new StringTokenizer(LineBuffer);
			String Returned[];
			int Counter = 0;
			Returned = new String[Parser.countTokens()];
			while (Parser.hasMoreTokens()) {
				Returned[Counter] = Parser.nextToken();
				Counter++;
			}
			// changed by Howard Choi
			// this.SetCriteria(Returned);
			return Returned;
		}
	}

	/*
	 * Parse through the passed parameters to set criteria
	 */
	public void SetCriteria(String args[]) {
		int i;
		for (i = 0; i < args.length; i++) {
			// OPTIONS
			if (args[i].equals("--dcn")) {
				// Set DeltCN option
				SetDeltCN = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("--enzyme")) {
				int EnzymeNumber = new Integer(args[i + 1]).intValue();
				i++;
				if (EnzymeNumber == 99) {
					EnzymeCutsAfter = "";
					EnzymeCutsBefore = "";
				} else if (EnzymeNumber == 1) {
					EnzymeCutsAfter = "KR";
					EnzymeCutsBefore = "";
				} else if (EnzymeNumber == 2) {
					EnzymeCutsAfter = "K";
					EnzymeCutsBefore = "";
				} else if (EnzymeNumber == 3) {
					EnzymeCutsAfter = "";
					EnzymeCutsBefore = "D";
				}
			} else if (args[i].equals("--before")) {
				EnzymeCutsBefore = args[i + 1];
				i++;
			} else if (args[i].equals("--after")) {
				EnzymeCutsAfter = args[i + 1];
				i++;
			} else if (args[i].equals("--pI")) {
				DisplayPI = true;
			} else if (args[i].equals("--KD")) {
				DisplayKD = true;
			} else if (args[i].equals("--BB")) {
				DisplayBB = true;
			} else if (args[i].equals("--HPLC")) {
				DisplayHPLC = true;
			} else if (args[i].equals("--seqpos")) {
				DisplaySeqPos = true;
			} else if (args[i].equals("--dm")) {
				DisplayDeltaMass = true;
			} else if (args[i].equals("--extra")) {
				ExtraColumns = true;
			} else if (args[i].equals("--iso")) {
				Isotopes = true;
			} else if (args[i].equals("--noiso")) {
				Isotopes = false;
			} else if (args[i].equals("--full")) {
				FullProteinName = true;
			}
			// UTILITIES
			else if (args[i].equals("--omssa")) {
				UseXCorrCriteria = false;
				DisplayXCorrCriteria = false;
				UseDeltCN = false;
				UseRawXCorr = true;
			} else if (args[i].equals("-n") || args[i].equals("--nofilter")) {
				// Run without criteria
				UseCriteria = false;
			} else if (args[i].equals("--NoBackground")) {
				System.out.println("The --NoBackground option is obsolete.");
			} else if (args[i].equals("--CGI")) {
				System.out.println("The --CGI option is obsolete.");
			} else if (args[i].equals("--Mascot")) {
				UseMascotOutput = true;
			} else if (args[i].equals("--copy")) {
				PrintCopyList = true;
			} else if (args[i].equals("--chroma")) {
				PrintChroma = true;
			} else if (args[i].equals("--mods")) {
				PrintMods = true;
			} else if (args[i].equals("--align")) {
				PrintAlignment = true;
			} else if (args[i].equals("--similar")) {
				PrintSimilarities = true;
			} else if (args[i].equals("--XML")) {
				PrintXML = true;
			} else if (args[i].equals("--DB")) {
				PrintDB = true;
			} else if (args[i].equals("--class")) {
				UseClassifications = true;
			} else if (args[i].equals("--aux")) {
				UseAuxInfo = true;
			} else if (args[i].equals("-.") || args[i].equals("--here")) {
				LookInCurrentDirectory = true;
			} else if (args[i].equals("--GUI")) {
				UseGUIInstead = true;
			} else if (args[i].equals("--BE")) {
				DisplayBirdsEye = true;
			}
			else if (args[i].equals("--compress")) {
				CompressDTAs = true;
			} else if (args[i].equals("--quiet")) {
				Quiet = true;
			} else if (args[i].equals("--brief")) {
				Brief = true;
			} else if (args[i].equals("--short")) {
				ShortForm = true;
			} else if (args[i].equals("--path")) {
				UseCustomPath = true;
				CustomPath = args[i + 1];
				i++;
			}

			// STATISTICS OPTIONS
			else if (args[i].equals("--nostats")) {
				UseStatistics = false;
				MaxPepFP = 1.0f;
				MinPepConf = 0.0f;
				MaxProtFP = 1.0f;
				MinProtConf = 0.0f;
				MinSingleChargeXCorr = 1.8f;
				MinDoubleChargeXCorr = 2.5f;
				MinTripleChargeXCorr = 3.5f;
				MinQuadrupleChargeXCorr = 4.0f;
				MinDeltCN = 0.08f;
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinPepsPerLocus = 2;
				MinTrypticPeps = 0;
			} else if (args[i].equals("--decoy")) {
				DecoyLabel = args[i + 1];
				i++;
			} else if (args[i].equals("--trueprot")) {
				UseTrueProtein = true;
				TrueProteinLabel = args[i + 1];
				i++;
			} else if (args[i].equals("--hidedecoy")) {
				HideDecoy = true;
			} else if (args[i].equals("--noemp")) {
				UseEmpirics = false;
			} else if (args[i].equals("--raw")) {
				UseRawXCorr = true;
			} else if (args[i].equals("--nochgstat")) {
				UseChargeState = false;
			} else if (args[i].equals("--modstat")) {
				UseModStats = true;
			} else if (args[i].equals("--trypstat")) {
				UseTrypticInfo = true;
			} else if (args[i].equals("--noxc")) {
				UseXCorr = false;
			} else if (args[i].equals("--nodcn")) {
				UseDeltCN = false;
			} else if (args[i].equals("--sp")) {
				UseSp = true;
			} else if (args[i].equals("--mass")) {
				UseMassDifference = true;
				DisplayDeltaMass = true;
			} else if (args[i].equals("--logmass")) {
				UseMassDifference = true;
				UseLogMass = true;
			} else if (args[i].equals("--logspr")) {
				UseLogSpRank = true;
			} else if (args[i].equals("--ionstat")) {
				UseIonStat = true;
			} else if (args[i].equals("--fit")) {
				// Set statistics level
				FitLevel = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("--cmax")) {
				// Set maximum statistics charge state
				MaxStatisticsCharge = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("--plot")) {
				// Set plotting option level
				PlotLevel = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("--fptol")) {
				FPTol = new Float(args[i + 1]).doubleValue();
				i++;
			}

			// SPECTRUM FILTERS
			else if (args[i].equals("-c")) {
				// Set minimum charge state
				MinChargeState = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-C")) {
				// Set max charge state
				MaxChargeState = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-1")) {
				// Set minimum single charge XCorr
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinSingleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-2")) {
				// Set minimum double charge XCorr
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinDoubleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-3")) {
				// Set minimum triple charge XCorr
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinTripleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-4")) {
				// Set minimum quadruple charge XCorr
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinQuadrupleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-d")) {
				// Set minimum DeltCN
				UseXCorrCriteria = true;
				DisplayXCorrCriteria = true;
				MinDeltCN = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--fp")) {
				// Set maximum peptide false discovery rate
				MaxPepFP = new Float(args[i + 1]).floatValue();
			//	System.out.println("\\t"+MaxPepFP);
				i++;
			} else if (args[i].equals("--pr")) {
				// Set minimum peptide probability
				MinPepConf = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-i")) {
				// Set minimum Ion Proportion
				MinIonProportion = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-y")) {
				// Include only trypics?
				RequireTryptic = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("--nm")) {
				// Set minimum number of modification sites
				MinimumModifications = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("--NM")) {
				// Set maximum number of modification sites
				MaximumModifications = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("--mc")) {
				// Set minimum number of missed tryptic sites
				MinimumMissedTryptic = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("--MC")) {
				// Set maximum number of missed tryptic sites
				MaximumMissedTryptic = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-a")) {
				// Permit IDs with multiple equivalent scores?
				PermitAmbiguous = new Boolean(args[i + 1]).booleanValue();
				i++;
			} else if (args[i].equals("-m")) {
				// Include only modified peptides
				HandleModified = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-v")) {
				// How handle validation for peptides?
				PeptideValidation = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-s")) {
				// Set maximum Sp Rank
				MaximumSp = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-S")) {
				// Set minimum Sp Score
				MinSpScore = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--mz")) {
				// Set minimum m/z for peptides
				UseMinMZ = true;
				MinMZ = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--MZ")) {
				// Set maximum m/z for peptides
				UseMaxMZ = true;
				MaxMZ = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-file")) {
				/* Require that spectra come from this file */
				UseFilenameFilter = true;
				FilenameFilter = args[i + 1];
				i++;
			} else if (args[i].equals("-Sic")) {
				/*
				 * Require that included peptides have sequences containing all
				 * of these residues somewhere before C terminus
				 */
				UseSequenceCriteria = true;
				SeqMustIncludeResidues = args[i + 1];
				i++;
			} else if (args[i].equals("-Sip")) {
				/*
				 * Require that included peptides have sequences containing
				 * specified pattern
				 */
				UseSequenceCriteria = true;
				SeqMustIncludePattern = args[i + 1];
				i++;
			} else if (args[i].equals("-Sec")) {
				/*
				 * Require that included peptides have sequences excluding all
				 * of the following residues before the c-terminus residue
				 */
				UseSequenceCriteria = true;
				SeqMustExcludeBeforeC = args[i + 1];
				i++;
			} else if (args[i].equals("-Stc")) {
				/*
				 * Require that included peptides have sequences ending in one
				 * of the following residues
				 */
				UseSequenceCriteria = true;
				SeqMustEndWithResidues = args[i + 1];
				i++;
			} else if (args[i].equals("-Stn")) {
				/*
				 * Require that the included peptides follow one of the
				 * following residues
				 */
				UseSequenceCriteria = true;
				SeqMustComeAfterResidues = args[i + 1];
				i++;
			} else if (args[i].equals("-Smn")) {
				/*
				 * Require that the included peptides be at least this long
				 */
				UseSequenceCriteria = true;
				SeqMinLength = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-Smx")) {
				/*
				 * Require that the included peptides be no longer than this
				 */
				UseSequenceCriteria = true;
				SeqMaxLength = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-Ser")) {
				/*
				 * If IDs are mixed complete and incomplete sequences, determine
				 * which to keep.
				 */
				// 0 Keep all sequences
				// 1 Keep only one-ended sequences
				// 2 Keep only two-ended sequences
				UseSequenceCriteria = true;
				SequenceCompleteness = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-Ngly")) {
				// 0 Filter out N-linked glycosylated peptides
				// 1 Keep only N-linked glycosylated peptides
				UseSequenceCriteria = true;
				NGlycosylation = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-Si2")) {
				UseSequenceCriteria = true;
				SeqMustIncludeTwo = args[i + 1];
				i++;
			} else if (args[i].equals("-Si1")) {
				UseSequenceCriteria = true;
				SeqMustIncludeOne = args[i + 1];
				i++;
			} else if (args[i].equals("-Si0")) {
				UseSequenceCriteria = true;
				SeqMustIncludeNone = args[i + 1];
				i++;
			}
			else if (args[i].equals("-SiN")) {
				UseSequenceCriteria = true;
				SeqMustIncludeEither = args[i + 1];
				i++;
			} else if (args[i].equals("-X1")) {
				/*
				 * Require that the included +1 DTAs score less than this cutoff
				 */
				UseMaxXCorrs = true;
				MaxSingleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-X2")) {
				/*
				 * Require that the included +2 DTAs score less than this cutoff
				 */
				UseMaxXCorrs = true;
				MaxDoubleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-X3")) {
				/*
				 * Require that the included +3 DTAs score less than this cutoff
				 */
				UseMaxXCorrs = true;
				MaxTripleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-X4")) {
				/*
				 * Require that the included +4 DTAs score less than this cutoff
				 */
				UseMaxXCorrs = true;
				MaxQuadrupleChargeXCorr = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-dm")) {
				/* Minimum delta mass */
				DisplayDeltaMass = true;
				UseMinDeltaMass = true;
				MinDeltaMassPPM = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-DM")) {
				/* Maximum delta mass */
				DisplayDeltaMass = true;
				UseMaxDeltaMass = true;
				MaxDeltaMassPPM = new Float(args[i + 1]).floatValue();
				i++;
			}

			// LOCUS FILTERS
			else if (args[i].equals("-V")) {
				// How handle validation for peptides?
				LocusValidation = new Byte(args[i + 1]).byteValue();
				i++;
			} else if (args[i].equals("-p")) {
				// Set minimum peptides per locus
				MinPepsPerLocus = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-N")) {
				// Set minimum peptides per locus
				MinSpectraPerLocus = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-tp")) {
				// Set minimum tryptic peptides per locus
				MinTrypticPeps = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-ty")) {
				// Set minimum protein tryptic status
				MinProtTryptic = new Byte(args[i + 1]).byteValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("-in")) {
				// Include subset proteins
				RemoveSubsets = false;
			} else if (args[i].equals("-t")) {
				// Don't purge on TIC
				PurgeDuplicateSequences = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("--FP")) {
				// Set maximum protein false discovery rate
				MaxProtFP = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-tfp")) {
				// Set maximum protein peptide false discovery rate
				MaxProtMinFP = new Float(args[i + 1]).floatValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("-tpr")) {
				// Set minimum protein peptide probability
				MinProtPepConf = new Float(args[i + 1]).floatValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("-tDM")) {
				// Set maximum protein PPM delta mass
				MaxProtDM = new Float(args[i + 1]).floatValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("-tS")) {
				// Set minimum protein Sp Score
				MinProtSpScore = new Float(args[i + 1]).floatValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("-tX")) {
				// Set minimum protein XCorr
				MinProtXCorr = new Float(args[i + 1]).floatValue();
				UseProteinFilters = true;
				i++;
			} else if (args[i].equals("--CL")) {
				// Set minimum protein probability
				MinProtConf = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-M")) {
				// Include only proteins with this many modified
				// peptides
				MinModPepsPerLocus = new Integer(args[i + 1]).intValue();
				i++;
			} else if (args[i].equals("-u")) {
				// Include only loci with at least one unique peptide
				IncludeOnlyUniques = true;
			} else if (args[i].equals("-e")) {
				NameExcludePattern = args[i + 1];
				i++;
			} else if (args[i].equals("-E")) {
				NameIncludePattern = args[i + 1];
				i++;
			} else if (args[i].equals("-l")) {
				DescripExcludePattern = args[i + 1].toUpperCase();
				i++;
			} else if (args[i].equals("-L")) {
				DescripIncludePattern = args[i + 1].toUpperCase();
				i++;
			} else if (args[i].equals("--mw")) {
				// Set minimum molwt for proteins
				UseMinMW = true;
				MinMW = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--MW")) {
				// Set maximum molwt for proteins
				UseMaxMW = true;
				MaxMW = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--Seq")) {
				// Set minimum sequence coverage
				MinSequenceCoverage = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("-StC")) {
				/*
				 * Require that protein sequence ends in one of the following
				 * residues
				 */
				ProteinMustEndWithResidues = args[i + 1];
				i++;
			} else if (args[i].equals("-StN")) {
				/*
				 * Require that protein sequence starts with one of the
				 * following residues
				 */
				ProteinMustStartWithResidues = args[i + 1];
				i++;
			} else if (args[i].equals("--pi")) {
				// Set minimum pI for proteins
				UseMinPI = true;
				MinPI = new Float(args[i + 1]).floatValue();
				i++;
			} else if (args[i].equals("--PI")) {
				// Set maximum molwt for proteins
				UseMaxPI = true;
				MaxPI = new Float(args[i + 1]).floatValue();
				i++;
			}
			else if(args[i].equals("-no-db") || args[i].equals("-noDB"))
			{
				noDB = true;
			}
			else if(args[i].equals("--readFasta"))
			{
				readFasta = true;
			}
			else if(args[i].startsWith("-javaagent:")|| args[i].startsWith("-Didea.launcher") || args[i].startsWith("-Dfile") || args[i].startsWith("com.intellij") || args[i].startsWith("DTASelect"))
			{
				//to allow Intellij to debug the program
			}
			else if(args[i].equals("-U"))
			{
				PrintOnlyUnique =true;
				IncludeOnlyUniques = true;
			}
			else if(args[i].equals("--ms"))
			{
				UseMultipleScore = true;
			} else if (args[i].equals("-ug")) {
				UseProteinGrouper = true;

			} else if(args[i].equals("-addPSM"))
			{
				addPSM = true;
			}
			else if(args[i].equals("-DMS"))
			{
				dms = true;
				DisplayDeltaMass = true;
				int index = i+1;
				int showCorrectedDmsArg = Integer.parseInt(args[index]);
				showCorrectedDmValue = showCorrectedDmsArg == 1;
				i++;

			}
			else if(args[i].equals("--printHTML"))
			{
				printHTML = true;
			}
			else if(args[i].startsWith("--xcorrFilterReadSqt"))
			{
				filterReadXorr = true;
				if(args[i].contains("="))
				{
					String limitStr = args[i].split("=")[1];
					double limit = Double.parseDouble(limitStr);
					minReadXorr = limit;
				}
				else {
					String limitStr = args[i + 1];
					double limit = Double.parseDouble(limitStr);
					minReadXorr = limit;
					i++;
				}
			}
			else if(args[i].startsWith("--filterMod="))
			{
				filterByModification = true;
				modFilterStr = args[i].split("=")[1];
			}
			else {
				System.out.println("I don't understand this option:  "
						+ args[i]);
				System.exit(0);
			}
		}
		// LocusGUI GUI = new LocusGUI();
	}

	/*
	 * Check the passed protein to see if its name and description match
	 * filters. Return true if it's okay, and return false if it's not.
	 */
	public boolean Permit(Protein TestSubject) {
		boolean Okay = true;
		// Handle -e, -E, -l, -L, and decoy options
		if (HideDecoy) {
			if (TestSubject.Locus.startsWith(DecoyLabel))
				Okay = false;
		}
		if (NameExcludePattern.length() > 0) {
			if (TestSubject.Locus.indexOf(NameExcludePattern) > -1)
				Okay = false;
		}
		if (ProteinMustStartWithResidues.length() > 0) {
			if (ProteinMustStartWithResidues.indexOf(TestSubject.NChar) == -1)
				Okay = false;
		}
		if (ProteinMustEndWithResidues.length() > 0) {
			if (ProteinMustEndWithResidues.indexOf(TestSubject.CChar) == -1)
				Okay = false;
		}
		if (NameIncludePattern.length() > 0) {
			if (TestSubject.Locus.indexOf(NameIncludePattern) == -1)
				Okay = false;
		}
		if (DescripExcludePattern.length() > 0) {
			if (TestSubject.Gene.toUpperCase().indexOf(DescripExcludePattern) > -1)
				Okay = false;
		}
		if (DescripIncludePattern.length() > 0) {
			if (TestSubject.Gene.toUpperCase().indexOf(DescripIncludePattern) == -1)
				Okay = false;
		}
		// Handle statistics options
		if (TestSubject.ProtConf < MinProtConf)
			Okay = false;
		if (TestSubject.ProtFP > MaxProtFP)
			Okay = false;

		// Handle --mw, --MW options
		if (UseMinMW) {
			if (MinMW > TestSubject.MolWt)
				Okay = false;
		}
		if (UseMaxMW) {
			if (MaxMW < TestSubject.MolWt)
				Okay = false;
		}
		// Handle --pi, --PI options
		if (UseMinPI) {
			if (MinPI > TestSubject.pI)
				Okay = false;
		}
		if (UseMaxPI) {
			if (MaxPI < TestSubject.pI)
				Okay = false;
		}
		// Handle validation info
		switch (LocusValidation) {
		case -1:
			if (TestSubject.Validated == 'N')
				Okay = true;
			else
				Okay = false;
			break;
		case 1:
			if (TestSubject.Validated == 'Y')
				Okay = true;
			else if (TestSubject.Validated == 'N')
				Okay = false;
			break;
		case 2:
			if (TestSubject.Validated == 'Y' || TestSubject.Validated == 'M')
				Okay = true;
			else if (TestSubject.Validated == 'N')
				Okay = false;
			break;
		case 3:
			if (TestSubject.Validated == 'Y' || TestSubject.Validated == 'M')
				Okay = true;
			else if (TestSubject.Validated == 'N'
					|| TestSubject.Validated == 'U')
				Okay = false;
			break;
		}
		return Okay;
	}

	public boolean Allow(DTAFile TestSubject, boolean fromDTAReader) {
		/*
		 * Evaluates an individual identification to determine if it fits
		 * selection criteria
		 */
		boolean PassesSequenceCriteria = true;
		boolean PassesMaxXCorrs = true;
		boolean PassesXCorrCriteria = true;
		boolean PassesMZ = true;
		boolean PassesDeltaMass = true;
		boolean PassesFilenameFilter = true;
		boolean PassesUnique = true;
		float ObservedMZ;
		boolean PassesOnMerits;
		String TSequence = TestSubject.TrimmedSequence();
		//System.out.println(">>>"+TestSubject.Sequence);
		String BetweenDots = TestSubject.BetweenDots();
		// Handle "Smn" option
		if (TSequence.length() < SeqMinLength) {
			PassesSequenceCriteria = false;
		}
		// Handle "Smx" option
		if (TSequence.length() > SeqMaxLength) {
			PassesSequenceCriteria = false;
		}
		if (UseSequenceCriteria) {
			int Counter;
			int HowMany;
			int InsideCounter;
			String ShortSeq = TSequence.substring(0, TSequence.length() - 1);
			// Handle "Sic" option
			for (Counter = 0; Counter < SeqMustIncludeResidues.length(); Counter++) {
				if (TSequence.indexOf(SeqMustIncludeResidues.charAt(Counter)) == -1)
					PassesSequenceCriteria = false;
			}
			// Handle "Sip" option
			if (BetweenDots.indexOf(SeqMustIncludePattern) == -1)
				PassesSequenceCriteria = false;
			// Handle "Sec" option
			for (Counter = 0; Counter < SeqMustExcludeBeforeC.length(); Counter++) {
				if (TSequence.indexOf(SeqMustExcludeBeforeC.charAt(Counter)) != -1)
					PassesSequenceCriteria = false;
			}
			// Handle "Stc" option
			if (SeqMustEndWithResidues.length() > 0) {
				if (SeqMustEndWithResidues.indexOf(TSequence.charAt(TSequence
						.length() - 1)) == -1)
					PassesSequenceCriteria = false;
			}
			// Handle "Stn" option
			if (SeqMustComeAfterResidues.length() > 0) {
				if (SeqMustComeAfterResidues.indexOf(TestSubject.Sequence
						.charAt(0)) == -1)
					PassesSequenceCriteria = false;
			}
			// Handle "Si2" option
			for (Counter = 0; Counter < SeqMustIncludeTwo.length(); Counter++) {
				HowMany = 0;
				for (InsideCounter = 0; InsideCounter < TSequence.length(); InsideCounter++) {
					if (TSequence.charAt(InsideCounter) == SeqMustIncludeTwo
							.charAt(Counter))
						HowMany++;
				}
				if ((HowMany < 2) || (HowMany > 2))
					PassesSequenceCriteria = false;
			}
			// Handle "Si1" option
			for (Counter = 0; Counter < SeqMustIncludeOne.length(); Counter++) {
				HowMany = 0;
				for (InsideCounter = 0; InsideCounter < TSequence.length(); InsideCounter++) {
					if (TSequence.charAt(InsideCounter) == SeqMustIncludeOne
							.charAt(Counter))
						HowMany++;
				}
				if ((HowMany < 1) || (HowMany > 1))
					PassesSequenceCriteria = false;
			}
			//handles "SiN" option
			int siNCount = 0;
			for (Counter = 0; Counter < SeqMustIncludeEither.length(); Counter++) {

				if (TSequence.contains(Character.toString(SeqMustIncludeEither.charAt(Counter)))) {
					siNCount++;
				}
			}
			if (SeqMustIncludeEither.length() > 0) {

				if (siNCount > 0) {
					PassesSequenceCriteria = true;
				} else {
					PassesSequenceCriteria = false;
				}
			}
			// Handle "Si0" option
			for (Counter = 0; Counter < SeqMustIncludeNone.length(); Counter++) {
				HowMany = 0;
				for (InsideCounter = 0; InsideCounter < TSequence.length(); InsideCounter++) {
					if (TSequence.charAt(InsideCounter) == SeqMustIncludeNone
							.charAt(Counter))
						HowMany++;
				}
				if (HowMany > 0)
					PassesSequenceCriteria = false;
			}
			// Handle "Ser" option
			if (SequenceCompleteness > 0) {
				boolean NComplete = BetweenDots.charAt(0) != '-';
				boolean CComplete = BetweenDots
						.charAt(BetweenDots.length() - 1) != '-';
				if (SequenceCompleteness == 1) {
					if ((NComplete && CComplete) || (!NComplete && !CComplete))
						PassesSequenceCriteria = false;
				} else if (SequenceCompleteness == 2) {
					if (!(NComplete && CComplete)) {
						PassesSequenceCriteria = false;
					}
				}
			}
			// Handle "Ngly" option
			if (NGlycosylation == 0 || NGlycosylation == 1) {
				boolean FoundNGly = false;

				for (Counter = 0; Counter < TSequence.length() - 2; Counter++) {
					if (TSequence.charAt(Counter) == 'N') {
						InsideCounter = Counter + 1;
						if (TSequence.charAt(InsideCounter) != 'P') {
							InsideCounter = Counter + 2;
							if (TSequence.charAt(InsideCounter) == 'S'
									|| TSequence.charAt(InsideCounter) == 'T') {
								FoundNGly = true;
							}
						}
					}
				}
				if (NGlycosylation == 1 && !FoundNGly) {
					PassesSequenceCriteria = false;
				} else if (NGlycosylation == 0 && FoundNGly) {
					PassesSequenceCriteria = false;
				}
			}
		}
		if (UseMaxXCorrs) {
			PassesMaxXCorrs = ((TestSubject.ChargeState == 1) && (TestSubject.XCorr < this.MaxSingleChargeXCorr))
					|| ((TestSubject.ChargeState == 2) && (TestSubject.XCorr < this.MaxDoubleChargeXCorr))
					|| ((TestSubject.ChargeState == 3) && (TestSubject.XCorr < this.MaxTripleChargeXCorr))
					|| ((TestSubject.ChargeState == 4) && (TestSubject.XCorr < this.MaxQuadrupleChargeXCorr));
		}
		if (UseXCorrCriteria) {
			PassesXCorrCriteria = (((TestSubject.ChargeState == 1) && (TestSubject.XCorr > this.MinSingleChargeXCorr))
					|| ((TestSubject.ChargeState == 2) && (TestSubject.XCorr > this.MinDoubleChargeXCorr))
					|| ((TestSubject.ChargeState == 3) && (TestSubject.XCorr > this.MinTripleChargeXCorr)) || ((TestSubject.ChargeState >= 4) && (TestSubject.XCorr > this.MinQuadrupleChargeXCorr)))
					&& (TestSubject.DeltCN > this.MinDeltCN);
		}
		ObservedMZ = (TestSubject.PrecursorMass - 1 + TestSubject.ChargeState)
				/ TestSubject.ChargeState;
		if (UseMinMZ) {
			PassesMZ = ObservedMZ > MinMZ;
		}
		if (UseMaxMZ) {
			PassesMZ &= ObservedMZ < MaxMZ;
		}
		if (!dms)
		{
			if (UseMinDeltaMass) {
				PassesDeltaMass = Math.abs(TestSubject.Adjusted_PPM_Offset) > MinDeltaMassPPM;
			}
			if (UseMaxDeltaMass) {
				PassesDeltaMass &= Math.abs(TestSubject.Adjusted_PPM_Offset) < MaxDeltaMassPPM;
			}
		}
		if (UseFilenameFilter) {
			if (TestSubject.RootFileName.indexOf(FilenameFilter) == -1) {
				PassesFilenameFilter = false;
			}
		}
		if(PrintOnlyUnique)
		{
			if(!TestSubject.UniqueToLocus)
			{
				PassesUnique  = false;
			}
		}
		PassesOnMerits = (PassesSequenceCriteria
				&& PassesMaxXCorrs
				&& PassesMZ
				&& PassesDeltaMass
				&& PassesXCorrCriteria
				&& PassesFilenameFilter
				&& PassesUnique
				&& (TestSubject.ChargeState >= this.MinChargeState)
				&& (TestSubject.ChargeState <= this.MaxChargeState)
				&& (((HandleModified == 0) && TestSubject.Modified)
				|| ((HandleModified == 2) && !TestSubject.Modified) || (HandleModified == 1))
				&& (!(TestSubject.Tryptic < RequireTryptic))
				&& (PermitAmbiguous || TestSubject.EquivSeq < 2)
				&& (TestSubject.PepConf >= this.MinPepConf)
				&& (TestSubject.PepFP <= this.MaxPepFP)
				&& (TestSubject.NMods >= this.MinimumModifications)
				&& (TestSubject.NMods <= this.MaximumModifications)
				&& (TestSubject.TrypticSites >= this.MinimumMissedTryptic)
				&& (TestSubject.TrypticSites <= this.MaximumMissedTryptic)
				&& (TestSubject.IonProportion >= this.MinIonProportion)
				&& (TestSubject.Sp < this.MaximumSp) && (TestSubject.SpScore > this.MinSpScore)  );
		switch (PeptideValidation) {
			case -1:
				if (TestSubject.Validated == 'N')
					return true;
				else
					return false;
			case 1:
				switch (TestSubject.Validated) {
					case 'N':
						return false;
					case 'Y':
						return true;
					case 'M':
					case 'U':
						return PassesOnMerits;
					default:
						System.out.println("Validation went awry for "
								+ TestSubject.FileName);
						return PassesOnMerits;
				}
			case 2:
				switch (TestSubject.Validated) {
					case 'N':
						return false;
					case 'Y':
					case 'M':
						return true;
					case 'U':
						return PassesOnMerits;
					default:
						System.out.println("Validation went awry for "
								+ TestSubject.FileName);
						return PassesOnMerits;
				}
			case 3:
				switch (TestSubject.Validated) {
					case 'N':
					case 'U':
						return false;
					case 'Y':
					case 'M':
						return true;
					default:
						System.out.println("Validation went awry for "
								+ TestSubject.FileName);
						return PassesOnMerits;
				}
			default:
				return PassesOnMerits;
		}
	}

	/*
	 * Apply the spectrum-specific filters. Return true if the identification
	 * meets the criteria or false if it doesn't.
	 */
	public boolean Allow(DTAFile TestSubject) {
		return Allow(TestSubject, false);
	}

	public boolean AllowShiftDM(DTAFile TestSubject) {
		boolean PassesDeltaMass = true;
		if (UseMinDeltaMass) {
			PassesDeltaMass = Math.abs(TestSubject.Shifted_PPM_Offset) > MinDeltaMassPPM;
		}
		if (UseMaxDeltaMass) {
			PassesDeltaMass &= Math.abs(TestSubject.Shifted_PPM_Offset) < MaxDeltaMassPPM;
		}
		return PassesDeltaMass;
	}
	/*
	 * Create a string announcing the current settings for each of the criteria.
	 */
	public String PrintCriteria(String RowStart, String RowEnd, String Separator) {
		String Purge = "false";
		String Tryptic = "Any";
		String Mods = "Include";
		String PepValid = "Ignore";
		String LocValid = "Ignore";
		if (PurgeDuplicateSequences == 1)
			Purge = "Salt step";
		else if (PurgeDuplicateSequences == 2)
			Purge = "XCorr";
		if (RequireTryptic == 1)
			Tryptic = "Half";
		else if (RequireTryptic == 2)
			Tryptic = "Full";
		if (HandleModified == 0)
			Mods = "Require";
		else if (HandleModified == 2)
			Mods = "Exclude";
		if (PeptideValidation == 1)
			PepValid = "+Y -N";
		else if (PeptideValidation == 2)
			PepValid = "+YM -N";
		else if (PeptideValidation == 3)
			PepValid = "+YM -NU";
		else if (PeptideValidation == -1)
			PepValid = "+N -YMU";
		if (LocusValidation == 1)
			LocValid = "+Y -N";
		else if (LocusValidation == 2)
			LocValid = "+YM -N";
		else if (LocusValidation == 3)
			LocValid = "+YM -NU";
		else if (LocusValidation == -1)
			LocValid = "+N -YMU";
		return (RowStart
				+ new Boolean(UseCriteria).toString()
				+ Separator
				+ "Use criteria"
				+ RowEnd
				+ RowStart
				+ new Float(MinPepConf).toString()
				+ Separator
				+ "Minimum peptide probability"
				+ RowEnd
				+ RowStart
				+ new Float(MaxPepFP).toString()
				+ Separator
				+ "Peptide global false discovery rate"
				+ RowEnd
				+ RowStart
				+ new Float(MinProtConf).toString()
				+ Separator
				+ "Minimum protein probability"
				+ RowEnd
				+ RowStart
				+ new Float(MaxProtFP).toString()
				+ Separator
				+ "Protein false discovery rate"
				+ RowEnd
				+ (DisplayXCorrCriteria ? (RowStart
						+ new Float(MinSingleChargeXCorr).toString()
						+ Separator + "Minimum +1 XCorr" + RowEnd + RowStart
						+ new Float(MinDoubleChargeXCorr).toString()
						+ Separator + "Minimum +2 XCorr" + RowEnd + RowStart
						+ new Float(MinTripleChargeXCorr).toString()
						+ Separator + "Minimum +3 XCorr" + RowEnd + RowStart
						+ new Float(MinQuadrupleChargeXCorr).toString()
						+ Separator + "Minimum +4 XCorr" + RowEnd + RowStart
						+ new Float(MinDeltCN).toString() + Separator
						+ "Minimum DeltCN" + RowEnd) : "")
				+ (UseMaxXCorrs ? (RowStart + MaxSingleChargeXCorr + Separator
						+ "Maximum +1 XCorr" + RowEnd + RowStart
						+ MaxDoubleChargeXCorr + Separator + "Maximum +2 XCorr"
						+ RowEnd + RowStart + MaxTripleChargeXCorr + Separator
						+ "Maximum +3 XCorr" + RowEnd + RowStart
						+ MaxQuadrupleChargeXCorr + Separator
						+ "Maximum +4 XCorr" + RowEnd) : "")
				+ RowStart
				+ new Integer(MinChargeState).toString()
				+ Separator
				+ "Minimum charge state"
				+ RowEnd
				+ RowStart
				+ new Integer(MaxChargeState).toString()
				+ Separator
				+ "Maximum charge state"
				+ RowEnd
				+ RowStart
				+ new Float(MinIonProportion).toString()
				+ Separator
				+ "Minimum ion proportion"
				+ RowEnd
				+ RowStart
				+ new Integer(MaximumSp).toString()
				+ Separator
				+ "Maximum Sp rank"
				+ RowEnd
				+ RowStart
				+ new Float(MinSpScore).toString()
				+ Separator
				+ "Minimum Sp score"
				+ RowEnd
				+ (UseMinDeltaMass ? (RowStart + MinDeltaMassPPM + Separator
						+ "Minimum mass deviation (ppm)" + RowEnd) : "")
				+ (UseMaxDeltaMass ? (RowStart + MaxDeltaMassPPM + Separator
						+ "Maximum mass deviation (ppm)" + RowEnd) : "")
				+ (DisplayDeltaMass ? (RowStart + Isotopes + Separator
						+ "Define delta mass with respect to nearest isotope" + RowEnd)
						: "")
				+ RowStart
				+ Mods
				+ Separator
				+ "Modified peptide inclusion"
				+ RowEnd
				+ RowStart
				+ Tryptic
				+ Separator
				+ "Tryptic status requirement"
				+ RowEnd
				+ RowStart
				+ PermitAmbiguous
				+ Separator
				+ "Multiple, ambiguous IDs allowed"
				+ RowEnd
				+ RowStart
				+ PepValid
				+ Separator
				+ "Peptide validation handling"
				+ RowEnd
				+ (UseSequenceCriteria ? (RowStart + SeqMustIncludeResidues
						+ Separator + "Sequences must include all of" + RowEnd
						+ RowStart + SeqMustIncludePattern + Separator
						+ "Sequences must include pattern" + RowEnd + RowStart
						+ SeqMustExcludeBeforeC + Separator
						+ "Sequences can contain none of" + RowEnd + RowStart
						+ SeqMustComeAfterResidues + Separator
						+ "Preceding residue must be one of" + RowEnd
						+ RowStart + SeqMustEndWithResidues + Separator
						+ "C terminal residue must be one of" + RowEnd
						+ RowStart + SeqMinLength + Separator
						+ "Minimum sequence length" + RowEnd + RowStart
						+ SeqMaxLength + Separator + "Maximum sequence length"
						+ RowEnd + RowStart + SequenceCompleteness + Separator
						+ "Sequence completeness required" + RowEnd) : "")
				+ (UseMinMZ ? (RowStart + MinMZ + Separator
						+ "Minimum peptide m/z" + RowEnd) : "")
				+ (UseMaxMZ ? (RowStart + MaxMZ + Separator
						+ "Maximum peptide m/z" + RowEnd) : "")
				+ RowStart
				+ Purge
				+ Separator
				+ "Purge duplicate peptides by protein"
				+ RowEnd
				+ RowStart
				+ new Boolean(IncludeOnlyUniques).toString()
				+ Separator
				+ "Include only loci with unique peptide"
				+ RowEnd
				+ RowStart
				+ new Boolean(RemoveSubsets).toString()
				+ Separator
				+ "Remove subset proteins"
				+ RowEnd
				+ RowStart
				+ LocValid
				+ Separator
				+ "Locus validation handling"
				+ RowEnd
				+ (NameExcludePattern.length() > 0 ? (RowStart
						+ NameExcludePattern + Separator
						+ "Exclude protein names matching" + RowEnd) : "")
				+ (NameIncludePattern.length() > 0 ? (RowStart
						+ NameIncludePattern + Separator
						+ "Include only protein names matching" + RowEnd) : "")
				+ (DescripExcludePattern.length() > 0 ? (RowStart
						+ DescripExcludePattern + Separator
						+ "Exclude protein descriptions matching" + RowEnd)
						: "")
				+ (DescripIncludePattern.length() > 0 ? (RowStart
						+ DescripIncludePattern + Separator
						+ "Include only protein descriptions matching" + RowEnd)
						: "")
				+ (UseMinMW ? (RowStart + MinMW + Separator
						+ "Minimum protein molecular weight" + RowEnd) : "")
				+ (UseMaxMW ? (RowStart + MaxMW + Separator
						+ "Maximum protein molecular weight" + RowEnd) : "")
				+ (UseMinPI ? (RowStart + MinPI + Separator
						+ "Minimum protein isoelectric point" + RowEnd) : "")
				+ (UseMaxPI ? (RowStart + MaxPI + Separator
						+ "Maximum protein isoelectric point" + RowEnd) : "")
				+ RowStart + new Integer(MinModPepsPerLocus).toString()
				+ Separator + "Minimum modified peptides per locus" + RowEnd
				+ RowStart + new Integer(MinPepsPerLocus).toString()
				+ Separator + "Minimum peptides per locus" + RowEnd);
	}

	/*
	 * Create a string announcing the current settings for each of the criteria.
	 * The formatting is XML. NOTE: THIS IS OUT-OF-DATE!
	 */
	public String PrintXMLCriteria() {
		String Purge = "false";
		String Tryptic = "Any";
		String Mods = "Include";
		if (PurgeDuplicateSequences == 1)
			Purge = "Salt step";
		else if (PurgeDuplicateSequences == 2)
			Purge = "XCorr";
		if (RequireTryptic == 1)
			Tryptic = "Half";
		else if (RequireTryptic == 2)
			Tryptic = "Full";
		if (HandleModified == 0)
			Mods = "Require";
		else if (HandleModified == 2)
			Mods = "Exclude";
		return ("\t<DTASelect:UseCriteria>"
				+ new Boolean(UseCriteria).toString()
				+ "</DTASelect:UseCriteria>\n" + "\t<DTASelect:Min1XCorr>"
				+ new Float(MinSingleChargeXCorr).toString()
				+ "</DTASelect:Min1XCorr>\n" + "\t<DTASelect:Min2XCorr>"
				+ new Float(MinDoubleChargeXCorr).toString()
				+ "</DTASelect:Min2XCorr>\n" + "\t<DTASelect:Min3XCorr>"
				+ new Float(MinTripleChargeXCorr).toString()
				+ "</DTASelect:Min3XCorr>\n" + "\t<DTASelect:Min4XCorr>"
				+ new Float(MinQuadrupleChargeXCorr).toString()
				+ "</DTASelect:Min4XCorr>\n" + "\t<DTASelect:Max1XCorr>"
				+ new Float(MaxSingleChargeXCorr).toString()
				+ "</DTASelect:Max1XCorr>\n" + "\t<DTASelect:Max2XCorr>"
				+ new Float(MaxDoubleChargeXCorr).toString()
				+ "</DTASelect:Max2XCorr>\n" + "\t<DTASelect:Max3XCorr>"
				+ new Float(MaxTripleChargeXCorr).toString()
				+ "</DTASelect:Max3XCorr>\n" + "\t<DTASelect:Max4XCorr>"
				+ new Float(MaxQuadrupleChargeXCorr).toString()
				+ "</DTASelect:Max4XCorr>\n" + "\t<DTASelect:MinDeltCN>"
				+ new Float(MinDeltCN).toString() + "</DTASelect:MinDeltCN>\n"
				+ "\t<DTASelect:MinPepConf>" + new Float(MinPepConf).toString()
				+ "</DTASelect:MinPepConf>\n" + "\t<DTASelect:MaxPepFP>"
				+ new Float(MaxPepFP).toString() + "</DTASelect:MaxPepFP>\n"
				+ "\t<DTASelect:MinProtConf>"
				+ new Float(MinProtConf).toString()
				+ "</DTASelect:MinProtConf>\n" + "\t<DTASelect:MaxProtFP>"
				+ new Float(MaxProtFP).toString() + "</DTASelect:MaxProtFP>\n"
				+ "\t<DTASelect:MinChargeState>"
				+ new Integer(MinChargeState).toString()
				+ "</DTASelect:MinChargeState>\n"
				+ "\t<DTASelect:MaxChargeState>"
				+ new Integer(MaxChargeState).toString()
				+ "</DTASelect:MaxChargeState>\n"
				+ "\t<DTASelect:MinIonProportion>"
				+ new Float(MinIonProportion).toString()
				+ "</DTASelect:MinIonProportion>\n" + "\t<DTASelect:MaxSpRank>"
				+ new Integer(MaximumSp).toString()
				+ "</DTASelect:MaxSpRank>\n" + "\t<DTASelect:SetDeltCN>"
				+ new Integer(SetDeltCN).toString()
				+ "</DTASelect:SetDeltCN>\n"
				+ "\t<DTASelect:ModifiedInclusion>" + Mods
				+ "</DTASelect:ModifiedInclusion>\n"
				+ "\t<DTASelect:TrypticInclusion>" + Tryptic
				+ "</DTASelect:TrypticInclusion>\n"
				+ "\t<DTASelect:SequenceIncRes>" + SeqMustIncludeResidues
				+ "</DTASelect:SequenceIncRes>\n"
				+ "\t<DTASelect:SequenceIncPat>" + SeqMustIncludePattern
				+ "</DTASelect:SequenceIncPat>\n"
				+ "\t<DTASelect:SequenceExcRes>" + SeqMustExcludeBeforeC
				+ "</DTASelect:SequenceExcRes>\n"
				+ "\t<DTASelect:SequenceFollows>" + SeqMustComeAfterResidues
				+ "</DTASelect:SequenceFollows>\n"
				+ "\t<DTASelect:SequenceEndsWith>" + SeqMustEndWithResidues
				+ "</DTASelect:SequenceEndsWith>\n"
				+ "\t<DTASelect:PurgeDuplicates>" + Purge
				+ "</DTASelect:PurgeDuplicates>\n"
				+ "\t<DTASelect:OnlyUniques>"
				+ new Boolean(IncludeOnlyUniques).toString()
				+ "</DTASelect:OnlyUniques>\n" + "\t<DTASelect:SubsetsRemoved>"
				+ new Boolean(RemoveSubsets).toString()
				+ "</DTASelect:SubsetsRemoved>\n" + "\t<DTASelect:ExcludeName>"
				+ NameExcludePattern + "</DTASelect:ExcludeName>\n"
				+ "\t<DTASelect:IncludeName>" + NameIncludePattern
				+ "</DTASelect:IncludeName>\n" + "\t<DTASelect:ExcludeDescrip>"
				+ DescripExcludePattern + "</DTASelect:ExcludeDescrip>\n"
				+ "\t<DTASelect:IncludeDescrip>" + DescripIncludePattern
				+ "</DTASelect:IncludeDescrip>\n"
				+ "\t<DTASelect:MinModPeptidesPerLocus>"
				+ new Integer(MinModPepsPerLocus).toString()
				+ "</DTASelect:MinModPeptidesPerLocus>\n"
				+ "\t<DTASelect:MinPeptidesPerLocus>"
				+ new Integer(MinPepsPerLocus).toString() + "</DTASelect:MinPeptidesPerLocus>\n");
	}

	class SpectrumGUI extends Frame {
		private Checkbox n = new Checkbox("Use Criteria", UseCriteria);
		private TextField c = new TextField(
				new Integer(MinChargeState).toString());
		private TextField C = new TextField(
				new Integer(MaxChargeState).toString());
		private TextField s = new TextField(new Integer(MaximumSp).toString());
		private TextField z1 = new TextField(
				new Float(MinSingleChargeXCorr).toString());
		private TextField z2 = new TextField(
				new Float(MinDoubleChargeXCorr).toString());
		private TextField z3 = new TextField(
				new Float(MinTripleChargeXCorr).toString());
		private TextField z4 = new TextField(
				new Float(MinQuadrupleChargeXCorr).toString());
		private TextField d = new TextField(new Float(MinDeltCN).toString());
		private TextField cl = new TextField(new Float(MinPepConf).toString());
		private TextField fp = new TextField(new Float(MaxPepFP).toString());
		private TextField i = new TextField(new Float(MinDeltCN).toString());
		private CheckboxGroup m = new CheckboxGroup();
		private Checkbox m0 = new Checkbox("Require peptides to be modified",
				m, false);
		private Checkbox m1 = new Checkbox(
				"Include peptides regardless of modification", m, true);
		private Checkbox m2 = new Checkbox("Exclude modified peptides", m,
				false);
		private CheckboxGroup y = new CheckboxGroup();
		private Checkbox y0 = new Checkbox(
				"Include peptides regardless of tryptic status", y, true);
		private Checkbox y1 = new Checkbox(
				"Include only half- or fully tryptic peptides", y, false);
		private Checkbox y2 = new Checkbox(
				"Include only fully tryptic peptides", y, false);

		public SpectrumGUI() {
			GridBagConstraints NormalConstraints = new GridBagConstraints();
			GridBagLayout Layout = new GridBagLayout();
			Label TempLabel;
			switch (RequireTryptic) {
			case 0:
				y.setSelectedCheckbox(y0);
				break;
			case 1:
				y.setSelectedCheckbox(y1);
				break;
			case 2:
				y.setSelectedCheckbox(y2);
				break;
			}
			this.setLayout(Layout);
			NormalConstraints.fill = GridBagConstraints.BOTH;
			NormalConstraints.gridx = 0;
			NormalConstraints.gridwidth = 2;
			TempLabel = new Label("Standard Individual Spectrum Filters");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			Layout.setConstraints(n, NormalConstraints);
			NormalConstraints.gridwidth = 1;
			Layout.setConstraints(z1, NormalConstraints);
			Layout.setConstraints(z2, NormalConstraints);
			Layout.setConstraints(z3, NormalConstraints);
			Layout.setConstraints(d, NormalConstraints);
			Layout.setConstraints(cl, NormalConstraints);
			Layout.setConstraints(i, NormalConstraints);
			Layout.setConstraints(s, NormalConstraints);
			Layout.setConstraints(c, NormalConstraints);
			Layout.setConstraints(C, NormalConstraints);
			add(n);
			add(z1);
			add(z2);
			add(z3);
			add(d);
			add(i);
			add(s);
			add(c);
			add(C);
			NormalConstraints.gridx = 1;
			TempLabel = new Label("Lowest +1 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest +2 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest +3 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest +4 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest DeltCN");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest peptide probability");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Peptide false discovery rate");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Lowest proportion of fragment ions");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Maximum rank by Sp");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Minimum charge state");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Maximum charge state");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			NormalConstraints.gridx = 0;
			NormalConstraints.gridwidth = 2;
			Layout.setConstraints(m0, NormalConstraints);
			Layout.setConstraints(m1, NormalConstraints);
			Layout.setConstraints(m2, NormalConstraints);
			Layout.setConstraints(y0, NormalConstraints);
			Layout.setConstraints(y1, NormalConstraints);
			Layout.setConstraints(y2, NormalConstraints);
			add(m0);
			add(m1);
			add(m2);
			add(y0);
			add(y1);
			add(y2);
			this.setTitle("DTASelect Spectrum Filters");
			this.setSize(350, 450);
			this.addWindowListener(new SelectCriteria.WindowEventAdapter());
			this.show();
		}
	}

	class ExtendedGUI extends Frame {
		private TextField Sic = new TextField(SeqMustIncludeResidues);
		private TextField Sip = new TextField(SeqMustIncludePattern);
		private TextField Sec = new TextField(SeqMustExcludeBeforeC);
		private TextField Stn = new TextField(SeqMustComeAfterResidues);
		private TextField Stc = new TextField(SeqMustEndWithResidues);
		private TextField X1 = new TextField(
				new Float(MaxSingleChargeXCorr).toString());
		private TextField X2 = new TextField(
				new Float(MaxDoubleChargeXCorr).toString());
		private TextField X3 = new TextField(
				new Float(MaxTripleChargeXCorr).toString());
		private TextField X4 = new TextField(
				new Float(MaxQuadrupleChargeXCorr).toString());

		public ExtendedGUI() {
			GridBagConstraints NormalConstraints = new GridBagConstraints();
			GridBagLayout Layout = new GridBagLayout();
			Label TempLabel;
			NormalConstraints.fill = GridBagConstraints.BOTH;
			NormalConstraints.gridx = 0;
			NormalConstraints.gridwidth = 2;
			this.setLayout(Layout);
			TempLabel = new Label("Sequence Filters");
			Layout.setConstraints(TempLabel, NormalConstraints);
			NormalConstraints.gridwidth = 1;
			Layout.setConstraints(Sic, NormalConstraints);
			Layout.setConstraints(Sip, NormalConstraints);
			Layout.setConstraints(Sec, NormalConstraints);
			Layout.setConstraints(Stn, NormalConstraints);
			Layout.setConstraints(Stc, NormalConstraints);
			Layout.setConstraints(X1, NormalConstraints);
			Layout.setConstraints(X2, NormalConstraints);
			Layout.setConstraints(X3, NormalConstraints);
			Layout.setConstraints(X4, NormalConstraints);
			add(TempLabel);
			add(Sic);
			add(Sip);
			add(Sec);
			add(Stn);
			add(Stc);
			add(X1);
			add(X2);
			add(X3);
			add(X4);
			NormalConstraints.gridx = 1;
			TempLabel = new Label(
					"Must include all these residues before C terminus");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Must include this pattern");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label(
					"Must exclude all these residues before C terminus");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Must be preceded by one of these residues");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label(
					"Must have one of these residues at C terminus");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Maximum +1 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Maximum +2 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Maximum +3 XCorr");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			this.setTitle("DTASelect Extended Filters");
			this.setSize(450, 350);
			this.addWindowListener(new SelectCriteria.WindowEventAdapter());
			this.show();
		}
	}

	class LocusGUI extends Frame {
		private TextField r = new TextField("100");
		private TextField p = new TextField(
				new Integer(MinPepsPerLocus).toString());
		private Checkbox u = new Checkbox(
				"Include only loci with uniquely matching peptides");
		private CheckboxGroup t = new CheckboxGroup();
		private Checkbox t0 = new Checkbox(
				"Do not purge duplicate spectra for each sequence", t, false);
		private Checkbox t1 = new Checkbox(
				"Purge duplicate spectra in each salt step", t, false);
		private Checkbox t2 = new Checkbox(
				"Purge duplicate spectra on basis of XCorr", t, true);
		private Checkbox compress = new Checkbox(
				"Compress .dta files to create IDX and SPM files", false);
		private Checkbox copy = new Checkbox("Create copy script", false);
		private Checkbox GUI = new Checkbox("Show GUI of results", true);
		private Checkbox dot = new Checkbox(
				"Run on .dta files in current directory only", false);

		public LocusGUI() {
			GridBagConstraints NormalConstraints = new GridBagConstraints();
			GridBagLayout Layout = new GridBagLayout();
			Label TempLabel = new Label("Locus Filters");
			switch (PurgeDuplicateSequences) {
			case 0:
				t.setSelectedCheckbox(t0);
				break;
			case 1:
				t.setSelectedCheckbox(t1);
				break;
			case 2:
				t.setSelectedCheckbox(t2);
				break;
			}
			this.setLayout(Layout);
			NormalConstraints.fill = GridBagConstraints.BOTH;
			NormalConstraints.gridx = 0;
			NormalConstraints.gridwidth = 2;
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			NormalConstraints.gridwidth = 1;
			Layout.setConstraints(r, NormalConstraints);
			Layout.setConstraints(p, NormalConstraints);
			NormalConstraints.gridx = 1;
			TempLabel = new Label(
					"Show all loci with peptides appearing this many times");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			TempLabel = new Label("Set minimum peptides per locus criterion");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			NormalConstraints.gridx = 0;
			NormalConstraints.gridwidth = 2;
			Layout.setConstraints(u, NormalConstraints);
			Layout.setConstraints(t0, NormalConstraints);
			Layout.setConstraints(t1, NormalConstraints);
			Layout.setConstraints(t2, NormalConstraints);
			add(r);
			add(p);
			add(u);
			add(t0);
			add(t1);
			add(t2);
			TempLabel = new Label("Utilities");
			Layout.setConstraints(TempLabel, NormalConstraints);
			add(TempLabel);
			Layout.setConstraints(compress, NormalConstraints);
			Layout.setConstraints(copy, NormalConstraints);
			Layout.setConstraints(GUI, NormalConstraints);
			Layout.setConstraints(dot, NormalConstraints);
			add(compress);
			add(copy);
			add(GUI);
			add(dot);
			this.setTitle("DTASelect Locus Filters and Utilities");
			this.setSize(400, 350);
			this.addWindowListener(new SelectCriteria.WindowEventAdapter());
			this.show();
		}
	}

	class WindowEventAdapter extends WindowAdapter {
		public void windowClosing(WindowEvent WhatHappened) {
			System.exit(0);
		}
	}
}
