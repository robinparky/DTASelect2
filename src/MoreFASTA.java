import java.io.*;
import java.util.*;

/* MoreFASTA
   David L. Tabb
   September 26, 2001

   * Creates database- and spreadsheet-compatible compilations of
   * FASTA databases.  Fields reported include locus name, length,
   * molecular weight, pI, and description string.
   */

public class MoreFASTA {

    public static void usage() {
	System.out.println("MoreFASTA");
	System.out.println("Creates summary table for a FASTA sequence database");
	System.out.println("Copyright 2001, David L. Tabb, Yates Lab, The Scripps Research Institute");
	System.out.println("usage: MoreFASTA filename1 filename2 filename3 ...");
    }

    public static void main(String args[]) throws IOException {
	int            Counter;
	MoreFASTA      App;
	//Look for help request
	for(Counter = 0; Counter < args.length; Counter++) {
	    if (args[Counter].equals("--help") || 
		args[Counter].equals("-h") ||
		args[Counter].equals("/?")) {
		usage();
		System.exit(0);
	    }
	}
	//Make CSV files for each enumerated database
	for(Counter= 0; Counter < args.length; Counter++) {
	    App = new MoreFASTA(args[Counter]);
	}
    }

    public MoreFASTA(String FileName) throws IOException {
	File                  DB = new File(FileName);
	FileReader            InputFileReader;
	BufferedReader        Incoming;
	String                LineBuffer="";
	StringTokenizer       Parser;
	StringBuffer          SequenceSB;
	String                Sequence;
	Protein               CurrentProtein = new Protein();
	File                    CurrentDirectory = new File(System.getProperty("user.dir"));
	File			OutputFile;
	FileWriter		OutputFileWriter;
	BufferedWriter		Outgoing;
	int                     Counter = 0;
	int                     ResCounts[] = new int[26];
	int                     LocCounts[] = new int[26];
	int                     Looper;
	int                     SeqLength;

	try {
	    System.out.println("Now processing " + DB.getAbsolutePath());
	    if (DB.canRead()) {
		OutputFile = new File(FileName + ".DB");
		OutputFileWriter = new FileWriter(OutputFile);
		Outgoing = new BufferedWriter(OutputFileWriter);
		InputFileReader = new FileReader(DB);
		Incoming = new BufferedReader(InputFileReader);
		LineBuffer = Incoming.readLine();
		Outgoing.write("LocusID\tLength\tMW\tpI\tDescription");
		for (Looper = 0; Looper < 26; Looper++) {
		    Outgoing.write("\t" + new Character((char)('A' + Looper)).toString());
		}
		Outgoing.write("\n");
		//Queue up the DB file to the first locus line
		while ( (LineBuffer.length() == 0) || (LineBuffer.charAt(0) != '>') )
		    LineBuffer = Incoming.readLine();
		// Until we hit the end of the file
		while ( (LineBuffer != null) ) {
		    Counter++;
		    Parser = new StringTokenizer(LineBuffer, " \t>");
		    if (!Parser.hasMoreTokens()) {
			LineBuffer = Incoming.readLine();
		    }
		    else {
			for (Looper = 0; Looper < 26; Looper++) {
			    LocCounts[Looper] = 0;
			}
			LineBuffer = Parser.nextToken();
			/* Deal with the fact that SEQUEST .output files
			 * list only first 40 characters of locus
			 * names. */
			if (LineBuffer.length() > 40) {
			    LineBuffer = LineBuffer.substring(0,40);
			}
			CurrentProtein.Locus = LineBuffer;
			/* Read gene name and sequence */
			CurrentProtein.Gene = "";
			while (Parser.hasMoreTokens()) {
			    CurrentProtein.Gene += Parser.nextToken() + " ";
			}
			SequenceSB = new StringBuffer();
			LineBuffer = Incoming.readLine();
			while ( (LineBuffer != null) && (
							 (LineBuffer.length() == 0) ||
							 (LineBuffer.charAt(0) != '>') ) ) {
			    SequenceSB.append(LineBuffer);
			    LineBuffer = Incoming.readLine();
			}
			// Strip out extraneous characters from sequence
			Sequence = DTAFile.JustLettersFrom(SequenceSB.toString());
			SeqLength = Sequence.length();
			for (Looper = 0; Looper < SeqLength; Looper++) {
			    ResCounts[Sequence.charAt(Looper) - 'A']++;
			    LocCounts[Sequence.charAt(Looper) - 'A']++;
			}
			CurrentProtein.SequenceLength = Sequence.length();
			/* Determine protein's sequence coverage and calculate
			 * molecular weight and pI.
			 */
			CurrentProtein.CalculateMWAndpI(Sequence);
			Outgoing.write(CurrentProtein.Locus + "\t" +
				       new Integer(CurrentProtein.SequenceLength) + "\t" +
				       new Float(CurrentProtein.MolWt).toString() + "\t" +
				       new Float(CurrentProtein.pI).toString() + "\t" +
				       CurrentProtein.Gene);
			for (Looper = 0; Looper < 26; Looper++) {
			    Outgoing.write("\t" + new Integer(LocCounts[Looper]).toString());
			}
			Outgoing.write("\n");
		    }
		}
		System.out.println("\tProcessed " + new
		    Integer(Counter).toString() + " proteins");
		Outgoing.write("\n\t\t\t\t");
		for (Looper = 0; Looper < 26; Looper++) {
		    Outgoing.write("\t" + new Integer(ResCounts[Looper]).toString());
		}
		Outgoing.flush();
		Outgoing.close();
	    }
	    else {
		System.out.println("Could not read " + FileName);
	    }
	}
	catch (IOException failure) {
	    System.out.println("Failed while summarizing database.");
	    System.out.println(failure);
	}
    }
}
