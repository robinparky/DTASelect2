import java.io.*;
import java.util.*;

/* DTASelect IniFile class
 Dave Tabb, Feb 21, 2001
 */

public class IniFile {
	String DTASpecDisplay = "http://localhost/cgi-shl/web_display.exe";
	String DTAOUTDisplay = "http://localhost/cgi-shl/web_showoutput.exe";
	String SQTDisplay = "http://localhost/cgi-shl/show.pl";
	String BasicSeqCov = "http://localhost/cgi-shl/web_retrieve.exe";
	String DepthSeqCov = "http://localhost/cgi-shl/SeqCov.pl";
	boolean UseDepthCGI = false;
	String BLAST = "http://www.ncbi.nlm.nih.gov:80/blast/Blast.cgi?PROGRAM=blastp&QUERY=";
	String BLASTArgs = "";
	String ProtValidation = "http://localhost/cgi-shl/EvalocusA";
	String ServerType = "Windows";
	String MascotPath = "c:\\inetpub\\mascot";
	int LocusLengthCutoff = 40;

	public void Initialize() {
		try {
			String ClassPath = System.getProperty("java.class.path", ".");
			StringTokenizer Parser = new StringTokenizer(ClassPath,
					System.getProperty("path.separator"));
			boolean FoundIni = false;
			File CurrentDirectory;
			File IniFile = null;
			String LineBuffer;
			FileReader InputFileReader;
			BufferedReader Incoming;

			while (Parser.hasMoreTokens() && !FoundIni) {

				// TODO

				 CurrentDirectory = new File(Parser.nextToken()); //for
				// deployment //TODO
			//	CurrentDirectory = new File("/home/diego/Desktop/jolene/");// local
																		// debug

				// TODO

				IniFile = new File(CurrentDirectory, "DTASelect.ini");
				FoundIni = IniFile.exists();
			}
			if ((FoundIni) && (IniFile.canRead())) {
				InputFileReader = new FileReader(IniFile);
				Incoming = new BufferedReader(InputFileReader);
				LineBuffer = Incoming.readLine();
				while (LineBuffer != null) {
					if (LineBuffer.startsWith("#")) {
						// Don't do anything on comments
					} else if (LineBuffer.length() == 0) {
						// Don't do anything on blank lines
					} else if (LineBuffer.startsWith("DTASpecDisplay")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						DTASpecDisplay = Parser.nextToken();
					} else if (LineBuffer.startsWith("DTAOUTDisplay")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						DTAOUTDisplay = Parser.nextToken();
					} else if (LineBuffer.startsWith("SQTDisplay")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						SQTDisplay = Parser.nextToken();
					} else if (LineBuffer.startsWith("BasicSeqCov")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						BasicSeqCov = Parser.nextToken();
					} else if (LineBuffer.startsWith("DepthSeqCov")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						DepthSeqCov = Parser.nextToken();
					} else if (LineBuffer.startsWith("UseDepthCGI")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						UseDepthCGI = new Boolean(Parser.nextToken())
								.booleanValue();
					} else if (LineBuffer.startsWith("ProtValidation")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						if (Parser.hasMoreTokens())
							ProtValidation = Parser.nextToken();
					} else if (LineBuffer.startsWith("BlastArgs")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						if (Parser.hasMoreTokens())
							BLASTArgs = Parser.nextToken();
					} else if (LineBuffer.startsWith("Blast")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						BLAST = Parser.nextToken();
					} else if (LineBuffer.startsWith("ServerType")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						ServerType = Parser.nextToken();
					} else if (LineBuffer.startsWith("LocusLengthCutoff")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						LocusLengthCutoff = new Integer(Parser.nextToken())
								.intValue();
					} else if (LineBuffer.startsWith("MascotPath")) {
						Parser = new StringTokenizer(LineBuffer);
						Parser.nextToken();
						if (Parser.hasMoreTokens())
							MascotPath = Parser.nextToken();
						while (Parser.hasMoreTokens()) {
							MascotPath = MascotPath + " " + Parser.nextToken();
						}
					} else {
						System.out
								.println("Didn't understand this line in DTASelect.ini:");
						System.out.println(LineBuffer);
						System.out.println();
					}
					LineBuffer = Incoming.readLine();
				}
			} else {
				if (FoundIni) {
					System.out.println("Couldn't read DTASelect.ini.");
				} else {
					System.out.println("Could not find DTASelect.ini.");
					System.out
							.println("Your classpath must include the directory where DTASelect is installed.");
				}
				// System.exit(0);
			}
		} catch (IOException failure) {
			System.out.println("Error reading DTASelect.ini.\n");
			System.exit(0);
		}
	}

}
