import java.io.*;
import java.util.*;

/**   Huffman class
      David Tabb, September 29, 2000
*/

public class Huffman {
    Integer       Intensity = null;
    float         Percentage;
    Huffman       LChild = null;
    Huffman       RChild = null;
    Huffman       Next = null;
    String        Path = "";

    public void Initialize() {
	int              TreeCount = 0;
	Huffman          Buffer;
	Huffman          Duplicate = new Huffman();
	Huffman          DRunner = Duplicate;
	try {
	    String           ClassPath = System.getProperty("java.class.path",".");
	    StringTokenizer  Path = new StringTokenizer(ClassPath,
							System.getProperty("path.separator"));
	    boolean          FoundIni = false;
	    File             CurrentDirectory;
	    File             IniFile = null;
	    FileReader       InputFileReader;
	    BufferedReader   Incoming;
	    String           LineBuffer;
	    StringTokenizer  Parser;
	    while (Path.hasMoreTokens() && !FoundIni) {
		CurrentDirectory = new File(Path.nextToken());
		IniFile = new File(CurrentDirectory,"Huffman.ini");
		FoundIni = IniFile.exists();
	    }
	    if ( (FoundIni) && (IniFile.canRead()) ){
		//CurrentDirectory = new File(Path.nextToken());
		InputFileReader = new FileReader(IniFile);
		Incoming = new BufferedReader(InputFileReader);
		LineBuffer = Incoming.readLine();
		while (LineBuffer != null) {
		    Parser = new StringTokenizer(LineBuffer);
		    if (Parser.hasMoreTokens()) {
			Buffer = new Huffman();
			TreeCount++;
			Buffer.Intensity = new Integer(Parser.nextToken());
			Buffer.Percentage = new Float(Parser.nextToken()).floatValue();
			Buffer.Next = this.Next;
			this.Next = Buffer;
			// Add DRunner copy, too
			DRunner.Next = new Huffman();
			DRunner = DRunner.Next;
			DRunner.Intensity = Buffer.Intensity;
		    }
		    LineBuffer = Incoming.readLine();
		}
	    }
	    else {
		System.out.println("Could not find Huffman.ini.");
		System.out.println("Your classpath must include the directory where DTASelect is installed.");
	    }
	}
	catch (IOException failure) {
	    System.out.println("Error reading Huffman.ini.\n");
	    System.out.println("Does your classpath include only the DTASelect directory?");
	    System.exit(0);
	}
	// Build Tree
	Huffman      Lowest;
	Huffman      NextLowest;
	while (TreeCount > 1) {
	    Lowest = new Huffman();
	    NextLowest = new Huffman();
	    Lowest.Percentage = 100;
	    NextLowest.Percentage = 100;
	    //Find the lowest two percentages in the tree
	    Buffer = this.Next;
	    while (Buffer != null) {
		if (Buffer.Percentage < Lowest.Percentage) {
		    NextLowest = Lowest;
		    Lowest = Buffer;
		}
		else if (Buffer.Percentage < NextLowest.Percentage) {
		    NextLowest = Buffer;
		}
		Buffer = Buffer.Next;
	    }
	    //Create a new tree in the list for the two we're pulling out
	    Buffer = new Huffman();
	    Buffer.LChild = Lowest;
	    Buffer.RChild = NextLowest;
	    Buffer.Percentage = Buffer.LChild.Percentage + Buffer.RChild.Percentage;
	    Buffer.Next = this.Next;
	    this.Next = Buffer;
	    //Remove the two lowest from the list
	    Buffer = this;
	    while (Buffer.Next != null) {
		if ( (Buffer.Next == Lowest) || (Buffer.Next == NextLowest) ) {
		    Buffer.Next = Buffer.Next.Next;
		}
		else
		    Buffer = Buffer.Next;
	    }
	    //Reduce the TreeCount appropriately
	    TreeCount--;
	}
	this.LChild = this.Next.LChild;
	this.RChild = this.Next.RChild;
	this.Intensity = this.Next.Intensity;
	this.Percentage = this.Next.Percentage;
	this.Next = Duplicate.Next;
	// Set up encoder
	this.DepthFirstMatchup(this, "");
	//this.PrintDepth();
    }

    private void DepthFirstMatchup(Huffman Root, String PathSoFar) {
	if (this.Intensity != null) {
	    Huffman  Runner = Root.Next;
	    int      KeyStored = this.Intensity.intValue();
	    while (Runner.Intensity.intValue() != KeyStored)
		Runner = Runner.Next;
	    Runner.Path = PathSoFar;
	}
	else {
	    LChild.DepthFirstMatchup(Root, PathSoFar + "0");
	    RChild.DepthFirstMatchup(Root, PathSoFar + "1");
	}
    }

    public Point DecodeIntensities(int PointCount, byte[] Code, float
				   MaxIntensity) {
	Point           Points = new Point();
	int             CharCounter;
	int             TwoPower = 256;
	char            CurrentByte = (char)Code[0];
	//	System.out.print(new Integer(CurrentByte).toString() + "\t");
	Huffman         HRunner;
	Point           PRunner = Points;
	for (CharCounter = 0; CharCounter < PointCount; CharCounter++) {
	    PRunner.Next = new Point();
	    PRunner = PRunner.Next;
	    PRunner.MOverZ = new Integer(CharCounter).floatValue();
	}
	PRunner = Points.Next;
	CharCounter = 0;
	while (PRunner != null) {
	    /* Each pass through this loop corresponds to a discovered
	     * intensity for a peak in the spectrum.  The while loop
	     * below uses up the appropriate number of bits from the
	     * Huffman code sequence.  Once the intensity has been
	     * decoded from the stream by the internal while loop,
	     * that intensity is written to the list of points to be
	     * returned by this function.  */
	    HRunner = this;
	    while (HRunner.Intensity == null) {
		/* Each pass through this loop handles one more bit of
		 * the Huffman code, heading left on 0 and right on
		 * 1.  When a node of the tree is reached that
		 * corresponds to a leaf of the tree, the loop
		 * exits. */
		if (TwoPower == 1) {
		    TwoPower = 128;
		    CharCounter++;
		    CurrentByte = (char)Code[CharCounter];
		    //		    System.out.print(new Integer(CurrentByte).toString() + "\t");
		}
		else TwoPower /=2;
		if ( (TwoPower & CurrentByte) > 0) {
		    HRunner = HRunner.RChild;
		}
		else {
		    HRunner = HRunner.LChild;
		}
	    }
	    PRunner.Intensity = (HRunner.Intensity.intValue() *
				 MaxIntensity) / 100f;
	    /*	    System.out.println(HRunner.Intensity.toString() + "\t" + new
		    Float(PRunner.Intensity).toString()); */
	    PRunner = PRunner.Next;
	}
	//	System.out.println();
	return Points;
    }

    public String EncodeIntensities(Point Points, float MaxIntensity) {
	StringBuffer     Code = new StringBuffer();
	Point            PRunner = Points.Next;
	int              TwoPower = 256;
	char             CurrentByte = 0;
	Huffman          HRunner;
	int              NormInt;
	int              PathLooper;
	int              PathLength;
	while (PRunner != null) {
	    // Normalize this intensity
	    NormInt = new Double(Math.ceil(100f * PRunner.Intensity / MaxIntensity)).intValue();
	    if (NormInt > 100)
		NormInt = 100;
	    //System.out.print(new Integer(NormInt).toString() + "\t");
	    // Move HRunner to acquire correct path
	    HRunner = this;
	    while ( (HRunner.Intensity == null) || (HRunner.Intensity.intValue() != NormInt) ) {
		HRunner = HRunner.Next;
	    }
	    PathLength = HRunner.Path.length();
	    for (PathLooper = 0; PathLooper < PathLength; PathLooper++) {
		if (TwoPower == 1) {
		    TwoPower = 128;
		    Code.append(CurrentByte);
		    /*		    System.out.print(new
				    Integer(CurrentByte).toString() +
				    * "\t"); */
		    CurrentByte = 0;
		}
		else TwoPower /= 2;
		if (HRunner.Path.charAt(PathLooper) == '1')
		    CurrentByte += TwoPower;
	    }
	    PRunner = PRunner.Next;
	}
	Code.append(CurrentByte);
	/*		    System.out.print(new
			    Integer(CurrentByte).toString() + "\t");*/
	//	System.out.println();
	return Code.toString();
    }

    public void PrintDepth() {
	System.out.println("Root\t" + new Float(Percentage).toString());
	LChild.PrintDepthRecursive("0");
	RChild.PrintDepthRecursive("1");
    }

    private void PrintDepthRecursive(String Path) {
	if (Intensity != null)
	    System.out.println(new Float(Percentage).toString()
		+ "\t" + Intensity.toString()
			       + "\t" + new Integer(Path.length()).toString()
				   + "\t" + Path);
	else {
	    System.out.println(new Float(Percentage).toString() + "\t\t" + Path);
	    LChild.PrintDepthRecursive(Path + "0");
	    RChild.PrintDepthRecursive(Path + "1");
	}
    }
}
