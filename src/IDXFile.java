import java.io.*;

/*  IDXFile Object
    For DTASelect
    Created by Dave Tabb
    Begun 10/03/2000

    * IDX File format:

    * Directory lines are followed by 0 or more Spectrum lines.

    * Directory lines consist of the following: Directory Name Length
    * (byte), Directory Name (byte[]), Filepointer to next directory
    * line (long), DTA count for this directory (short), Minimum low
    * scan number (short), Maximum low scan number (short).

    * Spectrum lines consist of the following: Low scan number
    * (short), High scan number (short), M/Z of precursor ion (short *
    * 10), Filepointer to spectrum location within SPM file (long)
    */

public class IDXFile {
    Subdir          Data = new Subdir();
    Subdir          CurrentDirectory = Data;
    final long      BytesPerSubDir = 15;
    final long      BytesPerDTAinIDX = 14;

    public void WriteToDisk() throws IOException {
	String              UserDirName = System.getProperty("user.dir");
	RandomAccessFile    RAHandle;
	File                UserDir;
	int                 DTACount;
	long                FilePos;
	UserDir = new File(UserDirName);
	RAHandle = new RandomAccessFile(new File(UserDir, "DTASelect.IDX"), "rw");
	RAHandle.writeByte(Protein.Version().length());
	RAHandle.writeBytes(Protein.Version());
	CurrentDirectory = Data.Next;
	while (CurrentDirectory != null) {
	    CurrentDirectory.DRunner = CurrentDirectory.DTAs;
	    DTACount = 0;
	    while (CurrentDirectory.DRunner.Next != null) {
		CurrentDirectory.DRunner = CurrentDirectory.DRunner.Next;
		DTACount++;
	    }
	    FilePos = RAHandle.getFilePointer();
	    RAHandle.writeByte(CurrentDirectory.Name.length());
	    RAHandle.writeBytes(CurrentDirectory.Name);
	    RAHandle.writeLong(FilePos + BytesPerSubDir + CurrentDirectory.Name.length() +
			       DTACount * BytesPerDTAinIDX);
	    RAHandle.writeShort(DTACount);
	    if (DTACount == 0) {
		RAHandle.writeShort(0);
		RAHandle.writeShort(0);
	    }
	    else {
		RAHandle.writeShort(CurrentDirectory.DTAs.Next.LowScan);
		RAHandle.writeShort(CurrentDirectory.DRunner.LowScan);
	    }
	    CurrentDirectory.DRunner = CurrentDirectory.DTAs.Next;
	    while (CurrentDirectory.DRunner != null) {
		RAHandle.writeShort(CurrentDirectory.DRunner.LowScan);
		RAHandle.writeShort(CurrentDirectory.DRunner.HighScan);
		RAHandle.writeShort(CurrentDirectory.DRunner.MOverZ);
		RAHandle.writeLong(CurrentDirectory.DRunner.SPMPosition);
		CurrentDirectory.DRunner = CurrentDirectory.DRunner.Next;
	    }
	    CurrentDirectory = CurrentDirectory.Next;
	}
    }

    public void PrintDir() throws IOException {
	String            UserDirName = System.getProperty("user.dir");
	RandomAccessFile  RAHandle;
	SPMFile           Spectra = new SPMFile();
	File              UserDir;
	long              FilePos;
	long              FPBuffer;
	int               Length;
	byte              Buffer[];
	char              CharBuffer;
	Spectrum          RetrievedSpectrum;
	Point             PRunner;
	UserDir = new File(UserDirName);
	RAHandle = new RandomAccessFile(new File(UserDir, "DTASelect.IDX"), "r");
	//Print header of file (version info)
	RAHandle.seek(0);
	Buffer = new byte[RAHandle.readByte()];
	RAHandle.read(Buffer);
	System.out.println(new String(Buffer));
	FilePos = RAHandle.getFilePointer();
        while (FilePos < RAHandle.length()) {
	    RAHandle.seek(FilePos);
	    Length = RAHandle.readByte();
	    Buffer = new byte[Length];
	    RAHandle.read(Buffer);
	    /* Fields listed for each subdirectory include:
	       Directory name
	       File pointer to next subdirectory
	       DTA file count for this subdirectory
	       Low scan number
	       High scan number
	    */	       
	    System.out.print(new String(Buffer) + "\t");
	    FilePos = RAHandle.readLong();
	    System.out.print(FilePos);
	    System.out.print("\t" + new Short(RAHandle.readShort()).toString());
	    System.out.print("\t" + new Short(RAHandle.readShort()).toString() + "\t");
	    System.out.println(RAHandle.readShort());
	    while (RAHandle.getFilePointer() < FilePos) {
		/* Fields listed for each spectrum include:
		   Low Scan Number
		   High Scan Number
		   M/Z of precursor ion ( x 10)
		   Position within SPM file
		*/
		System.out.print("\t" + new Short(RAHandle.readShort()).toString()
		    + "\t" + new Short(RAHandle.readShort()).toString()
			+ "\t" + new Short(RAHandle.readShort()).toString());
		FPBuffer = RAHandle.readLong();
		System.out.println("\t" + new Long(FPBuffer).toString());
		try {
		    RetrievedSpectrum = Spectra.GetSpectrumAt(FPBuffer);
		}
		catch (IOException failure) {
		    System.out.println("Failure while reading at " +
				       new Long(FPBuffer).toString() + ".");
		    System.out.println(failure);
		    System.out.println(failure.getMessage());
		    System.out.println(failure.getLocalizedMessage());
		    System.out.println(failure.toString());
		    System.exit(0);
		}
	    }
	}
	RAHandle.close();
    }

    public Spectrum GetSpectrum(String Subdir, short TargetLow, short TargetHigh) {
	String            UserDirName = System.getProperty("user.dir");
	RandomAccessFile  RAHandle;
	SPMFile           Spectra = new SPMFile();
	File              UserDir;
	long              FilePos;
	long              FPBuffer;
	long              StartOfScans;
	long              EndOfScans;
	int               Length;
	byte              Buffer[];
	char              CharBuffer;
	Point             PRunner;
	boolean           Match = false;
	short             DTACount = 0;
	short             LowScan = 0;
	short             HighScan = 0;
	Spectrum          Returned = null;
	float             PrecursorMOverZ;
	int               SearchLength = 0;
	UserDir = new File(UserDirName);
	try {
	    RAHandle = new RandomAccessFile(new File(UserDir, "DTASelect.IDX"), "r");
	    //Skip header of file (version info)
	    RAHandle.seek(0);
	    Buffer = new byte[RAHandle.readByte()];
	    RAHandle.read(Buffer);
	    FilePos = RAHandle.getFilePointer();
	    while ( (FilePos < RAHandle.length()) && (!Match) ) {
		//Seek through IDX file for subdirectory
		RAHandle.seek(FilePos);
		Length = RAHandle.readByte();
		Buffer = new byte[Length];
		RAHandle.read(Buffer);
		FilePos = RAHandle.readLong();
		if (new String(Buffer).equals(Subdir)) {
		    Match = true;
		    DTACount = RAHandle.readShort();
		    LowScan = RAHandle.readShort();
		    HighScan = RAHandle.readShort();
		}
	    }
	    if (Match) {
		float         ProjectedSkip = new Float(TargetLow - LowScan).floatValue() /
		    new Float(HighScan - LowScan).floatValue();
		long          FileOffset;
		short         CurrentLowScan;
		short         CurrentHighScan;
		if (ProjectedSkip > 1.0f)
		    ProjectedSkip = 1.0f;
		if (ProjectedSkip < 0f)
		    ProjectedSkip = 0f;
		EndOfScans = FilePos;
		StartOfScans = RAHandle.getFilePointer();
		FileOffset = new Float(ProjectedSkip * (DTACount-1)).longValue();
		FileOffset *= BytesPerDTAinIDX;
		FileOffset += RAHandle.getFilePointer();
		//Move to correct scan
		RAHandle.seek(FileOffset);
		CurrentLowScan = RAHandle.readShort();
		CurrentHighScan = RAHandle.readShort();
		SearchLength++;
		/* System.out.println("\n" + new Short(CurrentLowScan).toString() + "\t" +
		   new Short(CurrentHighScan).toString()); */
		while ( ( (CurrentLowScan < TargetLow) || (CurrentHighScan < TargetHigh) ) &&
			FileOffset < EndOfScans ) {
		    FileOffset += BytesPerDTAinIDX;
		    if (FileOffset < EndOfScans) {
			RAHandle.seek(FileOffset);
			CurrentLowScan = RAHandle.readShort();
			CurrentHighScan = RAHandle.readShort();
			/* System.out.println(new Short(CurrentLowScan).toString() + "\t" +
			   new Short(CurrentHighScan).toString()); */
			SearchLength++;
		    }
		}
		while ( ( (CurrentLowScan > TargetLow) || (CurrentHighScan > TargetHigh) ) &&
			FileOffset >= StartOfScans ) {
		    FileOffset -= BytesPerDTAinIDX;
		    if (FileOffset >= StartOfScans) {
			RAHandle.seek(FileOffset);
			CurrentLowScan = RAHandle.readShort();
			CurrentHighScan = RAHandle.readShort();
			/* System.out.println(new Short(CurrentLowScan).toString() + "\t" +
			   new Short(CurrentHighScan).toString()); */
			SearchLength++;
		    }
		}
		if ( (CurrentLowScan == TargetLow) && (CurrentHighScan == TargetHigh) ) {
		    /* System.out.println(new Float(ProjectedSkip).toString() + "\t" +
		       new Integer(SearchLength).toString()); */
		    PrecursorMOverZ = new Float(RAHandle.readShort()).floatValue() / 10f;
		    Returned = Spectra.GetSpectrumAt(RAHandle.readLong());
		    Returned.PrecursorMOverZ = PrecursorMOverZ;
		}
		else {
		    System.out.println("IDX: Spectrum not found at expected place " +
				       new Float(ProjectedSkip).toString() + "\t" +
				       new Short(TargetLow).toString() + "\t" +
				       new Short(TargetHigh).toString());
		}
	    }
	    else {
		System.out.println("IDX: Could not find subdirectory " + Subdir);
	    }
	    RAHandle.close();
	}
	catch (IOException failure) {
	    System.out.println("IDX: Failure while decompressing spectrum " +
			       new Short(TargetLow).toString() + " " +
			       new Short(TargetHigh).toString());
	}
	return Returned;
    }

    public void AddSubdir(String Name) {
	CurrentDirectory.Next = new Subdir();
	CurrentDirectory = CurrentDirectory.Next;
	CurrentDirectory.Name = Name;
    }

    public void AddDTAToIndex(short LowScan, short HighScan, short MOverZ, long Pointer) {
	CurrentDirectory.AddToIndex(LowScan, HighScan, MOverZ, Pointer);
    }

    class Subdir {
	String        Name;
	short         LowScan;
	short         HighScan;
	Subdir        Next;
	DTA           DTAs = new DTA();
	DTA           DRunner = DTAs;

	public void AddToIndex(short LowScan, short HighScan, short MOverZ, long Pointer) {
	    DRunner.Next = new DTA();
	    DRunner = DRunner.Next;
	    DRunner.LowScan = LowScan;
	    DRunner.HighScan = HighScan;
	    DRunner.MOverZ = MOverZ;
	    DRunner.SPMPosition = Pointer;
	}
    }

    class DTA {
	short         LowScan;
	short         HighScan;
	short         MOverZ;
	DTA           Next;
	long          SPMPosition;
    }
}
