import java.io.*;
import java.util.*;

/* SPMFile Object
 * For DTASelect
 * Created by Dave Tabb
 * Begun 10/03/2000

 * SPM File format:

 * File consists of spectrum records back-to-back.  Each spectrum
 * record consists of the following fields: Point count (short),
 * Intensity byte count (short), Maximum intensity (float), Intensity
 * huffman code (byte[]), M/Z stream(short[] * 10)
 */

public class SPMFile {
    DTA                      SpectraToAdd = new DTA();
    DTA                      STARunner = SpectraToAdd;
    Huffman                  Coder = new Huffman();
    long                     FilePos = 0;
    short                    LastAddedMOverZ = 0;

    public SPMFile() {
	Coder.Initialize();
    }

    public long AddSpectrumToQueue(File Subdir, String FileName) throws IOException {
	File                   SpectrumFile = new File(Subdir, FileName);
	BufferedReader         Incoming = new BufferedReader(new FileReader(SpectrumFile));
	String                 LineBuffer = Incoming.readLine();
	StringTokenizer        Parser = new StringTokenizer(LineBuffer);
	Point                  Points = new Point();
	Point                  PRunner = Points;
	int                    Looper;
	long                   Temp;
	float                  M;
	float                  Z;
	STARunner.Next = new DTA();
	STARunner = STARunner.Next;
	STARunner.PointCount = 0;
	STARunner.MaxInt = 0;
	M = new Float(Parser.nextToken()).floatValue();
	Z = new Float(Parser.nextToken()).floatValue();
	LastAddedMOverZ = new Integer(Math.round((M-1)/Z*10f)).shortValue();
	LineBuffer = Incoming.readLine();
	while (LineBuffer != null) {
	    Parser = new StringTokenizer(LineBuffer);
	    STARunner.PointCount++;
	    PRunner.Next = new Point();
	    PRunner = PRunner.Next;
	    PRunner.MOverZ = new Float(Parser.nextToken()).floatValue();
	    PRunner.Intensity = new Float(Parser.nextToken()).floatValue();
	    if (PRunner.Intensity > STARunner.MaxInt)
		STARunner.MaxInt = PRunner.Intensity;
	    LineBuffer = Incoming.readLine();
	}
	Incoming.close();
	/* System.out.println("\n" + FileName + "\t" +
	   new Float(STARunner.MaxInt).toString() + "\t" +
	   new Integer(STARunner.PointCount).toString()); */
	//Compress Intensity and M/Z information
	STARunner.IntensityStream = Coder.EncodeIntensities(Points,STARunner.MaxInt);
	STARunner.ByteCount = new Integer(STARunner.IntensityStream.length()).shortValue();
	STARunner.MOverZStream = new short[STARunner.PointCount];
	PRunner = Points.Next;
	for (Looper = 0; Looper < STARunner.PointCount; Looper++) {
	    STARunner.MOverZStream[Looper] =
		new Double(Math.round(PRunner.MOverZ*10f)).shortValue();
	    PRunner = PRunner.Next;
	}
	Temp = STARunner.RecordSize();
	FilePos += Temp;
	return (FilePos - Temp);
    }

    public void AddSpectraToFile(RandomAccessFile RAHandle) throws IOException {
	int    Counter;
	STARunner = SpectraToAdd.Next;
	while(STARunner != null) {
	    RAHandle.writeShort(STARunner.PointCount);
	    RAHandle.writeShort(STARunner.ByteCount);
	    RAHandle.writeFloat(STARunner.MaxInt);
	    RAHandle.writeBytes(STARunner.IntensityStream);
	    for (Counter = 0; Counter < STARunner.PointCount; Counter++)
		RAHandle.writeShort(STARunner.MOverZStream[Counter]);
	    STARunner = STARunner.Next;
	}
	SpectraToAdd = new DTA();
	STARunner = SpectraToAdd;
    }

    public Spectrum GetSpectrumAt(long FilePos) throws IOException{
	String              UserDirName = System.getProperty("user.dir");
	RandomAccessFile    RAHandle;
	File                UserDir;
	Point               PRunner;
	short               PointCount;
	short               IntByteCount;
	float               MaxInt;
	float               MaxMOverZ = 0f;
	byte                Buffer[];
	short               Looper;
	Spectrum            Returned = new Spectrum();
	UserDir = new File(UserDirName);
	RAHandle = new RandomAccessFile(new File(UserDir, "DTASelect.SPM"), "r");
	RAHandle.seek(FilePos);
	PointCount = RAHandle.readShort();
	IntByteCount = RAHandle.readShort();
	MaxInt = RAHandle.readFloat();
	Buffer = new byte[IntByteCount];
	RAHandle.read(Buffer);
	Returned.Points = Coder.DecodeIntensities(PointCount, Buffer, MaxInt);
	PRunner = Returned.Points.Next;
	for (Looper = 0; Looper < PointCount; Looper++) {
	    MaxMOverZ = new Short(RAHandle.readShort()).floatValue() / 10.0f;
	    //System.out.println(MaxMOverZ);
	    PRunner.MOverZ = MaxMOverZ;
	    PRunner = PRunner.Next;
	}
	Returned.LowMOverZ = Returned.Points.Next.MOverZ;
	Returned.HighMOverZ = MaxMOverZ;
	Returned.MaxIntensity = MaxInt;
	RAHandle.close();
	return Returned;
    }

    class DTA {
	short           PointCount;
	short           ByteCount;
	float           MaxInt;
	String          IntensityStream;
	short           MOverZStream[];
	DTA             Next;

	public long RecordSize() {
	    return (2 + 2 + 4 + ByteCount + (PointCount * 2));
	}
    }
}
