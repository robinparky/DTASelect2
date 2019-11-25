import java.io.*;

/* This class stores information about a class of proteins,
 * whether grouped by physical similarity, functional similarity,
 * or random choice. */
public class Classification {
    byte           Identifier = 127;
    String         Descriptor;
    int            NRCountObserved = 0;
    int            RCountObserved = 0;
    Classification Next;

    public Classification GetClass(byte Code) {
	Classification         CRunner = this.Next;
	while (CRunner != null && CRunner.Identifier != Code) {
	    CRunner = CRunner.Next;
	}
	return CRunner;
    }

    public void PrintDatabase(String RootName) {
	Classification          CRunner = this.Next;
	File                    CurrentDirectory = new File(System.getProperty("user.dir"));
	File		    OutputFile;
	FileWriter		    OutputFileWriter;
	BufferedWriter	    Outgoing;
	String                  tab = "\t";
	try {
	    CurrentDirectory = new File(CurrentDirectory.getCanonicalPath());
	    OutputFile = new File(CurrentDirectory, RootName + "-Classes.txt");
	    OutputFileWriter = new FileWriter(OutputFile);
	    Outgoing = new BufferedWriter(OutputFileWriter);
	    Outgoing.write("ClassID\tDescription\tNRCount\tRedundantCount\n");
	    while (CRunner != null) {
		Outgoing.write(CRunner.Identifier + tab +
			       CRunner.Descriptor + tab +
			       CRunner.NRCountObserved + tab +
			       CRunner.RCountObserved + "\n");
		CRunner = CRunner.Next;
	    }
	    Outgoing.flush();
	    Outgoing.close();
	}
	catch (IOException failure) {
	    System.out.println("Failed to write database files.");
	    System.out.println(failure);
	}
    }
}
