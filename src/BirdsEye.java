import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

public class BirdsEye extends Frame {
    ProtGrid         Visualization = new ProtGrid();
    BirdsEyeIni      Config = new BirdsEyeIni();

    public BirdsEye(DataSet pData) {
	Visualization.Data = pData;
	Config.Initialize();
	Visualization.Config = Config;
	this.setSize(800,600);
	this.setTitle(Config.Title);
	addWindowListener (new WindowEventAdapter());
	this.setLayout(new BorderLayout());
	this.add("Center", Visualization);
	this.setVisible(true);
    }

    public class ProtGrid extends Canvas {
	DataSet                 Data;
	BirdsEyeIni             Config;
	// Characteristics of display
	int              RowHeight;
	int              DotHeight;
	int              ColWidth;
	int              DotWidth;
	int              HMargin = 5;
	int              VMargin = 5;
	int              Header = 0;
	int              PostClassSpace = 5;
	// Dynamic characteristics of display
	Dimension               Bounds;
	int                     DotsPerRow;
	int                     FontHeight;
	FontMetrics             CurrentFontInfo;

	public void paint(Graphics GraphContext) {
	    Protein               PRunner = Data.LocusList.Next;
	    int                   CurrentPosInGroup = 0;
	    int                   CurrentPosInRow = 2;
	    int                   CurrentRowInGroup = 0;
	    int                   TopOfGroup = 0;
	    int                   ProtCount = 0;
	    int                   ClassCount = 0;
	    RowHeight = Config.DotSize+1;
	    DotHeight = Config.DotSize;
	    ColWidth = Config.DotSize+1;
	    DotWidth = Config.DotSize;
	    // Find out how big the canvas is
	    Bounds = getSize();
	    if (Config.ReverseColors)
		GraphContext.setColor(Color.white);
	    else
		GraphContext.setColor(Color.black);
	    GraphContext.fillRect(0,0, Bounds.width, Bounds.height);
	    GraphContext.setFont(new Font("SanSerif",Font.BOLD,Config.TitleFont));
	    CurrentFontInfo = GraphContext.getFontMetrics();
	    FontHeight = CurrentFontInfo.getHeight();
	    DotsPerRow = (Bounds.width - HMargin * 2) / ColWidth;
	    // Put Title onscreen
	    if (Config.ReverseColors)
		GraphContext.setColor(Color.black);
	    else
		GraphContext.setColor(Color.white);
	    if (Config.ShowTitle == true) {
		Header = FontHeight;
		GraphContext.drawString(Config.Title,
					(Bounds.width - CurrentFontInfo.stringWidth(Config.Title))/2,
					Header + VMargin);
	    }
	    GraphContext.setFont(new Font("SanSerif",Font.BOLD,Config.OtherFont));
	    CurrentFontInfo = GraphContext.getFontMetrics();
	    FontHeight = CurrentFontInfo.getHeight();
	    if (PRunner.Classification == 127) {
		// These proteins have not been classified
		TopOfGroup = TopOfGroup + RowHeight + PostClassSpace;
		while (PRunner != null) {
		    GraphContext.setColor(GetColorFor(PRunner.SequenceCoverage));
		    ProtCount++;
		    GraphContext.fillRect(HMargin +
					  (CurrentPosInGroup % DotsPerRow) * ColWidth,
					  Header + VMargin + TopOfGroup +
					  (CurrentRowInGroup * RowHeight),
					  DotWidth, DotHeight);
		    if (CurrentPosInRow % 11 == 0 && CurrentPosInRow != DotsPerRow + 1) {
			CurrentPosInGroup++;
			CurrentPosInRow++;
		    }
		    CurrentPosInGroup++;
		    CurrentPosInRow++;
		    if (CurrentPosInGroup / DotsPerRow > CurrentRowInGroup) {
			CurrentRowInGroup = CurrentPosInGroup / DotsPerRow;
			CurrentPosInRow = 2;
		    }
		    PRunner = PRunner.Next;
		}
		TopOfGroup += CurrentRowInGroup * RowHeight;
	    }
	    else {
		// These proteins have been classified
		int              LastClass = 127;
		String           ClassCountString = "";
		Classification   CurrentClass;
		while (PRunner != null) {
		    if (PRunner.Classification != LastClass) {
			// Move the Current drawing position down appropriately
			TopOfGroup += RowHeight * CurrentRowInGroup;
			CurrentRowInGroup = 0;
			CurrentPosInGroup = 0;
			CurrentPosInRow = 2;
			// Time to put in a class label
			if (Config.ReverseColors)
			    GraphContext.setColor(Color.black);
			else
			    GraphContext.setColor(Color.white);
			CurrentClass = Data.Classifieds.GetClass(PRunner.Classification);
			LastClass = PRunner.Classification;
			ClassCount = 0;
			if (Config.ClassStats) {
			    Protein  PRunner2 = PRunner;
			    while (PRunner2 != null &&
				   PRunner2.Classification == PRunner.Classification) {
				PRunner2 = PRunner2.Next;
				ClassCount++;
			    }
			    ClassCountString = ": " + new Integer(ClassCount).toString();
			}
			if (CurrentClass == null) {
			    GraphContext.drawString("Unclassified" + ClassCountString,
						    HMargin,
						    Header + VMargin + TopOfGroup +
						    FontHeight + Config.DotSize);
			}
			else {
			    GraphContext.drawString(CurrentClass.Descriptor + ClassCountString,
						    HMargin,
						    Header + VMargin + TopOfGroup +
						    FontHeight + Config.DotSize);
			}
			TopOfGroup = TopOfGroup + FontHeight + RowHeight + PostClassSpace;
		    }
		    GraphContext.setColor(GetColorFor(PRunner.SequenceCoverage));
		    ProtCount++;
		    ClassCount++;
		    GraphContext.fillRect(HMargin +
					  (CurrentPosInGroup % DotsPerRow) * ColWidth,
					  Header + VMargin + TopOfGroup +
					  (CurrentRowInGroup * RowHeight),
					  DotWidth, DotHeight);
		    if (CurrentPosInRow % 11 == 0 && CurrentPosInRow != DotsPerRow + 1) {
			CurrentPosInGroup++;
			CurrentPosInRow++;
		    }
		    CurrentPosInGroup++;
		    CurrentPosInRow++;
		    if (CurrentPosInGroup / DotsPerRow > CurrentRowInGroup) {
			CurrentRowInGroup = CurrentPosInGroup / DotsPerRow;
			CurrentPosInRow = 2;
		    }
		    PRunner = PRunner.Next;
		}
	    }
	    TopOfGroup += FontHeight;
	    if (Config.OverallStats) {
		if (Config.ReverseColors)
		    GraphContext.setColor(Color.black);
		else
		    GraphContext.setColor(Color.white);
		GraphContext.drawString("Total Proteins: " +
					new Integer(ProtCount).toString(),
					HMargin,
					Header + VMargin + TopOfGroup +
					FontHeight + Config.DotSize);
		TopOfGroup = TopOfGroup + FontHeight + RowHeight + PostClassSpace;
	    }
	    if (Config.Legend) {
		float       Looper;
		if (Config.ReverseColors)
		    GraphContext.setColor(Color.black);
		else
		    GraphContext.setColor(Color.white);
		GraphContext.drawString("Sequence Coverage Legend:", HMargin,
					Header + VMargin + TopOfGroup + FontHeight +
					Config.DotSize);
		TopOfGroup = TopOfGroup + FontHeight + RowHeight + PostClassSpace * 2;
		CurrentPosInGroup = HMargin;
		// Abuse this variable name a bit...
		CurrentPosInRow = (Bounds.width - 2 * HMargin) / 5;
		for (Looper = 100f; Looper > -1; Looper -= 25f) {
		    GraphContext.setColor(GetColorFor(Looper));
		    GraphContext.fillRect(CurrentPosInGroup,
					  Header + VMargin + TopOfGroup,
					  DotWidth, DotHeight);
		    if (Config.ReverseColors)
			GraphContext.setColor(Color.black);
		    else
			GraphContext.setColor(Color.white);
		    GraphContext.drawString(new Float(Looper).toString() + "%",
					    CurrentPosInGroup + ColWidth,
					    Header + VMargin + TopOfGroup + RowHeight);
		    CurrentPosInGroup += CurrentPosInRow;
		    //TopOfGroup = TopOfGroup + FontHeight + PostClassSpace;
		}
	    }
	}

	private Color GetColorFor(float Coverage) {
	    return new Color(Color.HSBtoRGB(0.8f - Coverage *
					    0.008f, 1.0f,.7f));
	}
    }
    
    // Class to listen for and to handle window events
    class WindowEventAdapter extends WindowAdapter {
	public void windowClosing(WindowEvent WhatHappened) {
	    // BirdsEye.this.dispose();
	    System.exit(0);
	}
    }

    // Class to hold the config for this view
    class BirdsEyeIni {
	public String        Title = "Bird's Eye Proteome View";
	public int           DotSize = 8;
	public boolean       Legend = true;
	public boolean       OverallStats = true;
	public boolean       ClassStats = false;
	public boolean       ReverseColors = false;
	public boolean       ShowTitle = true;
	public int           TitleFont = 18;
	public int           OtherFont = 14;

	public void Initialize() {
	    try {
		File          BEIni = new File (new File(System.getProperty("user.dir")),
						"BirdsEye.ini");
		if (BEIni.exists()) {
		    FileReader      InputFileReader = new FileReader(BEIni);
		    BufferedReader  Incoming = new BufferedReader(InputFileReader);
		    String          LineBuffer;
		    String          WholeLine;
		    StringTokenizer Parser;
		    WholeLine = Incoming.readLine();
		    while (WholeLine != null) {
			Parser = new StringTokenizer(WholeLine);
			if (Parser.hasMoreTokens()) {
			    LineBuffer = Parser.nextToken();
			    if (LineBuffer.equals("Title")) {
				StringBuffer  temp = new StringBuffer(Parser.nextToken());
				while (Parser.hasMoreTokens()) {
				    temp.append(" ");
				    temp.append(Parser.nextToken());
				}
				this.Title = temp.toString();
			    }
			    else if (LineBuffer.startsWith("#")) {
				// Don't do anything on comments
			    }
			    else if (LineBuffer.equals("Dot-Size")) {
				this.DotSize = new Integer(Parser.nextToken()).intValue();
			    }
			    else if (LineBuffer.equals("Show-Title")) {
				this.ShowTitle = new Boolean(Parser.nextToken()).booleanValue();
			    }
			    else if (LineBuffer.equals("Legend")) {
				this.Legend = new Boolean(Parser.nextToken()).booleanValue();
			    }
			    else if (LineBuffer.equals("Overall-Stats")) {
				this.OverallStats = new Boolean(Parser.nextToken()).booleanValue();
			    }
			    else if (LineBuffer.equals("Reverse-Colors")) {
				this.ReverseColors = new Boolean(Parser.nextToken()).booleanValue();
			    }
			    else if (LineBuffer.equals("Class-Stats")) {
				this.ClassStats = new Boolean(Parser.nextToken()).booleanValue();
			    }
			    else if (LineBuffer.equals("Title-Font-Size")) {
				this.TitleFont = new Integer(Parser.nextToken()).intValue();
			    }
			    else if (LineBuffer.equals("Other-Font-Size")) {
				this.OtherFont = new Integer(Parser.nextToken()).intValue();
			    }
			    else {
				System.out.println("Didn't understand this line in BirdsEye.ini");
				System.out.println(WholeLine);
				System.out.println();
			    }
			}
			WholeLine = Incoming.readLine();
		    }
		    // Close file
		    Incoming.close();
		}
	    }
	    catch (IOException failure) {
		System.out.println("Failure while reading BirdsEye.ini.");
		System.out.println(failure);
	    }
	}
    }

}
