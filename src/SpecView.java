import java.io.*;
import java.util.*;
import java.awt.*;

/**   SpecView class
      David Tabb, August 17, 2000
*/

/* SpecView is a JDK 1.1 component for viewing tandem mass spectra of
 * peptides as well as their interpretations.  It is designed to work
 * with DTASelect.  It associates each ion series with a particular
 * color; red represents b ions while blue represents y ions.  Lighter
 * shades of these colors represent +2 charged species.  The residues
 * of the sequence corresponding to each breakpoint are paired with
 * the peaks evidencing that breakpoint.
 */
public class SpecView extends Canvas {
    // DisplayedInterpretation holds the results taken from the .out
    // file associated with the displayed spectrum
    DTAFile         DisplayedInterpretation;
    boolean         PrefixPresent;
    boolean         SuffixPresent;
    // DisplayedSpectrum stores the peaks drawn from the displayed
    // spectrum (.dta file) along with summary information about range.
    Spectrum        DisplayedSpectrum;
    // If DisplayedSpectrum is null, sequence coverage for this
    // protein will instead be displayed.
    Protein         DisplayedCoverage;
    CoverageZone    DisplayedCZs;
    // DisplayedParams stores the relevant features of the
    // sequest.params file.
    ParamsFile      DisplayedParams;
    // DisplayedMatch is a list of peaks corresponding to particular
    // fragment ions.
    PointList       DisplayedMatch;
    PointList       DisplayedNeutrals;
    PointList       DisplayedLadder;
    boolean         ShowA = false;
    boolean         ShowB = true;
    boolean         ShowY = true;
    // Show0 determines whether star ions are shown
    boolean         Show0 = true;
    // Show1 determines whether +1 fragment ions are shown
    boolean         Show1 = true;
    // Show2 determines whether +2 fragment ions are shown
    boolean         Show2 = true;
    boolean         ShowPrecursorLosses = true;
    float           LoMZ = 0.0f;
    float           HiMZ = 0.0f;
    float           FragPercentTIC = 0.0f;
    float           FragIntensitySum = 0f;
    int             FragCount = 0;

    //Select series to highlight in spectrum
    public void SetVisibility(boolean pA,
			      boolean pB,
			      boolean pY,
			      boolean p0,
			      boolean p1,
			      boolean p2,
			      boolean pre) {
	ShowA = pA;
	ShowB = pB;
	ShowY = pY;
	Show0 = p0;
	Show1 = p1;
	Show2 = p2;
	ShowPrecursorLosses = pre;
    }

    //Show this spectrum.  Paint gets called when this is in place
    public void DisplaySpectrum(DTAFile PassedDTA, ParamsFile
				PassedParams, Spectrum PassedSpectrum) {
	if (PassedSpectrum == null)
	    System.out.println("Null spectrum passed to SpecView.DisplaySpectrum");
	else {
	    DisplayedInterpretation = PassedDTA;
	    DisplayedParams = PassedParams;
	    DisplayedSpectrum = PassedSpectrum;
	    this.LoMZ = DisplayedSpectrum.LowMOverZ;
	    this.HiMZ = DisplayedSpectrum.HighMOverZ;
	    DisplayedMatch = FindMatchingPeaks();
	    if (ShowPrecursorLosses) {
		DisplayedNeutrals = FindPrecursorLosses();
	    }
	    repaint();
	}
    }

    //Show the coverage for this protein.
    public void DisplayCoverage(Protein PassedProtein, CoverageZone PassedCZs) {
	if (PassedProtein == null)
	    System.out.println("Null protein passed to SpecView.DisplayCoverage");
	else {
	    DisplayedSpectrum = null;
	    DisplayedCoverage = PassedProtein;
	    DisplayedCZs = PassedCZs;
	    repaint();
	}
    }

    /*
     * This function was originally used to display spectra directly
     * from their associated .dta files.  It has fallen into disuse
     * since spectra can now be extracted from the SPM file.
     */
    /* 9/25/02 You know, once you write a feature, people will expect
     * it to work!  On with the show...
     */
    public void DisplayRaw(DTAFile NewFile, ParamsFile SEQUESTParams) {
	File       MS2File = new File(NewFile.Subdirectory + ".ms2");
	boolean    SuccessfulRead;
	DisplayedInterpretation = NewFile;
	DisplayedParams = SEQUESTParams;
	if (SEQUESTParams == null)
	    System.out.println("Null Params passed");
	// Set up DisplayedSpectrum
	if (MS2File.canRead()) {
	    SuccessfulRead = ReadMS2File();
	}
	else {
	    SuccessfulRead = ReadDTAFile();
	}
	if (SuccessfulRead) {
	    this.LoMZ = DisplayedSpectrum.LowMOverZ;
	    this.HiMZ = DisplayedSpectrum.HighMOverZ;
	    DisplayedMatch = FindMatchingPeaks();
	    this.FragIntensitySum = DisplayedSpectrum.IntensitySum();
	    this.FragCount = DisplayedSpectrum.PeakCount();
	    this.FragPercentTIC = DisplayedMatch.IntensitySum() /
		this.FragIntensitySum;
	    if (ShowPrecursorLosses) {
		DisplayedNeutrals = FindPrecursorLosses();
	    }
	}
	else {
	    this.DisplayedSpectrum = null;
	}
	repaint();
    }

    private boolean ReadMS2File() {
	try {
	    File              MS2File = new File(DisplayedInterpretation.Subdirectory +
						 ".ms2");
	    FileReader        MS2FileReader = new FileReader(MS2File);
	    BufferedReader    Incoming = new BufferedReader(MS2FileReader);
	    String            LineBuffer = Incoming.readLine();
	    StringTokenizer   Parser = new StringTokenizer(DisplayedInterpretation.FileName,
							   ".");
	    String            TargetScan;
	    String            CurrentScan = "";
	    boolean           FoundIt = false;
	    Point             PlaceHolder;
	    Parser.nextToken();
	    TargetScan = Parser.nextToken();
	    while ((!FoundIt) && (LineBuffer != null) ) {
		if (LineBuffer.startsWith("S\t")) {
		    Parser = new StringTokenizer(LineBuffer);
		    // Skip the S
		    Parser.nextToken();
		    CurrentScan = Parser.nextToken();
		    if (CurrentScan.equals(TargetScan)) {
			FoundIt = true;
		    }
		}
		LineBuffer = Incoming.readLine();
	    }
	    if (FoundIt) {
		DisplayedSpectrum = new Spectrum();
		PlaceHolder = DisplayedSpectrum.Points;
		DisplayedSpectrum.MaxIntensity = 0;
		DisplayedSpectrum.PrecursorMOverZ = (DisplayedInterpretation.PrecursorMass -
						     1.0f + DisplayedInterpretation.ChargeState) /
		    new Float(DisplayedInterpretation.ChargeState).floatValue();
		LineBuffer = Incoming.readLine();
		// There may be another charge state to which this spectrum is assigned
		if (LineBuffer.startsWith("Z\t")) {
		    LineBuffer = Incoming.readLine();
		}
		while ( (LineBuffer != null) &&
			(!LineBuffer.startsWith("S\t")) ) {
		    Parser = new StringTokenizer(LineBuffer);
		    PlaceHolder.Next = new Point();
		    PlaceHolder = PlaceHolder.Next;
		    PlaceHolder.MOverZ = new Float(Parser.nextToken()).floatValue();
		    PlaceHolder.Intensity = new Float(Parser.nextToken()).floatValue();
		    if (PlaceHolder.Intensity > DisplayedSpectrum.MaxIntensity) {
			DisplayedSpectrum.MaxIntensity = PlaceHolder.Intensity;
		    }
		    LineBuffer = Incoming.readLine();
		}
		DisplayedSpectrum.LowMOverZ = DisplayedSpectrum.Points.Next.MOverZ;
		DisplayedSpectrum.HighMOverZ = PlaceHolder.MOverZ;
		return true;
	    }
	    else {
		System.out.println("Could not find scan " + TargetScan + " in " +
				   MS2File.toString());
		return false;
	    }
	}
	catch (IOException failure) {
	    System.out.println("IO Error while reading " +
			       DisplayedInterpretation.Subdirectory +
			       ".ms2, spectrum " +
			       DisplayedInterpretation.FileName);
	    return false;
	}
    }

    /*
    * This function has fallen into disuse as well since spectra are
    * now extracted from the SPM file rather than their associated DTAs
    */
    /* 9/25/02 You know, once you write a feature, people will expect
     * it to work!  On with the show...
    */
    private boolean ReadDTAFile() {
	try {
	File                   CurrentDirectory = new File(System.getProperty("user.dir"));
	File                   SubDirectory = new File(CurrentDirectory, DisplayedInterpretation.Subdirectory);
	File                   DTAFile = new File(SubDirectory, DisplayedInterpretation.FileName + ".dta");
	FileReader             DTAFileReader = new FileReader(DTAFile);
	BufferedReader         Incoming = new BufferedReader(DTAFileReader);
	String                 LineBuffer = Incoming.readLine();
	StringTokenizer        Parser = new StringTokenizer(LineBuffer);
	Point                  PlaceHolder;
	int                    Counter = 0;
	// System.out.println("SpecView.ReadDTAFile() Called");
	DisplayedSpectrum = new Spectrum();
	PlaceHolder = DisplayedSpectrum.Points;
	DisplayedSpectrum.MaxIntensity = 0;
	DisplayedSpectrum.PrecursorMOverZ = (DisplayedInterpretation.PrecursorMass -
					     1.0f + DisplayedInterpretation.ChargeState) /
	    new Float(DisplayedInterpretation.ChargeState).floatValue();
	LineBuffer = Incoming.readLine();
	while (LineBuffer != null) {
	    Parser = new StringTokenizer(LineBuffer);
	    Counter++;
	    PlaceHolder.Next = new Point();
	    PlaceHolder = PlaceHolder.Next;
	    PlaceHolder.MOverZ = new Float(Parser.nextToken()).floatValue();
	    PlaceHolder.Intensity = new Float(Parser.nextToken()).floatValue();
	    if (PlaceHolder.Intensity > DisplayedSpectrum.MaxIntensity)
		DisplayedSpectrum.MaxIntensity = PlaceHolder.Intensity;
	    LineBuffer = Incoming.readLine();
	}
	DisplayedSpectrum.LowMOverZ = DisplayedSpectrum.Points.Next.MOverZ;
	DisplayedSpectrum.HighMOverZ = PlaceHolder.MOverZ;
	return true;
	}
	catch (IOException failure) {
	    System.out.println("IO Error while reading " +
			       DisplayedInterpretation.FileName +
			       ".dta");
	    return false;
	}
    }

    /* As always, this function is the most important of the
     * component.  Paint is called any time the visual display needs
     * to be restored.  If repaint() is called, paint() will surely
     * follow.
     */
    public void paint(Graphics GraphContext) {
	// Find out how big the canvas is
	Dimension                Bounds = getSize();
	// Variables determining look of spectrum
	// Labels and tick marks appear in margin!!!
	// Space for the y axis labels and tick marks:
	int                      MarginLeft = 15;
	// Space for the spectrum sequence inference
	int                      MarginAbove;
	// Space for the x axis labels and tick marks:
	int                      MarginBelow = Bounds.height - 15;
	// Space where labels and sequence can finish out.  Nothing
	// really drawn here:
	int                      MarginRight = Bounds.width - 15;
	// How many pixels wide is the spectrum?
	int                      SpectrumWidth = MarginRight - MarginLeft;
	// How many pixels tall is the spectrum?
	int                      SpectrumHeight;
	// How many pixels tall is the font?
	int                      FontHeight;
	// Gimme a handle to learn about the font
	FontMetrics              CurrentFontInfo;
	// Gimme a looper to roll through axes values
	float                    Counter;
	// Gimme a limit when generating axis values
	float                    Limit;
	// What's our ratio of pixels to m/z shown?
	float                    Density;
	// How far apart in m/z are the labels?
	float                    BigSep;
	// How far apart in m/z are the medium ticks?
	float                    MedSep;
	// How far apart in m/z are the tiny ticks?
	float                    TinySep;
	// How about a buffer for little bits of math?
	double                   DoubleBuffer;
	// Gimme a buffer to hold X coords when plotting things
	int                      XPos;
	// Gimme a buffer to hold Y coords when plotting things
	int                      YPos;
	// What's the first series not used in a sequence display at the top?
	int                      MaxSeries = 2;
	// If I need to roll through lists of peaks:
	Point                    Runner;
	// If I need to roll through lists of identified peaks:
	PointList                MatchRunner;
	// Is this series to be shown?
	boolean                  OkayToShow = true;
	GraphContext.setColor(Color.white);
	GraphContext.fillRect(0,0, Bounds.width, Bounds.height);
	CurrentFontInfo = GraphContext.getFontMetrics();
	FontHeight = CurrentFontInfo.getHeight();
	// Set aside different amounts of real estate at the top on
	// the basis of the precursor charge state.
	if (DisplayedInterpretation != null && DisplayedInterpretation.ChargeState > 2)
	    MarginAbove = FontHeight * 4 + FontHeight;
	else
	    MarginAbove = FontHeight * 2 + FontHeight;
	SpectrumHeight = MarginBelow - MarginAbove;
	GraphContext.setColor(Color.black);
	// If we really do have a spectrum to show...
	if (DisplayedSpectrum != null) {
	    //Set boundaries if they're not currently set
	    if (this.HiMZ == 0.0f) {
		this.LoMZ = DisplayedSpectrum.LowMOverZ;
		this.HiMZ = DisplayedSpectrum.HighMOverZ;
	    }
	    Density = (this.HiMZ - this.LoMZ) / SpectrumWidth;
	    if (Density > 0.75f) {
		BigSep = 100f;
		MedSep = 50f;
		TinySep = 25f;
	    }
	    else if (Density > 0.25f) {
		BigSep = 50f;
		MedSep = 10f;
		TinySep = 5f;
	    }
	    else {
		BigSep = 10f;
		MedSep = 5f;
		TinySep = 1f;
	    }
	    if (this.FragPercentTIC == 0.0f) {
		this.FragPercentTIC = DisplayedMatch.IntensitySum() /
		    DisplayedSpectrum.IntensitySum();
	    }
	    //Draw Axes
	    GraphContext.drawLine(MarginLeft, MarginAbove, MarginLeft, MarginBelow);
	    GraphContext.drawLine(MarginLeft, MarginBelow, MarginRight, MarginBelow);
	    //Plot X-axis tick marks
	    Counter = BigSep * new Double(Math.ceil(DisplayedSpectrum.LowMOverZ/BigSep)).floatValue();
	    Limit = DisplayedSpectrum.HighMOverZ;
	    while (Counter < Limit) {
		XPos = MarginLeft + getXPosFromMOverZ(Counter, SpectrumWidth);
		GraphContext.drawLine(XPos, MarginBelow, XPos, MarginBelow + 5);
		GraphContext.drawString(new
		    Integer(Math.round(Counter)).toString(),
					XPos, MarginBelow + FontHeight);
		Counter += BigSep;
	    }
	    Counter = MedSep * new Double(Math.ceil(DisplayedSpectrum.LowMOverZ/MedSep)).floatValue();
	    while (Counter < Limit) {
		XPos = MarginLeft + getXPosFromMOverZ(Counter, SpectrumWidth);
		GraphContext.drawLine(XPos, MarginBelow, XPos, MarginBelow + 3);
		Counter += MedSep;
	    }
	    Counter = TinySep * new Double(Math.ceil(DisplayedSpectrum.LowMOverZ/TinySep)).floatValue();
	    while (Counter < Limit) {
		XPos = MarginLeft + getXPosFromMOverZ(Counter, SpectrumWidth);
		GraphContext.drawLine(XPos, MarginBelow, XPos, MarginBelow + 2);
		Counter += TinySep;
	    }
	    //Plot Y-axis tick marks
	    Counter = 0f;
	    Limit = 1.1f;
	    while (Counter < Limit) {
		YPos = MarginAbove + getYPosFromIntensity(Counter *
			    DisplayedSpectrum.MaxIntensity, SpectrumHeight);
		GraphContext.drawLine(MarginLeft, YPos, MarginLeft - 5, YPos);
		Counter += .1f;
	    }
	    //Plot points
	    GraphContext.setColor(Color.gray);
	    Runner = DisplayedSpectrum.Points.Next;
	    while (Runner != null) {
		XPos = MarginLeft + getXPosFromMOverZ(Runner.MOverZ, SpectrumWidth);
		YPos = MarginAbove + getYPosFromIntensity(Runner.Intensity, SpectrumHeight);
		GraphContext.drawLine(XPos, YPos, XPos, MarginBelow);
		Runner = Runner.Next;
	    }
	    //Plot precursor ion position
	    GraphContext.setColor(GiveColorFor(6));
	    XPos = MarginLeft +
		getXPosFromMOverZ(DisplayedSpectrum.PrecursorMOverZ,
				  SpectrumWidth);
	    GraphContext.drawLine(XPos-1, MarginBelow, XPos-1, MarginBelow + 15);
	    GraphContext.drawLine(XPos, MarginBelow, XPos, MarginBelow + 15);
	    GraphContext.drawLine(XPos+1, MarginBelow, XPos+1, MarginBelow + 15);
	    if (ShowPrecursorLosses) {
		//Plot neutral losses from precursor
		MatchRunner = DisplayedNeutrals.Next;
		while (MatchRunner != null) {
		    GraphContext.setColor(GiveColorFor(MatchRunner.Series));
		    XPos = MarginLeft +
			getXPosFromMOverZ(MatchRunner.ThisPoint.MOverZ,
					  SpectrumWidth);
		    YPos = MarginAbove +
			getYPosFromIntensity(MatchRunner.ThisPoint.Intensity,
					     SpectrumHeight);
		    GraphContext.drawLine(XPos, YPos, XPos, MarginBelow);
		    GraphContext.drawString(MatchRunner.Identifier, XPos + 2, YPos - 2);
		    MatchRunner = MatchRunner.Next;
		}
	    }
	    //Label and color matching fragment ions
	    MatchRunner = DisplayedMatch.Next;
	    GraphContext.setFont(new Font("SanSerif",Font.BOLD,12));
	    this.FragPercentTIC = 0f;
	    while (MatchRunner != null) {
		GraphContext.setColor(GiveColorFor(MatchRunner.Series));
		if (MatchRunner.ThisPoint != null) {
		    switch (MatchRunner.Series) {
		    case 0:
			OkayToShow = ShowB && Show1;
			break;
		    case 1:
			OkayToShow = ShowY && Show1;
			break;
		    case 2:
			OkayToShow = ShowB && Show2;
			break;
		    case 3:
			OkayToShow = ShowY && Show2;
			break;
		    case 4:
			OkayToShow = ShowA && Show1;
			break;
		    case 5:
			OkayToShow = ShowA && Show2;
			break;
		    default:
			System.out.println("What's series " +
					   MatchRunner.Series + "?");
			System.exit(0);
		    }
		    if (OkayToShow) {
			this.FragPercentTIC += MatchRunner.ThisPoint.Intensity;
			XPos = MarginLeft +
			    getXPosFromMOverZ(MatchRunner.ThisPoint.MOverZ,
					      SpectrumWidth);
			YPos = MarginAbove +
			    getYPosFromIntensity(MatchRunner.ThisPoint.Intensity,
						 SpectrumHeight);
			GraphContext.drawLine(XPos, YPos, XPos, MarginBelow);
			GraphContext.drawString(MatchRunner.Identifier, XPos + 2, YPos - 2);
			if (Show0 && (MatchRunner.StarPoint != null) ) {
			    this.FragPercentTIC += MatchRunner.StarPoint.Intensity;
			    XPos = MarginLeft +
				getXPosFromMOverZ(MatchRunner.StarPoint.MOverZ,
						  SpectrumWidth);
			    YPos = MarginAbove +
				getYPosFromIntensity(MatchRunner.StarPoint.Intensity,
						     SpectrumHeight);
			    GraphContext.drawLine(XPos, YPos, XPos, MarginBelow);
			    GraphContext.drawString(MatchRunner.Identifier + "*", XPos + 2, YPos - 2);
			}
		    }
		}
		MatchRunner = MatchRunner.Next;
	    }
	    GraphContext.setFont(new Font("SanSerif",Font.PLAIN,12));
	    GraphContext.setColor(Color.black);
	    GraphContext.drawString("  Max:" + new
				    Float(DisplayedSpectrum.MaxIntensity).toString(),
				    MarginLeft, MarginAbove + FontHeight);
	    DoubleBuffer = Math.log(this.FragIntensitySum)*10f;
	    DoubleBuffer = Math.rint(DoubleBuffer);
	    GraphContext.drawString("  TIC:" +
				    Double.toString(DoubleBuffer/10f),
				    MarginLeft, MarginAbove + FontHeight * 2);
	    this.FragPercentTIC /= this.FragIntensitySum;
	    GraphContext.drawString("  Pks:" + Integer.toString(this.FragCount),
				    MarginLeft, MarginAbove + FontHeight * 3);
	    GraphContext.drawString("  Frags:" + new
				    Float(Protein.RoundTo(100* this.FragPercentTIC,1)).toString() +
				    "%",
				    MarginLeft, MarginAbove + FontHeight * 4);
	    //Plot sequence interpretations for each series
	    if (DisplayedInterpretation.ChargeState > 2)
		MaxSeries = 4;
	    MatchRunner = DisplayedMatch.Next;
	    while (MatchRunner != null) {
		if (MatchRunner.Series < MaxSeries) {
		    // Plot letter
		    if (MatchRunner.ThisPoint == null)
			GraphContext.setColor(Color.black);
		    else
			GraphContext.setColor(GiveColorFor(MatchRunner.Series));
		    XPos = MarginLeft + getXPosFromMOverZ(MatchRunner.MOverZ, SpectrumWidth) - 15;
		    if (MatchRunner.Modified)
			GraphContext.setFont(new Font("SanSerif",Font.ITALIC,12));
		    GraphContext.drawString(MatchRunner.Letter,
					    XPos, FontHeight + FontHeight * MatchRunner.Series);
		    if (MatchRunner.Modified)
			GraphContext.setFont(new Font("SanSerif",Font.PLAIN,12));
		    // Draw the little bracket pointing to the appropriate peak
		    XPos = MarginLeft + getXPosFromMOverZ(MatchRunner.MOverZ, SpectrumWidth);
		    GraphContext.drawLine(XPos - 5, 10 + FontHeight * MatchRunner.Series,
					  XPos, 10 + FontHeight * MatchRunner.Series);
		    GraphContext.drawLine(XPos, 10 + FontHeight * MatchRunner.Series,
					  XPos, 3 + FontHeight + FontHeight * MatchRunner.Series);
		}
		MatchRunner = MatchRunner.Next;
	    }

	    // Plot any ladders passed via the applet
	    if (DisplayedLadder != null) {
		// What was the series of the last peak in the ladder?
		int                      LastSeries = 73;
		boolean                  NewLadder = false;
		int                      LastX = 0;
		int                      LastY = 0;
		XPos = 0;
		YPos = 0;
		MatchRunner = DisplayedLadder.Next;
		while (MatchRunner != null) {
		    NewLadder = ( (MatchRunner.Series != LastSeries) ||
				  (MatchRunner.Series == 73) );
		    if (NewLadder) {
			GraphContext.setColor(GiveColorFor(MatchRunner.Series));
		    }
		    if (MatchRunner.ThisPoint != null) {
			XPos = MarginLeft +
			    getXPosFromMOverZ(MatchRunner.ThisPoint.MOverZ,
					      SpectrumWidth);
			YPos = MarginAbove +
			    getYPosFromIntensity(MatchRunner.ThisPoint.Intensity,
						 SpectrumHeight);
			GraphContext.drawLine(XPos-2, YPos-2, XPos+2, YPos+2);
			GraphContext.drawLine(XPos-2, YPos+2, XPos+2, YPos-2);
		    }
		    if (!NewLadder) {
			GraphContext.drawLine(LastX, LastY, XPos, YPos);
		    }
		    LastX = XPos;
		    LastY = YPos;
		    LastSeries = MatchRunner.Series;
		    MatchRunner = MatchRunner.Next;
		}
	    }
	}
	else if (DisplayedCoverage != null) {
	    /* Show sequence coverage for the selected protein.
	     * Print the protein name and description at the top.
	     * Generate CoverageZones for the protein, then render the
	     * coverage to fill the available space. */
	    GraphContext.setFont(new Font("Courier",Font.PLAIN,12));
	    CurrentFontInfo = GraphContext.getFontMetrics();
	    int            FontWidth = CurrentFontInfo.charWidth('W');
	    int            RowsOnScreen = Bounds.height /
		CurrentFontInfo.getHeight();
	    int            ColsOnScreen = Bounds.width / FontWidth;
	    int            DescripRows = DisplayedCoverage.Gene.length() /
		ColsOnScreen;
	    int            Looper = 0;
	    int            DrawingRow = 2;
	    int            EndIndex;
	    GraphContext.setColor(Color.black);
	    GraphContext.drawString(DisplayedCoverage.Locus,
				    0,FontHeight);
	    /* This section draws the gene description onscreen.  This
	     * part could definitely be improved.  It assumes all
	     * characters are as wide as a capital w.*/
	    EndIndex = DisplayedCoverage.Gene.length();
	    while(Looper < EndIndex) {
		GraphContext.drawString(DisplayedCoverage.Gene.substring(
			  Looper, Math.min(EndIndex,Looper +
					   ColsOnScreen) ),
					0, DrawingRow*FontHeight);
		Looper += ColsOnScreen;
		DrawingRow++;
	    }
	    int TableColumns = (ColsOnScreen -5) / 11;
	    /* Draw the horizontal and vertical position information
	     * onscreen.  */
	    for (Looper = 0; Looper < TableColumns;
		 Looper++) {
		GraphContext.drawString(new Integer(Looper * 10 +
						    1).toString(),
					FontWidth * (Looper *
					11 + 5),
					DrawingRow*FontHeight);
	    }
	    DrawingRow++;
	    for (Looper = 0; Looper < (RowsOnScreen - DrawingRow) &&
		     Looper * 10 * TableColumns <
		     DisplayedCoverage.SequenceLength;
		 Looper++) {
		GraphContext.drawString(new Integer(Looper * 10 * TableColumns).toString(),
						    0,
						    (DrawingRow+Looper)*FontHeight);
	    }
	    /* Draw in dashes representing sequence. Looper represents
	     * sequence location to draw - 1.  */
	    CoverageZone  CZRunner;
	    int           RowOffset;
	    int           ColumnOffset;
	    boolean       NotCovered;
	    for (Looper = 0; Looper <
		     DisplayedCoverage.SequenceLength; Looper++) {
		NotCovered = true;
		CZRunner = DisplayedCZs.Next;
		while (NotCovered && (CZRunner != null) ) {
		    if (Looper == CZRunner.Start) {
			Looper = CZRunner.Finish;
			NotCovered = false;
		    }
		    CZRunner = CZRunner.Next;
		}
		if (NotCovered) {
		    RowOffset = Looper / (10 * TableColumns);
		    ColumnOffset = Looper % (10 * TableColumns);
		    ColumnOffset += ColumnOffset / 10;
		    GraphContext.drawString("-", FontWidth *
					    (5 + ColumnOffset), FontHeight
					    * (DrawingRow + RowOffset));
		}
	    }
	    /* Draw in sequences of identified peptides.  */
	    GraphContext.setColor(Color.red);
	    char      Buffer;
	    CZRunner = DisplayedCZs.Next;
	    while (CZRunner != null) {
		for (Looper = CZRunner.Start; Looper <=
			 CZRunner.Finish; Looper++) {
		    RowOffset = (Looper) / (10 * TableColumns);
		    ColumnOffset = (Looper) % (10 * TableColumns);
		    ColumnOffset += ColumnOffset / 10;
		    Buffer = CZRunner.Sequence.charAt(Looper -
						      CZRunner.Start);
		    GraphContext.drawString(new Character(Buffer).toString(),
					    FontWidth *
					    (5 + ColumnOffset), FontHeight
					    * (DrawingRow + RowOffset));
		}
		CZRunner = CZRunner.Next;
	    }
	}
	else {
	    GraphContext.drawString("No spectrum to show", 20,20);
	}
    }

    // Given an MOverZ, determine where on screen it should show up.
    private int getXPosFromMOverZ(float MOverZ, int Width) {
	return new Float( ( (MOverZ - this.LoMZ) /
			    (this.HiMZ - this.LoMZ) ) * Width).intValue();
	/*
    	return new Float( ( (MOverZ - DisplayedSpectrum.LowMOverZ) /
			    (DisplayedSpectrum.HighMOverZ -
			     DisplayedSpectrum.LowMOverZ) ) * Width).intValue();
	*/
    }

    // Given an Intensity, determine where on screen it should show up.
    private int getYPosFromIntensity(float Intensity, int Height) {
	return new Float( ((DisplayedSpectrum.MaxIntensity -
			    Intensity) /
			   DisplayedSpectrum.MaxIntensity ) *
			  Height).intValue();
    }

    public String ReportFoundIons(int PassedSeries) {
	PointList      PLRunner = DisplayedMatch.Next;
	StringBuffer   Ions = new StringBuffer();
	while (PLRunner != null) {
	    if (PLRunner.Series == PassedSeries) {
		if (PLRunner.ThisPoint != null) {
		    Ions.append("0\n");
		    Ions.append(PLRunner.Number);
		    Ions.append("\t");
		    Ions.append(Protein.RoundTo(PLRunner.MOverZ, 1));
		    Ions.append("\t");
		    Ions.append(Protein.RoundTo(PLRunner.MOverZ - PLRunner.ThisPoint.MOverZ, 2));
		    Ions.append("\n");
		}
		else {
		    if ((PLRunner.MOverZ > DisplayedSpectrum.LowMOverZ) &&
			(PLRunner.MOverZ < DisplayedSpectrum.HighMOverZ)) {
			Ions.append("1\n");
			Ions.append(PLRunner.Number);
			Ions.append("\t");
			Ions.append(Protein.RoundTo(PLRunner.MOverZ, 1));
			Ions.append("\n");
		    }
		}
	    }
	    PLRunner = PLRunner.Next;
	}
	return Ions.toString();
    }

    // Given a precursor mass and a spectrum, determine which peaks
    // may correspond to neutral losses from the precursor.
    public PointList FindPrecursorLosses() {
	PointList     Matches = new PointList();
	PointList     MRunner = Matches;
	Point         PBuffer;
	PBuffer = DisplayedSpectrum.Find(DisplayedSpectrum.PrecursorMOverZ -
					 17f/DisplayedInterpretation.ChargeState);
	if (PBuffer != null) {
	    MRunner.Next = new PointList();
	    MRunner = MRunner.Next;
	    MRunner.ThisPoint = PBuffer;
	    MRunner.MOverZ = PBuffer.MOverZ;
	    MRunner.Identifier = "*";
	    MRunner.Series = 6;
	}
	PBuffer = DisplayedSpectrum.Find(DisplayedSpectrum.PrecursorMOverZ -
					 18f/DisplayedInterpretation.ChargeState);
	if (PBuffer != null) {
	    MRunner.Next = new PointList();
	    MRunner = MRunner.Next;
	    MRunner.ThisPoint = PBuffer;
	    MRunner.MOverZ = PBuffer.MOverZ;
	    MRunner.Identifier = "o";
	    MRunner.Series = 6;
	}
	PBuffer = DisplayedSpectrum.Find(DisplayedSpectrum.PrecursorMOverZ -
					 98f/DisplayedInterpretation.ChargeState);
	if (PBuffer != null) {
	    MRunner.Next = new PointList();
	    MRunner = MRunner.Next;
	    MRunner.ThisPoint = PBuffer;
	    MRunner.MOverZ = PBuffer.MOverZ;
	    MRunner.Identifier = "Phs";
	    MRunner.Series = 6;
	}
	PBuffer = DisplayedSpectrum.Find(DisplayedSpectrum.PrecursorMOverZ -
					 42f/DisplayedInterpretation.ChargeState);
	if (PBuffer != null) {
	    MRunner.Next = new PointList();
	    MRunner = MRunner.Next;
	    MRunner.ThisPoint = PBuffer;
	    MRunner.MOverZ = PBuffer.MOverZ;
	    MRunner.Identifier = "Act";
	    MRunner.Series = 6;
	}
	return Matches;
    }

    // Give the correct color for each series
    public static Color GiveColorFor(int Series) {
	switch (Series) {
	case 0:
	    return Color.red;
	case 1:
	    return Color.blue;
	case 2:
	    return Color.magenta;
	case 3:
	    return Color.cyan.darker();
	case 4:
	    return Color.green.darker();
	case 5:
	    return Color.green;
	case 6:
	    return Color.orange.darker();
	default:
	    return Color.black;
	}
    }

    public PointList FindMatchingPeaks() {
	String        AllSeq = DisplayedInterpretation.Sequence;
	String        Sequence = DisplayedInterpretation.TrimmedSequence();
	String        ModSequence = DisplayedInterpretation.SymbolModString();
	String        stSymbol;
	int           Looper;
	int           SeqLength = Sequence.length();
	float         MassAccumulator;
	PointList     Matches = new PointList();
	PointList     MatchesRunner = Matches;
	char          ThisSymbol;
	char          ThisMod;
	boolean       Modified;
	if (AllSeq.charAt(2) == '-') {
	    //If there's a '-' at the start of the sequence, it's just
	    //a suffix, not the full sequence.
	    MassAccumulator = 0f;
	    while (MatchesRunner.Next != null)
		MatchesRunner = MatchesRunner.Next;
	    // Find the y ions
	    for (Looper = SeqLength-1; Looper > -1; Looper--) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetCIonsFor(MassAccumulator, SeqLength - Looper, "",
					   stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	    // Find the b ions
	    // MassAccumulator now stores the sum of residue masses in this suffix
	    // Find the a and b ions for the first breakpoint
	    MassAccumulator = DisplayedInterpretation.PrecursorMass - MassAccumulator - 19;
	    MatchesRunner.Next = GetNIonsFor(MassAccumulator, -SeqLength, "L",
					     "-", false, DisplayedSpectrum);
	    while (MatchesRunner.Next != null)
		MatchesRunner = MatchesRunner.Next;
	    for (Looper = 0; Looper < SeqLength; Looper++) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetNIonsFor(MassAccumulator, 1-(SeqLength-Looper), "L",
						 stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	}
	else if (AllSeq.charAt(AllSeq.length() - 3) == '-') {
	    //If there's a '-' at the end of the sequence, it's just
	    //a prefix, not the full sequence.
	    MassAccumulator = 0f;
	    for (Looper = 0; Looper < SeqLength; Looper++) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetNIonsFor(MassAccumulator, Looper+1, "",
						 stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	    // MassAccumulator now stores the sum of residue masses in this suffix
	    // Find the y ions for the first breakpoint
	    MassAccumulator = DisplayedInterpretation.PrecursorMass - MassAccumulator - 19;
	    MatchesRunner.Next = GetCIonsFor(MassAccumulator, -SeqLength, "L",
					     "-", false, DisplayedSpectrum);
	    while (MatchesRunner.Next != null)
		MatchesRunner = MatchesRunner.Next;
	    for (Looper = SeqLength-1; Looper > -1; Looper--) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetCIonsFor(MassAccumulator, -Looper, "L",
						 stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	}
	else {
	    // The full sequence is here!  Even though fragment ions
	    // should occur only _between_ amino acids, we'll include
	    // the full sum, too.
	    MassAccumulator = 0f;
	    for (Looper = 0; Looper < SeqLength; Looper++) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetNIonsFor(MassAccumulator, Looper+1, "",
					   stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	    MassAccumulator = 0f;
	    for (Looper = SeqLength-1; Looper > -1; Looper--) {
		ThisSymbol = Sequence.charAt(Looper);
		stSymbol = new Character(ThisSymbol).toString();
		ThisMod = ModSequence.charAt(Looper);
		if (DisplayedParams.AvgTypeForFragmentIons) {
		    MassAccumulator += DisplayedParams.AvgMasses[ThisSymbol - 'A'];
		}
		else {
		    MassAccumulator += DisplayedParams.MonoMasses[ThisSymbol - 'A'];
		}
		Modified = (ThisMod != ' ');
		if (Modified) {
		    MassAccumulator += DisplayedParams.DiffMods.getMassShiftFor(ThisMod);
		}
		MatchesRunner.Next = GetCIonsFor(MassAccumulator, SeqLength - Looper, "",
					   stSymbol, Modified, DisplayedSpectrum);
		while (MatchesRunner.Next != null)
		    MatchesRunner = MatchesRunner.Next;
	    }
	}
	return Matches;
	}

    // Given a sequence inference, determine which peaks in spectrum
    // match it.
    public PointList bFindMatchingPeaks() {
	String        AllSeq = DisplayedInterpretation.Sequence;
	String        Sequence = DisplayedInterpretation.TrimmedSequence();
	String        ModString = DisplayedInterpretation.SymbolModString();
	int           Looper;
	int           LooperPlus1;
	int           SequenceLength = Sequence.length();
	float         NTermAccum = DisplayedParams.NPepMod + DisplayedParams.NProtMod;
	float         CTermAccum = DisplayedParams.CPepMod + DisplayedParams.CProtMod;
	PointList     Matches = new PointList();
	PointList     MatchRunner = Matches;
	char          CurrentN;
	char          CurrentC;
	if (AllSeq.charAt(2) == '-') {
	    //If there's a '-' at the start of the sequence, it's just
	    //a suffix, not the full sequence.
	    System.out.println("Suffix");
	}
	if (AllSeq.charAt(AllSeq.length() - 3) == '-') {
	    //If there's a '-' at the end of the sequence, it's just
	    //a prefix, not the full sequence.
	    System.out.println("Prefix");
	}
	for (Looper = 0; Looper < SequenceLength; Looper++) {
	    LooperPlus1 = Looper+1;
	    if (DisplayedParams.AvgTypeForFragmentIons) {
		NTermAccum += DisplayedParams.AvgMasses[Sequence.charAt(Looper) - 'A'];
		CTermAccum += DisplayedParams.AvgMasses[Sequence.charAt(SequenceLength - Looper - 1) - 'A'];
	    }
	    else {
		NTermAccum += DisplayedParams.MonoMasses[Sequence.charAt(Looper) - 'A'];
		CTermAccum += DisplayedParams.MonoMasses[Sequence.charAt(SequenceLength - Looper - 1) - 'A'];
	    }
	    CurrentN = ModString.charAt(Looper);
	    CurrentC = ModString.charAt(SequenceLength - Looper - 1);
	    if (CurrentN != ' ')
		NTermAccum += DisplayedParams.DiffMods.getMassShiftFor(CurrentN);
	    if (CurrentC != ' ')
		CTermAccum += DisplayedParams.DiffMods.getMassShiftFor(CurrentC);
	    //Find A ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find(NTermAccum - 27.0f),
					     DisplayedSpectrum.Find(NTermAccum - 44.0f),
					     NTermAccum - 27.0f,
					     "a" + new Integer(LooperPlus1).toString(),
					     4,
					     LooperPlus1,
					     new Character(Sequence.charAt(Looper)).toString(),
					     CurrentN != ' ');
	    MatchRunner = MatchRunner.Next;
	    //Find B ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find(NTermAccum + 1.0f),
					     DisplayedSpectrum.Find(NTermAccum - 16.0f),
					     NTermAccum + 1.0f,
					     "b" + new Integer(LooperPlus1).toString(),
					     0,
					     LooperPlus1,
					     new Character(Sequence.charAt(Looper)).toString(),
					     CurrentN != ' ');
	    MatchRunner = MatchRunner.Next;
	    //Find Y ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find(CTermAccum + 19.0f),
					     DisplayedSpectrum.Find(CTermAccum + 2.0f),
					     CTermAccum + 19.0f,
					     "y" + new Integer(LooperPlus1).toString(),
					     1,
					     LooperPlus1,
					     new Character(Sequence.charAt(SequenceLength - Looper - 1)).toString(),
					     CurrentC != ' ');
	    MatchRunner = MatchRunner.Next;
	    // Now repeat for doubly-charged fragment ions!
	    //Find A ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find((NTermAccum - 26.0f)/2.0f),
					     DisplayedSpectrum.Find((NTermAccum - 43.0f)/2.0f),
					     (NTermAccum - 26.0f)/2.0f,
					     "a" + new Integer(LooperPlus1).toString(),
					     5,
					     LooperPlus1,
					     new Character(Sequence.charAt(Looper)).toString(),
					     CurrentN != ' ');
	    MatchRunner = MatchRunner.Next;
	    //Find B ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find((NTermAccum + 2.0f)/2.0f),
					     DisplayedSpectrum.Find((NTermAccum - 15.0f)/2.0f),
					     (NTermAccum + 2.0f)/2.0f,
					     "b" + new Integer(LooperPlus1).toString(),
					     2,
					     LooperPlus1,
					     new Character(Sequence.charAt(Looper)).toString(),
					     CurrentN != ' ');
	    MatchRunner = MatchRunner.Next;
	    //Find Y ion
	    MatchRunner.Next = new PointList(DisplayedSpectrum.Find((CTermAccum + 20.0f)/2.0f),
					     DisplayedSpectrum.Find((CTermAccum + 3.0f)/2.0f),
					     (CTermAccum + 20.0f)/2.0f,
					     "y" + new Integer(LooperPlus1).toString(),
					     3,
					     LooperPlus1,
					     new Character(Sequence.charAt(SequenceLength - Looper - 1)).toString(),
					     CurrentC != ' ');
	}
	//Matches.DebugPrint();
	return Matches;
    }

    public PointList GetCIonsFor(float CResSum, int IonNumber, String ExtraLabel,
				 String Letter, boolean Modified, Spectrum Peaks) {
	PointList  Matches = new PointList();
	PointList  MatchRunner = Matches;
	String     stIonNumber = new Integer(IonNumber).toString();
	//Find Y ion
	MatchRunner.Next = new PointList(Peaks.Find(CResSum + 19.0f),
					 Peaks.Find(CResSum + 2.0f),
					 CResSum + 19.0f,
					 "y" + ExtraLabel + stIonNumber,
					 1,
					 IonNumber,
					 Letter,
					 Modified);
	MatchRunner = MatchRunner.Next;
	// Now repeat for doubly-charged fragment ions!
	//Find Y ion
	MatchRunner.Next = new PointList(Peaks.Find((CResSum + 20.0f)/2.0f),
					 Peaks.Find((CResSum + 3.0f)/2.0f),
					 (CResSum + 20.0f)/2.0f,
					 "y" + ExtraLabel + stIonNumber,
					 3,
					 IonNumber,
					 Letter,
					 Modified);
	return Matches.Next;
    }

    public PointList GetNIonsFor(float NResSum, int IonNumber, String ExtraLabel,
				 String Letter, boolean Modified, Spectrum Peaks) {
	PointList  Matches = new PointList();
	PointList  MatchRunner = Matches;
	String     stIonNumber = new Integer(IonNumber).toString();
	// Start with singly charged fragments...
	//Find A ion
	MatchRunner.Next = new PointList(Peaks.Find(NResSum - 27.0f),
					 Peaks.Find(NResSum - 44.0f),
					 NResSum - 27.0f,
					 "a" + ExtraLabel + stIonNumber,
					 4,
					 IonNumber,
					 Letter,
					 Modified);
	MatchRunner = MatchRunner.Next;
	//Find B ion
	MatchRunner.Next = new PointList(Peaks.Find(NResSum + 1.0f),
					 Peaks.Find(NResSum - 16.0f),
					 NResSum + 1.0f,
					 "b" + ExtraLabel + stIonNumber,
					 0,
					 IonNumber,
					 Letter,
					 Modified);
	MatchRunner = MatchRunner.Next;
	// Now repeat for doubly-charged fragment ions!
	//Find A ion
	MatchRunner.Next = new PointList(Peaks.Find((NResSum - 26.0f)/2.0f),
					 Peaks.Find((NResSum - 43.0f)/2.0f),
					 (NResSum - 26.0f)/2.0f,
					 "a" + ExtraLabel + stIonNumber,
					 5,
					 IonNumber,
					 Letter,
					 Modified);
	MatchRunner = MatchRunner.Next;
	//Find B ion
	MatchRunner.Next = new PointList(Peaks.Find((NResSum + 2.0f)/2.0f),
					 Peaks.Find((NResSum - 15.0f)/2.0f),
					 (NResSum + 2.0f)/2.0f,
					 "b" + ExtraLabel + stIonNumber,
					 2,
					 IonNumber,
					 Letter,
					 Modified);
	return Matches.Next;
    }

}
