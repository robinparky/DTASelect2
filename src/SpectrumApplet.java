import java.util.*;
import java.applet.*;
import java.awt.*;
import java.awt.event.*;

public class SpectrumApplet extends Applet
    implements ActionListener, ItemListener{
    SpecView          ShownSpec = new SpecView();
    TupleTable        ShownIonsList = new TupleTable();
    TextField         LoMZ = new TextField();
    TextField         HiMZ = new TextField();
    Panel             IonsList = new Panel();
    ParamsFile        ShownParams = new ParamsFile();
    Checkbox          ShowA = new Checkbox("a", false);
    Checkbox          ShowB = new Checkbox("b", true);
    Checkbox          ShowY = new Checkbox("y", true);
    Checkbox          Show1 = new Checkbox("+1s", true);
    Checkbox          Show2 = new Checkbox("+2s", true);
    Checkbox          Show0 = new Checkbox("*", true);
    boolean           ShowPrecursor = true;

    public void itemStateChanged(ItemEvent WhichItem) {
	ShownSpec.SetVisibility(ShowA.getState(),
				ShowB.getState(),
				ShowY.getState(),
				Show0.getState(),
				Show1.getState(),
				Show2.getState(),
				ShowPrecursor);
	ShownSpec.repaint();
	IonsList.setVisible(false);
	ShownIonsList.setVisible(0,false);
	ShownIonsList.setVisible(1,false);
	ShownIonsList.setVisible(2,false);
	ShownIonsList.setVisible(3,false);
	ShownIonsList.setVisible(10,false);
	ShownIonsList.setVisible(11,false);
	ShownIonsList.setVisible(12,false);
	ShownIonsList.setVisible(13,false);
	ShownIonsList.setVisible(70,false);
	ShownIonsList.setVisible(71,false);
	ShownIonsList.setVisible(72,false);
	ShownIonsList.setVisible(73,false);
	if (ShowY.getState() && Show1.getState()) {
	    ShownIonsList.setVisible(1, true);
	    ShownIonsList.setVisible(11, true);
	    ShownIonsList.setVisible(71, true);
	}
	else if (ShowB.getState() && Show1.getState()) {
	    ShownIonsList.setVisible(0, true);
	    ShownIonsList.setVisible(10, true);
	    ShownIonsList.setVisible(70, true);
	}
	else if (ShowY.getState() && Show2.getState()) {
	    ShownIonsList.setVisible(3, true);
	    ShownIonsList.setVisible(13, true);
	    ShownIonsList.setVisible(73, true);
	}
	else if (ShowB.getState() && Show2.getState()) {
	    ShownIonsList.setVisible(2, true);
	    ShownIonsList.setVisible(12, true);
	    ShownIonsList.setVisible(72, true);
	}
	else {
	    ShownIonsList.setVisible(1, true);
	    ShownIonsList.setVisible(11, true);
	    ShownIonsList.setVisible(71, true);
	}
	IonsList.setVisible(true);
    }

    public void actionPerformed(ActionEvent WhichButton) {
	String Selected = WhichButton.getActionCommand();
	try {
	    float      NewLow = new Float(LoMZ.getText()).floatValue();
	    float      NewHigh = new Float(HiMZ.getText()).floatValue();
	    if (NewHigh > NewLow) {
		ShownSpec.LoMZ = NewLow;
		ShownSpec.HiMZ = NewHigh;
	    }
	}
	catch (NumberFormatException failure) {
	    System.out.println(Selected);
	}
	ShownSpec.repaint();
    }

    public void init() {
	Dimension        Bounds = getSize();
	int              Counter = 1;
	String           MZBuffer;
	String           IntBuffer;
	String           ModBuffer;
	String           SymbolBuffer;
	Spectrum         PassedSpec = new Spectrum();
	DTAFile          PassedInterp = new DTAFile();
	Point            PRunner = PassedSpec.Points;
	PointList        PassedLadder = null;
	int              Z;
	int              IonMissing;
	float            MPlusH;
	float            MaxInt = 0;
	StringTokenizer  Parser;
	Panel            TopPanel = new Panel();
	Panel            SidePanel = new Panel();
	Panel            RangePanel = new Panel();
	Panel            IonPanel = new Panel();
	this.setVisible(false);
	this.setLayout(new BorderLayout());
	// Process passed DTA information
	MPlusH = new Float(this.getParameter("PreMPlusH")).floatValue();
	PassedInterp.PrecursorMass = MPlusH;
	Z = new Integer(this.getParameter("PreZ")).intValue();
	PassedSpec.PrecursorMOverZ = (MPlusH + Z - 1) / Z;
	PassedSpec.PrecursorZ = Z;
	MZBuffer = this.getParameter("MZ" + 
				   new Integer(Counter).toString());
	PassedSpec.LowMOverZ = new Float(MZBuffer).floatValue();
	while (MZBuffer != null) {
	    IntBuffer = this.getParameter("Int" +
					  new Integer(Counter).toString());
	    PRunner.Next = new Point();
	    PRunner = PRunner.Next;
	    PRunner.MOverZ = new Float(MZBuffer).floatValue();
	    PRunner.Intensity = new Float(IntBuffer).floatValue();
	    if (PRunner.Intensity > MaxInt)
		MaxInt = PRunner.Intensity;
	    Counter++;
	    MZBuffer = this.getParameter("MZ" + 
					 new Integer(Counter).toString());
	}
	PassedSpec.MaxIntensity = MaxInt;
	PassedSpec.HighMOverZ = PRunner.MOverZ;
	ModBuffer = this.getParameter("LoMZ");
	if (ModBuffer != null)
	    PassedSpec.LowMOverZ = new Float(ModBuffer).floatValue();
	ModBuffer = this.getParameter("HiMZ");
	if (ModBuffer != null)
	    PassedSpec.HighMOverZ = new Float(ModBuffer).floatValue();
	LoMZ.setText(new Float(PassedSpec.LowMOverZ).toString());
	HiMZ.setText(new Float(PassedSpec.HighMOverZ).toString());
	// Process passed display info from sequest.params
	if (this.getParameter("CPepMod") != null)
	    ShownParams.CPepMod = new Float(this.getParameter("CPepMod")).floatValue();
	if (this.getParameter("NPepMod") != null)
	    ShownParams.NPepMod = new Float(this.getParameter("NPepMod")).floatValue();
	if (this.getParameter("CProtMod") != null)
	    ShownParams.CProtMod = new Float(this.getParameter("CProtMod")).floatValue();
	if (this.getParameter("NProtMod") != null)
	    ShownParams.NProtMod = new Float(this.getParameter("NProtMod")).floatValue();
	if (this.getParameter("AvgForFrag") != null)
	    ShownParams.AvgTypeForFragmentIons = new Boolean(this.getParameter("AvgForFrag")).booleanValue();
	if (this.getParameter("AvgForParent") != null)
	    ShownParams.AvgTypeForParentIon = new Boolean(this.getParameter("AvgForParent")).booleanValue();
	if (this.getParameter("ShowPrecursorLosses") != null)
	    ShowPrecursor = new Boolean(this.getParameter("ShowPrecursorLosses")).booleanValue();
	// Process passed static modification information
	Counter = 1;
	ModBuffer = this.getParameter("SMM" +
				      new Integer(Counter).toString());
	while (ModBuffer != null) {
	    SymbolBuffer = this.getParameter("SMR" +
				      new Integer(Counter).toString());
	    ShownParams.AvgMasses[SymbolBuffer.charAt(0) - 65] = new Float(ModBuffer).floatValue();
	    ShownParams.MonoMasses[SymbolBuffer.charAt(0) - 65] = new Float(ModBuffer).floatValue();
	    Counter++;
	    ModBuffer = this.getParameter("SMM" +
					  new Integer(Counter).toString());
	}
	// Process passed diff modification information
	Counter = 1;
	ModBuffer = this.getParameter("DMM" +
				      new Integer(Counter).toString());
	while (ModBuffer != null) {
	    SymbolBuffer = this.getParameter("DMS" +
				      new Integer(Counter).toString());
	    ShownParams.addDiffMod(new Float(ModBuffer).floatValue(),
				   SymbolBuffer.charAt(0), "");
	    Counter++;
	    ModBuffer = this.getParameter("DMM" +
					  new Integer(Counter).toString());
	}
	ShownParams.ApplyStaticMods();
	// Process any passed ladder information
	ModBuffer = this.getParameter("LadderLength1");
	if (ModBuffer != null) {
	    int       CurrentLadderLength;
	    int       CurrentLadderSeries;
	    PointList PLRunner;
	    Point     CurrentPeak;
	    int       PeakCounter;
	    int       LadderCounter = 1;
	    PassedLadder = new PointList();
	    PLRunner = PassedLadder;
	    Counter = 0;
	    ModBuffer = this.getParameter("LadderLength" + LadderCounter);
	    while (ModBuffer != null) {
		CurrentLadderLength = new Integer(ModBuffer).intValue();
		ModBuffer = this.getParameter("LadderSeries" + LadderCounter);
		if (ModBuffer == null)
		    System.out.println("Invalid series for ladder " + LadderCounter);
		CurrentLadderSeries = new Integer(ModBuffer).intValue();
		for (PeakCounter = 1; PeakCounter < (CurrentLadderLength + 1); PeakCounter++) {
		    ModBuffer = this.getParameter("LadMZ" + (PeakCounter + Counter));
		    if (ModBuffer == null)
			System.out.println("Invalid MZ for ladder point " + (PeakCounter + Counter));
		    CurrentPeak = new Point();
		    CurrentPeak.MOverZ = new Float(ModBuffer).floatValue();
		    ModBuffer = this.getParameter("LadInt" + (PeakCounter + Counter));
		    if (ModBuffer == null)
			System.out.println("Invalid intensity for ladder point " + PeakCounter + Counter);
		    CurrentPeak.Intensity = new Float(ModBuffer).floatValue();
		    PLRunner.Next = new PointList();
		    PLRunner = PLRunner.Next;
		    PLRunner.ThisPoint = CurrentPeak;
		    PLRunner.Series = CurrentLadderSeries;
		}
		Counter += CurrentLadderLength;
		LadderCounter++;
		// Leave a marker pointlist between Ladders
		PLRunner.Next = new PointList();
		PLRunner = PLRunner.Next;
		PLRunner.Series = 73;
		ModBuffer = this.getParameter("LadderLength" + LadderCounter);
	    }
	    PassedLadder.DebugPrint();
	}
	// Get interpretation of this spectrum
	PassedInterp.Sequence = this.getParameter("MatchSeq");
	if (PassedInterp.Sequence == null)
	    PassedInterp.Sequence = "-.G.-";
	PassedInterp.ChargeState = new Byte(this.getParameter("PreZ")).byteValue();
	ShownSpec.DisplayedInterpretation = PassedInterp;
	ShownSpec.DisplayedParams = ShownParams;
	ShownSpec.DisplayedSpectrum = PassedSpec;
	ShownSpec.DisplayedMatch = ShownSpec.FindMatchingPeaks();
	ShownSpec.FragIntensitySum = PassedSpec.IntensitySum();
	ShownSpec.FragCount = PassedSpec.PeakCount();
	ShownSpec.DisplayedNeutrals = ShownSpec.FindPrecursorLosses();
	ShownSpec.DisplayedLadder = PassedLadder;
	Bounds = getSize();
	if (Z < 3) {
	    Show2.setState(false);
	}
	IonPanel.setLayout(new GridLayout(5,2));
	IonPanel.add(ShowA);
	IonPanel.add(Show0);
	IonPanel.add(ShowB);
	IonPanel.add(Show1);
	IonPanel.add(ShowY);
	IonPanel.add(Show2);
	IonPanel.add(new Label("Lo m/z"));
	IonPanel.add(new Label("Hi m/z"));
	IonPanel.add(LoMZ);
	IonPanel.add(HiMZ);
	ShowA.addItemListener(this);
	ShowB.addItemListener(this);
	ShowY.addItemListener(this);
	Show0.addItemListener(this);
	Show1.addItemListener(this);
	Show2.addItemListener(this);
	LoMZ.addActionListener(this);
	HiMZ.addActionListener(this);
	// Add types for ions that are not present
	ShownIonsList.addTupleType(0, Color.black,false,"0 3", "#\tCalcMZ", false, 0);
	ShownIonsList.addTupleType(1, Color.black,true,"0 3", "#\tCalcMZ", false, 0);
	ShownIonsList.addTupleType(2, Color.black,false,"0 3", "#\tCalcMZ", false, 0);
	ShownIonsList.addTupleType(3, Color.black,false,"0 3", "#\tCalcMZ", false, 0);
	// Add types for ions that are present
	ShownIonsList.addTupleType(10, ShownSpec.GiveColorFor(0),
				   false,"0 3 11", "#\tCalcMZ\tDiff", true, 0);
	ShownIonsList.addTupleType(11, ShownSpec.GiveColorFor(1),
				   true,"0 3 11", "#\tCalcMZ\tDiff", true, 0);
	ShownIonsList.addTupleType(12, ShownSpec.GiveColorFor(2),
				   false,"0 3 11","#\tCalcMZ\tDiff", true, 0);
	ShownIonsList.addTupleType(13, ShownSpec.GiveColorFor(3),
				   false,"0 3 11", "#\tCalcMZ\tDiff", true, 0);
	// Add headers so the user knows what series is shown
	ShownIonsList.addTupleType(70, ShownSpec.GiveColorFor(0),false,"0","B +1", true, 0);
	ShownIonsList.addTupleType(71, ShownSpec.GiveColorFor(1),true,"0","Y +1", true, 0);
	ShownIonsList.addTupleType(72, ShownSpec.GiveColorFor(2),false,"0","B +2", true, 0);
	ShownIonsList.addTupleType(73, ShownSpec.GiveColorFor(3),false,"0","Y +2", true, 0);
	// Add the ions from the B +1 series
	Parser = new StringTokenizer(ShownSpec.ReportFoundIons(0), "\n");
	while (Parser.hasMoreTokens()) {
	    // Determine if the ion is present
	    IonMissing = new Integer(Parser.nextToken()).intValue();
	    if (IonMissing == 0) {
		ShownIonsList.addTuple(10, Parser.nextToken(), true);
	    }
	    else {
		ShownIonsList.addTuple(0, Parser.nextToken(), true);
	    }
	}
	// Add the ions from the Y +1 series
	Parser = new StringTokenizer(ShownSpec.ReportFoundIons(1), "\n");
	while (Parser.hasMoreTokens()) {
	    // Determine if the ion is present
	    IonMissing = new Integer(Parser.nextToken()).intValue();
	    if (IonMissing == 0) {
		ShownIonsList.addTuple(11, Parser.nextToken(), true);
	    }
	    else {
		ShownIonsList.addTuple(1, Parser.nextToken(), true);
	    }
	}
	// Add the ions from the B +2 series
	Parser = new StringTokenizer(ShownSpec.ReportFoundIons(2), "\n");
	while (Parser.hasMoreTokens()) {
	    // Determine if the ion is present
	    IonMissing = new Integer(Parser.nextToken()).intValue();
	    if (IonMissing == 0) {
		ShownIonsList.addTuple(12, Parser.nextToken(), true);
	    }
	    else {
		ShownIonsList.addTuple(2, Parser.nextToken(), true);
	    }
	}
	// Add the ions from the Y +2 series
	Parser = new StringTokenizer(ShownSpec.ReportFoundIons(3), "\n");
	while (Parser.hasMoreTokens()) {
	    // Determine if the ion is present
	    IonMissing = new Integer(Parser.nextToken()).intValue();
	    if (IonMissing == 0) {
		ShownIonsList.addTuple(13, Parser.nextToken(), true);
	    }
	    else {
		ShownIonsList.addTuple(3, Parser.nextToken(), true);
	    }
	}
	SidePanel.setLayout(new BorderLayout());
	IonsList.setLayout(new BorderLayout());
	IonsList.add(ShownIonsList, "Center");
	SidePanel.add(IonsList, "Center");
	SidePanel.add(IonPanel, "North");
	ShownSpec.SetVisibility(ShowA.getState(),
				ShowB.getState(),
				ShowY.getState(),
				Show0.getState(),
				Show1.getState(),
				Show2.getState(),
				ShowPrecursor);
	ShownSpec.setSize(Bounds.width-150, Bounds.height);
	this.add(SidePanel, "East");
	this.add(ShownSpec, "Center");
	this.setVisible(true);
    }
}
