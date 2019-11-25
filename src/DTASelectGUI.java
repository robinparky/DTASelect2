import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

public class DTASelectGUI extends Frame implements ItemListener {
    TupleTable           Results = new TupleTable();
    SpecView             SpecWindow;
    Scrollbar            Scroller = new Scrollbar();
    Checkbox             ShowA = new Checkbox("a ions", false);
    Checkbox             ShowB = new Checkbox("b ions", true);
    Checkbox             ShowY = new Checkbox("y ions", true);
    Checkbox             Show0 = new Checkbox("* ions", true);
    Checkbox             Show1 = new Checkbox("+1 ions", true);
    Checkbox             Show2 = new Checkbox("+2 ions", true);
    DataSet              Data;
    StringTokenizer      Parser;
    String               DisplayedDTA = "";
    IDXFile              CompressedSpectra = null;
    int                  ProtCount = 0;
    int                  PepCount = 0;
    KeyboardHandler      KeyHearer = new KeyboardHandler();

    public DTASelectGUI(DataSet pData) {
	File        CurrentDirectory = new File(System.getProperty("user.dir"));
	File        ParticularFile = new File(CurrentDirectory, "DTASelect.IDX");
	this.Data = pData;
	if (ParticularFile.canRead()) {
	    ParticularFile = new File(CurrentDirectory, "DTASelect.SPM");
	    if (ParticularFile.canRead()) {
		CompressedSpectra = new IDXFile();
	    }
	}
	Protein             PRunner;
	Protein             InsidePRunner;
	DTAFile             DRunner;
	Panel               IonsPanel = new Panel();
	Panel               TopPanel = new Panel();
	int                 LocusTab = Data.LocusList.LongestName().length() + 1;
	int                 FilenameTab = Data.LocusList.LongestFileName().length() + 4;
	String              PeptideLabels;
	PeptideLabels = "Filename\t" + Data.hXCorr + "\t" +
	    Data.hDeltCN + "\tObsM+H+\t" + Data.hCalcPreMass +
	    "\tSpR\t" + Data.hSpScore + "\t%Ion\t#\tSequence";
	setSize(1000,700);
	setTitle("DTASelect");
	addWindowListener (new WindowEventAdapter());
	addKeyListener(KeyHearer);
	SpecWindow = new SpecView();
	SpecWindow.setSize(400,300);
	SpecWindow.addKeyListener(KeyHearer);
	IonsPanel.add(ShowA);
	IonsPanel.add(ShowB);
	IonsPanel.add(ShowY);
	IonsPanel.add(Show0);
	IonsPanel.add(Show1);
	IonsPanel.add(Show2);
	ShowA.addItemListener(this);
	ShowB.addItemListener(this);
	ShowY.addItemListener(this);
	Show0.addItemListener(this);
	Show1.addItemListener(this);
	Show2.addItemListener(this);
	TopPanel.setLayout(new BorderLayout());
	TopPanel.add("Center", SpecWindow);
	TopPanel.add("South", IonsPanel);
	this.setLayout(new BorderLayout());
	this.add("Center", Results);
	this.add("North", TopPanel);
	this.add("East", Scroller);
	Results.addScrollbar(Scroller);
	Results.addContainer(this);
	Results.addKeyListener(KeyHearer);
	Scroller.addAdjustmentListener(new ScrollHandler());
	Results.addTupleType(0,Color.black,true,
			     "0 10 20 30 40","Value\tProperty", false,1);
	Results.addTupleType(2,Color.blue,true,
			     "1 " + new Integer(FilenameTab).toString() + " " +
			     new Integer(FilenameTab + 8).toString() + " " +
			     new Integer(FilenameTab + 17).toString() + " " +
			     new Integer(FilenameTab + 26).toString() + " " +
			     new Integer(FilenameTab + 35).toString() + " " +
			     new Integer(FilenameTab + 40).toString() + " " +
			     new Integer(FilenameTab + 48).toString() + " " +
			     new Integer(FilenameTab + 56).toString() + " " +
			     new Integer(FilenameTab + 60).toString() + " 0",
			     PeptideLabels, true ,0);
	Results.addTupleType(1,Color.red,true,"0 " + new Integer(LocusTab).toString() + " " +
			     new Integer(LocusTab + 4).toString() + " " +
			     new Integer(LocusTab + 8).toString() + " " +
			     new Integer(LocusTab + 15).toString() + " " +
			     new Integer(LocusTab + 20).toString() + " " +
			     new Integer(LocusTab + 28).toString() + " " +
			     new Integer(LocusTab + 33).toString(),
			     "Locus\tDC\tSC\t%SC\tL\tMW\tpI\tDescriptive Name", true, 0);
	Parser = new StringTokenizer(Data.Cutoffs.PrintCriteria("","\n","\t"),"\n");
	while (Parser.hasMoreTokens())
	    Results.addTuple(0, Parser.nextToken(),true);
	Results.addTuple(0,"",true);
	PRunner = Data.LocusList.Next;
	ProtCount = 0;
	PepCount = 0;
	while (PRunner != null) {
	    Results.addTuple(1, PRunner.GetLocusTabbed(),true);
	    ProtCount++;
	    InsidePRunner = PRunner.IdenticalLoci;
	    while (InsidePRunner != null) {
		Results.addTuple(1, InsidePRunner.GetLocusTabbed(),true);
		InsidePRunner = InsidePRunner.IdenticalLoci;
	    }
	    DRunner = PRunner.DTAs.Next;
	    while (DRunner != null) {
		PepCount++;
		Results.addTuple(2, DRunner.GetTabDelimitedFields(),true);
		DRunner = DRunner.Next;
	    }
	    PRunner = PRunner.Next;
	}
	Results.addTuple(0, new Integer(ProtCount).toString() + "\tProteins", true);
	Results.addTuple(0, new Integer(PepCount).toString() + "\tPeptides", true);
	this.setVisible(true);
    }

    public void paint(Graphics GraphContext) {
	String       SelectedRow = Results.getSelectedPrimary();
	int          SelectedType = Results.getSelectedType();
	CoverageZone CZBuffer;
	if ( (SelectedType == 2) && (!SelectedRow.equals(DisplayedDTA)) ) {
	    DTAFile SelectedDTA = Data.LocusList.FindDTA(SelectedRow);
	    if (CompressedSpectra != null) {
		short              LowScan;
		short              HighScan;
		StringTokenizer    Parser = new StringTokenizer(SelectedDTA.FileName,".");
		Parser.nextToken();
		LowScan = new Short(Parser.nextToken()).shortValue();
		HighScan = new Short(Parser.nextToken()).shortValue();
		SpecWindow.DisplaySpectrum(SelectedDTA, Data.SequestParams,
			     CompressedSpectra.GetSpectrum(SelectedDTA.Subdirectory,
							   LowScan, HighScan));
	    }
	    else {
		SpecWindow.DisplayRaw(SelectedDTA, Data.SequestParams);
	    }
	    DisplayedDTA = SelectedRow;
	}
	else if ( (SelectedType == 1) && (!SelectedRow.equals(DisplayedDTA)) ) {
	    Protein SelectedProtein = Data.LocusList.FindProtein(SelectedRow);
	    if (SelectedProtein != null) {
		CZBuffer = SelectedProtein.GenerateZones();
		CZBuffer.MakeMinimal();
		SpecWindow.DisplayCoverage(SelectedProtein,
					   CZBuffer);
	    }
	    DisplayedDTA = SelectedRow;
	}
	Results.repaint();
	SpecWindow.repaint();
    }

    public void itemStateChanged(ItemEvent WhichItem) {
	SpecWindow.SetVisibility(ShowA.getState(),
				ShowB.getState(),
				ShowY.getState(),
				Show0.getState(),
				Show1.getState(),
				Show2.getState(),
				true);
	SpecWindow.repaint();
    }

    // Class to listen for and to handle window events
    class WindowEventAdapter extends WindowAdapter {
	public void windowClosing(WindowEvent WhatHappened) {
	    System.exit(0);
	}
    }

    class KeyboardHandler extends KeyAdapter {
	public void keyPressed(KeyEvent thing) {
	    switch (thing.getKeyCode()) {
	    case KeyEvent.VK_HOME:
		Scroller.setValue(0);
		Results.setSelectionByIndex(0);
		repaint();
		break;
	    case KeyEvent.VK_END:
		Scroller.setValue(Scroller.getMaximum() - Scroller.getVisibleAmount());
		Results.setSelectionByIndex(Results.TupleCount - 1);
		repaint();
		break;
	    case KeyEvent.VK_UP:
		Scroller.setValue(Scroller.getValue() - 1);
		Results.setSelectionByIndex(Results.getSelectionIndex() - 1);
		repaint();
		break;
	    case KeyEvent.VK_DOWN:
		Scroller.setValue(Scroller.getValue() + 1);
		Results.setSelectionByIndex(Results.getSelectionIndex() + 1);
		repaint();
		break;
	    case KeyEvent.VK_PAGE_UP:
		Scroller.setValue(Scroller.getValue() - Scroller.getVisibleAmount() + 1);
		Results.setSelectionByIndex(Results.getSelectionIndex() - Scroller.getVisibleAmount() + 1);
		repaint();
		break;
	    case KeyEvent.VK_PAGE_DOWN:
		Scroller.setValue(Scroller.getValue() + Scroller.getVisibleAmount() - 1);
		Results.setSelectionByIndex(Results.getSelectionIndex() + Scroller.getVisibleAmount() - 1);
		repaint();
		break;
	    }
	}
    }

    class ScrollHandler implements AdjustmentListener {
	public void adjustmentValueChanged(AdjustmentEvent thing) {
	    Results.repaint();
	}
    }
}
