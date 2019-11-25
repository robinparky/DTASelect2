import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/* SpecBrowser lets users view the spectra in MS2 files without regard
 * for their sequences.  Initiated 9/27/03 by David Tabb of ORNL. */

public class SpecBrowser extends Frame {
    SpecView             SpecWindow;
    TextField            StatusReport;
    TextField            ScanEntry;
    Scrollbar            Scroller;

    public SpecBrowser() {
	Panel            BottomBar = new Panel();
	setSize(1000,700);
	setTitle("SpecBrowser");
	this.setLayout(new BorderLayout());
	StatusReport = new TextField("SpecBrowser 0.1");
	StatusReport.setEditable(false);
	ScanEntry = new TextField();
	Scroller = new Scrollbar(Scrollbar.HORIZONTAL,
				 1, 1, 1, 10);
	BottomBar.setLayout(new GridLayout(1,3));
	BottomBar.add(StatusReport);
	BottomBar.add(Scroller);
	BottomBar.add(ScanEntry);
	this.add("South", BottomBar);
	SpecWindow = new SpecView();
	this.add("Center", SpecWindow);
	this.show();
    }

    public static void main(String args[]) {
	SpecBrowser         AppObject = new SpecBrowser();
    }
}
