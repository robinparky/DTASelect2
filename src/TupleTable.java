import	java.awt.*;
import  java.awt.event.*;
import	java.util.*;

/* TupleTable is a component designed for viewing and scrolling
   through heterogeneous data.  Each type of tuple can be shown with a
   different set of tab stops and color.  The component can be
   scrolled to the next tuple of a particular type with a particular
   value in the primary field.
   
   Initiated September 25, 2000
   Copyright (c) 2001, David L. Tabb, University of Washington

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA

   Contact David L. Tabb at dtabb@u.washington.edu for more
   information.  My current postal address follows:

   The Scripps Research Institute
   10550 North Torrey Pines Rd.
   SR11, Department of Cell Biology 
   La Jolla, CA 92037
   (858)784-8876
*/

public class TupleTable extends Canvas{
    Tuple		Tuples = new Tuple();
    //Store the tuple currently selected
    Tuple		SelectedTuple;
    //Store the top tuple onscreen
    Tuple		TopTuple;
    //Store types
    TupleType           Types = new TupleType();
    //Graphics concerns
    int                 FontHeight;
    //Vertical scrollbar
    Scrollbar           Scroller;
    //Values for scrollbar
    int                 TopIndex;
    int                 TupleCount = 0;
    Rectangle           Bounds;
    Font                CurrentFont = new Font("Courier",Font.PLAIN,12);
    FontMetrics         CurrentFontInfo = getFontMetrics(CurrentFont);
    Component           Holder = null;

    public TupleTable() {
	super();
	setFont(CurrentFont);
	addMouseListener(new MouseHandler());
    }

    public void addContainer(Component BigGuy) {
	Holder = BigGuy;
    }

    public void addScrollbar(Scrollbar handle) {
	if (handle == null) {
	    System.out.println("Null scrollbar passed to TupleTable");
	    System.exit(0);
	}
	Scroller = handle;
	Scroller.setMinimum(0);
    }

    public void paint(Graphics GraphContext) {
	Tuple                    Runner = Tuples.Next;
	TupleType                TRunner = Types.Next;
	StringTokenizer          FieldsParser;
	StringTokenizer          TabsParser;
	int                      Counter;
	int                      CurrentTab;
	int                      LeftMargin = 10;
	int                      VisibleHeaderCount;
	int                      RowsOnScreen;
	String                   CurrentField;
	Bounds = getBounds();
	GraphContext.setColor(Color.white);
	GraphContext.fillRect(0,0, Bounds.width, Bounds.height);
	VisibleHeaderCount = Types.getVisibleHeaderCount();
	RowsOnScreen = Bounds.height / CurrentFontInfo.getHeight() - VisibleHeaderCount;
	if (Scroller != null) {
	    Scroller.setVisibleAmount(RowsOnScreen);
	    Scroller.setBlockIncrement(RowsOnScreen - 1);
	    Scroller.setMaximum(TupleCount);
	}
	Counter = 0;
	//Display any visible headers
	while (TRunner != null) {
	    if (TRunner.Visible && TRunner.HeadVisible) {
		GraphContext.setColor(TRunner.Legend);
		FieldsParser = new StringTokenizer(TRunner.getHeader(),"\t");
		TabsParser = new StringTokenizer(TRunner.getTabStops()," ");
		while (TabsParser.hasMoreTokens() && FieldsParser.hasMoreTokens()) {
		    CurrentTab = new Integer(TabsParser.nextToken()).intValue();
		    CurrentField = FieldsParser.nextToken();
		    GraphContext.drawString(CurrentField,
					    LeftMargin + CurrentTab * CurrentFontInfo.charWidth('W'),
					    (Counter + 1) * CurrentFontInfo.getHeight()
					    );
		}
		Counter++;
	    }
	    TRunner = TRunner.Next;
	}
	GraphContext.setColor(Color.black);
	GraphContext.drawLine(0, VisibleHeaderCount * CurrentFontInfo.getHeight() + 2,
			      Bounds.width, VisibleHeaderCount * CurrentFontInfo.getHeight() + 2);
	//If there are any tuples to view
	if ( (Runner != null) && (Runner != Tuples) ) {
	    setCurrentPosition();
	    Runner = TopTuple;
	    Counter = 0;
	    while ( (Counter < RowsOnScreen) && (Runner != Tuples) ) {
		if (Runner.Visible && Runner.Type.Visible) {
		    if (Runner == SelectedTuple) {
			GraphContext.setColor(Color.yellow);
			GraphContext.fillRect(0, (Counter + VisibleHeaderCount) * CurrentFontInfo.getHeight()+2,
					      Bounds.width, CurrentFontInfo.getHeight());
		    }
		    GraphContext.setColor(Runner.Type.Legend);
		    FieldsParser = new StringTokenizer(Runner.getAllFields(),"\t");
		    TabsParser = new StringTokenizer(Runner.Type.getTabStops()," ");
		    while (TabsParser.hasMoreTokens() && FieldsParser.hasMoreTokens()) {
			CurrentTab = new Integer(TabsParser.nextToken()).intValue();
			CurrentField = FieldsParser.nextToken();
			GraphContext.drawString(CurrentField,
						LeftMargin + CurrentTab * CurrentFontInfo.charWidth('W'),
						(Counter + 1 + VisibleHeaderCount) * CurrentFontInfo.getHeight()
						);
		    }
		    Counter++;
		}
		Runner = Runner.Next;
	    }
	}
	GraphContext.setColor(Color.black);
	GraphContext.drawLine(0, VisibleHeaderCount * CurrentFontInfo.getHeight() + 2,
			      Bounds.width, VisibleHeaderCount * CurrentFontInfo.getHeight() + 2);
    }

    public int getSelectionIndex() {
	Tuple      TRunner = Tuples.Next;
	int        Counter = 0;
	if (SelectedTuple != null) {
	    while (TRunner != SelectedTuple) {
		TRunner = TRunner.Next;
		Counter++;
	    }
	    return Counter;
	}
	else return -1;
    }

    public void setSelectionByIndex(int NewIndex) {
	Tuple        TRunner = Tuples.Next;
	for (; NewIndex > 0; NewIndex--) {
	    if (TRunner.Next != Tuples)
		TRunner = TRunner.Next;
	}
	SelectedTuple = TRunner;
    }

    private void setCurrentPosition() {
	int         CurrentScrollPos;
	if (TopTuple == null) {
	    TopTuple = Tuples.Next;
	    TopIndex = 0;
	}
	if (Scroller != null) {
	    CurrentScrollPos = Scroller.getValue();
	    if (CurrentScrollPos > TopTuple.Index)
		while ( (CurrentScrollPos > TopTuple.Index) && (TopTuple.Next != Tuples) )
		    TopTuple = TopTuple.Next;
	    else if (CurrentScrollPos < TopTuple.Index)
		while ( (CurrentScrollPos < TopTuple.Index) && (TopTuple.Prev != Tuples) )
		    TopTuple = TopTuple.Prev;
	}
    }

    public void addTupleType(int Type, Color Legend, boolean Visible, String TabStops,
			     String Header, boolean HeadVisible, int PrimaryField) {
	TupleType          Buffer = new TupleType();
	Buffer.Next = Types.Next;
	Types.Next = Buffer;
	Buffer.Type = Type;
	Buffer.Legend = Legend;
	Buffer.Visible = Visible;
	Buffer.TabStops = TabStops;
	Buffer.Header = Header;
	Buffer.HeadVisible = HeadVisible;
	Buffer.PrimaryField = PrimaryField;
    }

    public void addTuple(int Type, String Values, boolean Visible) {
	//Adds the new item to the bottom end of the list
	Tuple            Buffer = Tuples.Prev;
	if ( (Tuples.Prev == null) || (Tuples.Next == null) ) {
	    Tuples.Next = new Tuple();
	    Tuples.Prev = Tuples.Next;
	    Buffer = Tuples.Next;
	    TopTuple = Tuples.Next;
	    SelectedTuple = Tuples.Next;
	}
	else {
	    Buffer = Tuples.Prev;
	    Tuples.Prev = new Tuple();
	    Tuples.Prev.Next = Tuples;
	    Tuples.Prev.Prev = Buffer;
	    Tuples.Prev.Prev.Next = Tuples.Prev;
	    Buffer = Tuples.Prev;
	}
	Buffer.Type = Types.getTupleType(Type);
	Buffer.Fields = Values;
	Buffer.Visible = Visible;
	Buffer.Index = TupleCount;
	if (Buffer.Type == null) {
	    System.out.println("Added tuple of nonexistent type");
	    System.exit(0);
	}
	TupleCount++;
    }

    public void setVisible(int Type, boolean Setting) {
	TupleType         CurrentType = Types.getTupleType(Type);
	if (CurrentType != null)
	    CurrentType.Visible = Setting;
	else
	    System.out.println("Can't set visibility for nonexistent type " + new Integer(Type).toString());
    }

    public void clear() {
	Tuples = new Tuple();
	SelectedTuple = null;
	TopTuple = null;
	TupleCount = 0;
    }

    public void clearTypes() {
	Types = new TupleType();
    }

    public int getSelectedType() {
	if (SelectedTuple == null)
	    return -1;
	else {
	    return SelectedTuple.Type.Type;
	}
    }

    //Return the primary field of the selected tuple
    public String getSelectedPrimary() {
	if (SelectedTuple == null)
	    return "";
	else {
	    return SelectedTuple.getPrimary();
	}
    }

    public void scrollToNext(int Type, String Target) {
	Tuple            Runner = TopTuple;
	boolean          Found = false;
	TupleType        ProbeType = Types.getTupleType(Type);
	//Proceed only if at least one tuple is in list!
	if (Runner != null) {
	    //Search the rest of the list
	    while ( (Runner != Tuples) && (!Found) ) {
		Found = ( (Runner.Type == ProbeType) && (Runner.getPrimary() == Target) );
		Runner = Runner.Next;
	    }
	    if (Found) {
		TopTuple = Runner;
		repaint();
	    }
	    else {
		//Search from the top
		Runner = Tuples.Next;
		while ( (Runner != Tuples) && (!Found) ) {
		    Found = ( (Runner.Type == ProbeType) && (Runner.getPrimary() == Target) );
		    Runner = Runner.Next;
		}
		if (Found) {
		    TopTuple = Runner;
		    repaint();
		}
		else {
		    System.out.println(Target + " not found");
		}
	    }
	}
    }

    class TupleType {
	int             Type;
	Color           Legend;
	boolean         Visible;
	String          TabStops;
	String          Header;
	boolean         HeadVisible;
	int             PrimaryField;
	TupleType       Next;

	//Return pointer to TupleType in this headered list which matches probe index
	public TupleType getTupleType(int Probe) {
	    TupleType              Runner = this.Next;
	    while ( (Runner != null) && (Runner.Type != Probe) ) {
		Runner = Runner.Next;
	    }
	    return Runner;
	}

	public int getVisibleHeaderCount() {
	    int            HCount = 0;
	    TupleType      Runner = this.Next;
	    while (Runner != null) {
		if (Runner.Visible && Runner.HeadVisible)
		    HCount++;
		Runner = Runner.Next;
	    }
	    return HCount;
	}

	public String getHeader() {
	    return Header;
	}

	public String getTabStops() {
	    return TabStops;
	}
    }

    class Tuple {
	TupleType	Type;
	String		Fields;
	boolean         Visible;
	int             Index;
	Tuple		Prev;
	Tuple		Next;

	public String getPrimary() {
	    return getField(Type.PrimaryField);
	}

	public String getField(int WhichField) {
	    String                 Extracted;
	    StringTokenizer        Parser = new StringTokenizer(Fields, "\t");
	    int                    Counter = 0;
	    while ( (Parser.hasMoreTokens()) && (Counter < WhichField) ) {
		Extracted = Parser.nextToken();
		Counter++;
	    }
	    if (Parser.hasMoreTokens())
		Extracted = Parser.nextToken();
	    else
		Extracted = "";
	    return Extracted;
	}

	public String getAllFields() {
	    return Fields;
	}
    }

    class MouseHandler implements MouseListener {
	/* MouseListener Interface requirements
	 */
	public void mouseClicked(MouseEvent WhatHappened) {
	    Bounds = getBounds();
	    int               VisibleHeaderCount = Types.getVisibleHeaderCount();
	    int               RowsBelowTop = WhatHappened.getY() /
		CurrentFontInfo.getHeight() - VisibleHeaderCount;
	    SelectedTuple = TopTuple;
	    while (RowsBelowTop > 0) {
		SelectedTuple = SelectedTuple.Next;
		RowsBelowTop--;
	    }
	    if (Holder != null)
		Holder.repaint();
	    else
		repaint();
	}

	public void mouseEntered(MouseEvent WhatHappened) {
	}

	public void mouseExited(MouseEvent WhatHappened) {
	}

	public void mousePressed(MouseEvent WhatHappened) {
	}

	public void mouseReleased(MouseEvent WhatHappened) {
	}
    }
}
