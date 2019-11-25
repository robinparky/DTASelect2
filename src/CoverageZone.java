/* Create a class to store information about the starting and
 * finishing indecies of a subsequence within a larger sequence.
 */
class CoverageZone {
    int            Start;
    int            Finish;
    int            Depth = 1;
    String         Sequence="";
    CoverageZone   Next;

    //INDIVIDUAL FUNCTIONS

    public CoverageZone() {
    }

    //Provide a constructor that sets field values.
    public CoverageZone(int Start, String PassedSequence) {
	this.Start = Start;
	this.Finish = Start + PassedSequence.length() - 1;
	//System.out.println(">>>"+Start +"\t"+Finish+"\t"+PassedSequence);
	this.Sequence = PassedSequence;
    }

    public CoverageZone MakeClone() {
	CoverageZone CZBuffer = new CoverageZone();
	CZBuffer.Start = this.Start;
	CZBuffer.Finish = this.Finish;
	CZBuffer.Depth = this.Depth;
	CZBuffer.Sequence = this.Sequence;
	return CZBuffer;
    }

    /* Does this CZ cover the same region as the parameter one?
     * The +2 on the comparisons is to ensure that peptides that
     * abut each other are combined.  One peptide, for example,
     * may end at position 6 while the next starts at position 7.
     * If 2 is added to 6, 8 > 7, so the two are combined. */
    private boolean CollidesWith(CoverageZone Other) {
	if (this.Start < Other.Start) {
	    return ((this.Finish+2) > Other.Start);
	}
	else return (this.Start < (Other.Finish+2));
    }

    private void DebugPrint() {
	System.out.println(this.Sequence + "\t" +
			   new Integer(this.Start).toString() +
			   "\t" + new Integer(this.Finish).toString() +
			   "\t" + new Integer(this.Depth).toString());
    }

    //Make the current CZ include the region of the passed one.
    private void Combine(CoverageZone Other) {
	boolean        OtherStartLower   = this.Start  > Other.Start;
	boolean        OtherFinishHigher = this.Finish <
	    Other.Finish;
	boolean        ThisFinishHigher = this.Finish >
	    Other.Finish;
	String         LeftSection;
	String         RightSection = "";
	if (OtherStartLower) {
	    LeftSection = Other.Sequence;
	}
	else {
	    LeftSection = this.Sequence;
	}
	if (OtherFinishHigher) {
	    if (!OtherStartLower)
		RightSection = Other.Sequence.substring(1+ (this.Finish - Other.Start));
	}
	else {
	    if(OtherStartLower && ThisFinishHigher)
		RightSection = this.Sequence.substring(1+ (Other.Finish - this.Start));
	}
	this.Sequence = LeftSection + RightSection;
	this.Start = Math.min(this.Start, Other.Start);
	this.Finish = Math.max(this.Finish, Other.Finish);
    }

    /* MergeWith clones the peptide sequences from the passed list
     * into this list.  This is useful if you have multiple protein
     * objects and want to assess their cumulative sequence
     * coverage.  */
    public void MergeWith(CoverageZone Other) {
	//Go to the end of this list
	CoverageZone      CZRunner = this;
	CoverageZone      OtherRunner = Other.Next;
	while (CZRunner.Next != null)
	    CZRunner = CZRunner.Next;
	//Stick clones of the other list here
	while (OtherRunner != null) {
	    CZRunner.Next = OtherRunner.MakeClone();
	    CZRunner = CZRunner.Next;
	    OtherRunner = OtherRunner.Next;
	}
    }

    //LIST FUNCTIONS

    // Header call for Quicksort function.  Sorts CoverageZones in list by depth
    public void	SortList() {
	//this is just an empty header; the real data starts at the next item
	if (this.Next != null)
	    this.Next = this.Next.Sort(null);
    }
    
    // Dave's quicksorter
    private CoverageZone Sort(CoverageZone Follower) {
	CoverageZone		ListAbove = null;
	CoverageZone		ListBelow = null;
	CoverageZone		PlaceHolder;
	CoverageZone		PlaceHolder2;
	PlaceHolder = this.Next;
	//Partition all remaining points of this linked list
	while (PlaceHolder != null) {
	    PlaceHolder2 = PlaceHolder.Next;
	    if (this.Depth > PlaceHolder.Depth) {
				//Move this item to list above this
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else if (this.Depth < PlaceHolder.Depth) {
				//Move this item to list below this point
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    else if (this.Start > PlaceHolder.Start) {
				//Move this item to list above this
		PlaceHolder.Next = ListAbove;
		ListAbove = PlaceHolder;
	    }
	    else {
				//Move this item to list below this point
		PlaceHolder.Next = ListBelow;
		ListBelow = PlaceHolder;
	    }
	    //Move to next item to be partitioned
	    PlaceHolder = PlaceHolder2;
	}
	if (ListBelow == null)
	    this.Next = Follower;
	else
	    this.Next = ListBelow.Sort(Follower);
	if (ListAbove == null)
	    return this;
	else
	    return ListAbove.Sort(this);
    }
    
    /* Provide sum of lengths of peptides in this list.  Used for
     * determining sequence coverage after MakeMinimal() has been
     * run. */
    public int SeqLength() {
	CoverageZone      CZRunner = this.Next;
	int               Sum = 0;
	while (CZRunner != null) {
	    Sum += CZRunner.Sequence.length();
	    CZRunner = CZRunner.Next;
	}
	return (Sum);
    }

    /* Compare each CoverageZone with every other exactly once for
     * this locus.  If two CZs cover the same region of sequence
     * or abut, we should combine them (make the first CZ span the
     * entire region and ditch the second).  If a CZ is expanded,
     * recompare it to everything following it in the list.  See
     * the bottom of this file for the CoverageZone class.  */
    public void MakeMinimal() {
	CoverageZone   OutsideRunner;
	CoverageZone   InsideRunner;
	OutsideRunner = this.Next;
	while ( (OutsideRunner != null) && (OutsideRunner.Next != null) ) {
	    InsideRunner = OutsideRunner;
	    while (InsideRunner.Next != null) {
		if (OutsideRunner.CollidesWith(InsideRunner.Next)) {
		    OutsideRunner.Combine(InsideRunner.Next);
		    InsideRunner.Next = InsideRunner.Next.Next;
		    InsideRunner = OutsideRunner;
		}
		else {
		    InsideRunner = InsideRunner.Next;
		}
	    }
	    OutsideRunner = OutsideRunner.Next;
	}
    }

    /*GetConsensusList returns a list of peptides from this set of
     * CoverageZones.  The list is in Consensus CGI format;
     * peptides are separated from each other by "+" symbols.
     * This function should be called on the null head of a
     * list. */
    public String GetConsensusList() {
	StringBuffer      Builder = new StringBuffer();
	CoverageZone      CZRunner = this.Next;
	int               Counter;
	int               CurrentLength;
	if (CZRunner != null) {
	    CurrentLength = CZRunner.Sequence.length();
	    for (Counter = 0; Counter < CurrentLength; Counter += 70) {
		Builder.append(CZRunner.Sequence.substring(Counter,
							   Math.min(CurrentLength, Counter+70)));
		Builder.append("+");
	    }
	    CZRunner = CZRunner.Next;
	}
	while (CZRunner != null) {
	    CurrentLength = CZRunner.Sequence.length();
	    for (Counter = 0; Counter < CurrentLength; Counter += 70) {
		Builder.append(CZRunner.Sequence.substring(Counter,
							   Math.min(CurrentLength, Counter+70)));
		Builder.append("+");
	    }
	    CZRunner = CZRunner.Next;
	}
	String        Returned = Builder.toString();
	if (Returned.length() == 0)
	    return "";
	else
	    return Returned.substring(0, Returned.length()-1);
    }

    /* GetSeqCovString returns a string describing the depth of
     * coverage throughout this protein's sequence.  The first number
     * in each reported pair is the position at which this region
     * starts, while the second is the length of the region.  Regions
     * before the first ampersand are single coverages, and the amount
     * of coverage increases one layer for each ampersand.  An
     * asterisk concludes the string. */
    public String GetSeqCovString() {
	StringBuffer       Builder = new StringBuffer();
	CoverageZone       CZRunner = this.Next;
	CoverageZone       ActiveList = null;
	CoverageZone       FinalList = new CoverageZone();
	CoverageZone       CZBuffer;
	int                CurrentPosition = 0;
	int                LastDepth = 0;
	int                CurrentDepth = 0;
	IntList            StartPositions = new IntList();
	IntList            ILStart = StartPositions;
	IntList            FirstUncovereds = new IntList();
	IntList            ILStop = FirstUncovereds;
	// Create Lists of start and stop positions for these CZs.
	while (CZRunner != null) {
	    ILStart.Next = new IntList();
	    ILStop.Next = new IntList();
	    ILStart = ILStart.Next;
	    ILStop = ILStop.Next;
	    ILStart.Data = CZRunner.Start + 1;
	    ILStop.Data = CZRunner.Finish + 2;
	    CZRunner = CZRunner.Next;
	}
	// Sort the start and stop positions
	StartPositions.SortList();
	FirstUncovereds.SortList();
	ILStart = StartPositions.Next;
	ILStop = FirstUncovereds.Next;
	/* For debugging
	System.out.println("Start");
	StartPositions.DebugPrint();
	System.out.println("Stop");
	FirstUncovereds.DebugPrint();
	*/
	while ( (ILStart != null) || (ILStop != null) ) {
	    // Move to the next earliest sequence position
	    if ( (ILStart != null) && (ILStart.Data < ILStop.Data) ) {
		CurrentPosition = ILStart.Data;
	    }
	    else {
		CurrentPosition = ILStop.Data;
	    }
	    // Assimilate any increases in sequence depth
	    while ( (ILStart != null) &&
		    (ILStart.Data == CurrentPosition) ) {
		CurrentDepth++;
		ILStart = ILStart.Next;
	    }
	    // Assimilate any decreases in sequence depth
	    while ( (ILStop != null) &&
		    (ILStop.Data == CurrentPosition) ) {
		CurrentDepth--;
		ILStop = ILStop.Next;
	    }
	    if (CurrentDepth > LastDepth) {
		// Create coverage zone objects for each level of depth
		while (CurrentDepth > LastDepth) {
		    CZBuffer = ActiveList;
		    ActiveList = new CoverageZone();
		    ActiveList.Next = CZBuffer;
		    ActiveList.Start = CurrentPosition;
		    LastDepth++;
		    ActiveList.Depth = LastDepth;
		}
	    }
	    else if (CurrentDepth < LastDepth) {
		// Conclude the appropriate number of coverage zone objects
		while (CurrentDepth < LastDepth) {
		    LastDepth--;
		    ActiveList.Finish = CurrentPosition;
		    CZBuffer = FinalList.Next;
		    FinalList.Next = ActiveList;
		    ActiveList = ActiveList.Next;
		    FinalList.Next.Next = CZBuffer;
		}
	    }
	    LastDepth = CurrentDepth;
	}
	FinalList.SortList();
	CZBuffer = FinalList.Next;
	LastDepth = 0;
	while ( (CZBuffer != null) && (CZBuffer.Depth < 6) ) {
	    if (CZBuffer.Depth > LastDepth) {
		Builder.append('&');
		LastDepth = CZBuffer.Depth;
	    }
	    else {
		Builder.append('^');
	    }
	    Builder.append(CZBuffer.Start);
	    Builder.append('+');
	    Builder.append(CZBuffer.Finish - CZBuffer.Start);
	    //	    CZBuffer.DebugPrint();
	    CZBuffer = CZBuffer.Next;
	}
	Builder.append('*');
	return Builder.toString();
    }

    class IntList {
	int       Data;
	IntList   Next;

	// LIST FUNCTIONS
	public int getLength() {
	    IntList    Runner = this.Next;
	    int        Counter = 0;
	    while (Runner != null) {
		Counter++;
		Runner = Runner.Next;
	    }
	    return Counter;
	}

	public void DebugPrint() {
	    IntList  Runner = this.Next;
	    while (Runner != null) {
		System.out.print(Runner.Data);
		System.out.print("\t");
		Runner = Runner.Next;
	    }
	    System.out.println();
	}

	// Header call for Quicksort function.  Sorts integers in list
	public void	SortList() {
	    //this is just an empty header; the real data starts at the next item
	    if (this.Next != null)
		this.Next = this.Next.Sort(null);
	}

	// Dave's quicksorter
	private IntList Sort(IntList Follower) {
	    IntList		ListAbove = null;
	    IntList		ListBelow = null;
	    IntList		PlaceHolder;
	    IntList		PlaceHolder2;
	    PlaceHolder = this.Next;
	    //Partition all remaining points of this linked list
	    while (PlaceHolder != null) {
		PlaceHolder2 = PlaceHolder.Next;
		if (this.Data > PlaceHolder.Data) {
				//Move this item to list above this
		    PlaceHolder.Next = ListAbove;
		    ListAbove = PlaceHolder;
		}
		else {
				//Move this item to list below this point
		    PlaceHolder.Next = ListBelow;
		    ListBelow = PlaceHolder;
		}
		//Move to next item to be partitioned
		PlaceHolder = PlaceHolder2;
	    }
	    if (ListBelow == null)
		this.Next = Follower;
	    else
		this.Next = ListBelow.Sort(Follower);
	    if (ListAbove == null)
		return this;
	    else
		return ListAbove.Sort(this);
	}
    }
}
