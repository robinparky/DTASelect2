public class PointList {
    Point        ThisPoint;
    Point        StarPoint;
    float        MOverZ;
    String       Identifier;
    int          Series;
    int          Number;
    String       Letter;
    boolean      Modified = false;
    PointList    Next;

    public PointList() {
    }

    public PointList(Point   pThisPoint,
		     Point   pStarPoint,
		     float   pMOverZ,
		     String  pIdentifier,
		     int     pSeries,
		     int     pNumber,
		     String  pLetter,
		     boolean pModified) {
	this.ThisPoint = pThisPoint;
	this.StarPoint = pStarPoint;
	this.MOverZ = pMOverZ;
	this.Identifier = pIdentifier;
	this.Series = pSeries;
	this.Number = pNumber;
	this.Letter = pLetter;
	this.Modified = pModified;
    }

    public float IntensitySum() {
	PointList     PLRunner = this.Next;
	float         Sum = 0f;
	while (PLRunner != null) {
	    if (PLRunner.ThisPoint != null) {
		Sum += PLRunner.ThisPoint.Intensity;
	    }
	    PLRunner = PLRunner.Next;
	}
	return Sum;
    }

    public void DebugPrint() {
	PointList   PLRunner = this.Next;
	while (PLRunner != null) {
	    System.out.print("PL\t" + PLRunner.MOverZ + "\t" +
			       PLRunner.Identifier + "\t" +
			       PLRunner.Series + "\t" +
			       PLRunner.Number + "\t" +
			       PLRunner.Letter + "\t" +
			       PLRunner.Modified);
	    if (PLRunner.ThisPoint != null) {
		System.out.println("\t" + PLRunner.ThisPoint.MOverZ + "\t" +
				   PLRunner.ThisPoint.Intensity);
	    }
	    else {
		System.out.println();
	    }
	    PLRunner = PLRunner.Next;
	}
    }
}
