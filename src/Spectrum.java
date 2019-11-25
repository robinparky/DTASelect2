/**        Spectrum class
	   David Tabb, October 12, 2000
*/
public class Spectrum {
    float              FragmentTolerance = 0.75f;
    Point              Points = new Point();
    float              LowMOverZ;
    float              HighMOverZ;
    float              MaxIntensity;
    float              PrecursorMOverZ;
    int                PrecursorZ;

    public float IntensitySum() {
	Point PRunner = this.Points.Next;
	float Sum = 0f;
	while (PRunner != null) {
	    Sum += PRunner.Intensity;
	    PRunner = PRunner.Next;
	}
	return Sum;
    }

    public int PeakCount() {
	Point PRunner = this.Points.Next;
	int Count = 0;
	while (PRunner != null) {
	    Count++;
	    PRunner = PRunner.Next;
	}
	return Count;
    }

    public Point Find(float MOverZ) {
	Point          PRunner = Points.Next;
	Point          BiggestPeak;
	float          EarlyEdge = MOverZ - FragmentTolerance;
	float          LateEdge = MOverZ + FragmentTolerance;
	//Advance to first possible peak
	while ((PRunner != null) && (PRunner.MOverZ < EarlyEdge))
	    PRunner = PRunner.Next;
	if ( (PRunner == null) || (PRunner.MOverZ > LateEdge) )
	    return null;
	else {
	    //Find tallest peak in this range
	    BiggestPeak = PRunner;
	    PRunner = PRunner.Next;
	    while ( (PRunner != null) && (PRunner.MOverZ < LateEdge) ) {
		if (PRunner.Intensity > BiggestPeak.Intensity)
		    BiggestPeak = PRunner;
		PRunner = PRunner.Next;
	    }
	    return BiggestPeak;
	}
    }
}
