/**   Point class
      David Tabb, August 17, 2000
    */

public class Point {
    float       MOverZ;
    float       Intensity;
    Point       Next;
    
    public void DebugPrint() {
	Point     PRunner = this.Next;
	System.out.println("Point.DebugPrint()");
	while(PRunner != null) {
	    System.out.println(new Float(PRunner.MOverZ).toString() +
			       "\t" + new Float(PRunner.Intensity).toString() );
	    PRunner = PRunner.Next;
	}
    }
}
