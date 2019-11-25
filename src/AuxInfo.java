/* AuxInfo object for DTASelect
 * David L. Tabb
 * April 29, 2002
 
 * Stores information describing add-in information for each protein.
 * Since all values are stored as floats, we note here whether
 * they're actually integers.
 */

public class AuxInfo {
    String        Descriptor;
    boolean       IsFloat = true;
    AuxInfo       Next;
}
