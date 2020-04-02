import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class IsoformUtils {

    public static final String ISOFORM_REGEX = "\\[([A-z]+)->([A-z]+)]";
    public static final Pattern ISOFORM_PATTERN = Pattern.compile(ISOFORM_REGEX);

    public static String getRegularSequence(String isoformSequence)
    {
        Matcher isoformMatcher = ISOFORM_PATTERN.matcher(isoformSequence);
        StringBuilder sb = new StringBuilder();
        int startLoc =0;
        while(isoformMatcher.find())
        {
            String prevSeq = isoformMatcher.group(1);
            sb.append(isoformSequence.substring(startLoc, isoformMatcher.start()));
            sb.append(prevSeq);
            startLoc = isoformMatcher.end();
        }
        sb.append(isoformSequence.substring(startLoc));
        return sb.toString();
    }


    public static void main(String [] args)
    {
        String test = "R.[R->C]RLLVA[AAFDA->C]K.M";
        String regeSequence = getRegularSequence(test);
        System.out.println(regeSequence);
    }


}
