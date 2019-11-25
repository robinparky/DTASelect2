package proteinGrouper;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/8/19.
 */

import java.util.Comparator;
//import java.util.List;

public class PeptideComparator implements Comparator {

    // how to comapare two PeptideComparaters
    private int sortBy = 1; // 1 sort by spectral count, 2 sort by XCorr, 3 sort by ZScore or Sp, 4 by scan number

    public PeptideComparator(int sortby) {
        this.sortBy = sortby;
    }
    public void setSortBy(int sortby) {
        sortBy = sortby;
    }


    public int compare(Object pep1, Object pep2) {

        Peptide p1 = (Peptide) pep1;
        Peptide p2 = (Peptide) pep2;
        double f = 0;  // difference

        switch (sortBy) {
            case 1:
                f = p1.getSpectralCount() - p2.getSpectralCount();
                break;
            case 2:
                f = p1.getXCorrValue() - p2.getXCorrValue();
                break;
            case 3:
                f = p1.getSpScoreValue() - p2.getSpScoreValue();
                break;
            case 4:
                f = p1.getScanNumber() - p2.getScanNumber();
                break;
        }



        if (f > 0) {
            return -1;
        } else if (f < 0) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * Return true if both the intensities and m2zs are the same,
     * otherwise return false
     */
    public boolean equals(Object o) {
        PeptideComparator c = (PeptideComparator)o;
        if(o == null) {
            return false;
        }
        return this.sortBy == c.sortBy;
    }
}


