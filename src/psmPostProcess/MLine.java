package psmPostProcess;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version 1.0
 */
public class MLine
{
    private String mLine;
    private ArrayList<String> lLineList;
    private String xcorrRank;
    private String spRank;
    
    private String calMZ;
    private String deltCN;
    
    private String xcorr;
    private double xcorrValue;
    private String sp;
    private int spRankInt;
    private int primaryRankInt;
    private double spVale;
    private String matchedIons;
    private String predictedIons;
    private String sequence;
    private String status;
    public MLine(String mLine)
    {
//	this(mLine.split("\t"));

                this.mLine = mLine;
        /*
                String[] arr = this.mLine.split("\t");
                this.xcorrRank = arr[1];
                this.spRank = arr[2];
                this.calMZ = arr[3];
                this.deltCN = arr[4];
                this.xcorr = arr[5];
                this.sp = arr[6];
                this.matchedIons = arr[7];
                this.predictedIons = arr[8];
                this.sequence = arr[9];
                this.status = arr[10];
        */
        try {
                lLineList = new ArrayList<String>();
                String [] arr = mLine.split("\t");
                this.xcorrRank = arr[1].trim();
                this.spRank = arr[2].trim();
                primaryRankInt = Integer.parseInt(xcorrRank);
                spRankInt = Integer.parseInt(spRank);
                this.calMZ = arr[3];
                this.deltCN = arr[4];
                this.xcorr = arr[5].trim();
                xcorrValue = Double.parseDouble(xcorr);
                this.sp = arr[6];
                //spValue = Double.parseDouble(sp);
                this.matchedIons = arr[7];
                this.predictedIons = arr[8];
                this.sequence = arr[9].trim();
                this.status = arr[10];
        } catch (Exception e) {

            System.err.println("Mline: " + mLine);
            e.printStackTrace();
        }
    }

    public String getFirstLLine() {
        return lLineList.get(0);
    }
    public int getPrimaryRankInt() {
        return primaryRankInt;
    }
    public int getSpRankInt() {
        return spRankInt;
    }
    public MLine(String[] arr)
    {
        lLineList = new ArrayList<String>();

        this.xcorrRank = arr[1].trim();
        this.spRank = arr[2].trim();
        primaryRankInt = Integer.parseInt(xcorrRank);
        spRankInt = Integer.parseInt(spRank);
        mLine = xcorrRank + "\t" + spRank;
        this.calMZ = arr[3];
        this.deltCN = arr[4];
        this.xcorr = arr[5].trim();
        xcorrValue = Double.parseDouble(xcorr);
        this.sp = arr[6];
        //spValue = Double.parseDouble(sp);
        this.matchedIons = arr[7];
        this.predictedIons = arr[8];
        this.sequence = arr[9].trim();
        this.status = arr[10];
    }

    // return true if at least one of the L lines starts with loc, 
    // return false otherwise
    public boolean contains(String loc) {
        for(String l : lLineList) {
            if(l.startsWith(loc)) {
                return true;
            }
        }
        return false;
    }
    public boolean isReverseHit() {
        for(String l : lLineList) {
            if(l.startsWith("Rever")) {
                return true;
            }
        }
        return false;
    }

    public String getMLine() {
        return mLine;
    }
    public double getXcorrValue() {
        return xcorrValue;
    }
    public boolean lLineStartsWith(String s) {
        for(String l : lLineList) {
            if(l.startsWith(s)) {
                return true;
            }
        }
        return false;
    }
    public boolean addLLine(String lline)
    {
//System.out.println("ind addLLine: " + lline);
        return lLineList.add( lline.split("\t")[1] );
    }

    public Iterator<String> getLLine()
    {
        return lLineList.iterator();
    }

    public void setMLine(String mLine)
    {
        this.mLine = mLine;
    }

    public void setXcorrRank(String xcorrRank)
    {
        this.xcorrRank = xcorrRank;
    }

    public void setSpRank(String spRank)
    {
        this.spRank = spRank;
    }

    public void setCalMZ(String calMZ)
    {
        this.calMZ = calMZ;
    }

    public void setDeltCN(String deltCN)
    {
        this.deltCN = deltCN;
    }

    public void setDeltCN(float deltCN)
    {
        this.deltCN = String.valueOf(deltCN);
    }

    public void setXcorr(String xcorr)
    {
        this.xcorr = xcorr;
    }

    public void setSp(String sp)
    {
        this.sp = sp;
    }

    public void setMatchedIons(String matchedIons)
    {
        this.matchedIons = matchedIons;
    }

    public void setPredictedIons(String predictedIons)
    {
        this.predictedIons = predictedIons;
    }

    public void setSequence(String sequence)
    {
        this.sequence = sequence;
    }

    public void setStatus(String status)
    {
        this.status = status;
    }

    public String getXcorrRank()
    {
        return xcorrRank;
    }

    public int getXcorrRankInt()
    {
        return Integer.parseInt(xcorrRank);
    }
    public String getSpRank()
    {
        return spRank;
    }

    public String getCalMZ()
    {
        return calMZ;
    }

    public String getDeltCN()
    {
        return deltCN;
    }
    public double getDeltaCnValue()
    {
        return Double.parseDouble(deltCN);
    }

    public String getXcorr()
    {
        return xcorr;
    }

    public String getSp()
    {
        return sp;
    }
    public double getSpValue()
    {
        return Double.parseDouble(sp);
    }

    public String getMatchedIons()
    {
        return matchedIons;
    }

    public String getPredictedIons()
    {
        return predictedIons;
    }

    public String getSequence()
    {
        return sequence;
    }

    public String getStatus()
    {
        return status;
    }

}
