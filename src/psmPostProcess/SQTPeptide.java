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
public class SQTPeptide
{
    private String sLine;
    private ArrayList<MLine> mLineList;
    private String loScan;
    private String hiScan;
    private String chargeState;
    private String timeToProcess;
    private String serverName;
    private String calMZ;
    private String totalIntensity;
    private String lowestSp;
    private String numSeq;
    private MLine topHit;
    public static final double MASSDIFFC12C13 = 1.003354826;
    public static final double HALFMASSDIFFC12C13 = 1.003354826/2;
    
    public SQTPeptide(String sLine)
    {
        mLineList = new ArrayList<MLine>();
        this.sLine = sLine;


        String[] arr = this.sLine.split("\t");
        this.loScan=arr[1].trim();
        this.hiScan=arr[2];
        this.chargeState=arr[3].trim();
        this.timeToProcess=arr[4].trim();
        this.serverName=arr[5];
        this.calMZ=arr[6];
        this.totalIntensity=arr[7];
        numSeq = arr[9].trim();
    }
    public double getDeltaMass() {
        return getDeltaMass(getTopHit());
    }
    // return delta mass in 
    public double getDeltaMass(MLine m) {
        double calMass = Double.parseDouble(m.getCalMZ());
        double observedMass = Double.parseDouble(calMZ);
        double deltaMass = observedMass - calMass;
        /*
        if(deltaMass > 0) {
            while(deltaMass > MASSDIFFC12C13) {
                deltaMass -= MASSDIFFC12C13;
            }
        } else if(deltaMass < 0) {
             
        }
        deltaMass = deltaMass/calMass*1000000;
        */
        return deltaMass;
    }
   
    public MLine getMLineStartsWith(String locStart) {
        for(MLine m : mLineList) {
            if(m.contains(locStart)) {
                return m;
            }
        }
        return null;
    }
    public int getNumMlines() {
        return mLineList.size();
    }
    public double getTScore() {
        double tScore = 0;
        if(mLineList.size() > 1) {
            for(MLine m : mLineList) {
                if(m.getXcorrRankInt() == 1) {
                    tScore = Double.parseDouble(m.getDeltCN());
                    break; 
                }
            }
        }
        return tScore;
    }
    public double getDeltaCn() {
        double deltaCn = 0;
        if(mLineList.size() > 0) {
            for(MLine m : mLineList) {
                if(m.getXcorrRankInt() == 2) {
                    deltaCn = Double.parseDouble(m.getDeltCN());
                    break; 
                }
            }
        }
        return deltaCn;
    }
    public MLine getSpRankOneHit() {
        int numM = mLineList.size();
        if(numM > 0) {
            return mLineList.get(numM-1);
        }
        return null;
    }
    public boolean isReverseHit() {
        if(getTopHit() == null) {
            return false;
        } 
        return getTopHit().isReverseHit();
    }
    public MLine getTopHit() {
        if(getNumMlines() > 0) {
            topHit = mLineList.get(0); 
        } else {
            System.out.println("numMline is 0");
            return null;
        }
        return topHit;
    }
    // return true if at least one of the L line starts with locStart
    // return false otherwise
    public boolean contains(String locStart) {
        for(MLine m : mLineList) {
           if(m.lLineStartsWith(locStart))  {
               return true;
           } 
        }
        return false;
    }
    public boolean topHitStartsWith(String locStart) {
        for(MLine m : mLineList) {
           if(Integer.parseInt(m.getXcorrRank()) == 1 && m.lLineStartsWith(locStart))  {
               return true;
           } 
        }
        return false;
    }
    public boolean addMLine(MLine mline)
    {
        return mLineList.add(mline);
    }

    public Iterator<MLine> getMLine()
    {
        return mLineList.iterator();
    }

    public void setSLine(String sLine)
    {
        this.sLine = sLine;
    }

    public void setLoScan(String loScan)
    {
        this.loScan = loScan;
    }

    public void setHiScan(String hiScan)
    {
        this.hiScan = hiScan;
    }

    public void setChargeState(String chargeState)
    {
        this.chargeState = chargeState;
    }

    public void setTimeToProcess(String timeToProcess)
    {
        this.timeToProcess = timeToProcess;
    }

    public void setServerName(String serverName)
    {
        this.serverName = serverName;
    }

    public void setCalMZ(String calMZ)
    {
        this.calMZ = calMZ;
    }

    public void setTotalIntensity(String totalIntensity)
    {
        this.totalIntensity = totalIntensity;
    }

    public void setLowestSp(String lowestSp)
    {
        this.lowestSp = lowestSp;
    }

    public void setNumSeq(String numSeq)
    {
        this.numSeq = numSeq;
    }

    public String getSLine()
    {
        return sLine;
    }

    public String getLoScan()
    {
        return loScan;
    }

    public String getHiScan()
    {
        return hiScan;
    }

    public String getChargeState()
    {
        return chargeState;
    }

    public int getChargeStateInt()
    {
        return Integer.parseInt(chargeState);
    }
    public String getTimeToProcess()
    {
        return timeToProcess;
    }

    public String getServerName()
    {
        return serverName;
    }

    public String getCalMZ()
    {
        return calMZ;
    }

    public String getTotalIntensity()
    {
        return totalIntensity;
    }

    public String getLowestSp()
    {
        return lowestSp;
    }

    public String getNumSeq()
    {
        return numSeq;
    }
}
