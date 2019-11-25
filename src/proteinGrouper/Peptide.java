package proteinGrouper;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/8/19.
 */
import java.util.*;

public class Peptide
{

    private boolean unique;
    private String fileName;
    private String xCorr;

    private double xCorrValue;
    private double deltCNValue;
    private double mhPlusValue;
    private double totalIntensityValue;
    private double calcMHplusValue;
    private int spRankValue;
    private double spScoreValue;
    private double zScoreValue;
    private double ionProportionValue;
    private int scanNumValue;
    private int chargeStateValue;
    private String kd;

    private String deltCN;
    private String mhPlus;
    private String calcMHplus;
    private String totalIntensity;
    private String spRank;
    private String spScore;
    private String zScore;
    private String ionProportion;
    private double ppm;
    private double retTime;
    private int redundancy;
    private String sequence;
    private String[] peptideLine;
    private String tmpStr;
    private String scanNum;
    private String conf;
    private String proteins;
    private String proteinDescription;
    private String chargeState;
    private  int uniqueIndex = -1;
    private  int scanNumIndex =-1;
    private  int xcorrIndex = -1;
    private  int dcnIndex = -1;
    private  int confIndex = -1;
    private  int mPlusHIndex = -1;
    private  int calcMassIndex = -1;
    private  int totalIntensityIndex = -1;
    private  int spRankIndex = -1;
    private  int spScoreIndex = -1;
    private  int ionProportionIndex = -1;
    private  int redundancyIndex = -1;
    private  int sequenceIndex = -1;
    private  int pIIndex = -1;
    private  int ppmIndex = -1;
    private  int retTimeIndex = -1;
    private  int zScoreIndex=-1;
    private boolean isDecoy=false;


    private int hashcode = -1;

    private List<Protein> proteinList = new ArrayList<Protein>();

    private boolean destroyedInCensusOut;

    public String peptideWholeLine;

    public int hashCode() {
        if(hashcode == -1) {
            // fileName contains the scan number
            hashcode = (getSequence() + fileName).hashCode();
            //hashcode = (fileName + scanNum).hashCode();
        }
        return hashcode;
    }
    public boolean equals(Object o) {
        Peptide p = (Peptide)o;
        return getSequence().equals(p.getSequence()) && fileName.equals(p.fileName);
    }

    // this function need to be changed if the format of getInfo() function changes
    public static String getInfoHeader() {
        return "Peptide\tXCorr\tDeltaCN\tMPlusH\tCalcMPlusH\tDeltaMass\tSpScore\tConfidence\tFileName";
    }
    public String getInfo() {
        StringBuffer sb = new StringBuffer(1000);
        sb.append(sequence + "\t" + xCorr + "\t" + deltCN);
        sb.append("\t" + mhPlus + "\t" + calcMHplus + "\t" + getDeltaMass());
        sb.append("\t" + spScore + "\t" + conf + "\t" + fileName);
        return sb.toString();
    }
    /* DTASelect 2.0 file */
    private void parseLine2() throws ArrayIndexOutOfBoundsException
    {
        if(scanNumIndex<0) System.out.println("---->ScanNumIndex"+scanNumIndex);
        scanNum = peptideLine[scanNumIndex];
        this.chargeState = scanNum.substring(scanNum.lastIndexOf(".")+1);
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);


        this.setUnique((peptideLine[uniqueIndex]).startsWith("*"));
        this.setFileName(peptideLine[scanNumIndex]);
        this.setXCorr(peptideLine[xcorrIndex]);
        this.setDeltCN(peptideLine[dcnIndex]);

        if(confIndex>0)
            this.setConf(peptideLine[confIndex]);
        this.setMhPlus(peptideLine[mPlusHIndex]);
        this.setCalcMHplus(peptideLine[calcMassIndex]);

        if(ppmIndex>=0)
            this.setPpm(peptideLine[ppmIndex]);

        if(retTimeIndex >=0) {
            this.setRetTime(peptideLine[retTimeIndex]);
        }

        this.setTotalIntensity(peptideLine[totalIntensityIndex]);
        this.setSpRank(peptideLine[spRankIndex]);
        if(spScoreIndex>0)
            this.setSpScore(peptideLine[spScoreIndex]);
        this.setIonProportion(peptideLine[ionProportionIndex]);
        this.setRedundancy(peptideLine[redundancyIndex]);
        this.setSequence(peptideLine[sequenceIndex]);

    }
    public  void clearIndex(){
        uniqueIndex = -1;
        scanNumIndex = -1;
        xcorrIndex = -1;
        dcnIndex = -1;
        confIndex = -1;
        mPlusHIndex = -1;
        calcMassIndex = -1;
        totalIntensityIndex = -1;
        spRankIndex = -1;
        spScoreIndex = -1;
        ionProportionIndex = -1;
        redundancyIndex = -1;
        sequenceIndex = -1;
        pIIndex = -1;
        ppmIndex = -1;

        zScoreIndex=-1;

    }
    public  void setFeatureIndex(String features) {
        clearIndex();
        //System.out.println("---indices>"+features);
        String [] contents = features.split("\t");
        for(int i = 0; i < contents.length; i++) {
            String s = contents[i].trim();
            uniqueIndex = s.startsWith("Uni")? i : uniqueIndex;
            scanNumIndex = s.startsWith("File")? i : scanNumIndex;
            xcorrIndex = s.startsWith("XC")? i : xcorrIndex;
            dcnIndex = (s.startsWith("DeltCN") || s.startsWith("DeltaCN"))? i : dcnIndex;
            confIndex = s.startsWith("Conf%")? i : confIndex;
            mPlusHIndex = s.startsWith("M")? i : mPlusHIndex;
            calcMassIndex = s.startsWith("CalcM")? i : calcMassIndex;
            totalIntensityIndex = s.startsWith("Total")? i : totalIntensityIndex;
            spRankIndex = s.startsWith("SpR")? i : spRankIndex;
            spScoreIndex = (s.startsWith("Prob")||s.startsWith("SpScore"))? i : spScoreIndex;
            ionProportionIndex = s.startsWith("IonP")? i : ionProportionIndex;
            redundancyIndex = s.startsWith("Red")? i : redundancyIndex;
            sequenceIndex = s.startsWith("Sequence")? i : sequenceIndex;
            pIIndex = s.startsWith("pI")? i : pIIndex;
            ppmIndex = s.startsWith("PPM")? i : ppmIndex;
            retTimeIndex = s.startsWith("RetTime")? i : retTimeIndex;
        }
//	if(scanNumIndex<0) System.out.println("-----ScanNumIndex>"+scanNumIndex);
    }

    public void copyFeatureIndex(Peptide peptide2)
    {

        uniqueIndex = peptide2.uniqueIndex;
        scanNumIndex = peptide2.scanNumIndex;
        xcorrIndex = peptide2.xcorrIndex;
        dcnIndex = peptide2.dcnIndex;
        confIndex = peptide2.confIndex;
        mPlusHIndex = peptide2.mPlusHIndex;
        calcMassIndex = peptide2.calcMassIndex;
        totalIntensityIndex = peptide2.totalIntensityIndex;
        spRankIndex =peptide2.spRankIndex;
        spScoreIndex = peptide2.spScoreIndex;
        ionProportionIndex = peptide2.ionProportionIndex;
        redundancyIndex = peptide2.redundancyIndex;
        sequenceIndex = peptide2.sequenceIndex;
        pIIndex = peptide2.pIIndex;
        ppmIndex = peptide2.ppmIndex;
        retTimeIndex = peptide2.retTimeIndex;

    }

    //For DTASelect version 2
    public Peptide() {

    }


    public Peptide(String peptideLine, boolean isV2, boolean isV1, Peptide peptide) {
        this(peptideLine.split("\t"), isV2, isV1,peptide);
        this.peptideWholeLine = peptideLine;
    }

    public Peptide(String[] peptideLine, boolean isV2, boolean isV1,Peptide peptide)
    {
        this.peptideLine = peptideLine;
        if(peptide!=null)copyFeatureIndex(peptide);
	    /*
		  ppmIndex = 7;
		  totalIntensityIndex = 8;
		  spRankIndex = 9;
		  spScoreIndex = 10;
		  ionProportionIndex = 12;
		  redundancyIndex = 13;
		  sequenceIndex = 14;
		  pIIndex = 11;
	     */
        if(isV1)
            parseLine();
        else
            parseLine2();
        //else
    }


    /* DTASelect file */
    private void parseLine() throws ArrayIndexOutOfBoundsException
    {
        scanNum = peptideLine[1];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        //this.setUnique(!"".equals(peptideLine[0]));
        this.setUnique((peptideLine[0]).startsWith("*"));
        this.setFileName(peptideLine[1]);
        this.setXCorr(peptideLine[2]);
        this.setDeltCN(peptideLine[3]);
        this.setMhPlus(peptideLine[4]);
        this.setCalcMHplus(peptideLine[5]);
        this.setTotalIntensity(peptideLine[6]);
        this.setSpRank(peptideLine[7]);
        this.setSpScore(peptideLine[8]);
        this.setIonProportion(peptideLine[9]);
        this.setRedundancy(peptideLine[10]);
        this.setSequence(peptideLine[11]);
    }
    public float getPi() {
        return Float.parseFloat(peptideLine[pIIndex]);
    }
    public float getDeltaMass() {
        return ppmIndex != -1 ? Float.parseFloat(peptideLine[ppmIndex]) : 1000;
    }
    public String getChargeState()
    {
        return this.fileName.substring( this.fileName.lastIndexOf(".") + 1 );
    }
    public int getChargeStateValue()
    {
        return Integer.parseInt(getChargeState());
    }
    public void setChargeState(String cs){
        chargeStateValue = Integer.parseInt(cs);
        this.chargeState=cs;
    }
    public String getLoScan()
    {
        tmpStr = this.fileName.substring( this.fileName.indexOf(".") +1 );
        return tmpStr.substring(0, tmpStr.indexOf(".") );
    }

    public double getPpm()
    {
        return this.ppm;
    }
    public String getFileName()
    {
        return fileName.substring(0, fileName.indexOf("."));
    }
    public String getFullFileName()
    {
        return fileName;
    }

    public String getXCorr()
    {
        return xCorr;
    }

    public double  getXCorrValue()
    {
        return Double.parseDouble(xCorr);
    }
    public String getDeltCN()
    {
        return deltCN;
    }
    public double  getDeltCnValue()
    {
        return Double.parseDouble(deltCN);
    }

    public String getMhPlus()
    {
        return mhPlus;
    }

    public String getSpRank()
    {
        return spRank;
    }

    public String getSpScore()
    {
        return spScore;
    }
    public String getZScore()
    {
        return zScore;
    }

    public double getSpScoreValue()
    {
        return Double.parseDouble(spScore);
    }

    public String getIonProportion()
    {
        return ionProportion;
    }

    public int getRedundancy()
    {
        return redundancy;
    }

    public int getSpectralCount()
    {
        return redundancy;
    }
    public String getSequence()
    {
        return sequence;
    }
    // return the peptide sequence without leading and tailing residues
    public String getMidSeq()
    {
        int lastindex = sequence.length() - 2;
        return sequence.substring(2, lastindex);
    }

    public boolean isUnique()
    {
        return unique;
    }

    public double getCalcMHplusValue()
    {
        return Double.parseDouble(calcMHplus);
    }
    public String getCalcMHplus()
    {
        return calcMHplus;
    }

    public String getTotalIntensity()
    {
        return totalIntensity;
    }

    public String getScanNum()
    {
        return scanNum;
    }
    public int getScanNumber()
    {
        return Integer.parseInt(scanNum);
    }

    public void setFileName(String fileName)
    {
        this.fileName = fileName;
    }
    public void setPpm(String ppm)
    {
        setPpm( Double.parseDouble(ppm) );
    }
    public void setPpm(double ppm)
    {
        this.ppm=ppm;
    }
    public void setXCorr(String xCorr)
    {
        this.xCorr = xCorr;
        xCorrValue = Double.parseDouble(xCorr);
    }

    public void setDeltCN(String deltCN)
    {
        deltCNValue= Double.parseDouble(deltCN);
        this.deltCN = deltCN;
    }

    public void setMhPlus(String mhPlus)
    {
        this.mhPlus = mhPlus;
        mhPlusValue = Double.parseDouble(mhPlus);
    }

    public void setSpRank(String spRank)
    {
        spRankValue = Integer.parseInt(spRank);
        this.spRank = spRank;
    }

    public void setSpScore(String spScore)
    {
        spScoreValue = Double.parseDouble(spScore);
        this.spScore = spScore;
    }
    public void setZScore(String zScore)
    {
        zScoreValue = Double.parseDouble(zScore);
        this.zScore = zScore;
    }

    public void setIonProportion(String ionProportion)
    {
        ionProportionValue = Double.parseDouble(ionProportion);
        this.ionProportion = ionProportion;
    }

    public void setRedundancy(String redundancy) {
        setRedundancy(Integer.parseInt(redundancy));
    }

    public void setRedundancy(int redundancy)
    {
        this.redundancy = redundancy;
    }

    public void setSequence(String sequence)
    {
        this.sequence = sequence;
    }

    public void setUnique(boolean unique)
    {
        this.unique = unique;
    }

    public void setCalcMHplus(String calcMHplus)
    {
        this.calcMHplus = calcMHplus;
    }

    public void setTotalIntensity(String totalIntensity)
    {
        totalIntensityValue = Double.parseDouble(totalIntensity);
        this.totalIntensity = totalIntensity;
    }

    public void setScanNum(String scanNum)
    {
        scanNumValue = Integer.parseInt(scanNum);
        this.scanNum = scanNum;
    }

    public String getConf() {
        return (null==conf)?"":conf;
    }
    public double  getConfValue() {
        if(conf == null) {
            return 0;
        }
        return Double.parseDouble(conf);
    }

    public void setConf(String conf) {
        this.conf = conf;
    }

    public String[] getPeptideLine()
    {
        return peptideLine;
    }
    public String[] getPeptideLineWithCleanPeptide()
    {
        String peptide = peptideLine[peptideLine.length-1];
        int start = peptide.indexOf('.');
        int end = peptide.lastIndexOf('.');
        String cleanPeptide = peptide.substring(start,end);

        return peptideLine;
    }

    public void setPeptideLine(String[] peptideLine) {
        this.peptideLine = peptideLine;
    }


    public void addProtein(Protein protein) {
        this.proteinList.add(protein);

    }

    public void setProtein(String protein)
    {
        this.proteins = protein;
    }

    public String getProteins()
    {
        return this.proteins;
    }


    public void setProteinDescription(String proteinDescription)
    {
        this.proteinDescription = proteinDescription;
    }

    public String getProteinDescription()
    {
        return this.proteinDescription;
    }



    public List<Protein> getProteinList() {
        return proteinList;
    }

    public void setProteinList(List<Protein> proteinList) {
        this.proteinList = proteinList;
    }

    public boolean isIsDecoy() {
        return isDecoy;
    }

    public void setIsDecoy(boolean isDecoy) {
        this.isDecoy = isDecoy;
    }

    public String getxCorr() {
        return xCorr;
    }

    public void setxCorr(String xCorr) {
        this.xCorr = xCorr;
    }


    public boolean isDestroyedInCensusOut() {
        return destroyedInCensusOut;
    }

    public void setDestroyedInCensusOut(boolean destroyedInCensusOut) {
        this.destroyedInCensusOut = destroyedInCensusOut;
    }

    public String getPeptideWholeLine() {
        return peptideWholeLine;
    }

    public void setPeptideWholeLine(String peptideWholeLine) {
        this.peptideWholeLine = peptideWholeLine;
    }

    public double getxCorrValue() {
        return xCorrValue;
    }

    public void setxCorrValue(double xCorrValue) {
        this.xCorrValue = xCorrValue;
    }

    public double getDeltCNValue() {
        return deltCNValue;
    }

    public void setDeltCNValue(double deltCNValue) {
        this.deltCNValue = deltCNValue;
    }

    public double getMhPlusValue() {
        return mhPlusValue;
    }

    public void setMhPlusValue(double mhPlusValue) {
        this.mhPlusValue = mhPlusValue;
    }

    public double getTotalIntensityValue() {
        return totalIntensityValue;
    }

    public void setTotalIntensityValue(double totalIntensityValue) {
        this.totalIntensityValue = totalIntensityValue;
    }

    public void setCalcMHplusValue(double calcMHplusValue) {
        this.calcMHplusValue = calcMHplusValue;
    }

    public int getSpRankValue() {
        return spRankValue;
    }

    public void setSpRankValue(int spRankValue) {
        this.spRankValue = spRankValue;
    }

    public void setSpScoreValue(double spScoreValue) {
        this.spScoreValue = spScoreValue;
    }

    public double getzScoreValue() {
        return zScoreValue;
    }

    public void setzScoreValue(double zScoreValue) {
        this.zScoreValue = zScoreValue;
    }

    public double getIonProportionValue() {
        return ionProportionValue;
    }

    public void setIonProportionValue(double ionProportionValue) {
        this.ionProportionValue = ionProportionValue;
    }

    public int getScanNumValue() {
        return scanNumValue;
    }

    public void setScanNumValue(int scanNumValue) {
        this.scanNumValue = scanNumValue;
    }

    public void setChargeStateValue(int chargeStateValue) {
        this.chargeStateValue = chargeStateValue;
    }

    public String getzScore() {
        return zScore;
    }

    public void setzScore(String zScore) {
        this.zScore = zScore;
    }

    public String getTmpStr() {
        return tmpStr;
    }

    public void setTmpStr(String tmpStr) {
        this.tmpStr = tmpStr;
    }

    public void setProteins(String proteins) {
        this.proteins = proteins;
    }

    public int getUniqueIndex() {
        return uniqueIndex;
    }

    public void setUniqueIndex(int uniqueIndex) {
        this.uniqueIndex = uniqueIndex;
    }

    public int getScanNumIndex() {
        return scanNumIndex;
    }

    public void setScanNumIndex(int scanNumIndex) {
        this.scanNumIndex = scanNumIndex;
    }

    public int getXcorrIndex() {
        return xcorrIndex;
    }

    public void setXcorrIndex(int xcorrIndex) {
        this.xcorrIndex = xcorrIndex;
    }

    public int getDcnIndex() {
        return dcnIndex;
    }

    public void setDcnIndex(int dcnIndex) {
        this.dcnIndex = dcnIndex;
    }

    public int getConfIndex() {
        return confIndex;
    }

    public void setConfIndex(int confIndex) {
        this.confIndex = confIndex;
    }

    public int getmPlusHIndex() {
        return mPlusHIndex;
    }

    public void setmPlusHIndex(int mPlusHIndex) {
        this.mPlusHIndex = mPlusHIndex;
    }

    public int getCalcMassIndex() {
        return calcMassIndex;
    }

    public void setCalcMassIndex(int calcMassIndex) {
        this.calcMassIndex = calcMassIndex;
    }

    public int getTotalIntensityIndex() {
        return totalIntensityIndex;
    }

    public void setTotalIntensityIndex(int totalIntensityIndex) {
        this.totalIntensityIndex = totalIntensityIndex;
    }

    public int getSpRankIndex() {
        return spRankIndex;
    }

    public void setSpRankIndex(int spRankIndex) {
        this.spRankIndex = spRankIndex;
    }

    public int getSpScoreIndex() {
        return spScoreIndex;
    }

    public void setSpScoreIndex(int spScoreIndex) {
        this.spScoreIndex = spScoreIndex;
    }

    public int getIonProportionIndex() {
        return ionProportionIndex;
    }

    public void setIonProportionIndex(int ionProportionIndex) {
        this.ionProportionIndex = ionProportionIndex;
    }

    public int getRedundancyIndex() {
        return redundancyIndex;
    }

    public void setRedundancyIndex(int redundancyIndex) {
        this.redundancyIndex = redundancyIndex;
    }

    public int getSequenceIndex() {
        return sequenceIndex;
    }

    public void setSequenceIndex(int sequenceIndex) {
        this.sequenceIndex = sequenceIndex;
    }

    public int getpIIndex() {
        return pIIndex;
    }

    public void setpIIndex(int pIIndex) {
        this.pIIndex = pIIndex;
    }

    public int getPpmIndex() {
        return ppmIndex;
    }

    public void setPpmIndex(int ppmIndex) {
        this.ppmIndex = ppmIndex;
    }

    public int getzScoreIndex() {
        return zScoreIndex;
    }

    public void setzScoreIndex(int zScoreIndex) {
        this.zScoreIndex = zScoreIndex;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public void setDecoy(boolean decoy) {
        isDecoy = decoy;
    }

    public int getHashcode() {
        return hashcode;
    }

    public void setHashcode(int hashcode) {
        this.hashcode = hashcode;
    }
    public double getRetTime() {
        return retTime;
    }

    public void setRetTime(double retTime) {
        this.retTime = retTime;
    }

    public void setRetTime(String retTime) {
        try {
            setRetTime(Double.parseDouble(retTime));
        } catch (Exception e) {
            // TODO: handle exception
        }
    }

    public int getRetTimeIndex() {
        return retTimeIndex;
    }
    public void setRetTimeIndex(int retTimeIndex) {
        this.retTimeIndex = retTimeIndex;
    }

    public String getKd() {
        return kd;
    }

    public void setKd(String kd) {
        this.kd = kd;
    }
}

