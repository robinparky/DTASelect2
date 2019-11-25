package proteinGrouper;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/8/19.
 */
import java.lang.String;
import java.util.*;


public class Protein
{
    private String locus;
    private String accessionDTASelectFormat;
    private String seqCount;//peptideNum
    private String spectrumCount;//SPEC_COUNT
    private String seqCoverage;
    private String length;
    private String molWt;
    private String pI;
    private String validation;
    private String description;
    private String heavySpec;
    private String lightSpec;
    private String listOfExpNames;
    private String listOfExpNamesWithLink;
    private String listOfSearchNames;
    private String listOfSeqCount;
    private String listOfSpecCount;
    private String listOfSeqCoverage;
    private double nsaf=0;
    private double empai=0;
    private String geneName = null;

    private boolean ampOnly = false;
    private ArrayList<Peptide> peptideList;

    public String getProteinLine() {
        return proteinLine;
    }

    public void setProteinLine(String proteinLine) {
        this.proteinLine = proteinLine;
    }

    private String proteinLine="";
    private boolean isProblematic = false;
    private ProteinGroup proteinGroup;

    private static int locusIndex = -1;
    private static int seqCountIndex = -1;
    private static int spectrumCountIndex = -1;
    private static int seqCoverageIndex = -1;
    private static int lengthIndex = -1;
    private static int molWtIndex = -1;
    private static int pIIndex = -1;
    private static int validationIndex = -1;
    private static int descriptionIndex = -1;
    private static int heavySpecIndex = -1;
    private static int lightSpecIndex = -1;
    private static int listOfExpNamesIndex = -1;
    private static int listOfSearchNamesIndex = -1;
    private static int listOfSeqCountIndex = -1;
    private static int listOfSpecCountIndex = -1;
    private static int listOfSeqCoverageIndex = -1;
    private static int nsafIndex = -1;
    private static int empaiIndex = -1;

    private String compositeRatio;
    private String reverseCompositeRatio;
    private String averageRatio;
    private String averageRatioReverse;
    private String standardDeviation;
    private String standardDeviationReverse;
    private String areaRatio;

    private String listOfAverageRatio;
    private String listOfStandardDeviation;
    private String listOfAreaRatio;
    private String sequence;

    public Protein(String proteinLine) throws ArrayIndexOutOfBoundsException
    {
        this( proteinLine.split("\t") );
        this.proteinLine = proteinLine;
    }



    public Protein(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        peptideList = new ArrayList<Peptide>();

        init(strArr);
    }
    public Protein()
    {
        peptideList = new ArrayList<Peptide>();


    }
    //The setProteinIndexess is created for the DTASelectFileterPeptideFilegenerarotr to set the indexes..
    public void setProteinIndexes(String[] strArr)
    {
        this.setLocus( strArr[0] );
        this.setSeqCount(strArr[1]);
        this.setSpectrumCount(strArr[2]);
        this.setSeqCoverage(strArr[3]);
        this.setLength(strArr[4]);
        this.setMolWt(strArr[5]);
        this.setPI(strArr[6]);
        this.setValidation(strArr[7] );
        this.setDescription(strArr[8]);

            /*
            //remove previous lines after testing following lines.  Robin

            if(locusIndex>=0) locus = strArr[locusIndex];
            if(seqCountIndex>=0) seqCount = strArr[seqCountIndex];
            if(spectrumCountIndex>=0) spectrumCount = strArr[spectrumCountIndex];
            if(seqCoverageIndex>=0) seqCoverage = strArr[seqCoverageIndex];
            if(lengthIndex>=0) length = strArr[lengthIndex];
            if(molWtIndex>=0) molWt = strArr[molWtIndex];
            if(heavySpecIndex>=0) heavySpec = strArr[heavySpecIndex];
            if(lightSpecIndex>=0) lightSpec = strArr[lightSpecIndex];
            if(listOfExpNamesIndex>=0) listOfExpNames = strArr[listOfExpNamesIndex];
            if(listOfSearchNamesIndex>=0) listOfSearchNames = strArr[listOfSearchNamesIndex];
            if(listOfSeqCountIndex>=0) listOfSeqCount = strArr[listOfSeqCountIndex];
            if(listOfSpecCountIndex>=0) listOfSpecCount = strArr[listOfSpecCountIndex];
            if(listOfSeqCoverageIndex>=0) listOfSeqCoverage = strArr[listOfSeqCoverageIndex];
            if(nsafIndex>=0) nsaf = strArr[nsafIndex];
            if(empaiIndex>=0) empai = strArr[empaiIndex];
            if(pIIndex>=0) pI = strArr[pIIndex];
            if(validationIndex>=0) validation = strArr[validationIndex]; //this.setValidation(strArr[7] );
            if(descriptionIndex>=0) description = strArr[descriptionIndex];
            */

    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }
    public void setPeptideList(ArrayList<Peptide> list) {
        peptideList = list;
    }
    public ArrayList<Peptide> getPeptideList() {
        return peptideList;
    }
    public void setProteinGroup(ProteinGroup pg) {
        proteinGroup = pg;
    }

    public ProteinGroup getProteinGroup() {
        return proteinGroup;
    }

    public boolean equals(Object o) {

        Protein p = (Protein) o;
        return p == null? false : getAccession().equals(p.getAccession());
    }
    public int hashCode() {
        return getAccession().hashCode();
    }

    public String getAccession() {
        return getLocus();
    }
    public String getAccessionWithoutVersion() {
        return Fasta.getAccessionWithNoVersion(getLocus());
    }


    public void setElement(String[] strArr)
    {

        init(strArr);
    }

    public boolean addPeptide(Peptide peptide)
    {
        return peptideList.add(peptide);
    }

    public void removePeptide(Peptide peptide)
    {
        peptideList.remove(peptide);
    }

    public Iterator<Peptide> getPeptides()
    {
        return peptideList.iterator();
    }


    private void init(String[] strArr) throws ArrayIndexOutOfBoundsException
    {

        try {

            this.setLocus( Fasta.getAccession(strArr[this.getLocusIndex()]) );
            this.setAccessionDTASelectFormat(strArr[this.getLocusIndex()]);
            this.setSeqCount(strArr[this.getSeqCountIndex()]);
            this.setSpectrumCount(strArr[this.getSpectrumCountIndex()]);
            this.setSeqCoverage(strArr[this.getSeqCoverageIndex()]);
            this.setLength(strArr[this.getLengthIndex()]);
            this.setMolWt(strArr[this.getMolWtIndex()]);
            this.setPI(strArr[this.getpIIndex()]);
            this.setValidation(strArr[this.getValidationIndex()]);
            this.setDescription(strArr[this.getDescriptionIndex()]);


            if(this.getNsafIndex()>=0) {
                if(!"NA".equals(strArr[this.getNsafIndex()]) && !"NA".equals(strArr[this.getEmpaiIndex()])) {
                    this.setNsaf(Double.parseDouble(strArr[this.getNsafIndex()]));
                    this.setEmpai(Double.parseDouble(strArr[this.getEmpaiIndex()]));
                }

            }
            if(strArr[this.getHeavySpecIndex()]!=null && strArr[this.getLightSpecIndex()]!=null){
                this.setHeavySpec(strArr[this.getHeavySpecIndex()]);
                this.setLightSpec(strArr[this.getLightSpecIndex()]);
            }





        }
        catch(ArrayIndexOutOfBoundsException ex) {
            isProblematic = true;
            //throw new ArrayIndexOutOfBoundsException("Error : Mal formed DTASelect-filter.txt file in protein line: " + strArr[0]);
        }
    }

    /*
    private void init(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        try {

            this.setLocus( Fasta.getAccession(strArr[0]) );
            this.setSeqCount(strArr[1]);
            this.setSpectrumCount(strArr[2]);
            this.setSeqCoverage(strArr[3]);
            this.setLength(strArr[4]);
            this.setMolWt(strArr[5]);
            this.setPI(strArr[6]);
            this.setValidation(strArr[7]);
            this.setDescription(strArr[8]);
		if(strArr[9]!=null && strArr[10]!=null){
			this.setHeavySpec(strArr[9]);
			this.setLightSpec(strArr[10]);
		}
        }
        catch(ArrayIndexOutOfBoundsException ex) {
            isProblematic = true;
            //throw new ArrayIndexOutOfBoundsException("Error : Mal formed DTASelect-filter.txt file in protein line: " + strArr[0]);
        }
    }*/

    public boolean isProblematic() {
        return isProblematic;
    }
    public String getLocus()
    {
        return locus;
    }

    public String getSeqCount()
    {
        return seqCount;
    }

    public String getSpectrumCount()
    {
        return spectrumCount;
    }

    public String getSeqCoverage()
    {
        return seqCoverage;
    }

    public String getLength()
    {
        return length;
    }

    public int getLengthValue()
    {
        return Integer.parseInt(length);
    }
    public String getMolWt()
    {
        return molWt;
    }

    public String getPI()
    {
        return pI;
    }

    public String getValidation()
    {
        return validation;
    }

    public String getDescription()
    {
        return description;
    }

    public void setLocus(String locus)
    {
        this.locus = locus;
    }

    public void setSeqCount(String seqCount)
    {
        this.seqCount = seqCount;
    }

    public void setSpectrumCount(String spectrumCount)
    {
        this.spectrumCount = spectrumCount;
    }

    public void setSeqCoverage(String seqCoverage)
    {
        this.seqCoverage = seqCoverage;
    }

    public void setLength(String length)
    {
        this.length = length;
    }

    public void setMolWt(String molWt)
    {
        this.molWt = molWt;
    }

    public void setPI(String pI)
    {
        this.pI = pI;
    }

    public void setValidation(String validation)
    {
        this.validation = validation;
    }

    public void setDescription(String description)
    {

        this.description = description;
        if(description.contains("Gene_Symbol=") ||description.contains("GN=") )
        {
            String[] temp = null;
            if(description.contains("Gene_Symbol"))
                temp = description.split("Gene_Symbol=");
            else if(description.contains("GN"))
                temp = description.split("GN=");


            String[] temp1 = temp[1].split(" ");

            if("-".equals(temp1[0])) {
                setGeneName("N/A");
                return;
            }

            String[] temp2 = temp1[0].split(";");

            setGeneName( removeRef(temp2[0]) );
        }
    }
    public void setHeavySpec(String hs){
        this.heavySpec = hs;
    }
    public void setLightSpec(String ls){
        this.lightSpec = ls;
    }
    public void setListOfSearchNames(String lsn){
        this.listOfSearchNames = lsn;
    }
    public void setListOfExpNames(String len){
        this.listOfExpNames = len;
    }
    public void setListOfSeqCount(String lsc){
        this.listOfSeqCount = lsc;
    }
    public void setListOfSeqCoverage(String lscv){
        this.listOfSeqCoverage = lscv;
    }
    public void setListOfSpecCount(String lspc){
        this.listOfSpecCount = lspc;
    }
    public String getHeavySpec(){
        return heavySpec;
    }
    public String getLightSpec(){
        return lightSpec;
    }
    public String getListOfSearchNames(){
        return listOfSearchNames;
    }
    public String getListOfExpNames(){
        return listOfExpNames;
    }
    public String getListOfSeqCount(){
        return listOfSeqCount;
    }
    public String getListOfSeqCoverage(){
        return listOfSeqCoverage;
    }
    public String getListOfSpecCount(){
        return listOfSpecCount;
    }
    public int getPeptideSize()
    {
        return peptideList.size();
    }
    public int getNumPeptides()
    {
        return peptideList.size();
    }



    public boolean getAmpOnly(){
        return ampOnly;
    }
    public void setAmpOnly(boolean ao){
        ampOnly = ao;
    }

    public static void setFeatureIndices(String features) {
        clearIndex();
        String [] contents = features.split("\t");
//Locus

        for(int i = 0; i < contents.length; i++) {
            String s = contents[i].trim();
            locusIndex = s.startsWith("Locus")? i :locusIndex;
            seqCountIndex = s.startsWith("Sequence Count")? i :seqCountIndex;
            spectrumCountIndex = s.startsWith("Spectrum Count")? i : spectrumCountIndex;
            seqCoverageIndex = s.startsWith("Sequence Coverage")? i :seqCoverageIndex;
            lengthIndex = s.startsWith("Length")? i :lengthIndex;
            molWtIndex = s.startsWith("MolWt")? i :molWtIndex;
            pIIndex = s.startsWith("pI")? i :pIIndex;
            validationIndex = s.startsWith("Validation Status")? i :validationIndex;
            descriptionIndex = s.startsWith("Descriptive Name")? i :descriptionIndex;
            heavySpecIndex = s.startsWith("HRedundancy")? i :heavySpecIndex;
            lightSpecIndex = s.startsWith("LRedundancy")? i :lightSpecIndex;
            nsafIndex = s.startsWith("NSAF")? i :nsafIndex;
            empaiIndex = s.startsWith("EMPAI")? i :empaiIndex;
            descriptionIndex = s.startsWith("Description")? i :descriptionIndex;
            //listOfSeqCountIndex = s.startsWith("")? i :listOfSeqCountIndex;
            //listOfSpecCountIndex = s.startsWith("")? i :listOfSpecCountIndex;
            //listOfSeqCoverageIndex = s.startsWith("")? i :listOfSeqCoverageIndex;

        }

        //   System.out.println("aaa");
    }

    public static void clearIndex(){
        seqCountIndex = -1;
        spectrumCountIndex = -1;
        seqCoverageIndex = -1;
        lengthIndex = -1;
        molWtIndex = -1;
        pIIndex = -1;
        validationIndex = -1;
        descriptionIndex = -1;
        heavySpecIndex = -1;
        lightSpecIndex = -1;
        listOfExpNamesIndex = -1;
        listOfSearchNamesIndex = -1;
        listOfSeqCountIndex = -1;
        listOfSpecCountIndex = -1;
        listOfSeqCoverageIndex = -1;
        nsafIndex = -1;
        empaiIndex = -1;
        descriptionIndex = -1;

    }

    public String getpI() {
        return pI;
    }

    public void setpI(String pI) {
        this.pI = pI;
    }

    public double getNsaf() {
        return nsaf;
    }

    public void setNsaf(double nsaf) {
        this.nsaf = nsaf;
    }

    public double getEmpai() {
        return empai;
    }

    public void setEmpai(double empai) {
        this.empai = empai;
    }

    public boolean isIsProblematic() {
        return isProblematic;
    }

    public void setIsProblematic(boolean isProblematic) {
        this.isProblematic = isProblematic;
    }

    public static int getLocusIndex() {
        return locusIndex;
    }

    public static void setLocusIndex(int locusIndex) {
        Protein.locusIndex = locusIndex;
    }

    public static int getSeqCountIndex() {
        return seqCountIndex;
    }

    public static void setSeqCountIndex(int seqCountIndex) {
        Protein.seqCountIndex = seqCountIndex;
    }

    public static int getSpectrumCountIndex() {
        return spectrumCountIndex;
    }

    public static void setSpectrumCountIndex(int spectrumCountIndex) {
        Protein.spectrumCountIndex = spectrumCountIndex;
    }

    public static int getSeqCoverageIndex() {
        return seqCoverageIndex;
    }

    public static void setSeqCoverageIndex(int seqCoverageIndex) {
        Protein.seqCoverageIndex = seqCoverageIndex;
    }

    public static int getLengthIndex() {
        return lengthIndex;
    }

    public static void setLengthIndex(int lengthIndex) {
        Protein.lengthIndex = lengthIndex;
    }

    public static int getMolWtIndex() {
        return molWtIndex;
    }

    public static void setMolWtIndex(int molWtIndex) {
        Protein.molWtIndex = molWtIndex;
    }

    public static int getpIIndex() {
        return pIIndex;
    }

    public static void setpIIndex(int pIIndex) {
        Protein.pIIndex = pIIndex;
    }

    public static int getValidationIndex() {
        return validationIndex;
    }

    public static void setValidationIndex(int validationIndex) {
        Protein.validationIndex = validationIndex;
    }

    public static int getDescriptionIndex() {
        return descriptionIndex;
    }

    public static void setDescriptionIndex(int descriptionIndex) {
        Protein.descriptionIndex = descriptionIndex;
    }

    public static int getHeavySpecIndex() {
        return heavySpecIndex;
    }

    public static void setHeavySpecIndex(int heavySpecIndex) {
        Protein.heavySpecIndex = heavySpecIndex;
    }

    public static int getLightSpecIndex() {
        return lightSpecIndex;
    }

    public static void setLightSpecIndex(int lightSpecIndex) {
        Protein.lightSpecIndex = lightSpecIndex;
    }

    public static int getListOfExpNamesIndex() {
        return listOfExpNamesIndex;
    }

    public static void setListOfExpNamesIndex(int listOfExpNamesIndex) {
        Protein.listOfExpNamesIndex = listOfExpNamesIndex;
    }

    public static int getListOfSearchNamesIndex() {
        return listOfSearchNamesIndex;
    }

    public static void setListOfSearchNamesIndex(int listOfSearchNamesIndex) {
        Protein.listOfSearchNamesIndex = listOfSearchNamesIndex;
    }

    public static int getListOfSeqCountIndex() {
        return listOfSeqCountIndex;
    }

    public static void setListOfSeqCountIndex(int listOfSeqCountIndex) {
        Protein.listOfSeqCountIndex = listOfSeqCountIndex;
    }

    public static int getListOfSpecCountIndex() {
        return listOfSpecCountIndex;
    }

    public static void setListOfSpecCountIndex(int listOfSpecCountIndex) {
        Protein.listOfSpecCountIndex = listOfSpecCountIndex;
    }

    public static int getListOfSeqCoverageIndex() {
        return listOfSeqCoverageIndex;
    }

    public static void setListOfSeqCoverageIndex(int listOfSeqCoverageIndex) {
        Protein.listOfSeqCoverageIndex = listOfSeqCoverageIndex;
    }

    public static int getNsafIndex() {
        return nsafIndex;
    }

    public static void setNsafIndex(int nsafIndex) {
        Protein.nsafIndex = nsafIndex;
    }

    public static int getEmpaiIndex() {
        return empaiIndex;
    }

    public static void setEmpaiIndex(int empaiIndex) {
        Protein.empaiIndex = empaiIndex;
    }
    public int checkNSAF()
    {
        if(nsafIndex == -1 )
        {
            return -1;
        }
        else
            return nsafIndex;
    }


    public static String removeRef(String gene) {

        int index = gene.indexOf("-Ref");

        if(index<=0) return gene;

        return gene.substring(0,index);

    }

    public String getListOfExpNamesWithLink() {
        return listOfExpNamesWithLink;
    }

    public void setListOfExpNamesWithLink(String listOfExpNamesWithLink) {
        this.listOfExpNamesWithLink = listOfExpNamesWithLink;
    }

    public String getAverageRatio() {
        return averageRatio;
    }

    public void setAverageRatio(String averageRatio) {
        this.averageRatio = averageRatio;
    }

    public String getStandardDeviation() {
        return standardDeviation;
    }

    public void setStandardDeviation(String standardDeviation) {
        this.standardDeviation = standardDeviation;
    }

    public String getAreaRatio() {
        return areaRatio;
    }

    public void setAreaRatio(String areaRatio) {
        this.areaRatio = areaRatio;
    }

    public String getListOfAverageRatio() {
        return listOfAverageRatio;
    }

    public void setListOfAverageRatio(String listOfAverageRatio) {
        this.listOfAverageRatio = listOfAverageRatio;
    }

    public String getListOfStandardDeviation() {
        return listOfStandardDeviation;
    }

    public void setListOfStandardDeviation(String listOfStandardDeviation) {
        this.listOfStandardDeviation = listOfStandardDeviation;
    }

    public String getListOfAreaRatio() {
        return listOfAreaRatio;
    }

    public void setListOfAreaRatio(String listOfAreaRatio) {
        this.listOfAreaRatio = listOfAreaRatio;
    }

    public String getCompositeRatio() {
        return compositeRatio;
    }

    public void setCompositeRatio(String compositeRatio) {
        this.compositeRatio = compositeRatio;
    }

    public String getReverseCompositeRatio() {
        return reverseCompositeRatio;
    }

    public void setReverseCompositeRatio(String reverseCompositeRatio) {
        this.reverseCompositeRatio = reverseCompositeRatio;
    }

    public String getAverageRatioReverse() {
        return averageRatioReverse;
    }

    public void setAverageRatioReverse(String averageRatioReverse) {
        this.averageRatioReverse = averageRatioReverse;
    }

    public String getStandardDeviationReverse() {
        return standardDeviationReverse;
    }

    public void setStandardDeviationReverse(String standardDeviationReverse) {
        this.standardDeviationReverse = standardDeviationReverse;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getAccessionDTASelectFormat() {
        return accessionDTASelectFormat;
    }

    public void setAccessionDTASelectFormat(String accessionDTASelectFormat) {
        this.accessionDTASelectFormat = accessionDTASelectFormat;
    }
}

