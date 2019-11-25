package proteinGrouper;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/8/19.
 */

import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.FileReader;
import java.io.RandomAccessFile;

public class DTASelectFilterReader
{
    private final long READ_FROM_THE_END = 500; //position from the end

    private String dbFileName;
    private String dbFilePathAndName;
    //    private InputStreamReader reader;
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private String criteria;
    private int unfilteredProteinNum;
    private int redundantProteinNum;
    private int nonRedundantProteinNum;
    private int unfilteredPeptideNum;
    private int redundantPeptideNum;
    private int nonRedundantPeptideNum;
    private double proteinFP;
    private double peptideFP;
    private double spectrumFP;
    private String fileType = null;
    private Peptide peptide = new Peptide();

    private boolean version2;  //DTASelect version
    private boolean version1;  //DTASelect version
    private boolean isMascot = false;
    //read total ipeptide number
    public int getTotalPeptideNumber() throws IOException
    {
        int totalPeptideCount=0;

        DTASelectFilterReader dtaReader = new DTASelectFilterReader(fileName);

        for (Iterator<Protein> itr1 = dtaReader.getProteins(); itr1.hasNext(); ) {
            Protein protein = itr1.next();

            totalPeptideCount += protein.getPeptideSize();
        }

        return totalPeptideCount;

    }

    public static void main(String args[]) throws IOException
    {
        //DTASelectFilterReader reader = new DTASelectFilterReader("/data/2/rpark/ip2_data//aslanian/frog/130831_2013_09_01_09_18659/search/projects2014_05_13_11_63282/DTASelect-filter.txt");
        DTASelectFilterReader reader =
                //new DTASelectFilterReader("/data/2/rpark/ip2_data/benstein/Mammalian_LKB1_AMPK_Interactome/20140824_flag_IP_S1FL_wt_FM_293T_2014_08_25_10_26833/search/projects2014_08_26_17_68592/phospho/DTASelect-filter.txt");
                //  new DTASelectFilterReader("/data/2/rpark/ip2_data/yrc/Mahjoub/MMEB_heat_08292014_2014_10_07_15_27782/search/projects2014_10_10_11_71308/DTASelect-filter.txt");

                //new DTASelectFilterReader("/data/2/rpark/ip2_data/yrc/Noriyuki/YFP_AFP3_3_2014_05_19_10_24980/search/projects2014_05_19_10_63655/phospho/DTASelect-filter.txt");
//	DTASelectFilterReader reader = new DTASelectFilterReader("/data/2/rpark/ip2_data/rpark/test4/tiny2_2010_09_29_11_2084/search/projects2013_09_27_14_50137/test/DTASelect-filter.txt");

                //new DTASelectFilterReader("/data/2/rpark/ip2_data//xmhan/Liwei/Mono_Tri_2015_07_02_13_33164/search/projects2015_07_02_14_82787//DTASelect-filter.txt");
                new DTASelectFilterReader("/home/rpark/titus/1311HumanVsCompil/compil/DTASelect-filter.txt0.5");


        //Peptide ipeptide;
        // ArrayList<Protein> aList = new ArrayList<Protein>();



        Iterator<Protein> pitr = reader.getProteins();
//			int matchNum = 0;
//
//

        int proteinGroup=0;
        boolean isUnitProt=false;

	/*
        for (Iterator<com.ipa.ip2.util.dtaselect.Protein> itr = pitr; itr.hasNext(); ) {
            com.ipa.ip2.util.dtaselect.Protein protein = itr.next();
//
            ProteinHit pHit = new ProteinHit();
            pHit.convertProteinHit(protein);
            System.out.println(""+protein.getPeptideList());

//            System.out.println("=======" + protein.getPeptideList().get(0).getSequenceIndex());

        }
        */


        //ProteinGroup pGroup = new ProteinGroup();
        List l = new ArrayList();
        Set s = new HashSet();
        for (Iterator<Protein> itr = pitr; itr.hasNext(); ) {
            Protein protein = itr.next();

//	    if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
//		continue




            //ArrayList<Protein> aList = new ArrayList<Protein>();
            HashSet set = new HashSet();

            s.add(protein.getLocus());
//	    s.setConfidence(Float.parseFloat(protein.getConfidence()));

            if(protein.getPeptideSize()<=0)
            {


                if(protein.getProteinLine().startsWith("sp")) {
                    isUnitProt = true;
                }
                //	pGroup.add(pHit);
            }
            else
            {


                if(protein.getProteinLine().startsWith("sp")) {
                    isUnitProt = true;
                }

                if(isUnitProt) proteinGroup++;

                l.add(s);
                s = new HashSet();

                isUnitProt=false;
//		aList.clear();
            }
        }



        System.out.println(proteinGroup);

        reader.close();
//	DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

//	System.out.println(reader.getRedundantProteinNum());
//	System.out.println(reader.getNonRedundantProteinNum());
//	System.out.println(reader.getProteinFP());

    }

    public DTASelectFilterReader(String fileName) throws IOException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);


        init();
    }

    /*
    public Iterator <Protein> getProteins(final InputStream is) throws IOException {
        return getProteins(new InputStreamReader(is));
    }
    */
    private void readSummary() throws IOException
    {
        RandomAccessFile file = null;

        try {
            file = new RandomAccessFile(fileName, "r");

            file.seek(file.length()-500);

            String eachLine;
            eachLine=file.readLine();

            while( (eachLine=file.readLine()) != null && !eachLine.startsWith("Unfiltered") );

            String[] arr = eachLine.split("\t");
            this.unfilteredProteinNum = Integer.parseInt(arr[1]);
            this.unfilteredPeptideNum = Integer.parseInt(arr[2]);

            eachLine=file.readLine(); //Redundant line of DTASelect-filter.txt
            arr = eachLine.split("\t");
            this.redundantProteinNum = Integer.parseInt(arr[1]);
            this.redundantPeptideNum = Integer.parseInt(arr[2]);

            eachLine=file.readLine(); //NonRedundant line of DTASelect-filter.txt
            if(null != eachLine) {
                arr = eachLine.split("\t");
                if(arr.length>3) {
                    this.nonRedundantProteinNum = Integer.parseInt(arr[1]);
                    this.nonRedundantPeptideNum = Integer.parseInt(arr[2]);
                }
            }

            while( (eachLine=file.readLine()) != null && !eachLine.startsWith("Forward FP") && !eachLine.startsWith("Forward FDR"));

            if(null != eachLine) {
                arr = eachLine.split("\t");
                if(arr.length>3) {
                    this.proteinFP = Double.parseDouble(arr[1]);

                    if(arr.length>2) {
                        this.peptideFP = Double.parseDouble(arr[2]);
                        this.spectrumFP = Double.parseDouble(arr[3]);
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("Error:" + e);
        } finally {

            if(null != file) file.close();
        }
    }

    private void init() throws IOException
    {
        readSummary();

        try {
            br = new BufferedReader(new FileReader(fileName));

            readHeader();

            //Move line to parameters.  We assume parameters start after carriage return
            while ((lastLine = br.readLine()) != null) {
                if(lastLine.startsWith("ProLuCID") || lastLine.startsWith("SEQUEST") || lastLine.startsWith("BlindPTM") || lastLine.startsWith("?") || lastLine.startsWith("MASCOT") || lastLine.startsWith("Blazmass") || lastLine.startsWith("Comet")) {
                    break;
                }
            }
            if(lastLine.equals("MASCOT"))
                isMascot = true;
            StringBuffer sb = new StringBuffer();
            sb.append(lastLine=br.readLine());  //Read this line, which can be either empty or not

            while (!((lastLine = br.readLine())).equals(""))
            {
                sb.append(lastLine);
                sb.append("\n");
            }

            criteria = sb.toString();

            //Move line to parse protein

            while (!(lastLine = br.readLine()).startsWith("Locus"));
            Protein.setFeatureIndices(lastLine);

            while (!(lastLine = br.readLine()).startsWith("Unique"));
            //add setFeatureIndex here

            peptide.setFeatureIndex(lastLine);

            lastLine = br.readLine();
        } catch (Exception e) {
            e.printStackTrace();
            try {   if(null != br) br.close();
                System.out.println("error: "  + e);
            }
            catch(IOException ie) { }
        }
    }

    private void readHeader() throws IOException
    {
        lastLine = br.readLine();

        if(lastLine.startsWith("DTASelect v2.0"))
            version2 = true;
        else if(lastLine.startsWith("DTASelect v1")) {
            version2 = false;
            version1 = true;

        }

        //whie (!(lastLine = br.readLine()).endsWith("fasta")); .startsWith("Unique"));
        for(int i=0; i<2; i++)
            lastLine = br.readLine();

        //Remove directory name
        this.dbFilePathAndName = lastLine;
        this.dbFileName = lastLine.substring(lastLine.lastIndexOf("/")+1);

        String dbFileUpperCase = this.dbFileName.toUpperCase();

        //if(this.dbFileName.startsWith("UniProt") || this.dbFileName.startsWith("uniprot"))
        if(dbFileUpperCase.startsWith("UNIPROT"))
            this.fileType = "UniProt";
        else if(this.dbFileName.startsWith("EBI"))
            this.fileType = "IPI";

    }
    public String getFileType() {
        return fileType;
    }

    public void setFileType(String fileType) {
        this.fileType = fileType;
    }
    public Iterator <Protein> getProteins() throws IOException {

        return new Iterator<Protein>() {
            private Protein protein;
            private Peptide ipeptide;

            public boolean hasNext() {
                return lastLine != null && !lastLine.startsWith("\tProteins\t");
            }

            public Protein next() {

                try {

                    String tempLastLine = lastLine;

                    protein = getProtein(lastLine.split("\t"));
                    protein.setProteinLine(tempLastLine);

                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private Protein getProtein(String[] strArr) throws IOException {

//                String[] strArr = lastLine.split("\t");
//		String[] peptideLine;
                //protein = new Protein(strArr);

                protein = new Protein();
                try {
                    protein.setElement(strArr);

                } catch (ArrayIndexOutOfBoundsException e ) {
                    System.out.println(e.getMessage());
                }
                /*
                The format of the DTASelect-filter.txt file ipeptide line
                 is the following:

                 - a star if the ipeptide is unique to that protein (optional)
                 - then an integer greater than 1 (optional)
                 - then one of the characters 'M', 'Y', or 'N' (optional)
                 - then tab (mandatory)
                 - then the rest of the fields...

                 Anything that comes until the first tab is optional. You can
                 have a line that starts with "\t", or "*\t", or "2\t", or
                 "*2\t", or "*2M\t", or "Y\t", or "2N\t", or "*Y\t", etc...

		Peptide line starts like (*)(int)(M||Y||N)\t
                 */

                /*
                 *  Some proteins does not have ipeptide lines, because those proteins
                 *  are assumed to have identical peptides as following protein has.
                 *
                 **/

           /*     if(protein.getLocus().equals("tr|E5RIR4|E5RIR4_HUMAN")){
                    System.out.println("");
                }*/
                while ( ((lastLine = br.readLine()) != null && !lastLine.equals(""))
                        && !lastLine.startsWith("\tProteins\t"))
                {
                    strArr = lastLine.split("\t");

                    // If Spectrum Count position does not have a decimal point,
                    // it is not a ipeptide line
                    if(strArr[2].indexOf(".")<=0)
                        break;

                    ipeptide = new Peptide(lastLine, version2, version1,peptide);
                    //ipeptide = new Peptide(strArr, version2, version1);
                    protein.addPeptide(ipeptide);
                }

                return protein;
            }
        };
    }



    public static String getFilterParams(String fname) {

        BufferedReader br = null;
        try {
            java.io.File f = new java.io.File(fname);
            if(!f.exists())
                return "NA";

            br = new BufferedReader(new FileReader(f));

            String lastLine = null;

            while ( null != (lastLine = br.readLine())) {
                if(lastLine.contains("SQT format"))
                    break;
            }

            lastLine  = br.readLine();

            br.close();

            return lastLine;


        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error: " + e);
            try {   if(null != br) br.close(); }
            catch(IOException ie) { }
            return "NA";

        }
    }


    public void close() throws IOException
    {
        br.close();
    }

    public void setDbFileName(String dbFileName)
    {
        this.dbFileName = dbFileName;
    }

    public void setUnfilteredProteinNum(int unfilteredProteinNum)
    {
        this.unfilteredProteinNum = unfilteredProteinNum;
    }

    public void setRedundantProteinNum(int redundantProteinNum)
    {
        this.redundantProteinNum = redundantProteinNum;
    }

    public void setNonRedundantProteinNum(int nonRedundantProteinNum)
    {
        this.nonRedundantProteinNum = nonRedundantProteinNum;
    }

    public void setUnfilteredPeptideNum(int unfilteredPeptideNum)
    {
        this.unfilteredPeptideNum = unfilteredPeptideNum;
    }

    public void setRedundantPeptideNum(int redundantPeptideNum)
    {
        this.redundantPeptideNum = redundantPeptideNum;
    }

    public void setNonRedundantPeptideNum(int nonRedundantPeptideNum)
    {
        this.nonRedundantPeptideNum = nonRedundantPeptideNum;
    }

    public String getDbFileName()
    {
        return dbFileName;
    }
    public String getDbFilePathAndName()
    {
        return dbFilePathAndName;
    }

    public String getCriteria()
    {
        return criteria;
    }

    public int getUnfilteredProteinNum()
    {
        return unfilteredProteinNum;
    }

    public int getRedundantProteinNum()
    {
        return redundantProteinNum;
    }

    public int getNonRedundantProteinNum()
    {
        return nonRedundantProteinNum;
    }

    public int getUnfilteredPeptideNum()
    {
        return unfilteredPeptideNum;
    }

    public int getRedundantPeptideNum()
    {
        return redundantPeptideNum;
    }

    public int getNonRedundantPeptideNum()
    {
        return nonRedundantPeptideNum;
    }

    public boolean isVersion2() {
        return version2;
    }

    public void setProteinFP(double proteinFP)
    {
        this.proteinFP = proteinFP;
    }

    public double getProteinFP() {
        return proteinFP;
    }

    public void setPeptideFP(double peptideFP)
    {
        this.peptideFP = peptideFP;
    }

    public double getPeptideFP() {
        return peptideFP;
    }
    public void setSpectrumFP(double spectrumFP)
    {
        this.spectrumFP = spectrumFP;
    }

    public double getSpectrumFP() {
        return spectrumFP;
    }
    public boolean getIsMascot(){
        return isMascot;
    }
    public void setIsMascot(boolean im){
        isMascot = im;
    }
}

