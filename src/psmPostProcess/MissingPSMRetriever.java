package psmPostProcess;

import proteinGrouper.DTASelectFilterReader;
import proteinGrouper.Peptide;
import proteinGrouper.Protein;
import sun.misc.FloatingDecimal;

import java.io.*;
import java.util.*;



/**
 * Created by Titus Jung titusj@scripps.edu on 8/22/19.
 */
public class MissingPSMRetriever {

    private boolean DisplayPI = false;
    private boolean DisplayKD = false;
    private boolean DisplayBB = false;
    private boolean DisplayHPLC = false;
    private boolean DisplayDeltaMass = false;
    private boolean extraColumns = false;
    private Map<String, List<PeptideContainer>> seqPSMMap = new HashMap<>();
    private Set<String> sequencesToLooksFor = new HashSet<>();
    private Set<String> existingPsm = new HashSet<>( );
    private Map<String, Integer> peptideRedundacyMap = new HashMap<>();
    private int psmCount =0;
    private boolean isotopes =true;


    public static void main(String [] args) throws IOException {
        String dtaselect = args[0];
        String searchDir = args[1];
        String output = args[2];
        boolean printPeptides = args[3].equals("1");

        MissingPSMRetriever missingPSMRetriever = new MissingPSMRetriever();
        missingPSMRetriever.collectMissedPSM(dtaselect, searchDir, output, printPeptides);
    }


    public Set<String> foundPSM(String dtaselect) throws IOException {
        DTASelectFilterReader filterReader = new DTASelectFilterReader(dtaselect);
        Set<String> foundPSM = new HashSet<>();
        for (Iterator<Protein> proteinItr = filterReader.getProteins(); proteinItr.hasNext(); )
        {
            Protein protein = proteinItr.next();
            for(Iterator<Peptide> peptideIterator = protein.getPeptides(); peptideIterator.hasNext();)
            {
                Peptide peptide = peptideIterator.next();
                String midSeq = getMiddleSequence(peptide.getSequence());
                String key = midSeq + peptide.getChargeState();
                peptide.getRedundancy();
                foundPSM.add(peptide.getFullFileName());
                sequencesToLooksFor.add(midSeq);
            }
        }
        return foundPSM;
    }

    private void generateNewPSM(String sqtFile, boolean printMode) throws IOException {
        SQTParser parser = new SQTParser(sqtFile);
        String fileName = sqtFile;
        int index =-1;
        if((index = sqtFile.lastIndexOf(File.separator))>=0)
        {
            fileName = sqtFile.substring(index+1).replace(".sqt","");
        }
        String cleanSeq = null;
        String xcorr= null;
        String mass= null;
        String calcMass= null;
        String totalIntensity = null;
        double ionProportion = 0;
        String spr = "1";
        String probScore = null;
        String oldSeq = null;
        for(Iterator<SQTPeptide> peptideIterator= parser.getSQTPeptide(); peptideIterator.hasNext();)
        {
            SQTPeptide peptide = peptideIterator.next();
            String key = fileName + "." + peptide.getHiScan() + "."+ peptide.getHiScan() + "." + peptide.getChargeState();
            boolean findMode = false;
            int mi =0;
            for(Iterator<MLine> mLineIterator = peptide.getMLine(); mLineIterator.hasNext(); )
            {

                MLine mLine = mLineIterator.next();
                if(mi == 0)
                {
                    String sequence  = getMiddleSequence(mLine.getSequence());
                    if(sequencesToLooksFor.contains(sequence))
                    {
                        psmCount++;
                        if(printMode)
                        {
                            findMode = true;
                            cleanSeq = cleanSequence(sequence);
                            oldSeq = sequence;
                            xcorr = mLine.getXcorr();
                            mass = peptide.getCalMZ();
                            calcMass = mLine.getCalMZ();
                            float tempCalcMass = Float.parseFloat(calcMass);
                            calcMass = Float.toString(tempCalcMass);
                            totalIntensity = peptide.getTotalIntensity();
                            double matched = Double.parseDouble(mLine.getMatchedIons());
                            double predictedIons = Double.parseDouble(mLine.getPredictedIons());
                            ionProportion = matched/predictedIons*100.0;

                            float tempProbScore = Float.parseFloat(mLine.getSp());
                            probScore = Float.toString(tempProbScore);
                        //    System.out.println(mLine.getMLine()+"\t"+key);
                        }
                    }
                }
                else if(findMode)
                {
                    String sequence = cleanSequence(getMiddleSequence(mLine.getSequence()));
                    if(!sequence.equals(cleanSeq) && ! existingPsm.contains(key))
                    {
                        String deltaCn = mLine.getDeltCN();
                        double prcMass = Double.parseDouble(mass);
                        double dCalcMass = Double.parseDouble(calcMass);
                        int cs = peptide.getChargeStateInt();
                        findMode = false;
                        double pprm = GetMassOffsets(prcMass,dCalcMass,cs);
                        PeptideContainer container = new PeptideContainer(oldSeq,
                                deltaCn, xcorr,mass,calcMass,totalIntensity, Float.toString(RoundTo((float)ionProportion,1))
                                , probScore, Float.toString(RoundTo((float)pprm,1)),key);

                        List<PeptideContainer> test =  seqPSMMap.getOrDefault(oldSeq, new ArrayList<>());
                        test.add(container);
                        seqPSMMap.put(oldSeq, test);
                    }
                }
                mi++;
            }

        }
    }

    public static String getMiddleSequence(String seq)
    {
        int indexA = seq.indexOf('.');
        int indexb = seq.lastIndexOf('.');
        return seq.substring(indexA+1,indexb);
    }

    public static String removeNumerics(String str) {
        if (str == null) {
            return null;
        }else {
            return str.replaceAll("[0-9|.]", "");
        }
    }

    public static String cleanSequence(String seq) {
        seq = removeNumerics(seq);
        seq = seq.replaceAll("\\(\\)", "");
        seq = seq.replaceAll("\\(-\\)", "");
        //seq = seq.replaceAll("\\(\\)", "\\*");
        seq = seq.replaceAll("\\*", "");
        seq = seq.replaceAll("#", "");
        seq = seq.replaceAll("@", "");
        Protein.getNsafIndex();

        return seq;
    }


    public void collectMissedPSM(String dtaselect, String searchDir, String output, boolean printPeptides)
            throws IOException {
        existingPsm =    foundPSM(dtaselect);
        File dir = new File(searchDir);
       String[] list =  dir.list(new FilenameFilter() {
            @Override
            public boolean accept(File file, String s) {
                return s.endsWith(".sqt");
            }
        });
       for(String s: list)
       {
           generateNewPSM(searchDir + File.separator+ s,printPeptides);
       }
       if(printPeptides)
       {
            updatePeptides(dtaselect,output);
       }
       else{
            updateFDROnly(dtaselect, output);
       }
    }


    public void updatePeptides(String dtaSelectInput, String output) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectInput));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        String line;
        String [] arr2;
        int peptideColNum = -1;
        int confIndex = -1;
        int ppmIndex = -1;
        int piIndex = -1;
        int kdIndex = -1;
        int bbIndex = -1;
        int hplcIndex = -1;
        int seqPosIndex = -1;
        int sequenceIndex = -1;
        int fileNameIndex = -1;
        String kd = null;
        String pi = null;
        String conf = null;
        String bb = null;
        String hplc = null;
        String seqPos = null;
        boolean findMode = false;
        String sequenceToPrint = null;
        double confDouble = Double.MAX_VALUE;
        String uniqueMarker = "";
        String oldSequence = null;
        int proteinLineLength = -1;

        while((line = br.readLine())!=null)
        {
            String [] arr = line.split("\t");

            if(line.startsWith("Forward matches\t"))
            {
                arr[3] = Integer.toString(psmCount);
                bw.append(arr[0]).append("\t").append(arr[1]).append("\t").append(arr[2]).append("\t").append(arr[3]);
                bw.newLine();
            }
            else if(line.startsWith("\tProteins\t") || (arr.length == proteinLineLength && findMode))
            {
                List<PeptideContainer> containerList = seqPSMMap.get(sequenceToPrint);
                for(PeptideContainer c: containerList)
                {
                    String row = c.generateRowString(uniqueMarker,kd,pi,bb,hplc,seqPos,conf,ppmIndex>0,oldSequence);
                    bw.write(row);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                findMode = false;

            }
            else if((line.startsWith("\t") || line.startsWith("*\t")) && (arr.length == peptideColNum ))
            {
                String filename = arr[fileNameIndex];
                String sequence = arr[sequenceIndex];
                String midSequence = getMiddleSequence(sequence);

                if(findMode && !midSequence.equals(sequenceToPrint))
                {
                    findMode = false;
                    List<PeptideContainer> containerList = seqPSMMap.get(sequenceToPrint);
                    for(PeptideContainer c: containerList)
                    {
                        String row = c.generateRowString(uniqueMarker,kd,pi,bb,hplc,seqPos,conf,ppmIndex>0,oldSequence);
                        bw.write(row);
                        bw.newLine();
                    }
                    confDouble = Double.MAX_VALUE;
                }

                if(seqPSMMap.containsKey(midSequence))
                {
                    findMode = true;
                    uniqueMarker = arr[0];
                    sequenceToPrint = midSequence;
                    kd = kdIndex != -1 ? arr[kdIndex] : "";
                    pi = piIndex != -1 ? arr[piIndex] : "";
                    bb = bbIndex != -1 ? arr[bbIndex] : "";
                    hplc = hplcIndex != -1 ? arr[hplcIndex] : "";
                    seqPos = seqPosIndex != -1 ? arr[seqPosIndex] : "";
                    String conftemp =  arr[confIndex];
                    double temp = Double.parseDouble(conftemp);
                    if(temp < confDouble)
                    {
                        confDouble = temp;
                        conf = conftemp;
                    }
                    oldSequence = sequence;
                }

                bw.write(line);
                bw.newLine();

            }
            else
            {
                if(line.startsWith("Unique\t"))
                {
                    peptideColNum = arr.length;
                    for(int i =0 ; i<arr.length; i++)
                    {
                        String colName = arr[i];
                        switch (colName)
                        {
                            case "Conf%":
                                confIndex = i;
                                break;
                            case "PPM":
                                ppmIndex =i;
                                break;
                            case "pI":
                                piIndex =i;
                                break;
                            case "KD":
                                kdIndex =i;
                                break;
                            case "BB":
                                bbIndex = i;
                                break;
                            case "HPLC":
                                hplcIndex = i;
                                break;
                            case "SeqPosition":
                                seqPosIndex = i;
                                break;
                            case "Sequence":
                                sequenceIndex =i;
                                break;
                            case "FileName":
                                fileNameIndex = i;
                                break;
                        }
                    }
                }
                if(line.startsWith("Locus"))
                {
                    proteinLineLength = arr.length;
                }
                bw.write(line);
                bw.newLine();

            }
        }
        br.close();
        bw.close();

    }


    public void updateFDROnly(String dtaSelectInput, String output) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectInput));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        String line;
        while((line = br.readLine())!=null)
        {
            if(line.startsWith("Forward matches\t"))
            {
                String [] arr = line.split("\t");
                arr[3] = Integer.toString(psmCount);
                bw.append(arr[0]).append("\t").append(arr[1]).append("\t").append(arr[2]).append("\t").append(arr[3]);
                bw.newLine();
            }
            else
            {
                bw.write(line);
                bw.newLine();
            }
        }
        br.close();
        bw.close();
    }

    public void setOptionDisplay(boolean displayBB, boolean displayHPLC, boolean displayKD, boolean displayPI)
    {
        this.DisplayBB = displayBB;
        this.DisplayHPLC = displayHPLC;
        this.DisplayKD = displayKD;
        this.DisplayPI = displayPI;
    }

    public float CalculatepI(String sequence) {
        String         Sequence = sequence;
        int            Length = Sequence.length();
        int            Looper;
        int            CountLys = 0;
        int            CountArg = 0;
        int            CountHis = 0;
        int            CountAsp = 0;
        int            CountGlu = 0;
        int            CountCys = 0;
        int            CountTyr = 0;
        float          CurrentPH = 7.0f;
        float          CurrentJump = 3.5f;
        float          CurrentCharge;
        float          LastCharge = 0;
        float          MWAccum = 0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                    case 'A':
                        MWAccum += 71.0;
                        break;
                    case 'C':
                        MWAccum += 103.0;
                        CountCys ++;
                        break;
                    case 'D':
                        MWAccum += 115.0;
                        CountAsp ++;
                        break;
                    case 'E':
                        MWAccum += 129.0;
                        CountGlu ++;
                        break;
                    case 'F':
                        MWAccum += 147.0;
                        break;
                    case 'G':
                        MWAccum += 57.0;
                        break;
                    case 'H':
                        MWAccum += 137.0;
                        CountHis ++;
                        break;
                    case 'I':
                        MWAccum += 113.0;
                        break;
                    case 'K':
                        MWAccum += 128.0;
                        CountLys ++;
                        break;
                    case 'L':
                        MWAccum += 113.0;
                        break;
                    case 'M':
                        MWAccum += 131.0;
                        break;
                    case 'N':
                        MWAccum += 114.0;
                        break;
                    case 'P':
                        MWAccum += 97.0;
                        break;
                    case 'Q':
                        MWAccum += 128.0;
                        break;
                    case 'R':
                        MWAccum += 156.0;
                        CountArg ++;
                        break;
                    case 'S':
                        MWAccum += 87.0;
                        break;
                    case 'T':
                        MWAccum += 101.0;
                        break;
                    case 'V':
                        MWAccum += 99.0;
                        break;
                    case 'W':
                        MWAccum += 186.0;
                        break;
                    case 'Y':
                        MWAccum += 176.0;
                        CountTyr ++;
                        break;
                }
            }
            /* Use a bracketing strategy to identify the isoelectric
             * point.  Calculate charge at pH of 7, and then move up
             * 3.5 or down 3.5 as necessary.  Make each successive
             * move up or down only half as large.  Keep going until
             * two successive charges reported match to one place past
             * the decimal point.
             */
            CurrentCharge = ChargeAtPH(CurrentPH, CountLys,
                    CountArg, CountHis,
                    CountAsp, CountGlu,
                    CountCys, CountTyr);
            while (RoundTo(CurrentCharge,1) != RoundTo(LastCharge,1)) {
                //		System.out.println("pH:\t" + new Float(CurrentPH).toString()
                //              + "\tCharge\t" + new Float(CurrentCharge).toString());
                if (CurrentCharge > 0)
                    CurrentPH += CurrentJump;
                else
                    CurrentPH -= CurrentJump;
                CurrentJump /= 2;
                LastCharge = CurrentCharge;
                CurrentCharge = ChargeAtPH(CurrentPH,
                        CountLys, CountArg,
                        CountHis, CountAsp,
                        CountGlu, CountCys,
                        CountTyr);
                if ( (CurrentPH > 14) || (CurrentPH < 0) ) {
                    System.out.println("pI can't be figured for ");
                    System.exit(0);
                }
            }
        }
        return CurrentPH;
    }

    public float CalculateKyteDoolittle(String sequence) {
        String         Sequence = sequence;
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                    case 'A':
                        MWAccum += 18.0;
                        break;
                    case 'C':
                        MWAccum += 25.0;
                        break;
                    case 'D':
                        MWAccum -= 35.0;
                        break;
                    case 'E':
                        MWAccum -= 35.0;
                        break;
                    case 'F':
                        MWAccum += 28.0;
                        break;
                    case 'G':
                        MWAccum -= 4.0;
                        break;
                    case 'H':
                        MWAccum -= 32.0;
                        break;
                    case 'I':
                        MWAccum += 45.0;
                        break;
                    case 'K':
                        MWAccum -= 39.0;
                        break;
                    case 'L':
                        MWAccum += 38.0;
                        break;
                    case 'M':
                        MWAccum += 19.0;
                        break;
                    case 'N':
                        MWAccum -= 35.0;
                        break;
                    case 'P':
                        MWAccum -= 16.0;
                        break;
                    case 'Q':
                        MWAccum -= 35.0;
                        break;
                    case 'R':
                        MWAccum -= 45.0;
                        break;
                    case 'S':
                        MWAccum -= 8.0;
                        break;
                    case 'T':
                        MWAccum -= 7.0;
                        break;
                    case 'V':
                        MWAccum += 42.0;
                        break;
                    case 'W':
                        MWAccum -= 9.0;
                        break;
                    case 'Y':
                        MWAccum -= 13.0;
                        break;
                }
            }
        }
        return MWAccum/10.0f;
    }


    public float CalculateBullBreese(String sequence) {
        String         Sequence = sequence;
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0.0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                    case 'A':
                        MWAccum += 610;
                        break;
                    case 'C':
                        MWAccum += 360;
                        break;
                    case 'D':
                        MWAccum += 610;
                        break;
                    case 'E':
                        MWAccum += 510;
                        break;
                    case 'F':
                        MWAccum += -1520;
                        break;
                    case 'G':
                        MWAccum += 810;
                        break;
                    case 'H':
                        MWAccum +=690;
                        break;
                    case 'I':
                        MWAccum += -1450;
                        break;
                    case 'K':
                        MWAccum += 460;
                        break;
                    case 'L':
                        MWAccum += -1650;
                        break;
                    case 'M':
                        MWAccum += -660;
                        break;
                    case 'N':
                        MWAccum += 890;
                        break;
                    case 'P':
                        MWAccum += -170;
                        break;
                    case 'Q':
                        MWAccum += 970;
                        break;
                    case 'R':
                        MWAccum += 690;
                        break;
                    case 'S':
                        MWAccum += 420;
                        break;
                    case 'T':
                        MWAccum += 290;
                        break;
                    case 'V':
                        MWAccum += -750;
                        break;
                    case 'W':
                        MWAccum += -1200;
                        break;
                    case 'Y':
                        MWAccum += -1430;
                        break;
                }
            }
        }
        return MWAccum;
    }


    public float CalculateHPLCpH34(String sequence) {
        String         Sequence = sequence;
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0.0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                    case 'A':
                        MWAccum += 42;
                        break;
                    case 'C':
                        MWAccum += 84;
                        break;
                    case 'D':
                        MWAccum += -51;
                        break;
                    case 'E':
                        MWAccum += -37;
                        break;
                    case 'F':
                        MWAccum += 174;
                        break;
                    case 'G':
                        MWAccum += 0;
                        break;
                    case 'H':
                        MWAccum += -228;
                        break;
                    case 'I':
                        MWAccum += 181;
                        break;
                    case 'K':
                        MWAccum += -203;
                        break;
                    case 'L':
                        MWAccum += 180;
                        break;
                    case 'M':
                        MWAccum += 118;
                        break;
                    case 'N':
                        MWAccum += -103;
                        break;
                    case 'P':
                        MWAccum += 86;
                        break;
                    case 'Q':
                        MWAccum += -96;
                        break;
                    case 'R':
                        MWAccum += -156;
                        break;
                    case 'S':
                        MWAccum += -64;
                        break;
                    case 'T':
                        MWAccum += -26;
                        break;
                    case 'V':
                        MWAccum += 134;
                        break;
                    case 'W':
                        MWAccum += 146;
                        break;
                    case 'Y':
                        MWAccum += 51;
                        break;
                }
            }
        }
        return MWAccum/100.0f;
    }

    public static float ChargeAtPH(float pH, int CountLys, int CountArg,
                                   int CountHis, int CountAsp, int CountGlu, int CountCys, int CountTyr) {
        // Start out accumulator with charge of termini
        float Accum = PercentPositive(pH, 8.0f) - PercentNegative(pH, 3.1f);
        Accum += CountLys * PercentPositive(pH, 10.0f);
        Accum += CountArg * PercentPositive(pH, 12.0f);
        Accum += CountHis * PercentPositive(pH, 6.5f);
        Accum -= CountAsp * PercentNegative(pH, 4.4f);
        Accum -= CountGlu * PercentNegative(pH, 4.4f);
        Accum -= CountCys * PercentNegative(pH, 8.5f);
        Accum -= CountTyr * PercentNegative(pH, 10.0f);
        return Accum;
    }

    public static float PercentPositive(float pH, float pK) {
        double ConcentrationRatio = Math.pow(10f, pK - pH);
        return new Double(ConcentrationRatio / (ConcentrationRatio + 1))
                .floatValue();
    }

    public static float PercentNegative(float pH, float pK) {
        double ConcentrationRatio = Math.pow(10, pH - pK);
        return new Double(ConcentrationRatio / (ConcentrationRatio + 1))
                .floatValue();
    }

    public static float RoundTo(float Value, int Places) {
        // Converts a value to a rounded value
        if (Places == 0) {
            return new Integer(Math.round(Value)).intValue();
        } else {
            double Multiplier = Math.pow(10, Places);
            return new Double(Math.rint(Value * Multiplier) / Multiplier)
                    .floatValue();
        }
    }

    public double  GetMassOffsets(double PrecursorMass, double CalcPreMass, int ChargeState) {
        int counter;
        double diffC12C13 = 1.003354826;
        double AbsoluteOffset = PrecursorMass - CalcPreMass;

        double Raw_PPM_Offset = AbsoluteOffset;
        if (isotopes) {
            for (counter = 1; counter < 2 * ChargeState + 1; counter++) {
                if (Math.abs(Raw_PPM_Offset) > (Math.abs(AbsoluteOffset - counter * diffC12C13)))
                    Raw_PPM_Offset = AbsoluteOffset - counter * diffC12C13;
            }
        }
        Raw_PPM_Offset = Raw_PPM_Offset * 1E6 / PrecursorMass;
        double Adjusted_PPM_Offset = Raw_PPM_Offset;
        return Adjusted_PPM_Offset;
    }


    public class PeptideContainer
    {
        public final String middleSeq;
        public final String deltaCn;
        public final String xcorr;
        public final String mass;
        public final String calcMass;
        public final String totalIntensity;
        public final String ionProportion;
        public final String probScore;
        public final String pprm;
        public final String fileNameKey;

        public PeptideContainer(String middleSeq, String deltaCn, String xcorr, String mass, String calcMass, String totalIntensity, String ionProportion, String probScore, String pprm, String fileNameKey) {
            this.middleSeq = middleSeq;
            this.deltaCn = deltaCn;
            this.xcorr = xcorr;
            this.mass = mass;
            this.calcMass = calcMass;
            this.totalIntensity = totalIntensity;
            this.ionProportion = ionProportion;
            this.probScore = probScore;
            this.pprm = pprm;
            this.fileNameKey = fileNameKey;
        }

        public String generateRowString(String uniqueString, String kd, String pi, String bb, String hplc, String seqPos,
                                        String conf, boolean printPPM, String sequence)
        {
            String ppm = pprm;
            if(!printPPM)
                ppm   ="";
            StringBuilder rowSB  = new StringBuilder();
            rowSB.append(uniqueString).append("\t");
            rowSB.append(fileNameKey).append("\t");
            rowSB.append(xcorr).append("\t");
            rowSB.append(deltaCn).append("\t");

            rowSB.append(conf).append("\t");
            rowSB.append(mass).append("\t");
            rowSB.append(calcMass).append("\t");
            if(ppm.length()>0)
                rowSB.append(ppm).append("\t");
            rowSB.append(totalIntensity).append("\t");
            rowSB.append("1").append("\t");
            rowSB.append(probScore).append("\t");
            if(pi.length()>0)
                rowSB.append(pi).append("\t");
            if(kd.length()>0)
                rowSB.append(kd).append("\t");
            if(bb.length()>0)
                rowSB.append(bb).append("\t");
            if(hplc.length()>0)
                rowSB.append(hplc).append("\t");
            rowSB.append(ionProportion).append("\t");
            rowSB.append("1").append("\t");
            rowSB.append(sequence);
            if(seqPos.length()>0)
                rowSB.append("\t").append(seqPos);
            return rowSB.toString();


        }
    }

    public boolean isIsotopes() {
        return isotopes;
    }

    public void setIsotopes(boolean isotopes) {
        this.isotopes = isotopes;
    }
}
