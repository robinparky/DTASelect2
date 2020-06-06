import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class DTAToSQT {

    public static class PSM
    {
        public final String psmStr;
        public final String sequence;
        public final String chargeState;
        public final String xcorr;
        public final String deltaCN;
        public final String peaksMatched;
        public final String numPeaksTotal;
        public final String scanNumber;
        public final String expMass;
        public final String theorMass;
        public final String totalIntensity;

        private Map<String, String> locusMap = new HashMap<>();

        public PSM(String[] lineArr,  String [] infoArr, Set<String> locusSet) {
            this.psmStr = lineArr[1];

            xcorr = lineArr[2];
            deltaCN = lineArr[3];
            expMass = lineArr[5];
            theorMass = lineArr[6];
            totalIntensity = lineArr[8];
            this.sequence = lineArr[lineArr.length-1];
            int hackPeaks = Integer.parseInt(lineArr[lineArr.length-2]);
            peaksMatched = Integer.toString(hackPeaks);
            numPeaksTotal = "100";
            scanNumber = infoArr[1];
            chargeState = infoArr[infoArr.length-1];
            for(String locus : locusSet)
            {
                this.locusMap.put(locus, sequence);
            }
        }

        public String ToSQTString()
        {
            StringBuilder sb = new StringBuilder();
            sb.append("S\t").append(scanNumber).append("\t").append(scanNumber).append("\t").append(chargeState);
            sb.append("\t").append(1000).append("\t").append("NOT_A_REAL_THREAD_1").append("\t").append(expMass)
                    .append("\t").append(totalIntensity);
            sb.append("\t").append(1000).append("\t").append(100).append("\t").append("100").append("\n");

            sb.append("M\t1\t1\t").append(theorMass).append("\t0.000\t").append(xcorr).append("\t1.000\t").append(peaksMatched);
            sb.append("\t").append(numPeaksTotal).append("\t").append(sequence).append("\tU\n");
            for(Map.Entry<String,String> entry: locusMap.entrySet())
            {
                sb.append("L\t").append(entry.getKey()).append("\t10\t").append(entry.getValue()).append("\n");
            }

            sb.append("M\t2\t2\t").append(theorMass).append("\t").append(deltaCN).append("\t").append(xcorr).append("\t1.000\t")
                    .append(numPeaksTotal);
            sb.append("\t").append(numPeaksTotal).append("\t").append(sequence).append("\tU\n");

            for(Map.Entry<String,String> entry: locusMap.entrySet())
            {
                sb.append("L\t").append(entry.getKey()).append("\t10\t").append(entry.getValue()).append("\n");
                break;
            }
            return sb.toString();
        }


        public Map<String, String> getLocusMap() {
            return locusMap;
        }

    }

    public static Map<String, List<PSM>> ReadDTASelectFilter(String path) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line;
        boolean readMode = false;
        Map<String, PSM> psmMap = new HashMap<>();
        Map<String, List<PSM>> filePSMMap  = new HashMap<>();
        Set<String> proteinLocusSet = new HashSet<>();
        boolean peptideMode = false;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("\tProtein")) {
                break;
            }
            if (readMode) {
                if(line.startsWith("\t") || line.startsWith("*"))
                {
                    peptideMode = true;
                    String [] arr = line.split("\t");
                    String infoStr= arr[1];
                    String[] infoArr = infoStr.split("\\.");
                    String seq = arr[arr.length-1];
                    String fileName = infoArr[0];
                    PSM psm = psmMap.get(infoStr);
                    if(psm !=null)
                    {
                        for(String locus: proteinLocusSet)
                        {
                            psm.getLocusMap().put(locus,seq);
                        }

                    }
                    else
                    {
                        psm =  new PSM(arr, infoArr, proteinLocusSet);
                        psmMap.put(infoStr,psm);
                        List<PSM> psmList = filePSMMap.getOrDefault(fileName, new ArrayList<>());
                        psmList.add(psm);
                        filePSMMap.put(fileName, psmList);
                    }

                }
                else
                {
                    if(peptideMode)
                    {
                        peptideMode = false;
                        proteinLocusSet = new HashSet<>();
                    }
                    String proteinLocus = line.substring(0, line.indexOf("\t"));
                    proteinLocusSet.add(proteinLocus);
                }
            } else if (line.startsWith("Unique\t")) {
                readMode = true;
            }
        }
        br.close();
        return filePSMMap;
    }

    public static void PrintPSM(String directory, String filename, List<PSM> psmList ) throws IOException {
        Path path = Paths.get(directory,filename+".sqt");
        BufferedWriter bw = new BufferedWriter(new FileWriter(path.toFile()));
        for(PSM psm :psmList)
        {
            bw.append(psm.ToSQTString());
            bw.newLine();
        }
        bw.close();
    }


    public static void main(String [] args) throws IOException {
        if(args.length <2 || args[0].equals("--help") || args[0].equals("-h"))
        {
            if(args.length <2 )
            {
                System.out.println("improper input");
            }
            System.out.println("java -cp /path/to/DTASelect2/ DTAToSqt /path/to/DTASelect-filter.txt /path/to/output/directory/");
            return;
        }
       String inputPath = args[0];
       String outputDirectory = args[1];
        Map<String, List<PSM>> fileNamePSMMap = ReadDTASelectFilter(inputPath);
       for(Map.Entry<String, List<PSM>> entry: fileNamePSMMap.entrySet())
       {
           PrintPSM(outputDirectory, entry.getKey(),entry.getValue());
       }
    }



}
