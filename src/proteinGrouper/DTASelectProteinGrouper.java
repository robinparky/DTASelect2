package proteinGrouper;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/8/19.
 */
import java.io.*;
import java.util.*;

/**
 * Created by Titus Jung titusj@scripps.edu on 11/1/18.
 */
public class DTASelectProteinGrouper {

    public static class ProteinContainer implements Comparable<ProteinContainer>
    {
        public final Protein protein;
        public ProteinContainer(Protein protein)
        {
            this.protein = protein;
        }
        public int hash()
        {
            return protein.getLocus().hashCode();
        }

        @Override
        public int compareTo(ProteinContainer proteinContainer) {
            if(proteinContainer==null)
                return -1;
            return this.protein.getLocus().compareTo(proteinContainer.protein.getLocus());
        }

        public boolean equals(Object object)
        {
            if(object !=null)
            {
                if(object instanceof ProteinContainer)
                {
                    ProteinContainer container = (ProteinContainer) object;
                    return compareTo(container)==0;
                }
            }
            return false;
        }
    }

    public static class PeptideContainer
    {
        public final Peptide peptide;

        public PeptideContainer(Peptide peptide) {
            this.peptide = peptide;
        }

        public int hash()
        {
            return peptide.getPeptideLine().hashCode();
        }
    }

    public static class ProteinGroup implements Comparable<ProteinGroup>
    {
        public final String locusKey;
        public final Set<ProteinContainer> containerList;

        public ProteinGroup( Set<ProteinContainer> containerList) {
            StringBuilder proteinGroupKeySB = new StringBuilder();
            for(ProteinContainer p: containerList)
            {
                proteinGroupKeySB.append(p.protein.getLocus());
            }
            locusKey = proteinGroupKeySB.toString();
            this.containerList = containerList;
        }

        public int hash()
        {
            return locusKey.hashCode();
        }


        @Override
        public int compareTo(ProteinGroup proteinGroup) {
            if(proteinGroup==null)
                return -1;
            return locusKey.compareTo(proteinGroup.locusKey);
        }

        public boolean equals(Object obj)
        {
            if(obj!=null)
            {
                if(obj instanceof ProteinGroup)
                {
                    ProteinGroup group2 = (ProteinGroup) obj;
                    return compareTo(group2)==0;
                }
            }
            return false;
        }

    }

    public static String convertArrToLine(String [] arr)
    {
        StringBuilder sb = new StringBuilder();
        for(String s: arr)
        {
            sb.append(s).append("\t");
        }
        return sb.toString();
    }


    public static void process(String input, String output) throws IOException {
        String dtaSelectFilterFile = input;

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaSelectFilterFile);
        String criteria = reader.getCriteria();
        String [] arr = criteria.split(" ");
        int pfilter = -1;
        for(int i=0; i<arr.length; i++)
        {
            if(arr[i].equals("-p"))
            {
                pfilter = Integer.parseInt(arr[i+1]);
            }
        }
        Map<ProteinContainer,List<Peptide>> proteinUniquePeptideMap = new HashMap<>();
        Map<String,Set<ProteinContainer>> peptideToProteinMap = new HashMap<>();
        Set<ProteinGroup> oldProteinSet = new HashSet<>();
        Set<String> oldProteinLocusSet = new HashSet<>();
        Map<String,String> slineMap = new HashMap<>();
        List<Protein> proteinList = new ArrayList<>();
        for(Iterator<Protein> p = reader.getProteins(); p.hasNext(); )
        {
            Protein protein = p.next();
            List<Peptide> uniquePeptideList = new ArrayList<>();
            List<Peptide> nonUniquePeptideList = new ArrayList<>();
            int numPeptides = 0;
            for(Iterator<Peptide> peptideItr = protein.getPeptides(); peptideItr.hasNext();)
            {
                Peptide peptide = peptideItr.next();
                if(peptide.isUnique())
                {
                    uniquePeptideList.add(peptide);
                }
                else
                {
                    nonUniquePeptideList.add(peptide);
                }
                numPeptides++;
            }
            proteinList.add(protein);

            if(numPeptides>0)
            {
                if(numPeptides>=pfilter)
                {
                    //System.out.println(protein.getLocus());
                    List<ProteinContainer> containerList = new ArrayList<>();
                    Set<ProteinContainer> containerSet = new TreeSet<>();
                    for(Protein prot: proteinList)
                    {
                        ProteinContainer container = new ProteinContainer(prot);
                        containerList.add(container);
                        containerSet.add(container);
                    }
                    ProteinGroup group = new ProteinGroup(containerSet);
                    oldProteinSet.add(group);
                    oldProteinLocusSet.add(group.locusKey);
                    for(Peptide peptide: nonUniquePeptideList)
                    {
                        String sline = convertArrToLine(peptide.getPeptideLine());
                        String fileScanCsStr = peptide.getFullFileName();
                        slineMap.put(fileScanCsStr, sline);
                        Set<ProteinContainer> proteinSet = peptideToProteinMap.getOrDefault(fileScanCsStr, new TreeSet<>());
                        proteinSet.addAll(containerSet);
                        peptideToProteinMap.put(fileScanCsStr,proteinSet);
                    }
                    for(ProteinContainer container: containerList)
                    {
                        proteinUniquePeptideMap.put(container,uniquePeptideList);
                    }
                }
                proteinList = new ArrayList<>();
            }
        }
        reader.close();


        Map<String,List<String>> locusToPeptideMap = new HashMap<>();
        Map<String,ProteinGroup> groupMap = new HashMap<>();
        for(Map.Entry<String,Set<ProteinContainer>> entry: peptideToProteinMap.entrySet())
        {
            String peptideLine = entry.getKey();
            Set<ProteinContainer> proteinSet = entry.getValue();
            if(proteinSet.size()>1)
            {
                ProteinGroup group = new ProteinGroup(proteinSet);
                if(!oldProteinLocusSet.contains(group.locusKey))
                {
                    List<String> peptideList = locusToPeptideMap.getOrDefault(group.locusKey,new ArrayList<>());
                    peptideList.add(peptideLine);
                    locusToPeptideMap.put(group.locusKey,peptideList);
                    groupMap.put(group.locusKey,group);
                }
            }
        }
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectFilterFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));

        boolean copyMode = true;
        String line;
        String endHeader = null;
        List<String> peptideLines = new ArrayList<>();
        List<String> proteinLines = new ArrayList<>();
        boolean proteinMode = false;
        while((line=br.readLine())!=null)
        {
            if(line.startsWith("\tProteins"))
            {
                if(!proteinMode)
                {
                    if(peptideLines.size()>0)
                    {
                        for(String s: proteinLines)
                        {
                            bw.write(s);
                            bw.newLine();
                        }
                        for(String s: peptideLines)
                        {
                            bw.write(s);
                            bw.newLine();
                        }
                    }
                    peptideLines = new ArrayList<>();
                    proteinLines = new ArrayList<>();
                }

                endHeader = line;
                break;
            }

            if(copyMode)
            {
                bw.write(line);
                bw.newLine();
            }
            else
            {

                if (line.startsWith("*")) {
                    peptideLines.add(line);
                    proteinMode = false;
                }
                else if(line.charAt(0)=='\t')
                {
                    proteinMode = false;
                }
                else
                {
                    if(!proteinMode)
                    {
                        if(peptideLines.size()>0)
                        {
                            for(String s: proteinLines)
                            {
                                bw.write(s);
                                bw.newLine();
                            }
                            for(String s: peptideLines)
                            {
                                bw.write(s);
                                bw.newLine();
                            }
                        }
                        peptideLines = new ArrayList<>();
                        proteinLines = new ArrayList<>();
                    }
                    proteinMode = true;
                    proteinLines.add(line);
                }
                //   bw.write(line);
                //bw.newLine();

            }
            if(line.startsWith("Unique\t"))
            {
                copyMode = false;
            }
        }
        for(Map.Entry<String,List<String>> entry: locusToPeptideMap.entrySet())
        {
            ProteinGroup group = groupMap.get(entry.getKey());
            List<String >  peptideLine = entry.getValue();
            if(peptideLine.size()>=pfilter)
            {
                for(ProteinContainer p : group.containerList) {
                    bw.write(p.protein.getProteinLine());
                    bw.newLine();
                }
                for(String s: peptideLine)
                {
                    String sline = slineMap.get(s);
                    bw.write("*");
                    bw.write(sline);
                    bw.newLine();

                }
            }

        }

        bw.write(endHeader);
        bw.newLine();
        while((line=br.readLine())!=null)
        {
            bw.write(line);
            bw.newLine();

        }

        bw.close();
        br.close();
    }




    public static void main(String [] args) throws IOException {
        process(args[0],args[1]);

    }
}

