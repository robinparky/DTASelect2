package psmPostProcess;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
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
 * @version $Id: SQTParser.java,v 1.1 2014/09/09 19:29:52 rpark Exp $
 */
  public class SQTParser
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;

    public static void main(String args[]) throws IOException
    {
      //String tempFile = "/data/2/rpark/ip2_data/rpark/microbiome_mongodb_test/merck_brian_2017_10_19_15_231415/spectra/re-search/local_run/QE-C_Ecoli0_5ug_Piercemix_1-2.sqt";
      String tempFile = "/data/2/rpark/ip2_data/rpark/microbiome_mongodb_test/merck_brian_2017_10_19_15_231415/spectra/re-search/QE-C_Ecoli0_5ug_Piercemix_1-2.sqt";
         //SQTParser reader = new SQTParser( args[0] );
      SQTParser reader = new SQTParser( tempFile );

         SQTPeptide peptide;
         MLine mLine;

 	 //System.out.println("start...11");
      System.out.println("XCORR\tCS\tSEQ\tdeltaCN\tsp\tspRank\ttheoMass\texpMass\tdeltaMass\tppm\ttarget\tprotein");
         for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
         {
             peptide = itr.next();

             //System.out.println("each====>>" + peptide.getSLine());


             //for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
           for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
             {
           //      mLine = peptide.getMLine().next(); //
               mLine = mItr.next();

//                 System.out.println("seq===>>" + mLine.getSequence().replace('#', "");
               //System.out.println("seq===>>" + mLine);


              // if(mLine.getXcorrRank().equals("1")) continue;

               boolean isTarget=false;
               double xcorr = Double.parseDouble(mLine.getXcorr());



                 System.out.print(xcorr+ "\t" + peptide.getChargeState() +
                 "\t" + mLine.getSequence() + "\t" + mLine.getDeltCN() + "\t" + mLine.getSp() + "\t" +
                   + mLine.getSpRankInt() +
                 "\t" + mLine.getCalMZ() + "\t" + peptide.getCalMZ()  + "\t" );

               double deltaMass = (Double.parseDouble(mLine.getCalMZ()) - Double.parseDouble(peptide.getCalMZ()));

               deltaMass = Math.abs(deltaMass);
               deltaMass = deltaMass - Math.floor(deltaMass);


               System.out.print(deltaMass + "\t");
               System.out.print(deltaMass*1000000/Double.parseDouble(peptide.getCalMZ()) + "\t");



               /*
               if(peptide.getChargeState().trim().equals("1")) {
                 if(xcorr<1.8) continue;

               } else if(peptide.getChargeState().trim().equals("2")) {
                 if(xcorr<2) continue;
               } else {
                 if(xcorr<3) continue;
               }
*/

               StringBuffer sb = new StringBuffer();
                 for(Iterator<String> lItr = mLine.getLLine(); lItr.hasNext(); )
                 {
                   String lLine = lItr.next();
                   if(!lLine.startsWith("Rev"))
                     isTarget=true;
                   sb.append(lLine).append("\t");
                     //System.out.print(lLine + "\t");
                 }

                 System.out.print(isTarget + "\t" + sb.toString());

                 System.out.println("");

                 /*
                 if(isTarget)
                   System.out.print("F\t");
               else
                   System.out.print("R\t");
*/
              // System.out.println(mLine.getXcorr() + "\t" + mLine.getDeltCN() + "\t" + peptide.getChargeState() + "\t" + peptide.getHiScan());


             }
         }
    }

    public SQTParser(String fileName) throws IOException, NullPointerException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        init();
    }

    private void init() throws IOException, NullPointerException
    {
        br = new BufferedReader(new FileReader(fileName));

        readHeader();
    }

    private void readHeader() throws IOException, NullPointerException
    {
        lastLine = br.readLine();
        while(lastLine != null && !lastLine.startsWith("S")) {
            lastLine = br.readLine();
        }
/*
        if(lastLine.startsWith("H\t"))
        {
            while ((lastLine = br.readLine()).startsWith("H\t")); //read head line
            while ((lastLine = br.readLine()).equals("")); //read empty line
        }
        else
        {
            while (!(lastLine = br.readLine()).startsWith("S\t"));
        }
//System.out.println("lastLine: " + lastLine);
*/
    }

    public Iterator<SQTPeptide> getSQTPeptide() throws IOException
    {
        return new Iterator<SQTPeptide> ()
        {
            private SQTPeptide peptide;
            private MLine mLine;

            public boolean hasNext()
            {
                return lastLine != null && lastLine.startsWith("S");
            }

            public SQTPeptide next()
            {
                try
		{
                    peptide = getPeptide();
		    while(!isValidPeptide)
		    {
			isValidPeptide = true;
			peptide = getPeptide();
		    }
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return peptide;
            }

            public void remove()
            {
                throw new UnsupportedOperationException("Not supported");
            }

            private SQTPeptide getPeptide() throws IOException
            {

//                String[] strArr = lastLine.split("\t");

		//System.out.println("last line==>>" + lastLine);

		while( (lastLine != null) && lastLine.equals("") && !lastLine.startsWith("S\t") )
			lastLine = br.readLine();
//System.out.println(lastLine);
                peptide = new SQTPeptide(lastLine);


//System.out.println("Reached here. lastLine: " + lastLine);
                java.util.ArrayList<MLine> list = new java.util.ArrayList<MLine>();

                while ( ( (lastLine = br.readLine()) != null && !lastLine.startsWith("S")))
                {
                    //sqt file has wrong deltaCN value.  So, grep the right value from it.
//System.out.println("Reached here too. lastLine: " + lastLine);
                    boolean isDeltaCN=false;
                    float deltaCN=0;
                    float temp;

                    list.clear();

                    while(lastLine!=null && lastLine.startsWith("M\t"))
                    {
			//System.out.println(lastLine);
			/*  If M line is messed up in sqt file, discard it */
//System.out.println("lastLine: " + lastLine);
			String[] arr = lastLine.split("\t");
			if(arr.length<10)
			{
				isValidPeptide = false;
				break;
			}

                        mLine = new MLine(arr);

                        if(!isDeltaCN)
                        {
                            temp = Float.parseFloat( mLine.getDeltCN() );

                            if(temp>0)
                            {
                                deltaCN=temp;
                                isDeltaCN=true;
                            }
                        }

                        while ( (lastLine = br.readLine()) != null &&  lastLine.startsWith("L\t") )
                        {
			//System.out.println("---->>" + lastLine);
                            mLine.addLLine(lastLine);
                        }

                        list.add(mLine);
                    }

                    for(int i=0;i<list.size();i++)
                    {
                        list.get(i).setDeltCN(deltaCN);
                        peptide.addMLine(list.get(i));
                    }

                    //peptide.addMLine(mLine);
                    if(lastLine==null || lastLine.startsWith("S\t"))
                       break;
                }

////System.out.println("lastLine: " + lastLine);
                return peptide;
            }
        };
    }

    public void close() throws IOException
    {
        br.close();
    }

}
