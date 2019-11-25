import java.io.*;

public class SingleDimensionScore {
    String	Name;
    int		NTotal, NDirect, NDecoy, NTrue, NValidated;
    int		NAboveCutoff;
    int		NBins;
    int[]	RankIndex;
    boolean[]	Positive;
    double      TPRateCutoff;
    double      ScoreCutoff, MaxScore, MinScore;
    double	ROCIntegral;
    double	BinWidth;
    double	DecoySigma, PositiveSigma, DecoyMean, PositiveMean;
    double[]	ROCCurve;
    double[]	Score;
    double[]	Probability;
    double[]	Confidence;
    double[]    TPRate;
    double[]	BinPoints;
    int[]	DirectHistogram;
    int[]	DecoyHistogram;
    int[]	PositiveHistogram;
    double[]	FitDecoyHistogram;
    double[]	FitPositiveHistogram;

    public SingleDimensionScore(String InputName, int InputNTotal, double InputScore[], boolean InputPositive[]) {
	int counter;
	double DecoySumSquare, PositiveSumSquare;

	Name = InputName;
	NTotal = InputNTotal;
	Score = new double[NTotal];
	Positive = new boolean[NTotal];
	RankIndex = new int[NTotal];
	NDirect = 0; NDecoy = 0;
	DecoyMean = 0; PositiveMean = 0; DecoySumSquare = 0; PositiveSumSquare = 0;
	for (counter = 0; counter < NTotal; counter++) {
	    Score[counter] = InputScore[counter];
	    Positive[counter] = InputPositive[counter];
	    if (Positive[counter]) {
		NDirect++;
		PositiveMean += Score[counter];
		PositiveSumSquare += (Score[counter] * Score[counter]);
	    }
	    else {
		NDecoy++;
		DecoyMean += Score[counter];
		DecoySumSquare += (Score[counter] * Score[counter]);
		PositiveMean -= Score[counter];
		PositiveSumSquare -= (Score[counter] * Score[counter]);
	    }
	}
	NTrue = NDirect - NDecoy;
	DecoyMean = DecoyMean/NDecoy;
	PositiveMean = PositiveMean/NTrue;
	DecoySigma = Math.sqrt((DecoySumSquare - NDecoy * DecoyMean * DecoyMean)/ (NDecoy - 1));
	PositiveSigma = Math.sqrt((PositiveSumSquare - NTrue * PositiveMean * PositiveMean)/ (NTrue - 1));
	NRUtils.QuickSortCollections(Score, RankIndex);
    }

    public void GetBins() {

	int counter, bincounter, rank;
	double TempScore;

	MaxScore = Score[NTotal - 1];
	MinScore = Score[0];
	BinWidth = 2 * (Score[3*NTotal/4] - Score[NTotal/4]) / Math.pow(NTotal, 1.0/3.0);
	TempScore = MinScore;
	bincounter = 1;
	while (TempScore + BinWidth/2 < MaxScore) {
	    TempScore = TempScore + BinWidth;
	    bincounter++;
	}
	NBins = bincounter;
	BinPoints = new double[NBins];
	DirectHistogram = new int[NBins];
        DecoyHistogram = new int[NBins];
	PositiveHistogram = new int[NBins];
        FitDecoyHistogram = new double[NBins];
        FitPositiveHistogram = new double[NBins];
	for (bincounter = 0; bincounter < NBins; bincounter++) {
	    BinPoints[bincounter] = MinScore + bincounter * BinWidth;
	    DirectHistogram[bincounter] = 0;
	    DecoyHistogram[bincounter] = 0;
	    FitDecoyHistogram[bincounter] = BinWidth * NDecoy * Math.exp(-Math.pow(BinPoints[bincounter] - DecoyMean, 2)/(2 * Math.pow(DecoySigma,2))) / (DecoySigma * Math.sqrt(2 * Math.PI));
	    FitPositiveHistogram[bincounter] = BinWidth * NTrue * Math.exp(-Math.pow(BinPoints[bincounter] - PositiveMean, 2)/(2 * Math.pow(PositiveSigma,2))) / (PositiveSigma * Math.sqrt(2 * Math.PI));
	}
	bincounter = 0;
	for (counter = 0; counter < NTotal; counter++) {
	    rank = RankIndex[counter];
	    while (Score[counter] > BinPoints[bincounter] + BinWidth/2) {
		bincounter++;
	    }
	    if (Positive[rank])
		DirectHistogram[bincounter]++;
	    else
		DecoyHistogram[bincounter]++;
	}
	for (bincounter = 0; bincounter < NBins; bincounter++) {
	    PositiveHistogram[bincounter] = DirectHistogram[bincounter] - DecoyHistogram[bincounter];
	}
    }

    public void GetROC() {

	int counter, rank, Nt, Nf;

	ROCCurve = new double[NDecoy];
	ROCIntegral = 0.0;
	Nt = 0; Nf = 0;
	for (counter = NTotal -1; counter > -1; counter--) {
	    rank = RankIndex[counter];
	    if (Positive[rank]) {
		Nt++;
	    }
	    else {
		Nf++;
		if (Nt < Nf)
		    ROCCurve[Nf - 1] = 0.0;
		else if (Nt - Nf > NTrue)
		    ROCCurve[Nf - 1] = 1.0;
		else
		    ROCCurve[Nf - 1] = 1.0 * (Nt - Nf) / NTrue;
		ROCIntegral += ROCCurve[Nf - 1];
	    }
	}
	ROCIntegral = ROCIntegral / NDecoy;
    }

    public void GetProbability() { 

        int counter, rank;
                                                                                              
        Probability = new double[NTotal];
	Confidence = new double[NTotal];
        for (counter = NTotal -1; counter > -1; counter--) {
            rank = RankIndex[counter];
            Confidence[rank] = (NTrue * Math.exp(-Math.pow(Score[counter] - PositiveMean, 2)/(2 * Math.pow(PositiveSigma,2))) / (PositiveSigma * Math.sqrt(2 * Math.PI))) /
                                   (NDecoy * Math.exp(-Math.pow(Score[counter] - DecoyMean, 2)/(2 * Math.pow(DecoySigma,2))) / (DecoySigma * Math.sqrt(2 * Math.PI)) + 
                                    NTrue * Math.exp(-Math.pow(Score[counter] - PositiveMean,
2)/(2 * Math.pow(PositiveSigma,2))) / (PositiveSigma * Math.sqrt(2 * Math.PI)));
	    Probability[rank] = (NTrue * Math.exp(-Math.pow(Score[counter] - PositiveMean,
2)/(2 * Math.pow(PositiveSigma,2))) / (PositiveSigma * Math.sqrt(2 * Math.PI))) /
                                (2 * NDecoy * Math.exp(-Math.pow(Score[counter] - DecoyMean, 2)/(2 * Math.pow(DecoySigma,2))) / (DecoySigma * Math.sqrt(2 * Math.PI)) + 
                                 NTrue * Math.exp(-Math.pow(Score[counter] - PositiveMean,
2)/(2 * Math.pow(PositiveSigma,2))) / (PositiveSigma * Math.sqrt(2 * Math.PI)));
        }
    }

    public void GetTPRate() {
	// Compute TPRate by counting elements
	int counter, rank, Nt, Nf;

	Nt = 0; Nf = 0;
	TPRate = new double[NTotal];
	for (counter = NTotal -1; counter > -1; counter--) {
	    rank = RankIndex[counter];
	    if (Positive[rank])
		Nt++;
	    else 
		Nf++;
            if (Nt > Nf)
                TPRate[rank] = 1.0 * (Nt - Nf) / Nt;
            else
                TPRate[rank] = 0.0;
	}

    }

    public void ChangeClassification(byte Classification[]) {

	int counter, rank, Nt;

	Nt = 0;
	for (counter = NTotal -1; counter > -1; counter--) {
	    rank = RankIndex[counter];
	    if (Positive[rank]) {
		Nt++;
		if (Nt > 2 * NTrue && Nt > NDirect/10 && Nt > 5)
		    Classification[rank] = 2;
	    }
	}
    }

    public void GetCutoffs(double TargetTPRateCutoff, boolean Validated[]) {
                                                                                
        int counter, rank, Nt, Nf, Nv;
        double CurrentTPRateCutoff;

	NValidated = 0;
	TPRateCutoff = TargetTPRateCutoff;
        Nt = 0; Nf = 0; Nv = 0; NAboveCutoff = 0;
        for (counter = NTotal -1; counter > -1; counter--) {
	    rank = RankIndex[counter];
            if (Positive[rank]) {
                Nt++;
		if (Validated[rank])
		    Nv++;
                CurrentTPRateCutoff = 1.0 * (Nt - Nf) / Nt;
                if (CurrentTPRateCutoff >= TPRateCutoff) {
                    NAboveCutoff = Nt;
		    ScoreCutoff = Score[counter];
		    NValidated = Nv;
                }
            }
            else {
                Nf++;
            }
        }
	this.DebugPrint();
    }

    public void DebugPrint() {
        System.out.println("Name is\t\t\t" + Name);
//	System.out.println("NTotal is\t\t" + NTotal);
//	System.out.println("NDirect is\t\t" + NDirect);
//	System.out.println("NDecoy is\t\t" + NDecoy);
//	System.out.println("TPRateCutoff is\t\t" + TPRateCutoff);
	System.out.println("ROCIntegral is\t\t" + ROCIntegral);
        System.out.println("ScoreCutoff is\t\t" + ScoreCutoff);
        System.out.println("NAboveCutoff is\t\t" + NAboveCutoff);
	System.out.println("NValidated is\t\t" + NValidated);
    }

    public void PlotROC() {
        File                    OutputFile;
        FileWriter              OutputFileWriter;
        BufferedWriter          Outgoing;
        String                  FileName;
	int counter;

	FileName = "ROC." + Name + ".dat";
	System.out.println("Writing ROC plot to " + FileName + "...");
        try {
            OutputFile = new File(FileName);
            OutputFileWriter = new FileWriter(OutputFile);
            Outgoing = new BufferedWriter(OutputFileWriter);
            for (counter = 0; counter < NDecoy; counter++) {
                Outgoing.write(new Double(1.0 * counter/NDecoy).toString() + "\t" + new Double(ROCCurve[counter]).toString() + "\n");
            }
            Outgoing.flush();
            Outgoing.close();
        }
        catch (IOException failure) {
            System.out.println("Something went wrong while writing " + FileName
+ ".");
        }
    }

    public void PlotHistograms() {
        File                    OutputFile;
        FileWriter              OutputFileWriter;
        BufferedWriter          Outgoing;
        String                  FileName;
        int counter;
                                                                                
        FileName = "Histograms." + Name + ".dat";
        System.out.println("Writing histogram plots to " + FileName + "...");
        try {
            OutputFile = new File(FileName);
            OutputFileWriter = new FileWriter(OutputFile);
            Outgoing = new BufferedWriter(OutputFileWriter);
            for (counter = 0; counter < NBins; counter++) {
                Outgoing.write(new Double(BinPoints[counter]).toString() + "\t" 
+ new Integer(DirectHistogram[counter]).toString() + "\t" 
+ new Integer(DecoyHistogram[counter]).toString() + "\t" 
+ new Integer(PositiveHistogram[counter]).toString() + "\t"
+ new Double(FitDecoyHistogram[counter]).toString() + "\t"
+ new Double(FitPositiveHistogram[counter]).toString() + "\n");
            }
            Outgoing.flush();
            Outgoing.close();
        }
        catch (IOException failure) {
            System.out.println("Something went wrong while writing " + FileName
+ ".");
        }
    }

    public void DebugPrintVerbose() {

	int counter;
	int rank;

	this.DebugPrint();
	System.out.println("Counter\tScore\tPositive\tTPRate\tConfidence\tProbability");
	for (counter = NTotal -1; counter > -1; counter--) {
	    rank = RankIndex[counter];
	    System.out.println(counter + "\t" + Score[counter] + "\t" + Positive[rank] + "\t" + TPRate[rank] + "\t" + Confidence[rank] + "\t" + Probability[rank]);
	}
    }
}

