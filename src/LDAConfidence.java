import java.io.*;
import java.util.*;

public class LDAConfidence {
	int ChargeState;
	boolean CalculateCoeffOnly;
	boolean ModifiedStatus;
	byte TrypticStatus;
	boolean UseModStats;
	boolean UseTrypticInfo;
	boolean UseChargeState, UseTrueProtein;
	String DecoyLabel, TrueProteinLabel;
	double[][] DataPoints;
	boolean[] Positive;
	boolean[] TrueProt;
	byte[] Classification;
	byte[] Validated;
	double[] Confidence;
	double[] Integral;
	double[] MDistance;
	double[] Probability;
	double[] CompositeScore;
	double[] BufferScore;
	double[] DeltaMass;
	double AverageMassDeviation = 0.0;
	int DeltaMassIndex = -1;
	String[] ParameterName;
	String FlagName;
	double[] LinearCoefficients;
	int[] Index;
	int NTotal, NValidated, NDirect, NReverse, NDim, Length;
	int SeqMinLength;
	int PlotLevel;
	double XCorr, LogNL, FPTol, determinant;
	boolean UseXCorr, UseDeltCN, UseSp, UseMassDifference, UseTertiaryScore, UseQuinary;
	boolean UseRawXCorr, UseLogMass, UseLogSpRank, UseIonStat;
	boolean isProtein = false;
	double Slope, Score, TempDeltaMass;

	public LDAConfidence(Protein LocusList, SelectCriteria Cutoffs) {
		Protein Runner;
		int counter;

		isProtein = true;

		ChargeState = 0;
		DecoyLabel = Cutoffs.DecoyLabel;
		// Determine NTotal
		NTotal = 0;
		Runner = LocusList;
		while (Runner.Next != null) {
			Runner = Runner.Next;
			NTotal++;
		}

		// Allocate arrays
		Confidence = new double[NTotal];
		Integral = new double[NTotal];
		MDistance = new double[NTotal];
		Probability = new double[NTotal];
		Positive = new boolean[NTotal];
		CompositeScore = new double[NTotal];
		Index = new int[NTotal];
		for (counter = 0; counter < NTotal; counter++)
			Index[counter] = counter;
		NDim = 1;
		DataPoints = new double[NDim][NTotal];
		UseTertiaryScore = Cutoffs.UseMultipleScore && Cutoffs.colCount>11;
		// Initialize DataPoints and Positive arrays
		Runner = LocusList;
		NDirect = 0;
		NReverse = 0;
		counter = 0;
		while (Runner.Next != null) {
			Runner = Runner.Next;
			DataPoints[0][counter] = Runner.GetRawProbs();
			if (Runner.Locus.startsWith(DecoyLabel)) {
				Positive[counter] = false;
				NReverse++;
			} else {
				Positive[counter] = true;
				NDirect++;
			}
			counter++;
		}
		this.PrintTotalProteins();
	}

	public LDAConfidence(OUTFile OutFileList, SelectCriteria Cutoffs,
			int Charge, boolean Modified, byte Tryptic, boolean Status) {
		int counter, dim_counter;
		OUTFile OUTRunner;
		String FileName;
		byte TempChargeState;

		CalculateCoeffOnly = Status;
		ChargeState = Charge;
		UseChargeState = Cutoffs.UseChargeState;
		ModifiedStatus = Modified;
		TrypticStatus = Tryptic;
		if (CalculateCoeffOnly) {
			UseModStats = false;
			UseTrypticInfo = false;
		} else {
			UseModStats = Cutoffs.UseModStats;
			UseTrypticInfo = Cutoffs.UseTrypticInfo;
		}
		FlagName = (UseChargeState ? new Integer(ChargeState).toString()
				: "all")
				+ "."
				+ (UseTrypticInfo ? new Byte(TrypticStatus).toString() : "all")
				+ "."
				+ (UseModStats ? new Boolean(ModifiedStatus).toString() : "all");
		DecoyLabel = Cutoffs.DecoyLabel;
		TrueProteinLabel = Cutoffs.TrueProteinLabel;
		FPTol = Cutoffs.FPTol;
		UseXCorr = Cutoffs.UseXCorr;
		UseDeltCN = Cutoffs.UseDeltCN;
		UseSp = Cutoffs.UseSp;
		UseMassDifference = Cutoffs.UseMassDifference;
		UseRawXCorr = Cutoffs.UseRawXCorr;
		UseLogSpRank = Cutoffs.UseLogSpRank;
		UseIonStat = Cutoffs.UseIonStat;
		UseLogMass = Cutoffs.UseLogMass;
		UseTrueProtein = Cutoffs.UseTrueProtein;
		SeqMinLength = Cutoffs.SeqMinLength;
		PlotLevel = Cutoffs.PlotLevel;

		// Determine NTotal
		counter = 0;
		NTotal = 0;
		FileName = "";
		OUTRunner = OutFileList;

		while (OUTRunner.OriginalNext != null) {
			OUTRunner = OUTRunner.OriginalNext;
			TempChargeState = OUTRunner.DTA.ChargeState;
			//System.out.println("Filename "+OUTRunner.DTA.FileName +"\t"+OUTRunner.DTA.ScanNumber);

			if (TempChargeState > Cutoffs.MaxStatisticsCharge)
				TempChargeState = Cutoffs.MaxStatisticsCharge;
			if ((!UseChargeState || TempChargeState == ChargeState)
					&& OUTRunner.DTA.TrimmedSequence().length() >= SeqMinLength
					&& !(OUTRunner.DTA.FileName.equals(FileName))
					&& (!UseModStats || OUTRunner.DTA.Modified == ModifiedStatus)
					&& (!UseTrypticInfo || OUTRunner.DTA.Tryptic == TrypticStatus)) {
				OUTRunner.OngoingStatisticsFlag = true;
				FileName = OUTRunner.DTA.FileName;
				NTotal++;
			}
		}

		// Allocate arrays
		Confidence = new double[NTotal];
		Integral = new double[NTotal];
		MDistance = new double[NTotal];
		BufferScore = new double[NTotal];
		CompositeScore = new double[NTotal];
		Positive = new boolean[NTotal];
		TrueProt = new boolean[NTotal];
		Probability = new double[NTotal];
		Classification = new byte[NTotal];
		Validated = new byte[NTotal];
		Index = new int[NTotal];


	//	UseQuinary = true;



		for (counter = 0; counter < NTotal; counter++)
			Index[counter] = counter;
		UseTertiaryScore = Cutoffs.UseMultipleScore && Cutoffs.colCount>11;
		// Statistics options
		NDim = 0;
		if (UseXCorr)
			NDim++;
		if (UseDeltCN)
			NDim++;
		if (UseMassDifference)
			NDim++;
		if (UseSp)
			NDim++;
		if (UseLogSpRank)
			NDim++;
		if (UseIonStat)
			NDim++;
		if (UseTertiaryScore)
			NDim++;
		if (UseQuinary)
			NDim+=2;
		if (NDim != 0) {
			DataPoints = new double[NDim][NTotal];
			LinearCoefficients = new double[NDim];
			for (dim_counter = 0; dim_counter < NDim; dim_counter++) {
				LinearCoefficients[dim_counter] = 0.0;
			}
			ParameterName = new String[NDim];
		} else {
			System.out
					.println("No statistics criteria selected! Use --nostats if you want that. Exiting...");
			System.exit(0);
		}

		// Initialize Positive array
		OUTRunner = OutFileList;
		counter = 0;
		NValidated = 0;
		NDirect = 0;
		NReverse = 0;
		while (OUTRunner.OriginalNext != null) {
			OUTRunner = OUTRunner.OriginalNext;
			if (OUTRunner.OngoingStatisticsFlag) {
				if (OUTRunner.Locus.startsWith(DecoyLabel)) {
					Positive[counter] = false;
					Classification[counter] = 1;
					NReverse++;
					TrueProt[counter] = false;
					Validated[counter] = 1;
				} else {
					Positive[counter] = true;
					Classification[counter] = 0;
					NDirect++;
					if (OUTRunner.Locus.startsWith(TrueProteinLabel)) {
						TrueProt[counter] = true;
						Validated[counter] = 0;
						NValidated++;
					} else {
						TrueProt[counter] = false;
						Validated[counter] = 2;
					}
				}
				counter++;
			}
		}

		dim_counter = 0;

		// Load MassDifference data
		if (UseMassDifference) {
			ParameterName[dim_counter] = (UseLogMass ? "LogDMass" : "DMass");
			DeltaMass = new double[NTotal];
			DeltaMassIndex = dim_counter;
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					if (UseLogMass)
						DataPoints[dim_counter][counter] = -Math.log(Math.max(
								Math.abs(OUTRunner.DTA.Adjusted_PPM_Offset),
								0.000001));
					else {
						DataPoints[dim_counter][counter] = -Math
								.abs(OUTRunner.DTA.Adjusted_PPM_Offset);
						DeltaMass[counter] = OUTRunner.DTA.Adjusted_PPM_Offset;
						if (UseTrueProtein && Validated[counter] == 0) {
							AverageMassDeviation += DeltaMass[counter];
						}
					}
					counter++;
				}
			}
			if (!UseLogMass && UseTrueProtein) {
				AverageMassDeviation = AverageMassDeviation / NValidated;
				System.out.println("\tMass Deviation is:\t"
						+ AverageMassDeviation + " PPM");
				for (counter = 0; counter < NTotal; counter++) {
					DeltaMass[counter] = DeltaMass[counter]
							- AverageMassDeviation;
					DataPoints[dim_counter][counter] = -Math
							.abs(DeltaMass[counter]);
				}
			}
			dim_counter++;
		}

		// Load XCorr data
		if (UseXCorr) {
			ParameterName[dim_counter] = (UseRawXCorr ? "XCorr" : "XCorrNorm");
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					XCorr = OUTRunner.DTA.XCorr;

					// Normalize XCorr data
					if (!UseRawXCorr) {
						TempChargeState = OUTRunner.DTA.ChargeState;
						Length = OUTRunner.DTA.TrimmedSequence().length();
						if (TempChargeState == 1) {
							LogNL = Math.log(1 * Math.min(10, Length));
						} else if (TempChargeState == 2) {
							LogNL = Math.log(2 * Math.min(15, Length));
						} else if (TempChargeState >= 3) {
							LogNL = Math.log(4 * Math.min(25, Length));
						}
						DataPoints[dim_counter][counter] = Math.log(Math.max(
								XCorr, 0.1)) / LogNL;
					} else {
						DataPoints[dim_counter][counter] = XCorr;
					}
					counter++;
				}
			}
			dim_counter++;
		}

		// Load DeltCN data
		if (UseDeltCN) {
			ParameterName[dim_counter] = "DeltCN";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = OUTRunner.DTA.DeltCN;
					counter++;
				}
			}
			dim_counter++;
		}

		// Load Sp data
		if (UseSp) {
			ParameterName[dim_counter] = "Sp";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = new Double(
							OUTRunner.DTA.SpScore).doubleValue();
					counter++;
				}
			}
			dim_counter++;
		}

		if (UseLogSpRank) {
			ParameterName[dim_counter] = "LogSpRank";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = -Math.log(new Double(
							OUTRunner.DTA.Sp).doubleValue());
					counter++;
				}
			}
			dim_counter++;
		}

		if (UseIonStat) {
			ParameterName[dim_counter] = "Ion";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = new Double(
							OUTRunner.DTA.IonProportion).doubleValue();
					counter++;
				}
			}
			dim_counter++;
		}

		if (UseTertiaryScore) {
			ParameterName[dim_counter] = "TertiaryScore";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = new Double(
							OUTRunner.DTA.tertiaryScore).doubleValue();
					counter++;
				}
			}
			dim_counter++;
		}
		if (UseQuinary) {
			ParameterName[dim_counter] = "PValue";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = new Double(
							OUTRunner.DTA.pvalue).doubleValue();
					counter++;
				}
			}
			dim_counter++;


			ParameterName[dim_counter] = "EValue";
			OUTRunner = OutFileList;
			counter = 0;
			while (OUTRunner.OriginalNext != null) {
				OUTRunner = OUTRunner.OriginalNext;
				if (OUTRunner.OngoingStatisticsFlag) {
					DataPoints[dim_counter][counter] = new Double(
							OUTRunner.DTA.evalue).doubleValue();
					counter++;
				}
			}
			dim_counter++;
		}
		this.PrintTotals();
	}

	public void ReportScores(Protein LocusList) {
		Protein Runner;
		int counter = -1;

		Runner = LocusList;
		while (Runner.Next != null) {
			Runner = Runner.Next;
			counter++;
			Runner.ProtFP = new Float(1.0 - Integral[counter]).floatValue();
			Runner.ProtConf = new Float(Confidence[counter]).floatValue();
		}
	}

	public void GetLinearCoefficients(double Coefficients[]) {
		int counter;
		for (counter = 0; counter < NDim; counter++) {
			LinearCoefficients[counter] = Coefficients[counter];
		}
	}

	public void SaveLinearCoefficients(double Coefficients[]) {
		int counter;
		for (counter = 0; counter < NDim; counter++) {
			Coefficients[counter] = LinearCoefficients[counter];
		}
	}

	public void CloseStatisticsFlag(OUTFile OutFileList) {
		OUTFile OUTRunner;

		OUTRunner = OutFileList;
		while (OUTRunner.OriginalNext != null) {
			OUTRunner = OUTRunner.OriginalNext;
			if (OUTRunner.OngoingStatisticsFlag) {
				OUTRunner.OngoingStatisticsFlag = false;
			}
		}
	}

	public void ReportScores(OUTFile OutFileList) {
		int counter = -1;
		OUTFile OUTRunner;
		String FileName = "";

		OUTRunner = OutFileList;
		while (OUTRunner.OriginalNext != null) {
			OUTRunner = OUTRunner.OriginalNext;
			if (OUTRunner.OngoingStatisticsFlag) {
				FileName = OUTRunner.DTA.FileName;
				counter++;
				OUTRunner.OngoingStatisticsFlag = false;
				OUTRunner.DTA.PepFP = new Float(1.0 - Integral[counter])
						.floatValue();
				OUTRunner.DTA.PepConf = new Float(Confidence[counter])
						.floatValue();
				OUTRunner.DTA.Probability = Probability[counter];
			} else {
				if (OUTRunner.DTA.FileName.equals(FileName)) {
					OUTRunner.DTA.PepFP = new Float(1.0 - Integral[counter])
							.floatValue();
					OUTRunner.DTA.PepConf = new Float(Confidence[counter])
							.floatValue();
					OUTRunner.DTA.Probability = Probability[counter];
				}
			}
		}
	}

	public void PrintTotals() {

		System.out.println("\tNumber of spectra:\t\t" + NTotal);
		System.out.println("\tMatching regular database:\t" + NDirect);
		System.out.println("\tMatching decoy database:\t" + NReverse);
		if (UseTrueProtein)
			System.out.println("\tMatching validated protein:\t" + NValidated);
	}

	public void PrintClassification() {
		int counter, Ngood, Nmaybe, Ndecoy;

		Ngood = 0;
		Nmaybe = 0;
		Ndecoy = 0;
		for (counter = 0; counter < NTotal; counter++) {
			if (Classification[counter] == 0)
				Ngood++;
			if (Classification[counter] == 1)
				Ndecoy++;
			if (Classification[counter] == 2)
				Nmaybe++;
		}
		System.out.print("Good: " + Ngood + "\t\t");
		System.out.print("Maybe: " + Nmaybe + "\t\t");
		System.out.print("Decoy: " + Ndecoy + "\n");
	}

	public void PrintTotalProteins() {

		System.out.println("\tNumber of proteins:\t" + NTotal);
		System.out.println("\tFrom regular database:\t" + NDirect);
		System.out.println("\tFrom decoy database:\t" + NReverse);
	}

	public void PlotDirectReverse(String FileName, int dim1, int dim2,
			boolean positive) {
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		int counter;

		System.out.println("Writing plotting information to " + FileName
				+ "...");
		try {
			OutputFile = new File(FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			for (counter = 0; counter < NTotal; counter++) {
				if (Positive[counter] == positive)
					Outgoing.write((DeltaMassIndex == dim1 ? new Double(
							DeltaMass[counter]).toString() : new Double(
							DataPoints[dim1][counter]).toString())
							+ "\t"
							+ (DeltaMassIndex == dim2 ? new Double(
									DeltaMass[counter]).toString()
									: new Double(DataPoints[dim2][counter])
											.toString()) + "\n");
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing " + FileName
					+ ".");
		}
	}

	public void PlotOutput() {
		File OutputFile;
		FileWriter OutputFileWriter;
		BufferedWriter Outgoing;
		String FileName;
		int counter;

		FileName = "Probability." + ChargeState + ".dat";
		try {
			OutputFile = new File(FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			for (counter = 0; counter < NTotal; counter++) {
				Outgoing.write(new Double(MDistance[counter]).toString() + "\t"
						+ new Double(Probability[Index[counter]]).toString()
						+ "\n");
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing " + FileName
					+ ".");
		}

		FileName = "FPRate." + ChargeState + ".dat";
		try {
			OutputFile = new File(FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			for (counter = 0; counter < NTotal; counter++) {
				if (Positive[Index[counter]]) {
					Outgoing.write(new Double(MDistance[counter]).toString()
							+ "\t"
							+ new Double(Integral[Index[counter]]).toString()
							+ "\n");
				}
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing " + FileName
					+ ".");
		}

		FileName = "Confidence." + ChargeState + ".dat";
		try {
			OutputFile = new File(FileName);
			OutputFileWriter = new FileWriter(OutputFile);
			Outgoing = new BufferedWriter(OutputFileWriter);
			for (counter = 0; counter < NTotal; counter++) {
				if (Positive[Index[counter]]) {
					Outgoing.write(new Double(MDistance[counter]).toString()
							+ "\t"
							+ new Double(Confidence[Index[counter]]).toString()
							+ "\n");
				}
			}
			Outgoing.flush();
			Outgoing.close();
		} catch (IOException failure) {
			System.out.println("Something went wrong while writing " + FileName
					+ ".");
		}
	}

	public void DebugPrintConf() {

		int counter;
		System.out
				.println("####################################################");
		this.PrintTotals();
		for (counter = 0; counter < NTotal; counter++) {
			System.out.println("RealDB " + Positive[counter] + "\tConfidence "
					+ Confidence[counter]);
		}
		System.out
				.println("----------------------------------------------------");
	}

	public void DebugPrintMD() {

		int counter;
		System.out
				.println("####################################################");
		this.PrintTotals();
		for (counter = 0; counter < NTotal; counter++) {
			System.out.println("RealDB " + Positive[counter] + "\tMDistance "
					+ MDistance[counter]);
		}
		System.out
				.println("----------------------------------------------------");
	}

	public void DebugPrintCompositeScore() {

		int counter;

		System.out.println("*****Regular order*****");
		for (counter = NTotal - 1; counter > -1; counter--) {
			System.out.println(counter + "\t" + CompositeScore[counter] + "\t"
					+ Positive[Index[counter]]);
		}
		System.out.println("*****End regular order*****");
	}

	public void ComputeConfidence() {
		int counter, dim_counter, dim1, dim2;
		SingleDimensionScore SingleScore;
		boolean NonZeroCoefficient = false;

		if (NTotal == 0) {
			System.out.println("Not data in this subset! Skipping...");
			return;
		}
		if (NDirect == 0 || NReverse == 0) {
			System.out
					.println("Not enough data for statistics on this set! Skipping...");
			return;
		}
		if (PlotLevel >= 1) {
			for (dim1 = 0; dim1 < NDim; dim1++) {
				for (dim2 = dim1 + 1; dim2 < NDim; dim2++) {
					PlotDirectReverse(ParameterName[dim1] + "."
							+ ParameterName[dim2] + ".Direct." + FlagName
							+ ".dat", dim1, dim2, true);
					PlotDirectReverse(ParameterName[dim1] + "."
							+ ParameterName[dim2] + ".Reverse." + FlagName
							+ ".dat", dim1, dim2, false);
				}
			}
		}
		// PrintClassification();
		for (dim_counter = 0; dim_counter < NDim; dim_counter++) {
			for (counter = 0; counter < NTotal; counter++) {
				BufferScore[counter] = DataPoints[dim_counter][counter];
			}
			SingleScore = new SingleDimensionScore(ParameterName[dim_counter]
					+ "." + FlagName, NTotal, BufferScore, Positive);
			if (PlotLevel >= 2) {
				SingleScore.GetROC();
				SingleScore.GetCutoffs(1.0 - FPTol, TrueProt);
				SingleScore.PlotROC();
			}
			if (PlotLevel >= 3) {
				SingleScore.GetBins();
				SingleScore.PlotHistograms();
			}
			SingleScore.ChangeClassification(Classification);
			// PrintClassification();
			SingleScore = null;
		}
		if (NDim > 1) {
			if (CalculateCoeffOnly) {
				if (UseTrueProtein)
					determinant = NRUtils.LDA(DataPoints, MDistance,
							LinearCoefficients, NDim, NTotal, Validated);
				else
					determinant = NRUtils.LDA(DataPoints, MDistance,
							LinearCoefficients, NDim, NTotal, Classification);
				if (determinant == 0.0) {
					System.out
							.println("Determinant is zero in inverse matrix calculation!");
					System.out
							.println("Skipping linear coefficients calculation!");
					return;
				}
			} else {
				NonZeroCoefficient = false;
				for (dim_counter = 0; dim_counter < NDim; dim_counter++) {
					if (LinearCoefficients[dim_counter] != 0.0) {
						NonZeroCoefficient = true;
					}
				}
				if (!NonZeroCoefficient) {
					System.out.println("Linear coefficients are zero!");
					System.out.println("Skipping statistics calculation!");
					System.out.println("Scores set to zero.");
					return;
				}
				for (counter = 0; counter < NTotal; counter++) {
					MDistance[counter] = 0.0;
					for (dim_counter = 0; dim_counter < NDim; dim_counter++) {
						MDistance[counter] += DataPoints[dim_counter][counter]
								* LinearCoefficients[dim_counter];
					}
				}
			}
			SingleScore = new SingleDimensionScore("Mahalanobis." + FlagName,
					NTotal, MDistance, Positive);
		} else {
			LinearCoefficients[0] = 1.0;
			for (counter = 0; counter < NTotal; counter++) {
				MDistance[counter] = BufferScore[counter];
			}
			SingleScore = new SingleDimensionScore(ParameterName[0] + "."
					+ FlagName, NTotal, MDistance, Positive);
		}
		if (PlotLevel >= 2) {
			SingleScore.GetROC();
			SingleScore.GetCutoffs(1.0 - FPTol, TrueProt);
			SingleScore.PlotROC();
		}
		if (PlotLevel >= 3) {
			SingleScore.GetBins();
			SingleScore.PlotHistograms();
		}
		SingleScore.GetTPRate();
		SingleScore.GetProbability();
		for (counter = 0; counter < NTotal; counter++) {
			Integral[counter] = SingleScore.TPRate[counter];
		}
		SingleScore = null;
	}
}
