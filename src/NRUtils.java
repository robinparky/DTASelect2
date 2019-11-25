import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.regex.*;

public class NRUtils {

	public static void PrintSymmetricMatrix(double a[][], int size) {
		int i, j;

		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				System.out.print(a[i][j] + "\t");
			}
			System.out.print("\n");
		}
	}

	public static double GaussJordan(double a[][], int n) {
		int indxc[] = new int[n];
		int indxr[] = new int[n];
		int ipiv[] = new int[n];
		int i, icol, irow, j, k, l, ll;
		double big, dum, pivinv, temp;

		icol = 0;
		irow = 0;
		for (j = 0; j < n; j++)
			ipiv[j] = 0;

		for (i = 0; i < n; i++) { /* Main loop over columns */
			big = 0.0;
			for (j = 0; j < n; j++) { /* Search for pivot element */
				if (ipiv[j] != 1.0) {
					for (k = 0; k < n; k++) {
						if (ipiv[k] == 0) {
							if (Math.abs(a[j][k]) >= big) {
								big = Math.abs(a[j][k]);
								irow = j;
								icol = k;
							}
						} else if (ipiv[k] > 1)
							return 0.0;
					}
				}
			}
			++(ipiv[icol]);

			/*
			 * We now have a pivot element, so we interchange rows (relabel
			 * indexes)
			 */

			if (irow != icol) {
				for (l = 0; l < n; l++) {
					temp = a[irow][l];
					a[irow][l] = a[icol][l];
					a[icol][l] = temp;
				}
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (a[icol][icol] == 0.0)
				return 0.0;
			pivinv = 1.0 / a[icol][icol];
			a[icol][icol] = 1.0;
			for (l = 0; l < n; l++)
				a[icol][l] *= pivinv;

			/* Reduce the rows */

			for (ll = 0; ll < n; ll++) {
				if (ll != icol) {
					dum = a[ll][icol];
					a[ll][icol] = 0.0;
					for (l = 0; l < n; l++)
						a[ll][l] -= a[icol][l] * dum;
				}
			}
		} /* End of main loop */

		/* Now unscramble the permutation of the columns */

		for (l = n - 1; l > -1; l--) {
			if (indxr[l] != indxc[l]) {
				for (k = 1; k < n; k++) { // (k=1; k<=n; k++) { I GO BACK TO <=
					temp = a[k][indxr[l]];
					a[k][indxr[l]] = a[k][indxc[l]];
					a[k][indxc[l]] = temp;
				}
			}
		}
		return 1.0;
	}

	public static double LDA(double data_points[][], double projections[],
			double linear_coefficients[], int N_dim, int N_points,
			byte positive[]) {
		double[][] Total = new double[N_dim][N_dim];
		double[][][] CovarianceMatrix = new double[2][N_dim][N_dim];
		double[] DiffVector = new double[N_dim];
		double[] Coefficients = new double[N_dim];
		double[][] GroupMean = new double[2][N_dim];
		int[] N_groups = new int[2];
		int dim_counter, point_counter, i, j, k;
		double determinant, coefficient;

		// Form group means
		for (i = 0; i < N_dim; i++) {
			for (j = 0; j < 2; j++) {
				GroupMean[j][i] = 0.0;
			}
		}

		N_groups[0] = 0;
		N_groups[1] = 0;
		for (point_counter = 0; point_counter < N_points; point_counter++) {
			i = positive[point_counter];
			if (i != 2) {
				N_groups[i]++;
				for (j = 0; j < N_dim; j++) {
					GroupMean[i][j] += data_points[j][point_counter];
				}
			}
		}
		for (i = 0; i < 2; i++) {
			if (N_groups[i] == 0)
				return 0.0;
			for (j = 0; j < N_dim; j++) {
				GroupMean[i][j] = GroupMean[i][j] / N_groups[i];
			}
		}

		// Form covariance matrix for each group

		for (k = 0; k < 2; k++) {
			for (i = 0; i < N_dim; i++) {
				for (j = 0; j < N_dim; j++) {
					CovarianceMatrix[k][i][j] = 0.0;
				}
			}
		}
		for (i = 0; i < N_dim; i++) {
			for (j = 0; j < N_dim; j++) {
				Total[i][j] = 0.0;
				for (point_counter = 0; point_counter < N_points; point_counter++) {
					k = positive[point_counter];
					if (k != 2)
						CovarianceMatrix[k][i][j] += (data_points[i][point_counter] - GroupMean[k][i])
								* (data_points[j][point_counter] - GroupMean[k][j]);
				}
				CovarianceMatrix[0][i][j] = CovarianceMatrix[0][i][j]
						/ N_groups[0];
				CovarianceMatrix[1][i][j] = CovarianceMatrix[1][i][j]
						/ N_groups[1];
				Total[i][j] += CovarianceMatrix[0][i][j]
						+ CovarianceMatrix[1][i][j];
			}
		}

		// Invert variance-covariance matrix
		// System.out.println("Initial matrix:");
		// NRUtils.PrintSymmetricMatrix(Total, N_dim);
		determinant = NRUtils.GaussJordan(Total, N_dim);
		// System.out.println("Inverted once:");
		// NRUtils.PrintSymmetricMatrix(Total, N_dim);
		if (determinant == 0.0) {
			return determinant;
		}

		// Form difference vector of group mean vectors
		for (j = 0; j < N_dim; j++) {
			DiffVector[j] = GroupMean[0][j] - GroupMean[1][j];
		}

		// Determine linear coefficients
		for (i = 0; i < N_dim; i++) {
			Coefficients[i] = 0;
			for (j = 0; j < N_dim; j++) {
				Coefficients[i] += Total[i][j] * DiffVector[j];
			}
		}
		coefficient = Coefficients[0];
		System.out.print("\tCoefficients: ");
		for (i = 0; i < N_dim; i++) {
			Coefficients[i] = Coefficients[i] / coefficient;
			linear_coefficients[i] = Coefficients[i];
			System.out.print(Coefficients[i] + "  ");
		}
		System.out.print("\n");

		// Determine projections
		for (point_counter = 0; point_counter < N_points; point_counter++) {
			projections[point_counter] = 0.0;
			for (i = 0; i < N_dim; i++) {
				projections[point_counter] += data_points[i][point_counter]
						* Coefficients[i];
			}
		}
		return determinant;
	}

	public static void QuickSortIndex(double a[], int Index[], int lo0, int hi0) {
		int lo = lo0;
		int hi = hi0;
                
               // System.out.println(a.length + " " + Index.length + " " + lo0 + " " + hi0);
                
		if (lo >= hi) {
			return;
		} else if (lo == hi - 1) {
			/*
			 * sort a two element list by swapping if necessary
			 */
			if (a[lo] > a[hi]) {
				double T = a[lo];
				int TIndex = Index[lo];
				a[lo] = a[hi];
				Index[lo] = Index[hi];
				a[hi] = T;
				Index[hi] = TIndex;
			}
			return;
		}

		/*
		 * Pick a pivot and move it out of the way
		 */
		double pivot = a[(lo + hi) / 2];
		int PIndex = Index[(lo + hi) / 2];
		a[(lo + hi) / 2] = a[hi];
		Index[(lo + hi) / 2] = Index[hi];
		a[hi] = pivot;
		Index[hi] = PIndex;

		while (lo < hi) {
			/*
			 * Search forward from a[lo] until an element is found that is
			 * greater than the pivot or lo >= hi
			 */
			while (a[lo] <= pivot && lo < hi) {
				lo++;
			}

			/*
			 * Search backward from a[hi] until element is found that is less
			 * than the pivot, or lo >= hi
			 */
			while (pivot <= a[hi] && lo < hi) {
				hi--;
			}

			/*
			 * Swap elements a[lo] and a[hi]
			 */
			if (lo < hi) {
				double T = a[lo];
				int TIndex = Index[lo];
				a[lo] = a[hi];
				Index[lo] = Index[hi];
				a[hi] = T;
				Index[hi] = TIndex;
			}
		}

		/*
		 * Put the median in the "center" of the list
		 */
		a[hi0] = a[hi];
		Index[hi0] = Index[hi];
		a[hi] = pivot;
		Index[hi] = PIndex;

		/*
		 * Recursive calls, elements a[lo0] to a[lo-1] are less than or equal to
		 * pivot, elements a[hi+1] to a[hi0] are greater than pivot.
		 */
                
                
		QuickSortIndex(a, Index, lo0, lo - 1);
		QuickSortIndex(a, Index, hi + 1, hi0);
	}

	public static void QuickSortIndex(double a[], int Index[]) {
		int counter;

		for (counter = 0; counter < Index.length; counter++)
			Index[counter] = counter;
		QuickSortIndex(a, Index, 0, a.length - 1);
	}

	public static class Container
	{
		public final double value;
		public final int index;

		public Container(double value, int index) {
			this.value = value;
			this.index = index;
		}

		public double getValue() {
			return value;
		}
	}

	public static void QuickSortCollections(double a[], int Index[]) {
		int counter;

		List<Container> containerList  = new ArrayList<>();
		for (counter = 0; counter < Index.length; counter++)
			containerList.add(new Container(a[counter],counter));

		containerList.sort(Comparator.comparingDouble(Container::getValue));

		for(int i=0; i<containerList.size(); i++)
		{
			a[i] = containerList.get(i).value;
			Index[i] = containerList.get(i).index;
		}



	}


	public static String getAccession(String accession) {
		Pattern pattern = Pattern.compile("(.*)\\d");

		// NCBI, IPI, or others such as UNIT_PROT, SGD, NCI
		// accession = getDefline().substring( getDefline().indexOf('>')+1 );
		// accession = getDefline();
		// There are many corruptted sqt file. Ignore it.
		try {
			if (accession.startsWith("gi") && accession.contains("|")) // NCBI
			{
				String[] arr = accession.split("\\|");

				if (arr.length >= 4
						&& ("sp".equals(arr[2]) || "ref".equals(arr[2]) || "gb"
								.equals(arr[2])))
					accession = arr[3];
				else {
					arr = accession.split(" ");
					accession = arr[1];
				}

				// Accession # should end with digit. If accession # does not
				// end with digit,
				// grap next string (We assume this next one ends with digit.)
				/*
				 * if( pattern.matcher(arr[3]).matches() ) accession = arr[3];
				 * else accession = arr[4].substring(0, arr[4].indexOf(" "));
				 */

			} else if (accession.startsWith("sp")) // OMSSA format
			{
				String arr[] = accession.split("\\|");
				accession = arr[1];
			} else if (accession.startsWith("Reverse_sp")) // Reverse OMSSA
															// format
			{
				String arr[] = accession.split("\\|");
				accession = "Reverse_" + arr[1];
			} else if (accession.startsWith("IPI")) // IPI
			{
				String arr[] = accession.split("\\|");
				String subArr[] = arr[0].split(":");
				if (subArr.length > 1)
					accession = subArr[1];
				else
					accession = subArr[0];
			} else if (accession.startsWith("Reverse_")
					|| accession.startsWith("contam")
					|| accession.startsWith("Contam")) // Reverse database
			{
				int space = accession.indexOf(" ");
				int tab = accession.indexOf("\t");
				if (space < 0)
					space = 40;
				if (tab < 0)
					tab = 40;
				int index = (tab > space) ? space : tab;
				int end;
				if (index <= 0 || index >= 40) // no space
				{
					int length = accession.length();
					end = (length > 40) ? 40 : length;
				} else
					// cut by the first space
					end = index;
				accession = accession.substring(0, end);
			} else // UNIT_PROT, NCI or SGD
			{
				int space = accession.indexOf(" ");
				int pipe = accession.indexOf("|");
				int tab = accession.indexOf("\t");
				if (space < 0)
					space = 40;
				if (pipe < 0)
					pipe = 40;
				if (tab < 0)
					tab = 40;
				int tempindex = (tab > space) ? space : tab;
				int index = (tempindex > pipe) ? pipe : tempindex;
				int end;
				if (index <= 0 || index >= 40) // no space or pipe
				{
					int length = accession.length();
					end = (length > 40) ? 40 : length;
				} else
					// cut by the first space or pipe
					end = index;
				accession = accession.substring(0, end);
			}
		} catch (Exception e) {
			// System.out.println("No Correct Accession found, but this will be handled by MSP system."
			// + accession + " " + e);
			int i = accession.indexOf(" ");
			if (i < 0)
				return accession;
			else
				return accession.substring(0, i);
		}
		return accession;
	}
}
