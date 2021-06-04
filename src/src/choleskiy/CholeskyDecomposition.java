package src.choleskiy;

import java.util.Arrays;

public class CholeskyDecomposition {

    public double[][] decompose(final double[][] A, final int dimension) {
        double[][] L = Arrays.stream(A).map(line -> new double[dimension]).toArray(double[][]::new);

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                if (i == j) {
                    L[i][j] = Math.sqrt(A[i][i] - sum);
                } else {
                    L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
                }
            }
        }

        return L;
    }
}
