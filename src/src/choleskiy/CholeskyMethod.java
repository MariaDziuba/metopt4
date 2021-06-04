package src.choleskiy;


import src.gauss.SLAEMethod;
import src.utils.MatrixUtil;

import java.util.Arrays;

public class CholeskyMethod implements SLAEMethod {

    private final double EPS;

    public CholeskyMethod(final double EPS) {
        this.EPS = EPS;
    }


    private void forwardStep(final double[][] L, final double[] b) {
        for (int i = 0; i < L.length; ++i) {
            for (int j = 0; j < i; ++j) {
                b[i] -= b[j] * L[i][j];
            }
            b[i] /= L[i][i];
        }
    }


    private void backSubstitution(final double[][] Lt, final double[] y) {
        for (int i = Lt.length - 1; i >= 0; --i) {
            for (int j = Lt.length - 1; j > i; --j) {
                y[i] -= y[j] * Lt[i][j];
            }
            y[i] /= Lt[i][i];
        }
    }


    @Override
    public double[] findSolution(final double[][] A, final double[] B) {
        CholeskyDecomposition cholesky = new CholeskyDecomposition();
        double[][] L = cholesky.decompose(A, A.length);

        double[][] Lt = Arrays.stream(L).map(line -> Arrays.copyOf(line, line.length)).toArray(double[][]::new);
        MatrixUtil.transposeMatrix(Lt);

        double[][] checkProduct = MatrixUtil.multiplyMatrices(L, Lt);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                if (Double.isNaN(checkProduct[i][j]) || Math.abs(A[i][j] - checkProduct[i][j]) > EPS) {
                    return null;
                }
            }
        }

        double[] answer = Arrays.copyOf(B, B.length);

        forwardStep(L, answer);
        backSubstitution(Lt, answer);

        return answer;
    }
}