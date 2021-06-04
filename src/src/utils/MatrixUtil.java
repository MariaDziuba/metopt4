package src.utils;

/**
 * Операции над матрицами
 */
public final class MatrixUtil {

    /**
     * Сложение двух векторов
     */
    public static double[] sumVectors(double[] vector1, double[] vector2) {
        double[] res = new double[vector1.length];

        for (int i = 0; i < vector1.length; ++i) {
            res[i] = vector1[i] + vector2[i];
        }
        return res;
    }

    /**
     * Сложение двух матриц
     */
    public static double[][] sumMatrices(double[][] matrix1, double[][] matrix2) {
        double[][] res = new double[matrix1.length][matrix1[0].length];

        for (int i = 0; i < matrix1.length; ++i) {
            for (int j = 0; j < matrix1.length; ++j) {
                res[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return res;
    }

    /**
     * Вычитание двух векторов
     */
    public static double[] subtractVectors(double[] vector1, double[] vector2) {
        double[] res = new double[vector1.length];

        for (int i = 0; i < vector1.length; ++i) {
            res[i] = vector1[i] - vector2[i];
        }
        return res;
    }

    /**
     * Вычитание двух матриц
     */
    public static double[][] subtractMatrices(double[][] matrix1, double[][] matrix2) {
        double[][] res = new double[matrix1.length][matrix1[0].length];
        for (int i = 0; i < matrix1.length; ++i) {
            for (int j = 0; j < matrix1.length; ++j) {
                res[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return res;
    }

    /**
     * Скалярное произведение векторов
     */
    public static double scalarProduct(double[] vector1, double[] vector2){
        double res = 0;
        for(int i = 0; i < vector1.length; i++){
            res += vector1[i]*vector2[i];
        }
        return res;
    }

    /**
     * Произведение вектора на число
     */
    public static double[] multiplyVectorOnScalar(double[] vector, double sc) {
        double[] res = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            res[i] = vector[i] * sc;
        }
        return res;
    }

    /**
     * Произведение вектора на вектор
     */
    public static double[][] multiplyVectors(double[] vector1, double[] vector2) {
        double[][] res = new double[vector1.length][vector1.length];
        for (int i = 0; i < vector1.length; i++) {
            for (int j = 0; j < vector1.length; j++) {
                res[i][j] += vector1[i] * vector2[j];
            }
        }
        return res;
    }

    /**
     * Произведение матрицы на число
     */
    public static double[][] multiplyMatrixOnScalar(double[][] matrix, double number) {
        double[][] res = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                res[i][j] += matrix[i][j] * number;
            }
        }
        return res;
    }

    /**
     * Произведение матрицы на вектор
     */
    public static double[] multiplyMatrixOnVector(double[][] matrix, double[] vector) {
        double[] res = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                res[i] += matrix[i][j] * vector[j];
            }
        }
        return res;
    }

    /**
     * Произведение двух матриц
     */
    public static double[][] multiplyMatrices(double[][] matrix1, double[][] matrix2) {
        double[][] resultMatrix = new double[matrix1.length][matrix1.length];
        for (int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1.length; j++) {
                for (int k = 0; k < matrix1.length; k++) {
                    resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return resultMatrix;
    }


    /**
     * Норма вектора
     */
    public static double norm(double[] vector) {
        double norm;
        double sum = 0;
        for (double v : vector) {
            sum += Math.pow(v, 2);
        }
        norm = Math.sqrt(sum);
        return norm;
    }

    /**
     * Транспонирование матрицы
     */
    public static void transposeMatrix(double[][] matrix) {
        if (matrix.length != matrix[0].length) {
            throw new IllegalArgumentException("Matrix should be quadratic");
        }
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i; j < matrix.length; ++j) {
                double temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
    }

    /**
     * Возвращает единичную матрицу
     */
    public static double[][] makeI(int length) {
        double[][] ans = new double[length][length];
        for (int i = 0; i < length; i++) {
            ans[i][i] = 1;
        }
        return ans;
    }
}
