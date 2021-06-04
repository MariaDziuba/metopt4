package src.newton;

import src.functions.NDimFunction;

import static src.utils.MatrixUtil.*;

/**
 * Метод Бройдена-Флетчера-Шено
 */
public class BFGS extends AbstractQuasiMethod {

    int countIterations = 0;
    /**
     * Текущее приближение (матрицу, близкую к гессиану)
     */
    double[][] C;
    double[] grad;

    public BFGS(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
    }

    @Override
    double[] iterate(double[] x) {
        int i = 0;
        while (norm(grad) > EPS && i < curX.length) {
            double[] p = multiplyVectorOnScalar(multiplyMatrixOnVector(C, grad), -1);
            double[] nextX = getNextX(func, curX, p);
            double[] nextGrad = func.getGrad(nextX);
            C = getNextC(C, subtractVectors(nextX, curX), subtractVectors(nextGrad, grad));
            curX = nextX;
            grad = nextGrad;
            i++;
        }
        countIterations += i;
        C = makeI(curX.length);
        grad = func.getGrad(curX);
        return curX;
    }

    @Override
    boolean cycleCondition() {
        return norm(grad) > EPS || countIterations == 0;
    }

    @Override
    void preCalc() {
        curX = x0;
        C = makeI(x0.length);
        grad = func.getGrad(x0);
    }

    /**
     * Возвращает следующее приближение
     * @param s - разность приближений
     */
    private double[][] getNextC(double[][] C, double[] s, double[] y) {
        double p = 1 / scalarProduct(y, s);
        double[][] nextC = subtractMatrices(makeI(C.length), multiplyMatrixOnScalar(multiplyVectors(s, y), p));
        nextC = multiplyMatrices(nextC, C);
        nextC = multiplyMatrices(nextC, subtractMatrices(makeI(C.length), multiplyMatrixOnScalar(multiplyVectors(y, s), p)));
        nextC = sumMatrices(nextC, multiplyMatrixOnScalar(multiplyVectors(s, s), p));
        return nextC;
    }
}