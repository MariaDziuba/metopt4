package src.newton;

import src.functions.NDimFunction;

import static src.utils.MatrixUtil.*;

public class PowellMethod extends AbstractQuasiMethod {
    double[][] C;
    double[] w;

    public PowellMethod(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
    }

    /**
     * Возвращает следующее приближение
     */
    private double[][] getNextC(double[][] C, double[] deltaX, double[] deltaW) {
        double k = 1 / scalarProduct(deltaW, deltaX);
        return subtractMatrices(C, multiplyMatrixOnScalar(multiplyVectors(deltaX, deltaX), k));
    }

    @Override
    double[] iterate(double[] x) {
        double[] p = multiplyMatrixOnVector(C, w);
        double[] nextX = getNextX(func, curX, p);
        double[] nextW = multiplyVectorOnScalar(func.getGrad(nextX), -1);
        // разность приближений
        double[] deltaX = subtractVectors(nextX, curX);
        // разность градиентов
        double[] deltaW = subtractVectors(nextW, w);
        C = getNextC(C, sumVectors(deltaX, multiplyMatrixOnVector(C, deltaW)), deltaW);
        curX = nextX;
        w = nextW;
        return curX;
    }

    @Override
    boolean cycleCondition() {
        return norm(w) > EPS;
    }

    @Override
    void preCalc() {
        curX = x0;
        C = makeI(x0.length);
        w = multiplyVectorOnScalar(func.getGrad(x0), -1);
    }
}