package src.newton;

import src.functions.NDimFunction;
import src.gauss.GaussSLAEMethod;

import static src.utils.MatrixUtil.*;

/**
 * Первый вариант метода Марквадрта
 */
public class MarquardtMethod1 extends AbstractMarquardt {

    public MarquardtMethod1(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps, new GaussSLAEMethod(), 2, 10000);
    }

    @Override
    double[] iterate(double[] x) {
        double[] antiGrad = multiplyVectorOnScalar(func.getGrad(curX), -1);
        double[][] hess = func.getHess(curX);
        double[] dir = slaemethod.findSolution(sumMatrices(hess, multiplyMatrixOnScalar(I, step)), antiGrad);
        double[] nextX = sumVectors(curX, dir);
        double fx = func.apply(curX);
        double fNext = func.apply(nextX);
        while (fNext > fx) {
            step *= beta;
            dir = slaemethod.findSolution(sumMatrices(hess, multiplyMatrixOnScalar(I, step)), antiGrad);
            nextX = sumVectors(curX, dir);
            fNext = func.apply(nextX);
        }

        step /= beta;
        curX = nextX;
        if (norm(dir) <= EPS) {
            breakFlag = true;
        }
        return curX;
    }
}
