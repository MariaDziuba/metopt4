package src.newton;

import src.functions.NDimFunction;
import src.gauss.GaussSLAEMethod;
import src.utils.MatrixUtil;

import static src.utils.MatrixUtil.*;

/**
 * Второй вариант метода Марквардта
 */
public class MarquardtMethod2 extends AbstractMarquardt {

    public MarquardtMethod2(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps, new GaussSLAEMethod(), 2, 0);
    }

    @Override
    double[] iterate(double[] x) {
        double[] antiGradient = MatrixUtil.multiplyVectorOnScalar(func.getGrad(curX), -1);
        double[][] hessian = func.getHess(curX);

        double[] dir;
        do {
            dir = slaemethod.findSolution(sumMatrices(hessian, multiplyMatrixOnScalar(I, step)), antiGradient);
            if (dir != null) {
                break;
            }
            step = Math.max(1, beta * step);
        } while (true);

        curX = sumVectors(curX, dir);
        step /= beta;

        if (norm(dir) <= EPS) {
            breakFlag = true;
        }

        return curX;
    }
}