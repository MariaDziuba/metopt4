package src.newton;

import src.functions.NDimFunction;
import src.gauss.SLAEMethod;
import src.gauss.GaussSLAEMethod;

import static src.utils.MatrixUtil.*;

/**
 * Классический метод Ньютона
 */
public class ClassicNewtonMethod extends AbstractMethod {

    SLAEMethod SLAEMethod;
    double[] prevX;
    double[] p;
    double diff;

    int countIterations = 0;

    public ClassicNewtonMethod(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
        this.SLAEMethod = new GaussSLAEMethod();
    }

    @Override
    double[] iterate(double[] x) {
        prevX = curX;
        p = SLAEMethod.findSolution(func.getHess(prevX), multiplyVectorOnScalar(func.getGrad(prevX), -1));
        curX = sumVectors(prevX, p);
        countIterations++;
        diff = norm(subtractVectors(curX, prevX));
        return curX;
    }

    @Override
    boolean cycleCondition() {
        return diff > EPS || countIterations == 0;
    }

    @Override
    void preCalc() {
        prevX = x0;
        p = SLAEMethod.findSolution(func.getHess(prevX), multiplyVectorOnScalar(func.getGrad(prevX), -1));
        curX = sumVectors(prevX, p);
    }
}