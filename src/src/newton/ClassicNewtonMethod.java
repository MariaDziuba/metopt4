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

    public ClassicNewtonMethod(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
        this.SLAEMethod = new GaussSLAEMethod();
    }

    @Override
    double[] iterate(double[] x) {
        prevX = curX;
        p = SLAEMethod.findSolution(func.getHess(prevX), multiplyVectorOnScalar(func.getGrad(prevX), -1));
        curX = sumVectors(prevX, p);
        return curX;
    }

    @Override
    boolean cycleCondition() {
        return norm(subtractVectors(curX, prevX)) < EPS || norm(p) < EPS;
    }

    @Override
    void preCalc() {
        prevX = x0;
        p = SLAEMethod.findSolution(func.getHess(prevX), multiplyVectorOnScalar(func.getGrad(prevX), -1));
        curX = sumVectors(prevX, p);
    }
}