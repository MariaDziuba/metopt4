package src.newton;

import src.functions.NDimFunction;
import src.functions.OneDimFunction;
import src.onedimmethods.BrentSearch;
import src.gauss.SLAEMethod;
import src.gauss.GaussSLAEMethod;
import src.utils.MatrixUtil;

import static src.utils.MatrixUtil.*;

/**
 * Метод Ньютона с направлением спуска
 */
public class NewtonMethodWithDescentDirection extends AbstractMethod {
    SLAEMethod SLAEMethod;
    double diff;
    int countIterations = 0;

    public NewtonMethodWithDescentDirection(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
        this.SLAEMethod = new GaussSLAEMethod();
    }

    @Override
    double[] iterate(double[] x) {
        double[] prevX = curX;
        double[] grad = func.getGrad(prevX);
        double[] antiGrad = multiplyVectorOnScalar(grad, -1);
        double[] p = SLAEMethod.findSolution(func.getHess(prevX), antiGrad);

        final double[] dir = MatrixUtil.scalarProduct(p, grad) >= 0 ? antiGrad : p;

        OneDimFunction f = alpha -> func.apply(MatrixUtil.sumVectors(prevX, multiplyVectorOnScalar(dir, alpha)));
        double alpha = new BrentSearch(f, 0, 10, EPS).findMin();
        curX = sumVectors(prevX, multiplyVectorOnScalar(dir, alpha));

        diff = MatrixUtil.norm(MatrixUtil.subtractVectors(curX, prevX));
        countIterations++;
        return curX;
    }

    @Override
    boolean cycleCondition() {
        return diff > EPS || countIterations == 0;
    }

    @Override
    void preCalc() {
        curX = x0;
    }
}