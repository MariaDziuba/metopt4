package src.newton;

import src.functions.NDimFunction;
import src.functions.OneDimFunction;
import src.onedimmethods.BrentSearch;
import src.onedimmethods.OneDimMethod;
import src.gauss.SLAEMethod;
import src.gauss.GaussSLAEMethod;

import static src.utils.MatrixUtil.*;

/**
 * Метод Ньютона с одномерным поиском
 */
public class LinarySearchNewtonMethod extends AbstractMethod {

    SLAEMethod SLAEMethod;
    double[] prevX;
    double[] p;

    public LinarySearchNewtonMethod(NDimFunction func, double[] x0, double eps) {
        super(func, x0, eps);
        this.SLAEMethod = new GaussSLAEMethod();
    }

    @Override
    double[] iterate(double[] x) {
        prevX = curX;

        p = SLAEMethod.findSolution(func.getHess(prevX), multiplyVectorOnScalar(func.getGrad(prevX), -1));

        double[] finalPrevX = prevX;
        double[] finalP = p;
        OneDimFunction fun = v -> func.apply(sumVectors(finalPrevX, multiplyVectorOnScalar(finalP, v)));
        OneDimMethod oneDimMethod = new BrentSearch(fun, -100.0, 100.0, EPS);//TODO: deal with borders
        double a = oneDimMethod.findMin();

        curX = sumVectors(prevX, multiplyVectorOnScalar(p, a));
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