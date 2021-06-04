package src.newton;


import src.functions.NDimFunction;
import src.functions.OneDimFunction;
import src.onedimmethods.GoldenRatio;
import src.onedimmethods.OneDimMethod;

import static src.utils.MatrixUtil.*;

/**
 * Абстрактный квазиньютоновский метод
 */
public abstract class AbstractQuasiMethod extends AbstractMethod {

    protected AbstractQuasiMethod(NDimFunction func, double[] x0, double eps) {
        super (func, x0, eps);
    }

    /**
     * Возвращает следующее приближение
     */
    protected double[] getNextX(NDimFunction func, double[] x0, double[] p) {
        double a = getLinearMinimum(func, x0, p);
        p = multiplyVectorOnScalar(p, a);
        return sumVectors(x0, p);
    }

    /**
     * Вычисление длины следующего шага одномерным поиском
     */
    protected double getLinearMinimum(NDimFunction func, double[] x, double[] p) {
        OneDimFunction f = a -> func.apply(sumVectors(x, multiplyVectorOnScalar(p, a)));
        OneDimMethod oneDimMethod = new GoldenRatio(f, 0, 10, EPS);
        return oneDimMethod.findMin();
    }
}
