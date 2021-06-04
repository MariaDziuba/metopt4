package src.newton;

import src.functions.NDimFunction;
import src.gauss.SLAEMethod;

import java.util.stream.IntStream;

import static src.utils.MatrixUtil.makeI;

/**
 * Абстрактный класс для метода Марквардта
 */
public abstract class AbstractMarquardt extends AbstractMethod {

    /**
     * Единичная матрица
     */
    double[][] I;
    boolean breakFlag = false;
    double step;

    protected final SLAEMethod slaemethod;

    /**
     * Коэффициент для изменения лямбды
     */
    protected final double beta;

    /**
     * Коэффициент, на который увеличивается диагональ гессиана
     */
    protected final double lambda;

    public AbstractMarquardt(final NDimFunction func, final double[] x0, final double eps,
                             final SLAEMethod slaemethod, final double beta, final double lambda) {
        super(func, x0, eps);
        this.slaemethod = slaemethod;
        this.beta = beta;
        this.lambda = lambda;
    }

    @Override
    boolean cycleCondition() {
        return !breakFlag;
    }

    @Override
    void preCalc() {
        I = makeI(x0.length);
        curX = x0;
        step = lambda;
    }
}