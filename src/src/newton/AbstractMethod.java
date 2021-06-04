package src.newton;

import src.functions.NDimFunction;

/**
 * Абстрактный ньютоновский или квазиньютоновский метод
 */
abstract class AbstractMethod implements Method {

    NDimFunction func;
    double[] x0;
    double[] curX;
    double EPS;

    public AbstractMethod(NDimFunction func, double[] x0, double EPS) {
        this.func = func;
        this.x0 = x0;
        this.EPS = EPS;
    }

    @Override
    public double[] findMinimum() {
        preCalc();
        while (cycleCondition()) {
            curX = iterate(curX);
        }
        return curX;
    }

    abstract double[] iterate(double[] x);

    abstract boolean cycleCondition();

    abstract void preCalc();
}
