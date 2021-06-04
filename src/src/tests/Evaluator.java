package src.tests;

import src.functions.NDimFunction;
import src.functions.QuadraticFunction;
import src.newton.*;

import java.util.Arrays;
import java.util.List;

public class Evaluator {
    protected final static double EPS = 0.000001;

    public void evaluate()  {
        double[] x0 = {10,
                10};

        NDimFunction func = new QuadraticFunction(new double[][]{{20, 20}, {20, 40}}, new double[]{10, 10}, 10);

        List<Method> methods = List.of(
                // ok
                new ClassicNewtonMethod(func, x0, EPS),
                // ok
                new LinarySearchNewtonMethod(func, x0, EPS),
                // ok
                new MarquardtMethod1(func, x0, EPS),
                // ok
                new NewtonMethodWithDescentDirection(func, x0, EPS),
                // ok
                new PowellMethod(func, x0, EPS),
                // ok
                new BFGS(func, x0, EPS),
                // ok
                new MarquardtMethod2(func, x0, EPS)
        );

        for (Method method : methods) {
            System.out.println(Arrays.toString(method.findMinimum()));
        }
    }

    public static void main(String[] args) {
        Evaluator evaluator = new Evaluator();
        evaluator.evaluate();
    }
}