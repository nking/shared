package algorithms.statistics;

import algorithms.misc.MiscMath0;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.util.Arrays;

public class LinearRegressionTest extends TestCase {

    public void testSimpleLinearRegression() throws NotConvergedException {

        /*
        https://en.m.wikipedia.org/wiki/Simple_linear_regression
        Weight (m), xi	1.47	1.50	1.52	1.55	1.57	1.60	1.63	1.65	1.68	1.70	1.73	1.75	1.78	1.80	1.83
        Mass (kg), yi	52.21	53.12	54.48	55.84	57.20	58.57	59.93	61.29	63.11	64.47	66.28	68.10	69.92	72.19	74.46
         */
        double[] x = new double[]{
                1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70,
                1.73, 1.75, 1.78, 1.80, 1.83};

        double[] y = new double[] {
                52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47,
                66.28, 68.10, 69.92, 72.19, 74.46};

        double[] moments = Util.caldc2DMomentsX2Y2(x, y);
        assertTrue(Math.abs(moments[0] - 24.76) < 1E-2);
        assertTrue(Math.abs(moments[1] - 931.17) < 1E-2);
        assertTrue(Math.abs(moments[2] - 41.0532) < 1E-2);
        assertTrue(Math.abs(moments[3] - 58498.5439) < 1E-2);
        assertTrue(Math.abs(moments[4] - 1548.2453) < 1E-2);

        /*
        y-intercept,
        slope of the line fit,
     * variance of the error (where error is y - yIntercept - slope*x),
     * variance of the y-intercept,
     and variance of the slope
         */
        double[] results = LinearRegression.simpleLinearRegression(x, y);

        double expYInter = -39.062;
        double expSlope = 61.272;
        double expVarError = 0.5762;
        double expVarYInter = 8.63185;
        double expVarSlope = 3.1539;

        assertTrue(Math.abs(expYInter - results[0]) < 1E-3);
        assertTrue(Math.abs(expSlope - results[1]) < 1E-3);
        assertTrue(Math.abs(expVarError - results[2]) < 1E-3);
        assertTrue(Math.abs(expVarYInter - results[3]) < 1E-3);
        assertTrue(Math.abs(expVarSlope - results[4]) < 1E-3);

        double pearsonCor = Covariance.correlationPearson(x, y);
        assertTrue(Math.abs(0.9946 - pearsonCor) < 1E-3);

        double[][] x2d = new double[x.length][1];
        for (int i = 0; i < x.length; ++i) {
            x2d[i][0] = x[i];
        }
        double[] b0 = LinearRegression.linearRegression(y, x2d);
        assertTrue(Math.abs(expYInter - b0[0]) < 1E-3);
        assertTrue(Math.abs(expSlope - b0[1]) < 1E-3);

    }
}
