package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import junit.framework.TestCase;

public class LeastSquaresTest extends TestCase {

    public void test0() throws Exception {
        /*//https://en.m.wikipedia.org/wiki/Ordinary_least_squares
        double[][] X = new double[][]{
                {1.47,1.50,1.52,1.55,1.57,1.60,1.63,1.65,1.68,1.70, 1.73,1.75,1.78,1.80,1.83}
        };
        X = MatrixUtil.transpose(X);
        double[] y = new double[]{
                52.21,53.12,54.48,55.84,57.20, 58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46
        };*/

        //https://www.math.umd.edu/~immortal/MATH401/book/ch_least_squares.pdf
        double[][] X = new double[][]{
                {1,2}, {1,1}, {1, -1}
        };
        double[] y = new double[]{6,4,1};

        double[] expBeta = new double[]{18./7, 23/14.};

        double[] betaEst = LeastSquares.ols(X, y);

        double tol = 1E-4;

        for (int i = 0; i < betaEst.length; ++i) {
            assertTrue(Math.abs(expBeta[i] - betaEst[i]) < tol);
        }

        int t = 2;
    }
}
