package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import no.uib.cipr.matrix.NotConvergedException;

import java.util.Arrays;

public class LeastSquares {

    /**
     * calculate the ordinary least squares fit to the observed data x
     * and the response y.
     * each row of x is a sample of p measuments (a.k.a. variables, features).
     * @param x a matrix of size [n X p] where n is the number of sawmple
     *          and p is the number of variables.
     * @param y an array of length n holding the response (dependent vars).
     * @return an array of length p giving the ordinary least squares fit
     * such that y_i = X_i^T * output + eps_i
     * where output is this output array,
     * X_i^t is a transposed row of X, and eps_i is the error term
     * for row i.
     */
    public static double[] ols(double[][] x, double[] y) throws NotConvergedException {

        double[][] aInv = MatrixUtil.pseudoinverseFullColumnRank(x);
        double[] betaEst = MatrixUtil.multiplyMatrixByColumnVector(aInv, y);
        return betaEst;
        /*
        yEst = X * betaEst = Projection * y
        and P = X * aInv
         */
    }

}
