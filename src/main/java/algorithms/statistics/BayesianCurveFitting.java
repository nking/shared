package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.NotConvergedException;

import java.util.Arrays;

public class BayesianCurveFitting {

    /**
     *
     * the likelihood and the prior are both Gaussian.
     * @param x training data x
     * @param t training data labels, a.k.a. targets
     * @param m the polynomial order found to fit the data (see ElastiNet for example).
     * @param alpha noise to add when generating predications for t
     * @param beta the precision, a.k.a. reciprocal of the variance.
     * @param predictionX x values for which to predict target t values.
     * @return
     */
    public static double[] predictUsingGaussianPrior(final double[] x, final double[] t, final int m, final double alpha, final double beta,
         double[] predictionX) throws NotConvergedException {

        int n = x.length;

        //[m X m]
        double[][] sInv = calcSInv(x, m, alpha, beta);
        //[m X m]
        double[][] s = MatrixUtil.pseudoinverseRankDeficient(sInv);
        //[n X 1]
        double[][] mean = calcMean(s, x, t, m, alpha, beta);
        //[n X n]
        double[][] cov = calcVariance(s, x, m, alpha, beta);

        /*
        Wasserman's eqn 2.10 from All of Statistics:
          f(x, mu, sigma) =
                1 / ( (2*p1)^(n/2) * |sigma|^(1/2) * exp( -(0.5)*(x-mu)^T * (sigma)^-1 * (x-mu) )
                where sigma is the covariance, mu is the mean,
                and |sigma| is the determinant of sigma

         from Wasserman's "All of Statistics", 2.43 Theorem:
        <pre>
           X = mu + sqrt(sigma)*Z where mu is the mean vector,
           sigma is the covariance, and Z is N(O, I) which is the
           unit standard normal distribution.
        </pre>
        */

        //from Wasserman's "All of Statistics", 2.43 Theorem:
        //N(0, I) ~ Σ^(−1/2) * (X−μ)
        double[][] sqrtCov = MatrixUtil.squareRoot(cov); //[n X n]
        double[][] invSqrtCov = MatrixUtil.pseudoinverseRankDeficient(sqrtCov);

        double[] meanV = MatrixUtil.transpose(mean)[0];
        double[] xMinusMean = MatrixUtil.subtract(x, meanV);

        double[] n0 = MatrixUtil.multiplyMatrixByColumnVector(invSqrtCov, xMinusMean);

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
     * implementing eqn 1.72 of Bishop's PRML.
     * @param x training data points
     * @param m order of polynomial originally used to fit the training data x,t.
     * @param alpha noise to add to the generated gaussian.
     * @param beta the precision, that is, the reciprocal of the variance.
     * @return size [x.length X x.length]
     */
    static double[][] calcVariance(final double[][] s, final double[] x, final int m, final double alpha, final double beta) throws NotConvergedException {

        /*
        create matrix phi(x)^T where each column is x and then the entries in it are taken to the power of the column number.
        result is size [x.length X m].
        */
        double[][] xMatrixT = createOrderMatrixTransposed(x, m);
        double[][] xMatrix = MatrixUtil.transpose(xMatrixT);
        double[][] tmp2 = MatrixUtil.multiply(xMatrixT, s);
        tmp2 = MatrixUtil.multiply(tmp2, xMatrix);

        // add 1/beta to tmp2
        double invBeta = 1./beta;
        int i;
        int j;
        for (i = 0; i < tmp2.length; ++i) {
            for (j = 0; j < tmp2[i].length; ++j) {
                tmp2[i][j] += beta;
            }
        }

        return tmp2;
    }

    /**
     * implementing eqn 1.72 of Bishop's PRML.
     * @param x training data points
     * @param t training data values
     * @param m order of polynomial originally used to fit the training data x,t.
     * @param alpha noise to add to the generated gaussian.
     * @param beta the precision, that is, the reciprocal of the variance.
     * @return size []
     */
    static double[][] calcMean(final double[][] s, final double[] x, final double[] t, final int m, final double alpha, final double beta) throws NotConvergedException {

        //[x.length X m]
        double[][] xMatrixT = createOrderMatrixTransposed(x, m);

        //[x.length x (m+1))] = [x.length x (m+1))] [(m+1) x (m+1)]
        double[][] out = MatrixUtil.multiply(xMatrixT, s);
        MatrixUtil.multiply(out, beta);

        double[] phiXn;

        double[] sum = new double[m + 1];
        for (int i = 0; i < x.length; ++i) {
            // phiXn is [(m+1) X 1]
            phiXn = generatePhiXn(m, x[i]);
            MatrixUtil.multiply(phiXn, t[i]);
            sum = MatrixUtil.add(phiXn, sum);
        }
        // [(m+1) X 1] in 2D array
        double[][] sumM = new double[m+1][];
        for (int i = 0; i < sum.length; ++i) {
            sumM[i] = new double[]{sum[i]};
        }

        out = MatrixUtil.multiply(out, sumM);

        return out;
    }

    /**
     * create a matrix where each column is x to the order j where j is the column number
     * and ranges from 0 to m, inclusive.
     * for example:
     <pre>
              x[0]    (x[0])^1    ...  (x[0])^m
              x[1]    (x[1])^1    ...  (x[1])^m
              ...
              x[n-1]  (x[n-1])^1  ...  (x[n-1])^m
     </pre>
     * @param x array of data points
     * @param m order of polynomial fit
     * @return a matrix where each column is x to the order j where j is the column number
     *      * and ranges from 0 to m, inclusive.
     */
    private static double[][] createOrderMatrixTransposed(double[] x, int m) {

        double[][] out = new double[x.length][];
        int i;
        int j;
        for (i = 0; i < x.length; ++i) {
            out[i] = new double[m + 1];
            out[i][0] = 1;//(x[i])^0;
            out[i][1] = x[i];
        }

        for (i = 0; i < x.length; ++i) {
            for (j = 2; j <= m; ++j) {
                out[i][j] = out[i][j - 1] * x[i];
            }
        }
        System.out.printf("out=\n%s\n", FormatArray.toString(out, "%.2e"));
        return out;
    }

    /**
     * implementing eqn 1.72 of Bishop's PRML.
     * @param x training data points
     * @param m order of polynomial originally used to fit the training data x,t.
     * @param alpha noise to add to the generated gaussian.
     * @param beta the precision, that is, the reciprocal of the variance.
     * @return
     */
    protected static double[][] calcSInv(final double[] x, final int m, final double alpha, final double beta) {

        double[][] sInv = MatrixUtil.createIdentityMatrix(m+1);
        MatrixUtil.multiply(sInv, alpha);

        double[] phiXn;
        double[][] tmp2;

        double[][] tmp = MatrixUtil.zeros(m+1, m+1);
        for (int i = 0; i < x.length; ++i) {
            // phiXn is [mX1]
            phiXn = generatePhiXn(m, x[i]);
            tmp2 = MatrixUtil.outerProduct(phiXn, phiXn);
            assert(tmp2.length == tmp.length);
            assert(tmp2[0].length == tmp[0].length);
            tmp = MatrixUtil.elementwiseAdd(tmp, tmp2);
        }
        MatrixUtil.multiply(tmp, beta);
        sInv = MatrixUtil.elementwiseAdd(sInv, tmp);

        return sInv;
    }

    /**
     * generate vector phiXn from data point xn up to order m where phiXn is
     *   ɸ_i(x) = x^i for i=0 to m.
     *   see comments under eqn 1.72 of Bishop's PRML.
     * @param m the largest order of the polynomial used to fit training data x, t.
     * @param xn point n in vector of training data x.
     * @return a vector of length m+1 composed of elements xn^0, xn^1, ... xn^m;
     */
    protected static double[] generatePhiXn(int m, double xn) {
        double[] phiXn = new double[m + 1];
        phiXn[0] = xn;
        for (int i = 1; i <= m; ++i) {
            phiXn[i] = phiXn[i - 1] * xn;
        }
        return phiXn;
    }
}
