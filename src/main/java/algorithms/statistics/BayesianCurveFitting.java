package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import no.uib.cipr.matrix.NotConvergedException;
import java.security.NoSuchAlgorithmException;
import java.util.logging.Logger;

/**
 * Class to preform Bayesian regression prediction of data points given a test data set.
 * <pre>
 *     Usage:
 *     The training data and test data are both transformed into a polynomial order feature matrix.
 *
 *     The training data are modeled using a Gaussian prior and likelihood.
 *
 *     to generate the training data polynomial feature matrix:
 *         double[][] phiX = BayesianCurveFitting.generatePhiX(xTrain, m);
 *
 *     to generate the test polynomial feature matrix:
 *         double[][] phiXTest = BayesianCurveFitting.generatePhiX(xTest, m);
 *
 *     the labels for the training data are double[] t;
 *
 *     perform the regression fit:
 *         ModelFit fit = BayesianCurveFitting.fit(phiX, t, alpha, beta);
 *
 *     predict the labels, given phiXTest:
 *         ModelPrediction prediction = BayesianCurveFitting.predict(fit, phiXTest);
  </pre>
 * The methods are adapted from github project PRML code
 * https://github.com/ctgk/PRML/blob/main/prml/linear/_bayesian_regression.py
 *
 * Their license is:
 *
 * MIT License
 *
 * Copyright (c) 2018 ctgk
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
public class BayesianCurveFitting {

    private static final Logger log = Logger.getLogger(BayesianCurveFitting.class.getSimpleName());

    /**
     * fit the training data
     * the likelihood and the prior are both Gaussian.
     * @param phiX training data feature matrix
     * @param t training data label vector
     * @param alpha the precision of the prior. the precision is the reciprocal of the variance.
     * @param beta the precision of the likelihood.  the precision is the reciprocal of the variance
     * @return the model fit data structures.
     */
    public static ModelFit fit(double[][] phiX, double[] t, final double alpha, final double beta) throws NotConvergedException {

        double[][] priorPrecision = MatrixUtil.createIdentityMatrix(phiX[0].length);
        MatrixUtil.multiply(priorPrecision, alpha);

        double[] priorMean = new double[phiX[0].length];

        return fit(phiX, t, alpha, beta, priorMean, priorPrecision);
    }

    /**
     * fit the training data.
     * the likelihood and the prior are both Gaussian.
     * @param phiX training data feature matrix
     * @param t training data label vector
     * @param alpha the precision of the prior. the precision is the reciprocal of the variance.
     * @param beta the precision of the likelihood.  the precision is the reciprocal of the variance
     * @param priorMean the prior mean
     * @param priorCov the prior covariance
     * @return the model fit data structures
     */
    public static ModelFit fit(double[][] phiX, double[] t, final double alpha, final double beta,
           final double[] priorMean, final double[][] priorCov) throws NotConvergedException {

        log.log(java.util.logging.Level.FINE, String.format("xTrain=phiX=\n%s",
                FormatArray.toString(phiX, "%.8e")));

        double[][] phiXT = MatrixUtil.transpose(phiX);

        //[m X m]
        double[][] sInv = calcSInv(priorCov, phiX, phiXT, alpha, beta);
        log.log(java.util.logging.Level.FINE, String.format("sInv=\n%s", FormatArray.toString(sInv, "%.4f")));

        // same result as calcSInv, difference in use of outer product
        //double[][] _sInv = _calcSInv(x, m, alpha, beta); where x is the original training data given to generator
        //log.log(java.util.logging.Level.FINE, String.format("_sInv=\n%s", FormatArray.toString(_sInv, "%.4f")));

        //[m X m]
        double[][] s = MatrixUtil.pseudoinverseRankDeficient(sInv);
        log.log(java.util.logging.Level.FINE, String.format("s=\n%s", FormatArray.toString(s, "%.4f")));

        // the mean is roughly similar to a slope term of yTrain/xTrain
        //[(m+1) X 1]
        double[] mean = calcMean(priorMean, priorCov, s, phiXT, t, alpha, beta);
        log.log(java.util.logging.Level.FINE, String.format("mean=\n%s", FormatArray.toString(mean, "%.4f")));

        //[m+1 X m+1]
        double[][] cov = s;

        log.log(java.util.logging.Level.FINE, String.format("covariance=\n%s", FormatArray.toString(cov, "%.4f")));

        log.log(java.util.logging.Level.FINE, String.format("phiX=\n%s", FormatArray.toString(phiX, "%.4f")));

        /*
         from Wasserman's "All of Statistics", 2.43 Theorem:
           X = mu + sqrt(sigma)*Z where mu is the mean vector,
           sigma is the covariance, and Z is N(O, I) which is the
           unit standard normal distribution.

        N(0, I) ~ Σ^(−1/2) * (X−μ)
        */

        ModelFit fit = new ModelFit();
        fit.phiX = phiX;
        fit.mean = mean;
        fit.cov = cov; // s
        fit.precision = sInv;
        fit.alpha = alpha;
        fit.beta = beta;

        return fit;
    }

    /**
     * predict the labels given the model and x test data.
     * @param fit model fit made from training data
     * @param phiXTest x value feature matrix for which to predict target t values.
     * @return
     */
    public static ModelPrediction predict(ModelFit fit, double[][] phiXTest) throws NotConvergedException {

        // [testX.length X (m+1)]  [m+1] = [testX.length]
        // use regression slope term:
        double[] y = MatrixUtil.multiplyMatrixByColumnVector(phiXTest, fit.mean);
        log.log(java.util.logging.Level.FINE, String.format("y=\n%s", FormatArray.toString(y, "%.4f")));

        //[testX.length X (m+1)] [(m+1)X(m+1)]  = [testX.length X (m+1)]
        double[][] ys1 = MatrixUtil.multiply(phiXTest, fit.cov);
        log.log(java.util.logging.Level.FINE, String.format("ys1=\n%s", FormatArray.toString(ys1, "%.4f")));
        //[testX.length X (m+1)] * [testX.length X (m+1)]  = [testX.length X (m+1)]
        ys1 = MatrixUtil.elementwiseMultiplication(ys1, phiXTest);
        log.log(java.util.logging.Level.FINE, String.format("ys1=\n%s", FormatArray.toString(ys1, "%.4f")));

        // propagation of errors:
        // sum ys1 along rows, add 1/beta, take sqrt:
        double[] yErr = new double[ys1.length];
        int j;
        for (int i = 0; i < ys1.length; ++i) {
            for (j = 0; j < ys1[i].length; ++j) {
                yErr[i] += ys1[i][j];
            }
            yErr[i] += (1./fit.beta);
            yErr[i] = Math.sqrt(yErr[i]);
        }
        log.log(java.util.logging.Level.FINE, String.format("yErr=\n%s", FormatArray.toString(yErr, "%.4f")));

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

        ModelPrediction model = new ModelPrediction();
        model.yFit = y;
        model.yErr = yErr;

        return model;
    }

    /**
     * predict values for test feature matrix phiXTest using randomly generated points
     * from a multivariate normal distribution based upon the model fit mean and covariance.
     * @param fit model fit made from training data
     * @param phiXTest x value feature matrix for which to predict target t values.
     * @return nSamples predicted from random samples of the model at the points from phiXTest.
     * the return matrix size is [nSamples X phiXTest.length] so that each row is a sample
     * generated for phiXTest.
     */
    public static double[][] predictRandomSample(ModelFit fit, double[][] phiXTest, int nSamples) throws NotConvergedException, NoSuchAlgorithmException {

        // [nSamples X fit.mean.length] = [nSamples X (M+1))]
        double[][] w = MultivariateNormalDistribution.sampleRandomlyFrom0(fit.mean, fit.cov, nSamples);

        // [phiXTest is [N2 X (M+1)]
        // [N2 X (M+1)] [(M+1) X nSamples]
        double[][] y = MatrixUtil.multiply(phiXTest, MatrixUtil.transpose(w));

        return MatrixUtil.transpose(y);
    }

    /**
     * predict values for test feature matrix phiXTest using randomly generated points
     * from a multivariate normal distribution based upon the model fit mean and covariance.
     * @param fit model fit made from training data
     * @param phiXTest x value feature matrix for which to predict target t values.
     * @return a sample predicted from random samples of the model at the points from phiXTest.
     */
    public static ModelPrediction predictRandomSample(ModelFit fit, double[][] phiXTest) throws NotConvergedException, NoSuchAlgorithmException {

        //  [(M+1))]
        double[] u = MultivariateNormalDistribution.sampleRandomlyFrom0(fit.mean, fit.cov);

        // [phiXTest is [N2 X (M+1)]
        // [N2 X (M+1)] [(M+1) X 1]
        double[] y = MatrixUtil.multiplyMatrixByColumnVector(phiXTest, u);

        //[testX.length X (m+1)] [(m+1)X(m+1)]  = [testX.length X (m+1)]
        double[][] ys1 = MatrixUtil.multiply(phiXTest, fit.cov);
        //[testX.length X (m+1)] * [testX.length X (m+1)]  = [testX.length X (m+1)]
        ys1 = MatrixUtil.elementwiseMultiplication(ys1, phiXTest);

        // propagation of errors:
        // sum ys1 along rows, add 1/beta, take sqrt:
        double[] yErr = new double[ys1.length];
        int j;
        for (int i = 0; i < ys1.length; ++i) {
            for (j = 0; j < ys1[i].length; ++j) {
                yErr[i] += ys1[i][j];
            }
            yErr[i] += (1./fit.beta);
            yErr[i] = Math.sqrt(yErr[i]);
        }
        log.log(java.util.logging.Level.FINE, String.format("yErr=\n%s", FormatArray.toString(yErr, "%.4f")));

        //from Wasserman's "All of Statistics", 2.43 Theorem:
        //N(0, I) ~ Σ^(−1/2) * (X−μ)

        ModelPrediction model = new ModelPrediction();
        model.yFit = y;
        model.yErr = yErr;

        return model;
    }

    /**
     * calculate the mean.
     * <pre>
     * the method is from eqn 3.5 of Bishop's PRML, following code in
     * https://github.com/ctgk/PRML/blob/main/prml/linear/_bayesian_regression.py
     * method fit().
     * </pre>
     *
     * @return size [(m+1)]
     */
    protected static double[] calcMean(final double[] priorMean, final double[][] priorPrecision,
            final double[][] s, final double[][] phiXT, final double[] t,  final double alpha, final double beta) throws NotConvergedException {

        /*
        solve for mean as x in a*x=b
        where
        a = sInv = precision_prev + beta * x_train.T @ x_train
                 = (alpha * I_(m+1)) + beta * (phiXT * phiX)
        b = precision_prev @ mean_prev + beta * x_train.T @ y_train
          = (alpha * I_(m+1)) * zero_(m+1) + beta * (phiXT * t)

        same as x = pinv(a) * b = s * b
        */

        double[] bV0 = MatrixUtil.multiplyMatrixByColumnVector(priorPrecision, priorMean);

        double[] bV = MatrixUtil.multiplyMatrixByColumnVector(phiXT, t);
        MatrixUtil.multiply(bV, beta);

        bV = MatrixUtil.add(bV0, bV);

        double[] mean = MatrixUtil.multiplyMatrixByColumnVector(s, bV);

        return mean;
    }

    /**
     * create a matrix generated from vector: φi(x) = xi for i = 0,...,M.
     <pre>
     ɸ(x)=  | x[0]    (x[0])^1    ...  (x[0])^m    |
            | x[1]    (x[1])^1    ...  (x[1])^m    |
            | ...                                  |
            | x[n-1]  (x[n-1])^1  ...  (x[n-1])^m  |
     </pre>
     * @param x array of data points
     * @param m order of polynomial fit
     * @return a matrix where each column is x to the order j where j is the column number
     * and ranges from 0 to m, inclusive.
     * output size is [N X (M+1)] where M is m and N is x.length.
     */
    public static double[][] generatePhiX(double[] x, int m) {

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
        return out;
    }

    /**
     * calculate S^-1=alpha*I + beta * (xT*x).
     * <pre>
     * This method
     * </pre>
     * @param x phi(x)
     * @param xT phi(x)^T
     * @param alpha noise to add to the generated gaussian.
     * @param beta the precision, that is, the reciprocal of the variance.
     * @return S^-1 = alpha * I + beta * (xT*x)
     */
    protected static double[][] calcSInv(final double[][] priorPrecision,
                                         final double[][] x, final double[][] xT, final double alpha, final double beta) {

        // eqn 1.72 from Bishop's PRML has an errata.
        //     corrected by the errata to: S^-1=alpha*I + beta * summation_from_n=1_to_N(phi(x_n) * phi(x_n)^T)
        //     the equation is corrected in Bishop's eqn 3.51 and in the
        // ctgk github PRML CODE : S^-1=alpha*I + beta * (xT*x)

        double[][] p1 = MatrixUtil.multiply(xT, x);
        MatrixUtil.multiply(p1, beta);

        double[][] sInv = MatrixUtil.elementwiseAdd(priorPrecision, p1);

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
    private static double[] generatePhiXn(int m, double xn) {
        double[] phiXn = new double[m + 1];
        phiXn[0] = 1;// xn^0=1
        for (int i = 1; i <= m; ++i) {
            phiXn[i] = phiXn[i - 1] * xn;
        }
        return phiXn;
    }

}
