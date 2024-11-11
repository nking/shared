package algorithms.statistics;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;
import thirdparty.scipy.optimization.ElasticNet;

import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class BayesianCurveFittingTest extends TestCase {

    private final Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public BayesianCurveFittingTest(String testName) {
        super(testName);
    }

    public void testSampleFromSine() throws NotConvergedException, NoSuchAlgorithmException {

        log.log(java.util.logging.Level.INFO, String.format("\ntestSampleFromSine"));

        // unit test made from Fig 1.17 of Bishop's PRML
        /*
        x = (0., 6./60, 13./60, 21./60, 27./60, 34./60, 40./60, 47./60, 54./60, 60./60)
        y = (3.5/9, 7.5/9, 9./9, 8.5/9, -1./9, -1.5/9, -8./9, -4./9, -5./9, 2.5.9)
        alpha = 5E-3
        beta = 11.1
        M = order of fit = 9 for this example
        changing to 8 to have a different dimension

        [junit] INFO: x avg and stdev = 0.5033, 0.3396
        [junit] INFO: t avg and stdev = 0.1278, 0.6657
                      1./0.1278 = 7.82
        [junit] INFO: cov=
        [junit] 0.1153, -0.1467
        [junit] -0.1467, 0.4431
        */
        double[] xTrain = new double[]{0., 6. / 60, 13. / 60, 21. / 60, 27. / 60, 34. / 60, 40. / 60, 47. / 60, 54. / 60, 60. / 60};
        double[] xTrainC = new double[]{(0. - 30.) / 60., (6. - 30.) / 60., (13. - 30.) / 60, (21. - 30.) / 60, (27. - 30.) / 60,
                (34. - 30.) / 60, (40. - 30.) / 60, (47. - 30.) / 60, (54. - 30.) / 60, (60. - 30.) / 60};

        double[] t = new double[]{3.5 / 9, 7.5 / 9, 9. / 9, 8.5 / 9, -1. / 9, -1.5 / 9, -8. / 9, -4. / 9, -5. / 9, 2.5 / 9};

        double[] xTest = new double[]{-0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0, 1.1};
        double[] xTestC = MatrixUtil.subtract(xTest, 0.5);

        final double alpha = 5e-3;
        final double beta = 11.1;
        final int m = 8;//9;

        double[][] phiX = BayesianCurveFitting.generatePolynomialPhiX(xTrain, m);
        double[][] phiXTest = BayesianCurveFitting.generatePolynomialPhiX(xTest, m);

        ModelFit fit = BayesianCurveFitting.fit(phiX, t, alpha, beta);

        ModelPrediction prediction = BayesianCurveFitting.predict(fit, phiXTest);

        double tol = 1E-3;
        double[] yExpected = new double[]{0.2343, 0.8219, 0.9126, 0.7776, 0.4582, 0.02793,
                -0.4061, -0.7102, -0.4186, 0.2561, 1.02878};
        double[] yErrExpected = new double[]{0.48093, 0.3475, 0.34986, 0.3489, 0.3433, 0.3435,
                0.3489, 0.3506, 0.3836, 0.4207, 1.6300};

        assertEquals(prediction.getYFit().length, yExpected.length);
        assertEquals(prediction.getYErr().length, yErrExpected.length);

        double diff;
        for (int i = 0; i < yExpected.length; ++i) {
            diff = prediction.getYFit()[i] - yExpected[i];
            assertTrue(Math.abs(diff) < tol);
            diff = prediction.getYErr()[i] - yErrExpected[i];
            assertTrue(Math.abs(diff) < tol);
        }

        // do the same for the zero-centered data.  note that they haven't been divided by standard deviation:
        double[][] phiXC = BayesianCurveFitting.generatePolynomialPhiX(xTrainC, m);
        double[][] phiXTestC = BayesianCurveFitting.generatePolynomialPhiX(xTestC, m);
        ModelFit fitC = BayesianCurveFitting.fit(phiXC, t, alpha, beta);
        ModelPrediction predictionC = BayesianCurveFitting.predict(fitC, phiXTestC);

        for (int i = 0; i < yExpected.length; ++i) {
            diff = predictionC.getYFit()[i] - yExpected[i];
            if (!(Math.abs(diff) < 2.5*prediction.getYErr()[i])) {
                System.out.printf("diff=%f, 2.5*err = %f\n", diff, 2.5*prediction.getYErr()[i]);
                System.out.flush();
            }
            assertTrue(Math.abs(diff) < 2.5*prediction.getYErr()[i]);
        }

	double sigF = 3.0;

        ModelPrediction predictionRandom = BayesianCurveFitting.predictRandomSample(fitC, phiXTestC);
        for (int i = 0; i < yExpected.length; ++i) {
            diff = predictionRandom.getYFit()[i] - yExpected[i];
            if (!(Math.abs(diff) < sigF*prediction.getYErr()[i])) {
                System.out.printf("diff=%f, (%.1f)*err = %f\n", diff, sigF, sigF*prediction.getYErr()[i]);
                System.out.flush();
            }
            assertTrue(Math.abs(diff) < sigF*prediction.getYErr()[i]);
        }

        // to see which of the (m+1) coefficients fit the data best, see the non-zero coefficients in:
        double alphaE =1e-15;//0.01
        double l1_ratio = 0.1;
        ElasticNet en = new ElasticNet(alphaE, l1_ratio);
        en.fit(phiX, t);
        double[] coef = en.getCoef();
        System.out.println("coef=" + Arrays.toString(coef));

        // for alpha=0.01, l1_ratio=0.7:
        //coef=[0.0, -0.7765111080212174, -1.5473265172779087, -0.0, -0.0, 0.0, 0.0, 0.5574769984932583, 1.071223968473462]
        //    [junit] -
        // alpha=1E-15, l1=0.1-0.9:
        //coef=[0.0, 5.788417822392037, -17.620929974823493, 5.174997856609265, 5.4221837073367585, 2.6105252550572944, 0.41117118941185976, -0.7737875928696337, -1.2109419254553988]

        /*
        https://towardsdatascience.com/ridge-lasso-and-elasticnet-regression-b1f9c00ea3a3
        LASSO (L1 regularization)
            regularization term penalizes absolute value of the coefficients
            sets irrelevant values to 0
            might remove too many features in your model
        Ridge regression (L2 regularization)
            penalizes the size (square of the magnitude) of the regression coefficients
            enforces the B (slope/partial slope) coefficients to be lower, but not 0
            does not remove irrelevant features, but minimizes their impact

         */

    }

    public void estSampleFrom() throws NotConvergedException, NoSuchAlgorithmException {

        // unit test from
        // https://jamesmccaffrey.wordpress.com/2017/11/03/example-of-calculating-a-covariance-matrix/

        double[][] a = new double[5][3];
        a[0] = new double[]{64.0,   580.0, 29.0};
        a[1] = new double[]{66.0,   570.0, 33.0};
        a[2] = new double[]{68.0,   590.0, 37.0};
        a[3] = new double[]{69.0,   660.0, 46.0};
        a[4] = new double[]{73.0,   600.0, 55.0};

        double[] m = new double[]{68.0,   600.0, 40.0};
        double[][] cov = new double[3][3];
        cov[0] = new double[]{11.50,    50.00,   34.75};
        cov[1] = new double[]{50.00,  1250.00,  205.00};
        cov[2] = new double[]{34.75,   205.00,  110.00};

        double[][] _cov = BruteForce.covariance(a);
        log.log(java.util.logging.Level.INFO, String.format(
                "_cov=\n%s", FormatArray.toString(_cov, "%.1f")));
        log.log(java.util.logging.Level.INFO, String.format(
                "cov=\n%s", FormatArray.toString(cov, "%.1f")));

        // compare Bayesian result to:
        //double[] x = MultivariateNormalDistribution.sampleFrom0(m, cov);

    }

}
