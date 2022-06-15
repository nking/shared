package algorithms.statistics;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.security.NoSuchAlgorithmException;
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
        double[] xTrain = new double[]{0., 6./60, 13./60, 21./60, 27./60, 34./60, 40./60, 47./60, 54./60, 60./60};
        double[] xTrainC = new double[]{(0.-30.)/60., (6.-30.)/60., (13.-30.)/60, (21.-30.)/60, (27.-30.)/60,
                (34.-30.)/60, (40.-30.)/60, (47.-30.)/60, (54.-30.)/60, (60.-30.)/60};

        double[] t = new double[]{3.5/9, 7.5/9, 9./9, 8.5/9, -1./9, -1.5/9, -8./9, -4./9, -5./9, 2.5/9};

        double[] xTest = new double[]{-0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0, 1.1};
        double[] xTestC = MatrixUtil.subtract(xTest, 0.5);

        final double alpha = 5e-3;
        final double beta = 11.1;
        final int m = 8;//9;

        double[][] phiX = BayesianCurveFitting.generatePhiX(xTrain, m);
        double[][] phiXTest = BayesianCurveFitting.generatePhiX(xTest, m);

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
