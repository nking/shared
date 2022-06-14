package algorithms.statistics;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

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

        [junit] INFO: x avg and stdev = 0.5033, 0.3396
        [junit] INFO: t avg and stdev = 0.1278, 0.6657
                      1./0.1278 = 7.82
        [junit] INFO: cov=
        [junit] 0.1153, -0.1467
        [junit] -0.1467, 0.4431
        */
        double[] x = new double[]{0., 6./60, 13./60, 21./60, 27./60, 34./60, 40./60, 47./60, 54./60, 60./60};
        double[] t = new double[]{3.5/9, 7.5/9, 9./9, 8.5/9, -1./9, -1.5/9, -8./9, -4./9, -5./9, 2.5/9};
        final double alpha = 5e-3;
        final double beta = 11.1;
        final int m = 9;

        double[][] a = new double[2][];
        a[0] = Arrays.copyOf(x, x.length);
        a[1] = Arrays.copyOf(t, t.length);
        a = MatrixUtil.transpose(a);
        double[][] covA = BruteForce.covariance(a);
        double[] xAvgStdev = MiscMath0.getAvgAndStDev(x);
        double[] tAvgStdev = MiscMath0.getAvgAndStDev(t);
        double[] meanV = new double[]{xAvgStdev[0], tAvgStdev[0]};
        double[] stdvV = new double[]{xAvgStdev[1], tAvgStdev[1]};
        double[] xtSample = MultivariateNormalDistribution.sampleRandomlyFrom0(meanV, covA);


        double[][] sInv = BayesianCurveFitting.calcSInv(x, m, alpha, beta);
        double[][] s = MatrixUtil.pseudoinverseRankDeficient(sInv);

        double[][] meanM = BayesianCurveFitting.calcMean(s, x, t, m, alpha, beta);

        double[][] varianceM = BayesianCurveFitting.calcVariance(s, x, m, alpha, beta);

        log.log(java.util.logging.Level.INFO, String.format(
                "x = %s", FormatArray.toString(x, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format(
                "t = %s", FormatArray.toString(t, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format(
                "x avg and stdev = %s", FormatArray.toString(xAvgStdev, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("t avg and stdev = %s",
                FormatArray.toString(tAvgStdev, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("brute force cov=\n%s", FormatArray.toString(covA, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("sInv=\n%s", FormatArray.toString(sInv, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("pseudoInv(sInv)=\n%s", FormatArray.toString(s, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("mean=\n%s", FormatArray.toString(meanM, "%.4f")));
        log.log(java.util.logging.Level.INFO, String.format("variance=\n%s", FormatArray.toString(varianceM, "%.4f")));


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
