package algorithms.statistics;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import algorithms.util.FormatArray;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.hash.TDoubleHashSet;
import junit.framework.TestCase;
import no.uib.cipr.matrix.*;

/**
 *
 * @author nichole
 */
public class MultivariateNormalDistributionTest extends TestCase {

    private final Logger log = Logger.getLogger(this.getClass().getSimpleName());
    private final Level level = Level.INFO;
    
    public MultivariateNormalDistributionTest(String testName) {
        super(testName);
    }
    
    // TODO: add positive definite test matrices from NIST:
    //    https://math.nist.gov/cgi-bin/matrix-query

    public void testCompare() throws NoSuchAlgorithmException, NotConvergedException {

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

        int n = m.length;

        double[] mean2 = new double[n];
        double[] stdev2 = new double[n];
        double[][] anorm = Standardization.standardUnitNormalization(a, mean2, stdev2);

        log.log(level, "mean=" + FormatArray.toString(mean2, "%.3f"));
        log.log(level, "stdev=" + FormatArray.toString(stdev2, "%.3f"));

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        log.log(level, "SEED=" + seed);
        rand.setSeed(seed);

        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, n);
        log.log(level, String.format("randomly chosen x from unit standard gaussian distr =%s\n",
                FormatArray.toString(u, "%.3f")));

        double[] x0 = MultivariateNormalDistribution._sampleFrom0(u, m, cov);
        double[] x1 = MultivariateNormalDistribution._sampleFrom1(u, m, cov);

        log.log(level, String.format("compare x0=%s\n", FormatArray.toString(x0, "%.3f")));
        log.log(level, String.format("compare x1=%s\n", FormatArray.toString(x1, "%.3f")));

        double diff;
        for (int i = 0; i < x0.length; ++i) {
            diff = Math.abs(x0[i] - x1[i]);
            // compare diff to stdev2[i] derived from distributions
            if (diff > 2.5*stdev2[i]) {
                log.log(level, String.format("difference between %f and %f is > 2.5*stdev where stdev=%.3f\n",
                        x0[i], x1[i], stdev2[i]));
            }
        }

    }

    public void testSampleFrom() throws NotConvergedException, NoSuchAlgorithmException {
        
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

        double[] x = MultivariateNormalDistribution.sampleRandomlyFrom0(m, cov);
        
        assertEquals(m.length, x.length);
        
        double[] eStdv = new double[]{
            Math.sqrt(11.5), Math.sqrt(1250), Math.sqrt(110)
        };
        
        int i;
        // create k vectors of x
        int k = 5;
        double[][] x2 = new double[k][x.length];
        for (i = 0; i < k; ++i) {
            x2[i] = MultivariateNormalDistribution.sampleRandomlyFrom0(m, cov);
            log.log(java.util.logging.Level.INFO, String.format( String.format(
                    "x2[%d]=%s", i, Arrays.toString(x2[i]))));
        }

        // the derived covariance is smaller than cov:
        double[][] cov2 = BruteForce.covariance(x2);
        for (i = 0; i < cov2.length; ++i) {
            log.log(java.util.logging.Level.INFO, String.format(
                    "cov2[%d]=%s", i, Arrays.toString(cov2[i])));
        }
        System.out.flush();
        
        double[][] x3 = new double[k][x.length];
        for (i = 0; i < k; ++i) {
            x3[i] = MultivariateNormalDistribution.sampleRandomlyFrom1(m, cov);
            log.log(java.util.logging.Level.INFO, String.format(
                    "x3[%d]=%s", i, Arrays.toString(x3[i])));
        }

        double[][] cov3 = BruteForce.covariance(x3);
        for (i = 0; i < cov3.length; ++i) {
            log.log(java.util.logging.Level.INFO, String.format(
                    "cov3[%d]=%s", i, Arrays.toString(cov3[i])));
        }

        // eigen values:
        SVD svd0 = SVD.factorize(new DenseMatrix(cov));
        double[] s0 = svd0.getS();
        
        // A * v = lambda * v
        EVD evd0 = EVD.factorize(new DenseMatrix(cov));
        DenseMatrix eVL0 = evd0.getLeftEigenvectors();
        DenseMatrix eVR0 = evd0.getRightEigenvectors();
        
        //principal directions from the sample mean are the directions of eigenvectors of C
        
        // left eigen vector of A*A^T:
        
        SVD svd2 = SVD.factorize(new DenseMatrix(cov2));
        double[] s2 = svd2.getS();
        
        EVD evd2 = EVD.factorize(new DenseMatrix(cov2));
        DenseMatrix eVL2 = evd2.getLeftEigenvectors();
        DenseMatrix eVR2 = evd2.getRightEigenvectors();

        log.log(java.util.logging.Level.INFO, String.format(
                "eVL0=%s", i, eVL0.toString()));
        log.log(java.util.logging.Level.INFO, String.format(
                "eVR0=%s", i, eVL2.toString()));
        log.log(java.util.logging.Level.INFO, String.format(
                "eVL2=%s", i, eVR0.toString()));
        log.log(java.util.logging.Level.INFO, String.format(
                "eVR2=%s", i, eVR2.toString()));

        // to examine the conditioning, following example in:
        // https://blogs.mathworks.com/cleve/2018/08/20/reviving-wilsons-matrix/
        //xact = ones(4,1);
        DenseMatrix xact = new DenseMatrix(3, 1);
        for (i = 0; i < 3; ++i) {
            xact.set(i, 0, 1.);
        }
        //b = A * xact
        DenseMatrix A = new DenseMatrix(a);
        DenseMatrix b = MatrixUtil.multiply(A, xact);
        log.log(java.util.logging.Level.INFO, String.format(
                "b=" + b.toString()));
        //x = A\b
        DenseMatrix _x = new DenseMatrix(A.numColumns(), b.numColumns());
        A.solve(b, _x);
        // expect this to be xact:
        log.log(java.util.logging.Level.INFO, String.format(
                "x = A\\b = " + _x.toString()));
        
        //perturb the righthand side a little:
        double[][] b_prime = MatrixUtil.convertToRowMajor(b);
        b_prime[0][0] += 0.01;
        b_prime[1][0] += -0.01;
        b_prime[2][0] += 0.01;
        
        //x_prime = A\b_prime
        DenseMatrix _x_prime = new DenseMatrix(A.numColumns(), b.numColumns());
        A.solve(b, _x_prime);
        // expect this to be xact:
        log.log(java.util.logging.Level.INFO, String.format(
                "x_prime = A\\b_prime = " + _x_prime.toString()));
    }

}
