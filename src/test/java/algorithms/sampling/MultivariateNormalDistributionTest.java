package algorithms.sampling;

import algorithms.correlation.BruteForce;
import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 *
 * @author nichole
 */
public class MultivariateNormalDistributionTest extends TestCase {
    
    public MultivariateNormalDistributionTest(String testName) {
        super(testName);
    }
    
    // TODO: add positive definite test matrices from NIST:
    //    https://math.nist.gov/cgi-bin/matrix-query

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
                
        double[] x = MultivariateNormalDistribution.sampleFrom0(m, cov);
        
        assertEquals(m.length, x.length);
        
        double[] eStdv = new double[]{
            Math.sqrt(11.5), Math.sqrt(1250), Math.sqrt(110)
        };
        
        double diff, chk;
        int i, j;
        for (i = 0; i < x.length; ++i) {
        }
        
        // create k vectors of x
        int k = 5;
        double[][] x2 = new double[k][x.length];
        for (i = 0; i < k; ++i) {
            x2[i] = MultivariateNormalDistribution.sampleFrom0(m, cov);
            System.out.printf("x2[%d]=%s\n", i, Arrays.toString(x2[i]));
        }
        System.out.println();
        
        // the derived covariance is smaller than cov:
        double[][] cov2 = BruteForce.covariance(x2);
        for (i = 0; i < cov2.length; ++i) {
            System.out.printf("cov2[%d]=%s\n", i, Arrays.toString(cov2[i]));
        }
        System.out.flush();
        
        double[][] x3 = new double[k][x.length];
        for (i = 0; i < k; ++i) {
            x3[i] = MultivariateNormalDistribution.sampleFrom1(m, cov);
            System.out.printf("x3[%d]=%s\n", i, Arrays.toString(x3[i]));
        }
        System.out.println();
        System.out.flush();
        
        double[][] cov3 = BruteForce.covariance(x3);
        for (i = 0; i < cov3.length; ++i) {
            System.out.printf("cov3[%d]=%s\n", i, Arrays.toString(cov3[i]));
        }
        System.out.flush();
        
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
        
        System.out.printf("eVL0=%s\n", i, eVL0.toString());
        System.out.printf("eVR0=%s\n", i, eVL2.toString());
        System.out.printf("eVL2=%s\n", i, eVR0.toString());
        System.out.printf("eVR2=%s\n", i, eVR2.toString());
        System.out.flush();
        
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
        System.out.println("b=" + b.toString());
        //x = A\b
        DenseMatrix _x = new DenseMatrix(A.numColumns(), b.numColumns());
        A.solve(b, _x);
        // expect this to be xact:
        System.out.println("x = A\\b = " + _x.toString());
        
        //perturb the rigth-hand side a little:
        double[][] b_prime = Matrices.getArray(b);
        b_prime[0][0] += 0.01;
        b_prime[1][0] += -0.01;
        b_prime[2][0] += 0.01;
        
        //x_prime = A\b_prime
        DenseMatrix _x_prime = new DenseMatrix(A.numColumns(), b.numColumns());
        A.solve(b, _x_prime);
        // expect this to be xact:
        System.out.println("x_prime = A\\b_prime = " + _x_prime.toString());
    }
    
}
