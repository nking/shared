package algorithms.correlation;

import algorithms.matrix.MatrixUtil;
import java.util.Arrays;

/**
 * implemented from 
        
 * useful for Mahalanobis distance among many things.
 * 
 * @author nichole
 */
public class BruteForce {
    
    /**
     * calculate the correlation matrix for a using a brute force method
     * @param a an mxn matrix of data with the dimensions being columns
     * and the datum number being rows.
     * <pre>
     * e.g.  a[0] = new double[]{10, 100, 1000}
     *       a[1] = new double[]{ 9, 110, 900}
     * </pre>
     * @return the correlation matrix as a double array of size [a[0].length][a[0].length]
     */
    public static double[][] correlation(double[][] a) {
        
        double eps= 1.e-15;
        
        // cor_i_j = cov_i_j / (sqrt(var_i)*sqrt(var_j))
        
        double[][] cov = covariance(a);
        
        double[][] cor = new double[cov.length][cov[0].length];
        int i, j;
        for (i = 0; i < cov.length; ++i) {
            cor[i] = new double[cov[i].length];
        }
        double si, sj;
        for (i = 0; i < cor.length; ++i) {
            si = (cov[i][i] > eps) ? Math.sqrt(cov[i][i]) : 0;
            for (j = i; j < cor[i].length; ++j) {
                sj = (cov[j][j] > eps) ? Math.sqrt(cov[j][j]) : 0;
                if (si > eps && sj > eps) {
                    cor[i][j] = cov[i][j]/(si*sj);
                    if (i != j) {
                        cor[j][i] = cor[i][j];
                    }
                }
            }
        }
        
        return cor;
    }
    
    /**
     * calculate the covariance matrix for a using a brute force method.
     * The covariance matrix is also known as auto-covariance matrix, 
     * dispersion matrix, variance matrix, and the variance–covariance matrix.
     * @param a an mxn matrix of data with the dimensions being columns
     * and the datum number being rows.
     * <pre>
     * e.g.  a[0] = new double[]{10, 100, 1000}
     *       a[1] = new double[]{ 9, 110, 900}
     * </pre>
     *          Note, aside from unit tests, a quick comparison with python numpy cov function produces same result.
     *
     * @return the covariance matrix as a double array of size [a[0].length][a[0].length]
     */
    public static double[][] covariance(double[][] a) {
        
        int nRows = a.length;
        int nCols = a[0].length;
        
        int i, j;
        
        // mean of each column:
        double[] mean = MatrixUtil.mean(a);
        
        double[][] diffs = new double[nRows][];
        for (i = 0; i < nRows; ++i) {
            diffs[i] = new double[nCols];
            for (j = 0; j < nCols; ++j) {
                diffs[i][j] = (a[i][j] - mean[j]);
            }
        }
        
        /*System.out.printf("bf stand. means=%s%n", Arrays.toString(mean));
        
        System.out.flush();
        System.out.printf("bf stand. diffs=%n");
        for ( i = 0; i < diffs.length; ++i) {
            for ( j = 0; j < diffs[i].length; ++j) {
                System.out.printf("%11.3e  ", diffs[i][j]);
            }
            System.out.printf("%n");
        }
        System.out.flush();
        */
        
        double[][] cov = new double[nCols][];
        for (i = 0; i < nCols; ++i) {
            cov[i] = new double[nCols];
        }
     
        double sum;
        int ii;
        for (i = 0; i < nRows; ++i) {
            for (j = i; j < nCols; ++j) {
                sum = 0;
                // multiply diffs[*][i] by diffs[*][j]
                for (ii = 0; ii < nRows; ++ii) {
                    sum += (diffs[ii][i] * diffs[ii][j]);
                }
                sum /= ((double)(nRows-1.));
                cov[i][j] = sum;
                if (i != j) {
                    cov[j][i] = sum;
                }
            }
        }
        
        /*System.out.printf("bf stand. cov=%n");
        for ( i = 0; i < cov.length; ++i) {
            for ( j = 0; j < cov[i].length; ++j) {
                System.out.printf("%11.3e  ", cov[i][j]);
            }
            System.out.printf("%n");
        }
        System.out.flush();*/
        
        return cov;
    }
    
    /**
     * calculate the covariance matrix for a using simplest brute force method
     * Cov_a = (1/n) a * a^T for a being standardized to subtract the mean of
     * each type from the measurement of same type.
     * The covariance matrix is also known as auto-covariance matrix, 
     * dispersion matrix, variance matrix, and the variance–covariance matrix.
     * @param a an mxn matrix of data with the dimensions being rows
     * and the datum number being columns.  NOTE that a must be "zero mean" centered.
     * <pre>
     * e.g.  a[0] = new double[]{10,  9}
     *       a[1] = new double[]{99,  100}
     *       a[2] = new double[]{1000, 900}
     * </pre>
     * @return the covariance matrix as a double array of size [a[0].length][a[0].length]
     */
    public static double[][] covariance2(double[][] a) {
        
        double[][] cov = MatrixUtil.createATransposedTimesA(a);

        int i, j;
        int n = a[0].length;
        double invN2 = 1./((double)n - 1.);
        invN2 *= invN2;
        MatrixUtil.multiply(cov, invN2);

        /*
        System.out.printf("bf cov=%n");
        for ( i = 0; i < cov.length; ++i) {
            for ( j = 0; j < cov[i].length; ++j) {
                System.out.printf("%11.3e  ", cov[i][j]);
            }
            System.out.printf("%n");
        }
        System.out.flush();
        */

        return cov;
    }
   
}
