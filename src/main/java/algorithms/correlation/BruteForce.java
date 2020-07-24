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
    
    public static double[][] correlation(double[][] a) {
        
        /*
       can implement as pairs in the matrix:
           aCov= covariance(a)
           covXY = covariance(x,y) result is in aCov.
           covXX  = covariance(x,x) result is in aCov.
           covXX  = covariance(y,y) result is in aCov.
       
           
       
         */
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * calculate the covariance matrix for a using a brute force method
     * @param a
     * @return 
     */
    public static double[][] covariance(double[][] a) {
        
        int nRows = a.length;
        int nCols = a[0].length;
        
        int i, j;
        
        double[] mean = new double[nCols];
        double sum;
        for (i = 0; i < nCols; ++i) {
            sum = 0;
            for (j = 0; j < nRows; ++j) {
                sum += (a[j][i]);
            }
            mean[i] = sum/(double)nRows;
        }
        
        double[][] diffs = new double[nRows][];
        for (i = 0; i < nRows; ++i) {
            diffs[i] = new double[nCols];
            for (j = 0; j < nCols; ++j) {
                diffs[i][j] = (a[i][j] - mean[j]);
            }
        }
        
        System.out.printf("bf stand. means=%s\n", Arrays.toString(mean));
        
        System.out.flush();
        System.out.printf("bf stand. diffs=\n");
        for ( i = 0; i < diffs.length; ++i) {
            for ( j = 0; j < diffs[i].length; ++j) {
                System.out.printf("%11.3e  ", diffs[i][j]);
            }
            System.out.printf("\n");
        }
        System.out.flush();
        
        double[][] cov = new double[nCols][];
        for (i = 0; i < nCols; ++i) {
            cov[i] = new double[nCols];
        }
     
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
        
        System.out.printf("bf stand. cov=\n");
        for ( i = 0; i < cov.length; ++i) {
            for ( j = 0; j < cov[i].length; ++j) {
                System.out.printf("%11.3e  ", cov[i][j]);
            }
            System.out.printf("\n");
        }
        System.out.flush();
        
        return cov;
    }
   
}
