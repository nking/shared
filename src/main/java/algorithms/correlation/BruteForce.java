package algorithms.correlation;

/**
 *
 * @author nichole
 */
public class BruteForce {
    
    public static double[][] correlation(double[][] a) {
        
        /*
       can implement as pairs in the matrix:
           aCov= covariance(a)
           covXY = covariance(x,y) result in aCov.
           covXX  = covariance(x,x) result in aCov.
           covXX  = covariance(y,y) result in aCov.
       if abs(covXX * covYY)<1e-10
           cor = 0;
       else
           cor = sign(covXY) * sqrt( abs(covXY) / sqrt(covXX * covYY) );
       end
        
        also, exploration of correlation:
        https://numbersandcode.com/some-interesting-observations-with-distance-correlation-coefficients
        
         */
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * calculate the covariance matrix for a
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
        
        return cov;
    }
    
}
