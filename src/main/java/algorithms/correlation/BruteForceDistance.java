package algorithms.correlation;

import algorithms.matrix.MatrixUtil;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * empirical distance correlation matrix.
          https://arxiv.org/pdf/0803.4101.pdf
          MEASURING AND TESTING DEPENDENCE BY CORRELATION OF DISTANCES
          Szekely, Rizzo, and Bakirov
          2007, Vol. 35, No. 6, 2769–2794
            dist. correlation = sqrt(  (1/n^2)*sAB /((1/n^2)sqrt(sAA*sBB) )
                              = sqrt(  sAB /(sqrt(sAA*sBB) )
        
 * useful for Mahalanobis distance among many things.
 * 
 * @author nichole
 */
public class BruteForceDistance {
    
    private static final Level LEVEL = Level.FINE;
    private static final Logger log;
    static {
        log = Logger.getLogger(BruteForceDistance.class.getSimpleName());
        log.setLevel(LEVEL);
    }
        
    /**
     * euclidean distance matrix
     * (ported from https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
     * 
     @param a
     @return 
     */
    public static double[][] distanceMatrix(double[][] a) {
        int n = a.length;
        int m = a[0].length;
        
        double[][] d = new double[n][n];
        
        int i, j, col;
        double diff;
        for (i = 0; i < n; ++i) {
            d[i] = new double[n];
            for (j = (i+1); j < n; ++j) {
                // for each column, subtract cell in row i from row j, square it, and add to total
                for (col = 0; col < m; ++col) {
                    diff = a[i][col] - a[j][col];
                    d[i][j] += (diff * diff);
                }
                d[i][j] = Math.sqrt(d[i][j]);
                d[j][i] = d[i][j];
            }
        }
        
        return d;
    }
 
    private static double[] _meanPerColumn(double[][] a) {
        
        double[] m = new double[a[0].length];
        
        for (int j = 0; j < a[0].length; ++j) {
            for (int i = 0; i < a.length; ++i) {
                m[j] += a[i][j];
            }
            m[j] /= (double)(a.length);
        }
        
        return m;
    }
    
    private static double[] _meanPerRow(double[][] a) {
        
        double[] m = new double[a.length];
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                m[i] += a[i][j];
            }
            m[i] /= (double)(a[i].length);
        }
        
        return m;
    }
    
    private static double _grandMean(double[][] a) {
        
        double m = 0;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                m += a[i][j];
            }
        }
        m /= (double)(a.length * a[0].length);
        
        return m;
    }
    
    private static void _normalizeDistanceMatrix(double[][] a) {
        
        double[] colMeans = _meanPerColumn(a);
        double[] rowMeans = _meanPerRow(a);
        double grandMean = _grandMean(a);
        
        if (log.isLoggable(LEVEL)) {
            StringBuilder sb = new StringBuilder();
            sb.append("col means=").append(Arrays.toString(colMeans)).append("\n");
            sb.append("row means=").append(Arrays.toString(rowMeans)).append("\n");
            sb.append("grand mean=").append(Double.toString(grandMean)).append("\n");
            log.log(Level.INFO, sb.toString());
        }
        
        // double center the matrix by 
        //    subtracting colMeans
        //    subtracting rowMeans
        //    add grandMeans
        double m;
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                m = - colMeans[j] - rowMeans[i] + grandMean;
                a[i][j] += m;
            }
        }
    }
    
    /**
     * calculate the empirical distance correlation matrix.
          https://arxiv.org/pdf/0803.4101.pdf
          MEASURING AND TESTING DEPENDENCE BY CORRELATION OF DISTANCES
          Szekely, Rizzo, and Bakirov
          2007, Vol. 35, No. 6, 2769–2794
            dist. correlation = sqrt(  (1/n^2)*sAB /((1/n^2)sqrt(sAA*sBB) )
                              = sqrt(  sAB /(sqrt(sAA*sBB) )
     * also see:
     * https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
     * https://numbersandcode.com/some-interesting-observations-with-distance-correlation-coefficients
     @param X data in double array format with number of rows n being same as
     * nRows for Y, though the number of dimensions, that is nCols can be different.
     @param Y data in double array format with number of rows n being same as
     * nRows for X, though the number of dimensions, that is nCols can be different.
     @return returns distance correlation and intermediate data.
     * (i) If E(|X|p + |Y|q) .lt. infinity, then 0 .lte. R .lte. 1, 
     *     and R(X, Y ) = 0 if and only if X and Y are independent.
       (ii) 0 .lte. R_n .lte. 1.
       (iii) If R_n(X, Y) = 1, then there exist a vector a, a nonzero real number
             b and an orthogonal matrix C such that Y = a + bXC.
     */
    public static DCOV correlation1(double[][] X, double[][] Y) {
        
        int n = X.length;
        if (Y.length != n) {
            throw new IllegalArgumentException("X and Y must have same number of rows");
        }
        
        X = MatrixUtil.copy(X);
        Y = MatrixUtil.copy(Y);
        
        // these may have negative values because of the double-centering normalization:
        double[][] dMX = distanceMatrix(X);
        double[][] dMY = distanceMatrix(Y);

        _normalizeDistanceMatrix(dMX);
        _normalizeDistanceMatrix(dMY);
       
        // === empirical distance correlation: ==
        // https://arxiv.org/pdf/0803.4101.pdf
        // MEASURING AND TESTING DEPENDENCE BY CORRELATION OF DISTANCES
        // Szekely, Rizzo, and Bakirov
        // 2007, Vol. 35, No. 6, 2769–2794
        //   dist. correlation = sqrt(  (1/n^2)*sXY /((1/n^2)sqrt(sXX*sYY) )
        //                     = sqrt(  sXY /(sqrt(sXX*sYY) )
        
        double sXY = 0;
        double sXX = 0;
        double sYY = 0;
                
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < dMX[i].length; ++j) {
                sXY += (dMX[i][j] * dMY[i][j]);
                sXX += (dMX[i][j] * dMX[i][j]);
                sYY += (dMY[i][j] * dMY[i][j]);
            }
        }
        
        double invN = 1./((double)n-1.);
        
        double dCovSq = sXY;
        double dVarXSq = sXX;
        double dVarYSq = sYY;
                
        //emp. dist. correlation =  sXY /(sqrt(sXX*sYY) )
        double corSq = (dVarXSq > 0 && dVarYSq > 0) ? dCovSq/Math.sqrt(dVarXSq * dVarYSq) : 0.;
        
        /*
        to test indepence:
           The Chi-Squared Test of Distance Correlation
           Shen & Vogelstein 2019 
           https://www.groundai.com/project/the-chi-square-test-of-distance-correlation/2
           testing independence using distance correlation now runs in linear time complexity
        
        */
        // T = nrow(DX)*V2XY;   /* test statistic p. 2783. Reject indep when T>=z */
        // z-score =( M - population mean)/(population stdev) using sample mean and stdev in place of population
        //
        //    
        
        if (log.isLoggable(LEVEL)) {
            StringBuilder sb = new StringBuilder();
            sb.append("n=").append(n).append("\n");
            sb.append("dVarXSq=").append(dVarXSq).append("\n");
            sb.append("dVarYSq=").append(dVarYSq).append("\n");
            sb.append("dCovXYSq=").append(dCovSq).append("\n");
            sb.append("dCorSq=").append(corSq).append("\n");
            log.log(Level.INFO, sb.toString());
        }
        
        DCOV dc = new DCOV();
        dc.corSq = corSq;
        dc.dCovSq = dCovSq;
        dc.dVarXSq = dVarXSq;
        dc.dVarYSq = dVarYSq;
        
        return dc;
    }
    
    /**
     *
     */
    public static class DCOV {
        double dCovSq;
        double dVarXSq;
        double dVarYSq;                
        double corSq;
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("corSq=").append(corSq).append("\n");
            sb.append("dCovSq: ").append(dCovSq).append("\n");
            sb.append("dVarXSq: ").append(dVarXSq).append("\n");
            sb.append("dVarYSq: ").append(dVarYSq).append("\n");
            return sb.toString();
        }
    }
}
