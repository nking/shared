package algorithms.statistics;

import algorithms.QuickSort;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class LinearRegression {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
    TODO: consider implementing Siegel Repeated Median estimator in this class too.

     * calculate the Theil-Sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
       https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
       In non-parametric statistics, the Theil–Sen estimator is a method for 
       robustly fitting a line to sample points in the plane (simple linear 
       regression) by choosing the median of the slopes of all lines through pairs of points. 
       ...  This estimator can be computed efficiently, and is insensitive to outliers.
     It can be significantly more accurate than non-robust simple linear regression
     (least squares) for skewed and heteroskedastic data, and competes well against
     least squares even for normally distributed data in terms of statistical power.[10]
     It has been called "the most popular nonparametric technique for estimating a linear trend".
     @param x
     @param y
     @return array holding y-intercept and slope
     */
    public float[] calculateTheilSenEstimatorParams(int[] x, int[] y) {
        
        int n = x.length;
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (float)(y[j] - y[i])/((float)x[j] - x[i]);
                count++;
            }
        }
        
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        float median;
        if ((idx & 1) == 0 && idx > 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.fine("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
        
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        float yIntercept = y[medianIdx] - median * x[medianIdx];
        
        //the estimation of yIntercept needs to be improved:
        // TODO: correct this to calculate yIntercept from  median of yi − mxi
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sum = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sum += (y[j] - median * x[j]);
            }
            yIntercept = sum/((float)(2*np + 1));
        }
        
        return new float[]{yIntercept, median};
    }

    /**
     * estimate the y-intercept and slope for the given x,y dataset
     * The estimates are closed form.  r.t.c. is O(n).
     * For a method more robust to outliers, see the Theil-Sen methods in this class.
     <pre>
     https://en.m.wikipedia.org/wiki/Simple_linear_regression
     </pre>
     The r.t.c. is O(n).
     * @param x array of x data points.
     * @param y array of y data points
     * @return array of {
     * y-intercept,
     * slope of the line fit,
     * variance of the error (where error is y - yIntercept - slope*x),
     * variance of the y-intercept, and variance of the slope
     */
    public static double[] simpleLinearRegression(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;
        //x    y    x^2    y^2    xy   = moments
        //0    1    2      3       4   = indexes in output
        double[] moments = Util.caldc2DMomentsX2Y2(x, y);

        // beta is slope
        double slopeEst = (n * moments[4] - moments[0] * moments[1]) /
                (n * moments[2] - moments[0] * moments[0]);

        // alpha is y intercept
        double yInterEst = (1./n) * moments[1] - (slopeEst * moments[0]/ n);

        double varError = (1/(n*(n-2.))) *
                (n * moments[3] - moments[1]*moments[1] -
                ((slopeEst*slopeEst)*(n * moments[2] - moments[0] * moments[0])));

        double varSlopeEst = (n * varError) / (n * moments[2] - moments[0]*moments[0]);

        double varYInterEst = varSlopeEst * moments[2]/n;

        /*
        {y-intercept, slope of the line fit,
        variance of the error (where error is y - yIntercept - slope*x),
        variance of the y-intercept, and variance of the slope
         */
        return new double[]{yInterEst, slopeEst, varError, varYInterEst, varSlopeEst};
    }

    /**
     * calculate simple linear regression using a design matrix approach.
     <pre>
     https://en.m.wikipedia.org/wiki/Design_matrix#Simple_linear_regression
     </pre>
     The r.t.c. is O(n^3), dominated by the pseudoinverse.
     * @param x matrix of size n x p where n is the number of samples and p is the number
     *          of variables (a.k.a. features).
     * @return an array of the intercept and slopes.  the array length withh be p + 1.
     */
    public static double[] linearRegression(double[] y, double[][] x) throws NotConvergedException {
        int n = y.length;
        if (x.length != n) {
            throw new IllegalArgumentException("x and y must be same lengths");
        }

        /*
        solves minimization problem:
           argmin of sum_{i=1 to n}( (y_i - alpha - beta*x_i))^2 )

              can find min estimates in the optimization eqn by setting derivs to 0

           alpha_est = y_mean - (beta_est * x_mean)

           beta_est = sum_{i=1 to n}( (x_i - x_mean) * (y_i - y_mean) )
                       / sum_{i=1 to n}( (x_i - x_mean)^2)

           those are efficient if x_mean and y_mean are known.

           else a design matrix approach, similar to the design matrix for polynomial regressions.
         */
        /*
        y = X * b + eps
        (y - eps) = X*b
        X^T*(y - eps) = X^T*X*b
        (X^T*X)*^-1 * X^T*(y - eps) = b
         */
        int p = x[0].length;
        double[][] X = new double[n][p + 1];
        for (int i = 0; i < n; ++i) {
            X[i][0] = 1;
            System.arraycopy(x[i], 0, X[i], 1, p);
        }

        double[][] a = MatrixUtil.pseudoinverseFullColumnRank(X);
        double[] b = MatrixUtil.multiplyMatrixByColumnVector(a, y);

        return b;
    }

    /**
     * estimate the y-intercept and slope for the given x,y dataset where the means
     * of the datasets are already known.
     * The estimates are closed form.  r.t.c. is O(n).
     * For a method more robust to outliers, see the Theil-Sen methods in this class.
     <pre>
     https://en.m.wikipedia.org/wiki/Simple_linear_regression
     </pre>
     The r.t.c. is O(n).
     * @param x array of x data points.
     * @param y array of y data points
     * @param meanX known mean of the x data
     * @param meanY known mean of the y data
     * @return array of y-intercept and slope of the line fit
     */
    public static double[] lineFit(double[] x, double[] y, double meanX, double meanY) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;;
        double dx, dy;
        double a = 0;
        double b = 0;
        for (int i = 0; i < n; ++i) {
            dx = x[i] - meanX;
            dy = y[i] - meanY;
            a += (dx * dy);
            b += (dx * dx);
        }
        // beta estimate:
        double slopeEst = a/b;
        // alpha estimate:
        double yInterceptEst = meanY - (slopeEst * meanX);
        return new double[]{yInterceptEst, slopeEst};
    }

    /**
     * make a plot of the linear regression of arrays x and y.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     @param x
     @param y 
     @return file name
     */
    public String plotTheLinearRegression(int[] x, int[] y) {
                        
        float[] tsbParams = calculateTheilSenEstimatorParams(x, y);
        
        float yIntercept = tsbParams[0];
        
        float slope = tsbParams[1];
        
        /*
        plot dx, dy
        and plot a line generated from the yIntercept and median: yIntercept − median*x_i
        */        
        int xMin = MiscMath0.findMin(x);
        int xMax = MiscMath0.findMax(x);
        int len = xMax - xMin + 1;
        int[] tsbX = new int[len];
        int[] tsbY = new int[len];
        int count = 0;
        for (int xCoord = xMin; xCoord <= xMax; xCoord++) {
            float yCoord = yIntercept + slope * (float)xCoord;
            tsbX[count] = xCoord;
            tsbY[count] = Math.round(yCoord);
            count++;
        }
        
        int yMin = MiscMath0.findMin(y);
        int yMax = MiscMath0.findMax(y);
       
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                x, y, 
                tsbX, tsbY,
                "X vs Y and thiel sen beta linear regression line");

            return plotter.writeFile();
            
        } catch(IOException e) {
            
            log.severe("ERROR while trying to write plot: " + e.getMessage());
        }
        return "";
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     @param x
     @param y
     @return 
     */
    public float[] calculateTheilSenEstimatorParams(float[] x, float[] y) {
        
        int n = x.length;
        
        if (n > 46340) {
            throw new IllegalArgumentException("x and y lengths must be "
                + "less than 46340 for indexing an array of size length*lnegth");
        }
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (y[j] - y[i])/(x[j] - x[i]);
                count++;
            }
        }
        
        if (count == 0) {
            // this can happen for vertical lines
            return new float[]{Float.NaN, Float.MAX_VALUE};
        }
        
        float median;
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        if ((idx & 1) == 0 && idx > 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.fine("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values 
        //     y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
        
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        float yIntercept = y[medianIdx] - median * x[medianIdx];
        
        //the estimation of yIntercept needs to be improved:
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sum = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sum += (y[j] - median * x[j]);
            }
            yIntercept = sum/((float)(2*np + 1));
        }
        
        return new float[]{yIntercept, median};
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     @param x
     @param y
     @return 
     */
    public double[] calculateTheilSenEstimatorParams(double[] x, double[] y) {
        
        int n = x.length;
        
        if (n > 46340) {
            throw new IllegalArgumentException("x and y lengths must be "
                + "less than 46340 for indexing an array of size length*lnegth");
        }
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        double[] s = new double[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (y[j] - y[i])/(x[j] - x[i]);
                count++;
            }
        }
        
        if (count == 0) {
            // this can happen for vertical lines
            return new double[]{Double.NaN, Double.MAX_VALUE};
        }
        
        double median;
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        if ((idx & 1) == 0 && idx > 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.fine("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values 
        //     y[i] − median * x[i]
        double[] s2 = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        int[] idxs = MiscSorter.mergeSortIncreasing(s2);
        
        //QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
        
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        double yIntercept = y[idxs[medianIdx]] - median * x[idxs[medianIdx]];
        
        //the estimation of yIntercept needs to be improved:
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            double sum = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sum += (y[idxs[j]] - median * x[idxs[j]]);
            }
            yIntercept = sum/(2.*np + 1.);
        }
        
        return new double[]{yIntercept, median};
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     @param x
     @param y
     @return 
     */
    public float[] calculateTheilSenEstimatorMedian(float[] x, float[] y) {
        
        int n = x.length;
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (y[j] - y[i])/(x[j] - x[i]);
                count++;
            }
        }
        
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        float median;
        if ((idx & 1) == 0 && idx > 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.fine("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
       
        float xMedian = x[medianIdx];
        float yMedian = y[medianIdx];
        // improve the vlue over several points
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sumX = 0;
            float sumY = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sumX += x[j];
                sumY += y[j];
            }
            xMedian = sumX/((float)(2*np + 1));
            yMedian = sumY/((float)(2*np + 1));
        }
        
        return new float[]{xMedian, yMedian};
    }
    
    /**
     * make a plot of the linear regression of arrays x and y.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     @param x
     @param y 
     @param xMin 
     @param xMax 
     @param yMin 
     @param yMax 
     @return plot file name
     */
    public String plotTheLinearRegression(float[] x, float[] y, int xMin, int xMax, int yMin, int yMax) {
                        
        float[] tsbParams = calculateTheilSenEstimatorParams(x, y);
        
        float yIntercept = tsbParams[0];
        
        float slope = tsbParams[1];
        
        /*
        plot dx, dy
        and plot a line generated from the yIntercept and median: yIntercept − median*x_i
        */    
        int len = xMax - xMin + 1;
        float[] tsbX = new float[len];
        float[] tsbY = new float[len];
        int count = 0;
        for (int xCoord = xMin; xCoord <= xMax; xCoord++) {
            float yCoord = yIntercept + slope * (float)xCoord;
            tsbX[count] = xCoord;
            tsbY[count] = yCoord;
            count++;
        }
        
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                x, y, 
                tsbX, tsbY,
                "X vs Y and thiel sen beta linear regression line");

            return plotter.writeFile();
            
        } catch(IOException e) {
            
            log.severe("ERROR while trying to write plot: " + e.getMessage());
        }
        
        return "";
    }
    
}
