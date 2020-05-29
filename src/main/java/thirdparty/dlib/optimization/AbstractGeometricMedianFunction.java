package thirdparty.dlib.optimization;

import algorithms.util.IFunction;
import java.util.Arrays;

/**
 * objective function implementations for use with search algorithms using IFunction
 * as input for the geometric-median.
 * 
   calculating the geometric median:
   definition: the point which minimizes the sum of the euclidean distance of that 
     point to all other points in the set.
     the sum is a convex function (i.e. local search will work).
     Unfortunately, no algorithms are closed form, that is no algorithms have a
     finite number of computational operations.
     The geometric median is a rotation and translation invariant estimator that 
     achieves the optimal breakdown point of 0.5, i.e. it is a good estimator
     even when up to half of the input data is arbitrarily corrupted.
     (https://dl.acm.org/doi/pdf/10.1145/2897518.2897647)
    
    It's the Fermat-Weber problem.
    given observed data points obs = (x_i, y_i, ...)
    want to solve for X=(x_geo_median, y_geo_median, ...).
    X = arg min of || X-obs ||_2, that is, the X which minimizes the
    sum of the differences where _2 is notation for using L2 (euclidean) distances.
    
      <pre>
      f = summation_i=1_n( || X - obs_i || )/n     
               where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_0 = (0.5/n) * ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(-1/2) * (-2*(X_0-obs_i_0))
              = (-1./n) * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_1 = (-1./n) * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)          
      </pre>
    
     The weighted geometric median is known as the Weber problem.
     
     IRLS, can use the Weiszfeld algorithm for geometric median:
         Two steps:
           (1) w_i = 1/[max(|| X-z_i ||_2, delta)]
           (2) X = (summation_i=1_n(w_i*z_i))/(summation_i=1_n(w_i))
                       where (2) is derived from setting the deriv to zero for:
                           X = minimum X in (summation_i=1_n(w_i*(z_i-X_i)^2))
         And a strategy for optimization:
            fix X{y} and minimize f with respect to X{x},
            fix X{x} and minimize f with respect to X{y},
            and repeat:
              for i ← 1, 2, . . .
                 X{y}_i+1 = min_(X{y}) f(X{x}_i+1, X{y}) . Optimize X{y} with X{x} fixed
                 X{x}_i+1 = min_(X{x}) f(X, X{y}_i) . Optimize X{x} with X{y} fixed
            
   TODO: consider including observational errors.  Those would affect the starting point 
   estimate and the weights.
   
 * @author nichole
 */
public abstract class AbstractGeometricMedianFunction implements IFunction {
        
    public static double eps = 1e-17;
    
    protected abstract int getNDimensions();
    
    protected abstract double[] getObs();
    
    /**
     * calculate centroid of data for each dimension. the centroid can be used
     * as a possible starting point for the search.
     *
     * @return an array of the centroids of the observations for each dimension
     */
     double[] calculateCentroid() {

        final int nDimensions = getNDimensions();
        final double[] obs = getObs();

        double[] sum = new double[nDimensions];

        int nData = (int) (obs.length / nDimensions);

        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                sum[d] += obs[j];
            }
        }
        for (d = 0; d < nDimensions; ++d) {
            sum[d] /= (double) nData;
        }

        return sum;
    }
     
    double[] calculateDifferences(double[] geoMedian) {

        final int nDimensions = getNDimensions();

        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length must == nDimensions");
        }

        final double[] obs = getObs();
        double[] diff = new double[obs.length];
        int nData = (int) (obs.length / nDimensions);

        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                diff[j] = geoMedian[d] - obs[j];
            }
        }
        return diff;
    }
    
    /**
     * calculate the sum of square differences for each point.
     * Note: the result is length nData, not obs.length or 1.
     * @param diffs the differences between each point and the current
     * geometric-median coordinates.
     * @return for each point, the sum of differences.
     */
    double[] calculateSSDPerPoint(double[] diffs) {
        
        final int nDimensions = getNDimensions();
        final double[] obs = getObs();
        int nData = (int) (obs.length / nDimensions);
        
        if (diffs.length != obs.length) {
            throw new IllegalArgumentException("diffs length must == obs.length");
        }
        
        double[] ssd = new double[nData];
        
        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                ssd[i] += (diffs[j] * diffs[j]);
            }
        }
        return ssd;
    }
    
    public static String toString(double[] a) {
        StringBuilder sb = new StringBuilder("[");
        for (int i=0;i<a.length;i++) {
            sb.append(String.format("%.3e", a[i]));
            if (i < a.length - 1) {
                sb.append(", ");
            }
        }
        sb.append("]");
        return sb.toString();
    }
}
