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
    
    public abstract int getNDimensions();
    
    public abstract double[] getObs();
    
    public abstract double getFDEps();
    
    public int getNData() {
        return getObs().length/getNDimensions();
    }
    
    /**
     * calculate centroid of data for each dimension. the centroid can be used
     * as a possible starting point for the search.
     *
     * @return an array of the centroids of the observations for each dimension
     */
     public double[] calculateCentroid() {

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

    /**
     * calculate the difference of each observation by dimension with the
     * geometric median by dimension.  The resulting array is ordered in the
     * same manner as the instance observations, that is
     * the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3, observations={x0, y0, z0, x1, y1, z1, ...
     *     x_(nPoints-1), y_(nPoints-1), z_(nPoints-1)}.
     * @param geoMedian
     * @return array of geometric median - observations.  the array is ordered
     * in the same manner as obs and is the same length.
     */     
    public double[] calculateDifferences(final double[] geoMedian) {

        final int nDimensions = getNDimensions();

        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length must == nDimensions");
        }

        final double[] obs = getObs();
        double[] diff = new double[obs.length];
        int nData = (int) (obs.length / nDimensions);

        int i, j, d;
        double a, b;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                a = geoMedian[d];
                b = obs[j];
                diff[j] = a - b;
            }
        }
        return diff;
    }
    
    /**
     * calculate the difference of each observation by dimension with the
     * geometric median by dimension.  The resulting array is ordered in the
     * same manner as the instance observations, that is
     * the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3, observations={x0, y0, z0, x1, y1, z1, ...
     *     x_(nPoints-1), y_(nPoints-1), z_(nPoints-1)}.
     * @param geoMedian
     * @param output output array of length obs.length to hold results
     * @return array of geometric median - observations.  the array is ordered
     * in the same manner as obs and is the same length.
     */
    void calculateDifferences(final double[] geoMedian, double[] output) {
        double[] diffs = calculateDifferences(geoMedian);
        System.arraycopy(diffs, 0, output, 0, diffs.length);
    }
    
    /**
     * given observed the geometric median X=(x_geo_median, y_geo_median, ...) 
     * and instance data points obs = (x_i, y_i, ...) evaluate 
     * || X-obs ||_2, that is, the X which minimizes the sum of the differences, 
     * where _2 is notation for using L2 distances.
     *
     * @param geoMedian input variable holding coordinates of current
     * estimate of geometric median.
     * @return evaluation of the objective, summation_i=1_n(|| geoMedian - obs_i
     * ||^2)/n
     */
    double evaluateGeoMedian(double[] geoMedian) {
        
        int nDimensions = geoMedian.length;

        //double[] geoMedian0 = Arrays.copyOf(geoMedian, geoMedian.length);

        double[] diffs = calculateDifferences(geoMedian);
        
        int nData = (int) (diffs.length / nDimensions);
        
        int i, j, d;
        double dist;
        double sum = 0;
        for (i = 0; i < nData; ++i) {
            dist = 0;
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                dist += (diffs[j] * diffs[j]);
            }
            dist = Math.sqrt(dist);
            sum += dist;
        }
        //sum /= (double)nData;
        
        return sum;
    }
    
    /**
     * given observed the geometric median X=(x_geo_median, y_geo_median, ...) 
     * and instance data points obs = (x_i, y_i, ...) evaluate 
     * || X-obs ||_2, that is, the X which minimizes the sum of the differences, 
     * where _2 is notation for using L2 distances.
     *
     * @param geoMedian input variable holding coordinates of current
     * estimate of geometric median.
     * @return evaluation of the objective, summation_i=1_n(|| geoMedian - obs_i
     * ||^2)/n
     */
    public double[] evaluateGeoMedianPerDimension(double[] geoMedian) {
        
        int nDimensions = geoMedian.length;
        
        double[] out = new double[nDimensions];

        //double[] geoMedian0 = Arrays.copyOf(geoMedian, geoMedian.length);

        double[] diffs = calculateDifferences(geoMedian);
        
        int nData = (int) (diffs.length / nDimensions);
        
        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                out[d] += (diffs[j] * diffs[j]);
            }
        }
        for (d = 0; d < nDimensions; ++d) {
            out[d] /= (double)nData;
        }
        
        return out;
    }
    
    /**
     * calculate the sum of square differences for each point.
     * Note: the result is length nData, not obs.length or 1.
     * @param diffs the differences between each point and the current
     * geometric-median coordinates.
     * @return for each point, the sum of differences.
     */
    public double[] calculateSSDPerPoint(double[] diffs) {
        
        final int nDimensions = getNDimensions();
        final double[] obs = getObs();
        int nData = (obs.length / nDimensions);
        
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
    
    /**
    adapted from dlib optimization.h
    Copyright (C) 2008  Davis E. King (davis@dlib.net)
    License: Boost Software License   See LICENSE.txt for the full license.
    */
    double[] finiteDifference(double[] coeffs) {

        //System.out.println("a1  x.size=" + coeffs.length);

        int n = coeffs.length;
        
        //TDO: set this be a fraction of the diagonal of the bounds of obs
        final double fds = 0.001;//getFDEps();

        double[] der = new double[n];
        double[] e = Arrays.copyOf(coeffs, n);

        for (int i = 0; i < n; ++i) {
            final double old_val = e[i];
            e[i] += fds;
            final double delta_plus = f(e);
            e[i] = old_val - fds;
            final double delta_minus = f(e);

            // finite difference:  this is the approx jacobian
            der[i] = (delta_plus - delta_minus)/(2.*fds); 

            //NOTE: newton's method would continue with:
            // x_(i+1) = x_i - (delta_plus/der(i))

            // and finally restore the old value of this element
            e[i] = old_val;
        }

        return der;
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
