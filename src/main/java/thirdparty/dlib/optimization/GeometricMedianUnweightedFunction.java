package thirdparty.dlib.optimization;

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
    X = arg min of summation over all points ( || X-obs ||_2 ), that is, the X which minimizes the
    sum of the differences where _2 is notation for using L2 (euclidean) distances.
    
      <pre>
      f = summation_i=1_n( || X - obs_i || )/n     
               where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_0 = (0.5/n) * ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(-1/2) * (2*(X_0-obs_i_0))
              = (1./n) * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_1 = (1./n) * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)     
      </pre>
      
      Caveat is that the algorithm can fail within a stop strategy tolerance of
      the centroid because the derivative becomes 0.
      
      Either the algorithm has to be alert for that and adapt to it (by checking
      for centroid and testing points around it outside of the tolerance, presumably),
      OR you should use the weighted geometric-median.
      
      NOTE that the initial point for the algorithm should be within bounds of the 
      data.
      
 * @author nichole
 */
public class GeometricMedianUnweightedFunction extends AbstractGeometricMedianFunction {
 
    /**
     * the number of dimensions present in the observations.  e.g. 2 for x and y axes.
     */
    final int nDimensions;
    
    /**
     * the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3, observations={x0, y0, z0, x1, y1, z1, ...
     *     x_(nPoints-1), y_(nPoints-1), z_(nPoints-1)}.
     */
    final double[] obs;
    
    /**
     * a rough number that isn't properly normalized for use in the finiteDifference
     * method.
     */
    final double fDEps = 1e3*eps;

    /**
     * class holding the objective and derivative of the geometric-median for use
     * with search algorithms such as the LBFGs.
     * 
     * @param observations the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3, observations={x0, y0, z0, x1, y1, z1, ...
     *     x_(nPoints-1), y_(nPoints-1), z_(nPoints-1)}.
     * @param numberOfDimensions the number of data dimensions present in
     * observations array.
     */
    public GeometricMedianUnweightedFunction(double[] observations, int numberOfDimensions) {
        if (numberOfDimensions < 1) {
            throw new IllegalArgumentException("numberOfDimensions must be > 0");
        }
        this.nDimensions = numberOfDimensions;
        this.obs = Arrays.copyOf(observations, observations.length);

        int n = obs.length;

        if (n < 1) {
            throw new IllegalArgumentException("observations must have length > 0");
        }

        if ((n % nDimensions) != 0) {
            throw new IllegalArgumentException("observations.length must be a multiple of "
                    + "numberOfDimensions");
        }
    }

    /**
     * given observed data points obs = (x_i, y_i, ...) want to solve for
     * X=(x_geo_median, y_geo_median, ...), that is X = arg min of || X-obs
     * ||_2, that is, the X which minimizes the sum of the differences, where _2
     * is notation for using L2 distances.
     *
     * the objective is calculated by setting the derivative of the cost function to 0.
       <pre>
       The cost function is:
          C(X) = summation_i=1_n( || X - obs_i || ) / n
      
          where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)  
      </pre>
     *
     * @param geoMedian input variable holding coordinates of current
     * estimate of geometric median.
     * @return evaluation of the objective, summation_i=1_n(|| geoMedian - obs_i
     * ||^2)/n
     */
    @Override
    public double f(double[] geoMedian) {

        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length should equal nDimensions");
        }

        return evaluateGeoMedian(geoMedian);
    }

    /**
     *
     * d/dX of the objective C(X) = summation_i=1_n(|| X-obs_i ||^2)/n where the
     * observed data points are obs = (x_i, y_i, ...) and X=(x_geo_median,
     * y_geo_median, ...)
     
       <pre>
       f = summation_i=1_n( || X - obs_i || )/n
           where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_0 = (0.5/n) * ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(-1/2) * (2*(X_0-obs_i_0))
               = (1./n) * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_1 = (1./n) * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       ...
       
       </pre>
     
      @param geoMedian coordinates of current estimate of geometric median
      @return evaluation of the derivative
     */
    @Override
    public double[] der(double[] geoMedian) {

        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length should equal nDimensions");
        }
        
        if (false) {
            return finiteDifference(geoMedian);
        }

        //double[] geoMedian0 = Arrays.copyOf(geoMedian, geoMedian.length);

        double[] diffs = calculateDifferences(geoMedian);
        double[] ssdPerPoint = calculateSSDPerPoint(diffs);
        
        int nData = (int) (obs.length / nDimensions);
        
        assert(ssdPerPoint.length == nData);
        
        double s = 0;
        int i;
        for (i = 0; i < ssdPerPoint.length; ++i) {
            s += Math.sqrt(ssdPerPoint[i]);
        }
        
        // to avoid divide by 0:
        s += eps;
        
        //df/dX_0 = (1./n) * (X_0 - obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2)
        double[] dfDX = new double[nDimensions];
        
        // NOTE: can see that the unweighted algorithm will not make progress when
        //    the current geometric-median estimate is the centroid.
        int j, d;
        for (i = 0; i < nData; ++i) {
            
            // skip a point if it is the same as the median
            boolean skip = true;
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                if (Math.abs(geoMedian[d] - obs[j]) > eps) {
                    skip = false; 
                    break;
                }
            }
            if (skip) {
                continue;
            }
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                dfDX[d] += (geoMedian[d] - obs[j]);
            }
        }
        
        for (d = 0; d < nDimensions; ++d) {
        //    dfDX[d] /= (nData * s);
        }
        
        return dfDX;
    }
    
    /*
    df/dX_0 = (1./n) * (X_0 - obs_i_0) * ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(-1/2)
    
    d/dx of df/dX_0 = (1./n)*( (-1/2)* ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(-3/2) ) ) * 2 * (X_0-obs_i_0)
        = (-1/n)*(X_0-obs_i_0)/ ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(3/2) ) )
    
    */
    
    public int getNDimensions() {
        return nDimensions;
    }
    
    public double[] getObs() {
        return obs;
    }
    
    public double getFDEps() {
        return fDEps;
    }
}
