package thirdparty.dlib.optimization;

import java.util.Arrays;

/**
 * objective function implementations for use with search algorithms using
 * IFunction as input for the geometric-median.
 *
 * calculating the geometric median: definition: the point which minimizes the
 * sum of the euclidean distance of that point to all other points in the set.
 * the sum is a convex function (i.e. local search will work). Unfortunately, no
 * algorithms are closed form, that is no algorithms have a finite number of
 * computational operations. The geometric median is a rotation and translation
 * invariant estimator that achieves the optimal breakdown point of 0.5, i.e. it
 * is a good estimator even when up to half of the input data is arbitrarily
 * corrupted. 
 * 
 * (https://dl.acm.org/doi/pdf/10.1145/2897518.2897647)
 * and "Numerical Algorithms" by Solomon https://people.csail.mit.edu/jsolomon/share/book/numerical_book.pdf
 *
 * given observed data points obs = (x_i, y_i,
 * ...) want to solve for X=(x_geo_median, y_geo_median, ...). X = arg min of ||
 * X-obs ||_2, that is, the X which minimizes the sum of the differences where
 * _2 is notation for using L2 (euclidean) distances.
 *
 * The weighted geometric median is known as the Weber problem.
 *
 * <pre>
 *   f = summation_i=1_n( w_i * || X - obs_i || )/n
 *      where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
 * </pre>
 *
 *  -- iteratively re-weighted least squares (IRLS)
        given observed points Z_I=(x_i, y_i)
        and wanting to solve for X=(x_median, y_median)
        let weight w_i = 1/( || X-z_i || )
        note that the geometric median usually uses equal weights w_i=1/n

        the cost function to minimize C_w(X) = summation_i=1_n( w_i*|| X-z_i ||^2 )
        set deriv to zero
        then iterated updates of estimates of X in time steps t:
           X_(t+1) = summation_i=1_n( w_t*z_i ) / summation_i=1_n( w_t^2 )

        for L1 version given by Ostresh:
           X_(t+1) = X_t - alpha*summation_i=1_n( w_t*z_i / || X-z_i || ) /
                       summation_i=1_n( w_t / || X-z_i || )

        the gradient of the Riemannian sum-of-distances function is given by
           grad( f(X) ) = -summation_i=1_n( w_i*Log(X_i) / (X-z_i))

           using a steepest descent iteration with step size alpha:
               X_(t+1) = Exp_(X_t)( alpha * v_t)
           where
               v_t = summation_i=1_n( w_i*Log(X_i) / (X-z_i)) /
                       summation_i=1_n( w_i / (X-z_i))
           authors found alpha=1 leads to convergence

      *IRLS, can use the Weiszfeld algorithm for geometric median:
         Two steps:
           (1) w_i = 1/[max(|| X-obs_i ||_2, delta)]
           (2) X = (summation_i=1_n(w_i*obs_i))/(summation_i=1_n(w_i))
                       where (2) is derived from setting the deriv to zero for:
                                 X = minimum X in (summation_i=1_n(w_i*(obs_i-X_i)^2))
         And a strategy for optimization:
            fix X{y} and minimize f with respect to X{x},
            fix X{x} and minimize f with respect to X{y},
            and repeat:
              for i ← 1, 2, . . .
                 X{y}_i+1 = min_(X{y}) f(X{x}_i+1, X{y}) . Optimize X{y} with X{x} fixed
                 X{x}_i+1 = min_(X{x}) f(X, X{y}_i) . Optimize X{x} with X{y} fixed

*  NOTE that the initial point for the algorithm should be within bounds of the 
      data.
      
 * TODO: consider including observational errors. Those would affect the
 * starting point estimate and the weights.
 *
 * @author nichole
 */
public class GeometricMedianWeightedFunction extends AbstractGeometricMedianFunction {

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
    final double fDEps = 1e-3;//1e3*eps;
    
    /**
    the weights of each point summed in quadrature over dimension.
    * array length is obs.length/nDimensions
     */
    final double[] w;
   
    /**
     * geomedian[dimension] - obs_dimension_point.
     * has length obs.length.
     */
    final double[] diffs;
    
    /**
     * a rough number used to avoid divide by zero
     * method.
     */
    final double fds = 1e3*eps;
    
    final double invfds = 1./fds;
    
    /**
     *
     * @param observations the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3 observations={x0, y0, z0, x1, y1, z1, ...
     * x_(nPoints-1), y_(nPoints-1), z_(nPoints-1),].
     * @param numberOfDimensions the number of data dimensions present in
     * observations array.
     */
    public GeometricMedianWeightedFunction(double[] observations, int numberOfDimensions) {
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
        
        int nData = (int) (obs.length / nDimensions);
        this.w = new double[nData];
        
        this.diffs = new double[obs.length];
    }
    
    /**
     * given observed data points obs = (x_i, y_i, ...) want to solve for
     * X=(x_geo_median, y_geo_median, ...), that is X = arg min of || X-obs
     * ||_2, that is, the X which minimizes the sum of the differences, where _2
     * is notation for using L2 distances.
     *
     * the objective is calculated by setting the derivative of the weighted
     * cost function to 0.
     * <pre>
     * The weighted cost function is:
     *    C_w(X) = min_X in summation_i=1_n( w_i*|| X-z_i ||^2 )
     *
     * The weight is
     *    w_i = 1/( || X - obs_i || ), and adding eps to denom to avoid divide by 0..
     *    NOTE that the geometric median usually uses equal weights w_i=1/nPoints.
     *
     * setting deriv of C_w(X) to zero:
     *    X = (summation_i=1_n(w_i*obs_i))/(summation_i=1_n(w_i))
     *
     * Then iterated updates of estimates of X in time steps t:
     * X_(t+1) = summation_i=1_n( w_t*obs_i ) / summation_i=1_n( w_t^2 )
     * </pre>
     *
     * @param geoMedian input variable holding coordinates of current estimate
     * of geometric median.
     * @return evaluation of the objective, summation_i=1_n(|| geoMedian - obs_i
     * ||^2)/n
     */
    @Override
    public double f(double[] geoMedian) {

        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length should equal nDimensions");
        }
        
        int nData = (int) (obs.length / nDimensions);

        calculateDifferences(geoMedian, diffs);
        
        Arrays.fill(w, 0.);
        
        int i, j, d;
        double ds, ds0;
        double fSum = 0;
        //w_i = 1/[max(|| X-z_i ||_2, delta)]
        for (i = 0; i < nData; ++i) {
            ds = 0;
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;        
                ds += diffs[j]*diffs[j];                
            }
            ds0 = ds;
            ds = Math.sqrt(ds);
            w[i] = 1./(ds + eps);
            fSum += (w[i] * ds0);
        }         
                
        return fSum;
    }
    
    /**
     *
     * <pre>
       f = summation_i=1_n( w_i * || X - obs_i || )/n
           where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_0 = w_i * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_1 = w_i * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       ...
       
       </pre>
     * @param geoMedian coordinates of current estimate of geometric median
     * @return evaluation of the derivative
     */
    @Override
    public double[] der(double[] geoMedian) {

        if (true) {
            return finiteDifference(geoMedian);
        }
        
        int nData = (int) (obs.length / nDimensions);
        
        int i, j, d;
        
        
        double[] deltaWeiszfeldUpdate = new double[geoMedian.length];
        
        // a hack to deltaX with assumption of small steps, by
        //    calculating the updated geometric-median and subtracting that
        //    from the geometric-median given to the method.
        
        // update the geometric median X = (summation_i=1_n(w_i*z_i))/(summation_i=1_n(w_i))
        double wSum = 0;
        for (i = 0; i < w.length; ++i) {
            wSum += w[i];
        }
        
        assert(wSum <= (1 + fDEps));
        
        double[] geoMedian1 = new double[geoMedian.length];
        for (d = 0; d < nDimensions; ++d) {
            for (i = 0; i < nData; ++i) {                
                j = i * nDimensions + d;        
                geoMedian1[d] += (w[i] * obs[j]);
            }
            geoMedian1[d] /= wSum;
            
            deltaWeiszfeldUpdate[d] = geoMedian1[d] - geoMedian[d];
        }

        return deltaWeiszfeldUpdate;      
    }

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
