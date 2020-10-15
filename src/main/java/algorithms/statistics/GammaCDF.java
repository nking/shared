package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class GammaCDF {
    
    public static double cdf(double x, double shape, double scale) {
        
        double numer = 
            thirdparty.smile.math.special.Gamma.regularizedIncompleteGamma(shape, x/scale);
        
        double denom = thirdparty.smile.math.special.Gamma.gamma(shape);
        
        return numer/denom;
    }
    
    public static double inverseCdf(double shape, double scale, double alpha) {
        // TODO: will use cdf within iteration such as Newton's method
        //   to find p = cdf(x, shape, scale) which is close to alpha
        throw new UnsupportedOperationException("not yet implemented");
    }
}
