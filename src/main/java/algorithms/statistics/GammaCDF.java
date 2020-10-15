package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class GammaCDF {
    
    public static double cdf(double x, double shape, double scale) {
     
        // from https://www.mathworks.com/help/stats/gamcdf.html
        
        double a = shape;
        double b = scale;
        
        if (x < 0) {
            throw new IllegalArgumentException("x cannot be a negative number");
        }
        if (a <= 0) {
            throw new IllegalArgumentException("a must be a positive number");
        }
        if (b <= 0) {
            throw new IllegalArgumentException("b must be a positive number");
        }
        
        return thirdparty.smile.math.special.Gamma.regularizedIncompleteGamma(a, x/b);
    }
    
    public static double inverseCdf(double shape, double scale, double alpha) {
        // TODO: will use cdf within iteration such as Newton's method
        //   to find when p = cdf(x, shape, scale) is approx alpha
        
        // the parameters to use to produce maximum "x" seems to be when shape=x
        
        throw new UnsupportedOperationException("not yet implemented");
    }
}
