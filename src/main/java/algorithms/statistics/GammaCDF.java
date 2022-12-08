package algorithms.statistics;

/**
 *
 *  note, can find method of moments for Gamma and many other functions here:
 *  https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/07%3A_Point_Estimation/7.02%3A_The_Method_of_Moments
 *
 * @author nichole
 */
public class GammaCDF {
    
    /**
     *
     @param x
     @param shape
     @param scale
     @return
     */
    public static double cdf(double x, double shape, double scale) {
     
        // from https://www.mathworks.com/help/stats/gamcdf.html
                
        if (x < 0) {
            throw new IllegalArgumentException("x cannot be a negative number");
        }
        if (shape <= 0) {
            throw new IllegalArgumentException("shape must be a positive number");
        }
        if (scale <= 0) {
            throw new IllegalArgumentException("scale must be a positive number");
        }
        
        return thirdparty.smile.math.special.Gamma.regularizedIncompleteGamma(shape, x/scale);
    }
    
    /**
     *
     @param shape
     @param scale
     @param alpha
     @return
     */
    public static double inverseCdf(double shape, double scale, double alpha) {
     
        if (alpha < 0 || alpha > 0.999) {
            throw new IllegalArgumentException("alpha must be [0, 1] and has an artificial upper limit here of 0.999 for search");
        }
        if (shape <= 0) {
            throw new IllegalArgumentException("shape must be a positive number");
        }
        if (scale <= 0) {
            throw new IllegalArgumentException("scale must be a positive number");
        }
        
        if (shape > 1000) {
            throw new IllegalArgumentException("currently using an artificial maximum allowed value of 1000 for shape");
        }
        if (scale > 1000) {
            throw new IllegalArgumentException("currently using an artificial maximum allowed value of 1000 for scale");
        }
        
        // use binary search
    
        final double mean = shape * scale;
        // variance = shape * scale*scale
        final double variance = mean * scale;
        
       
        final int nMaxIter = 1000;
        final double tolP = 1.e-3;
        final double tolX = 1.e-3;
        
        double xLow = 0;
        double xHigh = mean + 3.2*variance;                
        double x = (xLow + xHigh)/2;
        
        // may need to change how this is estimated.
        double delta = Math.min(xHigh - x, x - xLow);
        delta /= 2;
        
        double diffP = Double.POSITIVE_INFINITY;
        
        double p0, pm, pp;
        
        int nIter = 0;
        while (nIter < nMaxIter && xHigh > xLow) {
            
            x = (xLow + xHigh)/2;
            
            delta = Math.min(xHigh - x, x - xLow);
            delta /= 2;
            
            p0 = cdf(x, shape, scale);
            diffP = p0 - alpha;
            if (Math.abs(diffP) < tolP) {
                return x;
            }
                        
            /*System.out.printf("%d) alpha=%.4f p=%.4f diff=%.4f   x=%.4f delta=%.4f\n", 
                nIter, alpha, p0, diffP, x, delta);
            System.out.flush();*/
            
            if (diffP > 0) {
                if (Math.abs(xHigh - x) < tolX) {
                    xHigh -= delta;
                } else {
                    xHigh = x;
                }
            } else if (diffP < 0) {
                if (Math.abs(xLow - x) < tolX) {
                    xLow += delta;
                } else {
                    xLow = x;
                }
            }
            
            nIter++;
        }
        double minX = x;
        p0 = cdf(x, shape, scale);
        diffP = p0 - alpha;
        delta = Math.min(xHigh - x, x - xLow);
        delta /= 2;
        
        if ((x + delta) < xHigh) {
            if (Math.abs(cdf(x + delta, shape, scale) - alpha) < diffP) {
                minX = x + delta;
            }
        }
        if ((x - delta) >= 0) {
            if ( Math.abs(cdf(x - delta, shape, scale) - alpha) < diffP) {
                minX = x - delta;
            }
        }
        
        return minX;
    }
}
