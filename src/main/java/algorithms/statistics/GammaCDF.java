package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class GammaCDF {
    
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
    
    public static double inverseCdf(double shape, double scale, double alpha) {
        
        double mean = shape * scale;
        double variance = mean * scale;
        
        // so range of "x" is appprox 0 to mean + 3*variance
        
        //Gauss-Newton's:
        //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
        
        //Gauss-Newton's with a factor (alpha) to decrease step size:
        //x_{t+1} = x_{t} - alpha*(f(x_{t}) / f'(x_{t}))
        
        final int nMaxIter = 100;
        final double tol = 1.e-4;
        
        //double xLow = 0;
        //double xHigh = mean + 3.*variance;                
        double x = mean;
        // TODO: x = (alpha > 0.5) ? (xLow + xHigh)/2 : mean;
        
        double delta = x/4.;
        
        double diffP = Double.POSITIVE_INFINITY;
        
        double p, pm, pp;
        double t1, t2;
        
        // approximated as finite difference
        double deriv;
        
        int nIter = 0;
        while (diffP > tol && nIter < nMaxIter) {
            
            p = cdf(x, shape, scale);
            
            diffP = Math.abs(alpha - p);
            
            System.out.printf("%d) alpha=%.4f p=%.4f diff=%.4f   x=%.4f delta=%.4f\n", 
                nIter, alpha, p, diffP, x, delta);
            System.out.flush();
            
            // finite difference:
            final double old_val = x;
            x += delta;
            pp = cdf(x, shape, scale);
            x = old_val - delta;
            pm = cdf(x, shape, scale);
            
    //TODO: fix error here:
    
            //newton's method step = f(x_{t}) / f'(x_{t}) where f'(x_{t} =  (pp - pm)/(2*delta)
            deriv = (pp - pm)/(2*delta);
            
            t1 = ((alpha - p)/deriv);
            
            x = old_val + t1;
            
            t2 = Math.abs(old_val - x)/2;
            
            delta = Math.min(Math.abs(t1), t2);
            
            if ((x - delta) < 0) {
                delta = x/2.;
            }
                   
            nIter++;
        }
        
        return x;
    }
    
}
