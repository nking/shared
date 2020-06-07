package algorithms.optimization;

import java.util.Arrays;
import thirdparty.dlib.optimization.AbstractGeometricMedianFunction;

/**
 *
 * @author nichole
 */
public class GeometricMedian {
 
    /**
     * use first derivative to find minimum sum of function defining geometric-median.
     * newton's method has been altered to use back-tracking to previous solution
     * if function evaluation becomes larger, then proceeds forward again
     * with a smaller step.
     * 
     * @param function
     * @param init
     * @return 
     */
    public double newtonsMethodAltered(AbstractGeometricMedianFunction function,
        double[] init) {
        
        //TODO: for nPoints==nDimensions, is the centroid always a g.m.?
        
        //Newton's:
        //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
        
        //Newton's with a factor to decrease step size:
        //x_{t+1} = x_{t} - alpha*(f(x_{t}) / f'(x_{t}))
        
        double alpha = 1.0;
        double alphaFactor = 0.2;//0.33;
        
        // iterate until derivative is 0 or essentially no change of estimate
        int nDimensions = function.getNDimensions();
        int nData = function.getNData();
        double fds = 0.001;//TODO set this as a fraction of the data ranges
        double eps = 1.e-17;
        
        if (init.length != nDimensions) {
            throw new IllegalArgumentException("init.length must equal function.getNDimensions()");
        }
                
        double[] geoMedian1 = Arrays.copyOf(init, init.length);
        double[] prevGeoMedian1;
        double fEval = function.f(geoMedian1);
        double prevFEval;
        double[] fDerEval = function.der(geoMedian1);
        double[] fEval2 = function.evaluateGeoMedianPerDimension(geoMedian1);
        double[] prevFEval2;
        
        int i, j, d, cg;
        double r2;
        int c = 0;
        int backtracked = 0;
                
        while (true) {
            prevFEval = fEval;
            prevFEval2 = Arrays.copyOf(fEval2, fEval2.length);
            prevGeoMedian1 = Arrays.copyOf(geoMedian1, geoMedian1.length);
            
            for (i = 0; i < nDimensions; ++i) {
                //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
                r2 = alpha * (fEval2[i]/(double)nData) / (fDerEval[i] + eps);
                geoMedian1[i] -= r2;
            }
            
            fEval = function.f(geoMedian1);
            fDerEval = function.der(geoMedian1);
            fEval2 = function.evaluateGeoMedianPerDimension(geoMedian1);
            
            System.out.printf("f=%.3e  der=%s\n  => gm=(%s)\n  prevf=%.3e %s\n  prevgm=(%s)\n",
                fEval, AbstractGeometricMedianFunction.toString(fDerEval),
                AbstractGeometricMedianFunction.toString(geoMedian1),
                prevFEval, AbstractGeometricMedianFunction.toString(prevFEval2),
                AbstractGeometricMedianFunction.toString(prevGeoMedian1)
                );
            System.out.flush();
            
            // backtrack and reduce step size.
            //  following description of backtracking in Newton's method from
            //  "Numerical Recipes in C", Section 9.7
            //  by William H. Press, Saul A. Teukolsky, William T. Vetterling and Brian P. Flannery
            if ((fEval > prevFEval) || (Math.abs(fEval - prevFEval) <= fds))  {
                c++;
                if (backtracked > 1) {
                    // return this solution or previous
                    System.arraycopy(prevGeoMedian1, 0, init, 0, init.length);
                    return prevFEval;
                }
                backtracked++;
                fEval = prevFEval;
                fEval2 = Arrays.copyOf(prevFEval2, fEval2.length);
                geoMedian1 = Arrays.copyOf(prevGeoMedian1, geoMedian1.length);
                alpha *= alphaFactor;
                continue;
            } else {
                backtracked = 0;
            }
            
            // stopping strategy
            cg = 0;
            for (i = 0; i < nDimensions; ++i) {
                if (Math.abs(fDerEval[i]) <= fds) {
                    cg++;
                }
            }
            if (cg == nDimensions) {
                System.arraycopy(geoMedian1, 0, init, 0, init.length);
                return fEval;
            }
            c++;
        }
    }
    
}
