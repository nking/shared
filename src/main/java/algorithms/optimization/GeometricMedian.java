package algorithms.optimization;

import algorithms.misc.MiscMath0;
import java.util.Arrays;
import thirdparty.dlib.optimization.AbstractGeometricMedianFunction;

/**
 *
 * @author nichole
 */
public class GeometricMedian {
   
    /**
     * use first derivative to find minimum sum of function defining geometric-median.
     * newton's method has been altered to use back-tracking with a deceasing step size.
     * 
     * NOTE: if nDimensions == nDataPoints, the init is replaced with the centroid.
     * euclidean centroid.
     * @param function
     * @param init
     * @return 
     */
    public double newtonsMethod2(AbstractGeometricMedianFunction function,
        double[] init) {
        
        int nDimensions = function.getNDimensions();
        int nData = function.getNData();
        
        // assuming for now, that this is true:
        if (nData == nDimensions) {
            init = function.calculateCentroid();
        }
        
        //NOTE: if observed data given to function has been standardized,
        //   then it is "zero centered", that is, has a mean of 0.
        //   The derivative of the geometric-median will be 0 for points
        //   equal to the centroid whether standardized or not.  
        //   The centroid is a true geometric-median for
        //   some data, but for those that it is not, alternative methods of
        //   estimating the gradient at that point are needed.
        //   ... added logic to try finite difference if the derivatives are
        //   zero.
        
        //Newton's:
        //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
        
        //Newton's with a factor (alpha) to decrease step size:
        //x_{t+1} = x_{t} - alpha*(f(x_{t}) / f'(x_{t}))
        
        double alpha = 1.0;
        double alphaFactor = 0.2;//0.33;
        
        // iterate until derivative is 0 or essentially no change of estimate
        
        //TODO set this as a fraction of the data ranges if data is not standardized
        double fds = 0.001;
        
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
        boolean usedFiniteDifference;
        double[] fd;
        
        int i, j, d, cg;
        double r2;
        int c = 0;
        int backtracked = 0;
                
        while (true) {
            prevFEval = fEval;
            prevFEval2 = Arrays.copyOf(fEval2, fEval2.length);
            prevGeoMedian1 = Arrays.copyOf(geoMedian1, geoMedian1.length);
            
            // re-calculate any values of fDerEval that are 0, using finite difference:
            usedFiniteDifference = false;
            System.out.printf("[gm=(%s)  der=(%s)]\n",
                AbstractGeometricMedianFunction.toString(geoMedian1),
                AbstractGeometricMedianFunction.toString(fDerEval));
            System.out.flush();
            for (i = 0; i < nDimensions; ++i) {
                if (Math.abs(fDerEval[i]) < eps) {
                    // TODO: create a per-dimension method for finite difference and use it here as alt
                    fd = function.finiteDifference(geoMedian1);
                    System.arraycopy(fd, 0, fDerEval, 0, nDimensions);
                    usedFiniteDifference = true;
                    break;
                }
            }
            if (usedFiniteDifference) {
                System.out.printf("  [finite difference=(%s)]\n",
                    AbstractGeometricMedianFunction.toString(fDerEval));
                System.out.flush();
            }
            
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
