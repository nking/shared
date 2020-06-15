package algorithms.optimization;

import java.util.Arrays;
import thirdparty.dlib.optimization.AbstractGeometricMedianFunction;
import thirdparty.dlib.optimization.GeometricMedianUnweightedFunction;
import thirdparty.dlib.optimization.GeometricMedianWeightedFunction;

/**
 * the geometric median:
   a.k.a. spatial median.
   the weighted version is a.k.a. L1-median (though it uses the euclidean distance)
   and the multivariate L1 -median (L1 -MM).
   
   definition: the point which minimizes the sum of the euclidean distance of that
     point to all other points in the set.
     the sum is a convex function (i.e. local search will work).
     Unfortunately, no algorithms are closed form, that is no algorithms have a
     finite number of computational operations.
     The geometric median is a rotation and translation invariant estimator that
     achieves the optimal breakdown point of 0.5, i.e. it is a good estimator
     even when up to half of the input data is arbitrarily corrupted.
     (https://dl.acm.org/doi/pdf/10.1145/2897518.2897647)
 * 
 * https://feb.kuleuven.be/public/u0017833/PDF-FILES/l1medianR2.pdf
 * NOTE: the Vardi-Zhang 2000 update method with the author's implemented non-linear
 * optimization, appears to have a runtime of O(N log_2 N)
 * 
 * 
 * @author nichole
 */
public class GeometricMedian {
   
    /**
     * use first derivative to find minimum sum of function defining geometric-median.
     * Newton's method has been altered to use back-tracking with a deceasing step size.
     * The first derivative of the geometric-median is not zero so can use
     * first order Householder methods, but cannot use 2nd order Householder methods
     * because the 2nd derivative is 0.
     * <pre>
       example for nDimensions = 2:
       f = summation_i=1_n( || X - obs_i || )/n
           where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2)
       df/dX_0 = d/dx( (1/n)*summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2) ))
               = (1/n)*(1/2)*2*(X_0-obs_i_0)*(1)
                  / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2) )
               = (1./n)*(X_0-obs_i_0) / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2) )
       df/dX_1 = (1./n)*(X_1-obs_i_1) / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2) )

       d/dX_0 of df/dX_0 = (1./n) * (1) * (-1/2)*summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(-3/2) )
                           *2*(X_0-obs_i_0)*(1)
                         = (-1./n)*(X_0-obs_i_0) / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(3/2) )
       
       Hessian d/dX of d/dX where p is nDimensions and I is the identity matrix of size pxp:
                            ( (    I_p      )   ( (X - obs_i)*(X - obs_i)^T )
           = summation_i=1_n( (-------------) - ( --------------------------)
                            ( (||X - obs_i||)   (      ||X - obs_i||^3      )
         ...

     * </pre>
     * 
     * NOTE: if nDimensions == nDataPoints, the init is replaced with the
     * euclidean centroid.  These are "colinear" cases extended to more than
     * 2 dimensions.  There are many solutions for these cases, but the 
     * euclidean centroid is chosen here.
     * 
     * NOTE, algorithm currently may fail to minimize if arrives at demand points.
     * 
     * @param function
     * @param init
     * @return 
     */
    public double newtonsMethod2(AbstractGeometricMedianFunction function,
        double[] init) {
        
        int nDimensions = function.getNDimensions();
        int nData = function.getNData();
        
        if (init.length != nDimensions) {
            throw new IllegalArgumentException("init.length must equal nDimensions");
        }
        
        //TODO: handle this condition for unweighted spatial median quickly, 
        // no iteration needed:
        // "If more than n/2 observations are concentrated in one point, say y, 
        // the solution of the L1-median is y"
        /*if (function instanceof GeometricMedianUnweightedFunction) {
            double[] pt = occursAsMoreThanHalf(function);
            if (pt != null) {
                System.arraycopy(pt, 0, init, 0, init.length);
                return function.f(init);
            }
        }*/
        
        // assuming for now, that this is true:
        if (nData == nDimensions) {
            double[] cen = function.calculateCentroid();
            System.arraycopy(cen, 0, init, 0, nDimensions);
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
            System.out.printf("[gm=(%s)  f=%.3e der=(%s)]\n",
                AbstractGeometricMedianFunction.toString(prevGeoMedian1),
                prevFEval,
                AbstractGeometricMedianFunction.toString(fDerEval));
            System.out.flush();
           
            for (i = 0; i < nDimensions; ++i) {
                if (Math.abs(fDerEval[i]) < eps) {
                    // finite difference, secant method, broyden's method:
                    fd = function.finiteDifference(prevGeoMedian1);
                    System.arraycopy(fd, 0, fDerEval, 0, nDimensions);
                    usedFiniteDifference = true;
                    break;
                }
            }
            if (usedFiniteDifference) {
                System.out.printf("  [finite difference=(%s)]\n",
                    AbstractGeometricMedianFunction.toString(fDerEval));
                System.out.flush();
                
                // stopping strategy
                cg = 0;
                for (i = 0; i < nDimensions; ++i) {
                    if (Math.abs(fDerEval[i]) <= fds) {
                        cg++;
                    }
                }
                if (cg == nDimensions) {
                    System.arraycopy(prevGeoMedian1, 0, init, 0, init.length);
                    return prevFEval;
                }
            }
            
            // first derivative is used here, but not the hessian.
            // could use the hessian to refine the data step size.
            // hessian: d^2f(x)/dx_i dx_j
            //    x_{t+1} = x_{t} - (1/H(x_{t})) * f'(x_{t})
            //       with alpha as a line search step size again.
            // NOTE that if the equation were not convex, could get additional
            //    information about the point where gradient becomes 0 using
            //    the hessian (pt being min, max, or saddle point)
            
            for (i = 0; i < nDimensions; ++i) {
                //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
                if (Math.abs(fEval2[i]) < eps) {
                    r2 = 0;
                } else {
                    r2 = alpha * (fEval2[i]/(double)nData) / (fDerEval[i] + eps);
                }
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
    
    /**
     * if init is a root of the function, and error is errorTolerance or smaller,
       then the function derivative evaluated at the root plus and at the root minus
       the error will have opposing signs.
     * @param function
     * @param init
     * @return 
     */
    public boolean verify(AbstractGeometricMedianFunction function,
        double[] init, double errorTolerance) {
        
        double d = errorTolerance;
        double[] minus = Arrays.copyOf(init, init.length);
        double[] plus = Arrays.copyOf(init, init.length);
        for (int i = 0; i < minus.length; ++i) {
            minus[i] -= d;
            plus[i] += d;
        }
        double[] fMinus = function.der(minus);
        double[] fPlus = function.der(plus);
        for (int i = 0; i < minus.length; ++i) {
            if ((fMinus[i] < 0 && fPlus[i] > 0) || (fMinus[i] > 0 && fPlus[i] < 0)) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * run newtonsMethod2 to completion and use Vardi Zhang (2000) algorithm
     * to update argmin X if it's a point in obs.
     *
     * Vardi & Zhang 2000:
     * "The multivariate L1-median and associated data depth",
       Yehuda Vardi and Cun-Hui Zhang
     *     https://www.pnas.org/content/pnas/97/4/1423.full.pdf
     * 
     * @param function
     * @param geoMedian
     * @return 
     */
    public double newtonsThenVardiZhang(GeometricMedianWeightedFunction function,
        double[] geoMedian) {
        
        // scale these by the data weights?
        double fds0 = 1.e-4;
        double fds = 1.e-3;
        
        boolean[] checks = new boolean[2];
        int checkTol = 0;
        double f0, diffF;
        double f1 = Double.MAX_VALUE;
        int[] isMedian = new int[function.getNData()];
        double[] prevGeoMedian = new double[geoMedian.length];
        int nIter = 0;
        int nIterMax = 100;
        while (nIter < nIterMax) {
            
            f0 = newtonsMethod2(function, geoMedian);

            System.out.printf("%d) newton's f=%.7e gm=%s\n", nIter, f0, 
                AbstractGeometricMedianFunction.toString(geoMedian));
            
            function.isMedian(geoMedian, isMedian);

            System.arraycopy(geoMedian, 0, prevGeoMedian, 0, geoMedian.length);
            
            f1 = function.vardiZhang(isMedian, geoMedian, checks);
        
            // this is negative while minimizing
            diffF = f1 - f0;
            if (!(diffF <= 0.) ) {
                // NOTE: if this happens, I may have a bug in vardiZhang
                f1 = f0;
                System.arraycopy(prevGeoMedian, 0, geoMedian, 0, geoMedian.length);
                break;
            }
            
            System.out.printf("  vz f=%.7e gm=%s checks=%s\n", f1, 
                AbstractGeometricMedianFunction.toString(geoMedian),
                Arrays.toString(checks));
            System.out.flush();
            
            // stopping criteria:
            if (checks[0] && checks[1]) {
                break;
            }
            
            // stopping criteria:
            if (Math.abs(diffF) < fds0) {
                for (checkTol = 0; checkTol < geoMedian.length; ++checkTol) {
                    if (Math.abs(prevGeoMedian[checkTol] - geoMedian[checkTol])
                        > fds) {
                        break;
                    }
                }
                if (checkTol == geoMedian.length) {
                    break;
                }
            }
            
            nIter++;
        }

        return f1;
    }

    private double[] occursAsMoreThanHalf(AbstractGeometricMedianFunction function) {
        double[] obs = function.getObs();
        int nDimensions = function.getNDimensions();
        int nData = function.getNData();
        
        // store points and frequency
        // TODO: create a data structure like PaiInt to hold a double array
        //   with edited equals for a tolerance of equality
        //   (presumably 1e-17) and a hashcode to return same hash for
        //   similar entries truncated to tolerance.
        int i, j, d;
        double a, b;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                
            }
        }
        throw new UnsupportedOperationException("not currently implemented");
    }

    
}
