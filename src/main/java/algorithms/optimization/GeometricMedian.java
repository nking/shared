package algorithms.optimization;

import java.util.Arrays;
import thirdparty.dlib.optimization.AbstractGeometricMedianFunction;
import thirdparty.dlib.optimization.GeometricMedianWeightedFunction;

/**
 * the geometric median:
   a.k.a. spatial median, 1-median, Euclidean minisum point, and Torricelli point.
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
     https://dl.acm.org/doi/pdf/10.1145/2897518.2897647
     
     From  https://www.pnas.org/content/pnas/97/4/1423.full.pdf
     "...the problem of minimizing the weighted sum of the
    Euclidean distances from m points, in Real space of dimension d (R^d). 
    In industrial applications, this is known as the optimal location problem 
    of Weber (2). In statistics, the solution of this optimization problem
    is the spatial median or L1-MM, considered by Brown (3) and
    Small (4). As noted by Kuhn (5), the problem goes back to
    Fermat in the early seventeenth century and was generalized
    to the current form by Simpson in his Doctrine and Application of Fluxions (6). 
    In the nineteenth century, Steiner made
    significant contributions to this problem and its extensions (cf.
    Courant and Robbins; ref. 7). Thus, the problem is known as
    the Fermat–Weber location problem and also as the Euclidean–
    Steiner problem."
     
 * 
 * https://feb.kuleuven.be/public/u0017833/PDF-FILES/l1medianR2.pdf
 * NOTE: the Vardi-Zhang 2000 update method with the author's implemented non-linear
 * optimization, appears to have a runtime of O(N log_2 N)
 * 
 * NOTE: that there is research on handling sparse matrices in other applications
 * such as bundle adjustment where a sparse version of Levenberg-Marquardt is used.
 * see papers by  Lourakis and Argyros, etc.
 * http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
 * https://www.researchgate.net/profile/Antonis_Argyros/publication/221111908_Is_Levenberg-Marquardt_the_Most_Efficient_Optimization_Algorithm_for_Implementing_Bundle_Adjustment/links/00b7d51c7d377ba56e000000/Is-Levenberg-Marquardt-the-Most-Efficient-Optimization-Algorithm-for-Implementing-Bundle-Adjustment.pdf
 
 * @author nichole
 */
public class GeometricMedian {
   
    /**
     * use first derivative to find minimum sum of function defining geometric-median.
     * Gauss-Newton's method has been altered to use back-tracking with a deceasing step size.
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
     * One should use instead, the method newtonsThenVardiZhang which checks for 
     * that and updates the solution if needed.
     * 
     * @param function
     *  @param init input output variable holding the estimates for the
     * geometric median in all dimensions.
     * @return the minimum of the sum of the squared sum of thedifferences of
     * the observed points from the geometric median.
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
        
        //TODO: rewrite this search one day using more rigorous methods for
        //    implementation of finite difference within newtons
        
        //Gauss-Newton's:
        //x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
        
        //Gauss-Newton's with a factor (alpha) to decrease step size:
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
            // NOTE that if the equation were not convex, would also pursue additional
            //    information about the point where gradient becomes 0 using
            //    the hessian (pt being min, max, or saddle point).
            // critical point:
            //     min for hessian positive definite (all the eigenvalues are positive)
            //     max for hessian negative definit
            //     saddle for hessian indefinite
            
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
     * this methods uses newtonsMethod2 to completion and then uses the 
     * Vardi Zhang (2000) algorithm to update argmin X if it's a point in obs.
     *
     * Vardi and Zhang 2000:
     * "The multivariate L1-median and associated data depth",
       Yehuda Vardi and Cun-Hui Zhang
     *     https://www.pnas.org/content/pnas/97/4/1423.full.pdf
     * 
     * @param function
     * @param geoMedian input output variable holding the estimates for the
     * geometric median in all dimensions.
     * @return the minimum of the sum of the squared sum of the differences of
     * the observed points from the geometric median.
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
                System.err.println("check for an error in impl of method " +
                    function.getClass().getSimpleName() + ".vardiZhang");
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
    
    /**
     * The Levenberg-Marquardt method is similar to the Gauss-Newton method, but adds
     * a term for invertibility and L1 regularization to 
     * f' as the denominator of the update.
     * The added term is λ_t*I_nxn.  
     * It's a damping parameter, for the iteration function that can be changed 
     based on the result in each step.
     When λ → 0 , the added term vanishes and the technique reverts to Newton’s 
     method.   
     When λ becomes large, the scheme becomes the gradient descent method. 
     This improves the robustness of the algorithm when the initial values are 
     far from the final minimum.
     
     Marquardt suggested starting with a small value for λ. 
     If the results of the previous step improves the objective function, 
     x_{t+1} is incremented by the update and the the value of λ is decreased 
     (say by a factor of 2) and the method continues. 
     If the method (unfortunately) increases the objective function, the step 
     is discarded and λ is increased.
     
     Marquardt also replaced the identity matrix I with the diagonal matrix 
     consisting of the elements of (DF(x_{t}))^T*(DF(x_{t})) to create
     varied step sizes with direction (small gradient gets larger step).
     
     Improved Levenberg-Marquardt by Jia borrows some ideas from the simulated 
     annealing method (SA).   It interprets slow cooling as a slow decrease in 
     the probability of accepting worse solutions as it explores the solution 
     space.
     
     The ILM is provided a method that optimizer has a chance to get out of the 
     local minimum, but there is no guarantee that the global minimum can be 
     reached at the end.
     
       <pre>
       Gauss-Newton:
           x_{t+1} = x_{t} - f(x_{t}) / f'(x_{t})
        
           let F(x_{t}) = stacks of f(x_{t})_i into a column vector,
           and DF is the jacobian of F
            
           x_{t+1} = x_{t} - ( (DF(x_{t}))^T*(DF(x_{t})) )^-1 * (DF(x_{t}))^T*F(x_{t})
           
      Levenberg-Marquardt:
           x_{t+1} = x_{t} - ( (DF(x_{t}))^T*(DF(x_{t})) + lambda_t*I_nxn )^-1 * (DF(x_{t}))^T*F(x_{t})
           
      Improved Levenberg-Marquard (ILM)t
           see pg 199 og Jia thesis (can start at pg 208-ish, update is on pg 220).
           
     </pre>
     
        material is from the book "Numerical Algorithms" by Justin Solomon,
        and from publications of
        "Fitting a Parametric Model to a Cloud of Points Via Optimization Methods"
        by Pengcheng Jia
        https://surface.syr.edu/cgi/viewcontent.cgi?article=1673&context=etd
        
     * @param function
     * @param init
     * @return 
     */
    /*public double levenbergMarquardt(AbstractGeometricMedianFunction function,
        double[] init) {
        
    }*/

   /*
    NOTE: this method can only be used on 2 dimensional data.  
    NOTE: unit standardization is performed internally and de-normalization
    if needed, so no need to pre-process the data for that.
    
    an iterative linear least squares.
    
    <pre> 
    F(M) = Summation_over_i( || obs_i - M || )
        goal: find M which minimizes F(M)
    example with n=2 dimensions and objective 
          
    
    Full Rank (includes OVER-DETERMINED, m .gt. n, rank==n):
        A * x = b
        A is mxn
        x is m length
        b is n length
           has no exact solutions.
           minimize the error b - A*x to find x.
           set deriv to 0.
    
        x = (A^T A)^−1 A^T b where (A^T A)^−1 A^T (a m×n matrix) is the pseudo-inverse
    
    Rank Deficient (includes UNDER-DETERMINED, n .lt. m):
        A * x = b
        A is mxn
        x is m length
        b is n length
           has infinitely many solutions.
           columns of A are independent.
           find smallest by minimizing x subject
           to constraint b = A*x.
           adding a lagrange multiplier then set deriv to 0.
    
        x = A^T*(A A^T)^−1 * b where (A^T A)^−1 A^T (a m×n matrix) is the pseudo-inverse
            where A^T*(A A^T)^−1 (a m×n matrix) is the pseudo-inverse
    
        the pseudo-inverse is calculated from the SVD:
            
        A = U*S*V^T from SVD  
            U is mxn orthonormal columns
            S is nxn with non-negative singular values.  rank is number of non-zero entries
            V is  nxn
        x_LS = summation_over_rank((u_i^T * b / s_i) * v_i)
    
        note, the pseudo-inverse was V*R*U^T where R is 1/diagonal of S
    
    
    RANK==NDATA
        r == m  and  r == n  square an invertible   A*x=b   has 1 soln
        columns of A are independent.
        use Full-Rank solution, but if I is invertible, can use the inverse of A instead.
    </pre>
    @param init input variable holding coordinates of current estimate of
    geometric median.
    @return evaluation of the objective, summation_i=1_n(||geoMedian - obs_i||^2)/n
    */
    //private double leastSquares(AbstractGeometricMedianFunction function,
    //    double[] init) throws NotConvergedException {
    
        /*        
        following Strang, Linear Algebra for the full rank case.
        can solve in many ways.
        
        Example for unweighted geometric median:
           Input: (-20, 48), (-20, -48), (20, 0), (59, 0)
           Geometric Median = (20,0) with minimum distance = 163.964
        
        By Algebra:
        =============
            ||A*x-b||^2 isn't solvable, but using the projection is:
               ||A*x-p||^2 + ||e||^2 where e is small
        
        let A = point matrix
        -20   48
        -20  -48
         20    0
         59    0
        
        b would have to be a vector of the calculated distance of each point from
          the current geometric-median.
        
        solving for x
        
        presumably, one could use x as a factor on a property of the points
        to update the current geometric-median.
        
        let a = point matrix    calculated b=
        -20   48             66.06814663663572
        -20  -48             61.554853586049575
         20    0             3.605551275463989
         59    0             37.12142238654117
        
        a^T*a =
      -20  -20  20  59 ]  * [ -20   48 =[ -20*-20-20*-20+20*20+59*59  -20*48-20*-48+0+0 ] = [4681 0   ]
       48  -48   0   0 ]    [ -20  -48  [  48*-20-48*-20*0+0+0        48*48-48*-48+0+0      [-960 4608]
                            [  20    0
                            [  59    0

      a^T*b =
      -20  -20  20  59 ]  * [66.068  = [-20*66.068 + -20*61.555 + 20*3.606 + 59.*37.121 = [-290.26
       48  -48   0   0 ]    [61.555    [ 48*66.068 + -48*61.555 +  0*3.606 +  0.*37.121   [216.624
                            [3.606
                            [37.121

       a^T*a*x = a^T*b
      ( x = (a^T/(a^T*a)) * b )
      pseudoinverse(a^T*a) =  4608  960         -290.26  = [-0.0523
                              0     4681         216.624    0.0470
                              --------------
                              4681*4608 + 0

                 ? -0.0523*x_mean = -0.509925
                    0.0470*y_mean = 0   <--- can never use the found x to update when the mean is 0
                    standard deviation would presumably be better:
                    stdev = (37.9726, 39.3446)
                    update of median_0 = -0.0523*37.973 = -1.99
                    update of median_1 = 0.0470*39.345   = 1.85
        
        
        ==============
        By Calculus
        ==============
        df/dx_0 =  2*(-20*x_0 + 48*x_1 - 66.068)*(-20) + 2*(-20*x_0 - 48*x_1 - 61.555)*(-20) + 2*(20*x_0 - 3.606)*(20) +2*(59*x_0 - 37.121)*(59)
           = x_0*(9362.0) + x_1*(0) + 580.402
           = x_0*(9362.0) + 580.402 = 0  x_0==> -0.061
        
        df/dx_1 =  2*(-20*x_0 + 48*x_1 - 66.068)*(48) + 2*(-20*x_0 - 48*x_1 - 61.555)*(-48)
           =  x_0*(2.*-20*48 + 2.*-20*-48) +x_1*(2.*48.*48 + 2.*-48*-48) + 2.48.*-66.068 + 2.*-61.555*-48.
           =  x_0*(0) +x_1*(9216) -433.248
           =  x_1*(9216) -433.248 = 0 ==> x_1 = 0.047
        
        so the values of x_0 and x_1 would give same suggested updates of median.
        
        Then another iteration would be needed, and the cycle repeated until
        stopping conditions were met.
        
        */        
        
    //}
}
