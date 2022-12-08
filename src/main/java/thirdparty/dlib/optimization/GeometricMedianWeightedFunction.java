package thirdparty.dlib.optimization;

import java.util.Arrays;

/**
 * objective function implementations for use with search algorithms using IFunction
 * as input for the geometric-median.
 * 
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
     
      <pre>
      f = summation_i=1_n( || X - obs_i || )/n     
               where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_0 = (0.5/n) * ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(-1/2) * (2*(X_0-obs_i_0))
              = (1./n) * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)
      df/dX_1 = (1./n) * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 ...)^(1/2)  
      
       d/dX_0 of df/dX_0 = (1./n) * (1) * (-1/2)*summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(-3/2) )
                           *2*(X_0-obs_i_0)*(1)
                         = (-1./n)*(X_0-obs_i_0) / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(3/2) )
       d/dX_1 of df/dX_0 = 0
       Hessian d/dX of d/dX where p is nDimensions and I is the identity matrix of size pxp:
                                ( (    I_p      )   ( (X - obs_i)*(X - obs_i)^T )
               = summation_i=1_n( (-------------) - ( --------------------------)
                                ( (||X - obs_i||)   (      ||X - obs_i||^3      )
      </pre>
      
      For the weighted geometric-median, there are 2 papers will use (combined with an improved Newton's method):
     The Vardi-Zhang 2000 algorithm:
       "The multivariate L1-median and associated data depth", Yehuda Vardi† and Cun-Hui Zhang
       https://www.pnas.org/content/pnas/97/4/1423.full.pdf
     This review which compares several algorithms, including Vardi-Zhang 2000:
       "A comparison of algorithms for the multivariate L1-median",
       Heinrich Fritz, Peter Filzmoser, Christophe Croux
       https://feb.kuleuven.be/public/u0017833/PDF-FILES/l1medianR2.pdf
       (eqns (7)-(11) are Vardi-Zhang)

    <pre>
     let X be the geometric-median
     let n be the number of observations obs_i
     let M be the L1-median:
        M is a function of observed data obs_i and multiplicites eta_i
           where eta_i are used to make weights.
        M = argmin X of C(X)
           where C(X) = summation_i_1_to_n( eta_i*d_i(X) )
           where d_i is the euclidean distance between a point obs_i and X in all dimensions.
     ==> X=M iff T(X) .eq. X iff r(X) .lte. eta(X)􏰌􏰅
           where eta(X) = eta_i if X=obs_i, i.e. geometric-median is a point in the set obs
                        else eta(X) = 0;
           where r(X) = ||R̃(X)||
           where R̃(X) = summation_over_obs_except_X( eta_i*(obs_i-X)/||obs_i-X|| )
           where T(X) = (1-(eta(X)/r(X))) * T̃(X) + X*math.min(1, (eta(X)/r(X)))
             where 0/0 = 0 in the computation of eta(X)/r(X),
           where T̃(X) = summation_over_obs_except_X( eta_i*(obs_i)/||obs_i-X|| )
                         / summation_over_obs_except_X( eta_i/||obs_i-X|| )
         NOTE: after iterative algorithm solves M, the test
               X=M iff T(X) .eq. X iff r(X) .lte. eta(X) is performed.
    </pre>
     
  NOTE: consider using Standardization normalization  on
  the data before using this class then use Standardization de-normalization
  on the resulting geometric-median afterwards and run the evaluation on that
  result.
  
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
     * a rough number for use in the finiteDifference and for a tolerance.
     */
    final double fDEps = 0.001;
    
    /**
     * eta_i = factor per point given by user, acting as a multiplicity of the
     * point essentially.
     */
    final double[] eta;
   
    /**
     *
     @param observations the observations from all dimensions in format of
     * point_0 in all dimensions, followed by point_1 in all dimensions, etc.
     * e.g. for numberOfDimensions=3 observations={x0, y0, z0, x1, y1, z1, ...
     * x_(nPoints-1), y_(nPoints-1), z_(nPoints-1),].
     @param numberOfDimensions the number of data dimensions present in
     * observations array.
     @param multiplicities factor per point given by user, acting as a 
     * multiplicity of each point essentially.  it's a weight vector.
       <pre>
     * the weight w_i as defined by m_i/d_i
     * </pre>
     */
    public GeometricMedianWeightedFunction(double[] observations, 
        int numberOfDimensions, double[] multiplicities) {
        if (numberOfDimensions < 1) {
            throw new IllegalArgumentException("numberOfDimensions must be > 0");
        }
        this.nDimensions = numberOfDimensions;
        this.obs = Arrays.copyOf(observations, observations.length);
        
        int n = obs.length;
        int nData = n / nDimensions;

        if (n < 1) {
            throw new IllegalArgumentException("observations must have length > 0");
        }

        if ((n % nDimensions) != 0) {
            throw new IllegalArgumentException("observations.length must be a multiple of "
                + "numberOfDimensions");
        }
        
        if (multiplicities.length != nData) {
            throw new IllegalArgumentException("multiplicities.length must "
                + "equal nData");
        }
        
        this.eta = Arrays.copyOf(multiplicities, nData);
    }
    
    /**
     * 
     * <pre>
     * The weighted cost function is:
     *   C_w(X) = summation_i_1_to_n( eta_i*d_i(X) )
           where d_i is the euclidean distance between a point obs_i and X in all dimensions.
     *
     @param geoMedian input variable holding coordinates of current estimate
     * of geometric median.
     @return evaluation of the cost function for the given geometric-median
     */
    @Override
    public double f(double[] geoMedian) {
  
        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length should equal nDimensions");
        }
        
        int nData = (obs.length / nDimensions);
        
        //geoMedian - obs_i
        double[] diffs = calculateDifferences(geoMedian);
                
        int i, j, d;
        double dist;
        double sum = 0;
        for (i = 0; i < nData; ++i) {
            dist = 0;
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                dist += (diffs[j] * diffs[j]);
            }
            dist = Math.sqrt(dist);
            sum += eta[i] * dist;
        }
        
        return sum;                
    }
    
    /**
     *
       <pre>
       f = summation_i=1_n( eta_i * || X - obs_i || )/n
           where || X - obs_i ||_2 is ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_0 = eta_i * (X_0-obs_i_0) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
       df/dX_1 = eta_i * (X_1-obs_i_1) / ( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 + ...)^(1/2)
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
        
        //geoMedian - obs_i
        double[] diffs = calculateDifferences(geoMedian);
        
        double[] ssdPerPoint = calculateSSDPerPoint(diffs);
        
        int nData = (obs.length / nDimensions);
        
        assert(ssdPerPoint.length == nData);
        
        double s = 0;
        int i;
        for (i = 0; i < ssdPerPoint.length; ++i) {
            s += Math.sqrt(ssdPerPoint[i]);
        }
        
        // to avoid divide by 0:
        s += eps;
        
        //df/dX_0 = (1./n)*(X_0-obs_i_0) / summation_i=1_n( (X_0-obs_i_0)^2 + (X_1-obs_i_1)^2 )^(1/2) )
        double[] dfDX = new double[nDimensions];
        
        // NOTE: can see that the unweighted algorithm will not make progress when
        //    the current geometric-median estimate is the centroid.
        int j, d;
        for (i = 0; i < nData; ++i) {            
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                dfDX[d] += eta[i] * (geoMedian[d] - obs[j]);
            }
        }
        
        for (d = 0; d < nDimensions; ++d) {
            dfDX[d] /= s;
        }
        
        return dfDX;
    }
    
    /**
     * evaluating
     * <pre>
     * 
     * </pre>
     * 
     * 
     @param isMedian array of length super.getNData() indicating whether an 
     * observation is equal to the given geomedian in all dimensions.
     @param geoMedian array of length getNDimensions() holding the estimate
     * of the geometric median.
     @param checks output array of size 2 holding the conditions to check after the algorithm
     * has completed.  checks[0 = (T(X)==X); checks[1]=(X.lte.eta(X))
     @return 
     */
    public double vardiZhang(int[] isMedian, double[] geoMedian, boolean[] checks) {
        
        if (geoMedian.length != nDimensions) {
            throw new IllegalArgumentException("geoMedian length should equal nDimensions");
        }
        if (checks.length != 2) {
            throw new IllegalArgumentException("checks.length must == 2");
        }
        
        int nData = (obs.length / nDimensions);
        
        /*
         let X be the geometric-median
         let n be the number of observations obs_i
         let M be the L1-median:
            M is a function of observed data obs_i and multiplicites eta_i
               where eta_i are used to make weights.
            M = argmin X of C(X)
               where C(X) = summation_i_1_to_n( eta_i*d_i(X) )
               where d_i is the euclidean distance between a point obs_i and X in all dimensions.
         ==> X=M iff T(X)=X iff r(X)<=eta(X)􏰌􏰅
               where eta(X) = eta_i if X=obs_i, i.e. geometric-median is a point in the set obs
                            else eta(X) = 0;
               where r(X) = ||R̃(X)||
               where R̃(X) = summation_over_obs_except_X( eta_i*(obs_i-X)/||obs_i-X|| )
               where T(X) = (1-(eta(X)/r(X))) * T̃(X) + X*math.min(1, (eta(X)/r(X)))
                 where 0/0 = 0 in the computation of eta(X)/r(X),
               where T̃(X) = summation_over_obs_except_X( eta_i*(obs_i)/||obs_i-X|| )
                             / summation_over_obs_except_X( eta_i/||obs_i-X|| )
                            (and was derived by setting derivatives to zero)
         NOTE: after iterative algorithm solves M, the test
               X=M iff T(X)=X iff r(X)<=eta(X) is performed.
        */
                     
        //geoMedian - obs_i for each dimension for each point
        double[] diffs = calculateDifferences(geoMedian);
        assert(diffs.length == obs.length);
        
        // length is nData.  sum of squared diffeences
        double[] ssdPerPoint = calculateSSDPerPoint(diffs);        
        assert(ssdPerPoint.length == nData);
        
        int i, j, d;
        double s;
        double[] t1Numer = new double[nDimensions];
        double[] t1Denom = new double[nDimensions];
        double[] etaMu = new double[nData];
        double[] rMu = new double[nDimensions];
        for (i = 0; i < nData; ++i) {
            if (isMedian[i] == 1) {
                //etaMu[i] = 1;
                etaMu[i] = this.eta[i];
                continue;
            }
            s = Math.sqrt(ssdPerPoint[i]);
            if (s < 1.e-17) {
                continue;
            }
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                t1Numer[d] += (this.eta[i]*obs[j]/s);
                t1Denom[d] += (this.eta[i]/s);
                rMu[d] += (this.eta[i]*diffs[j]/s);
            }
        }
        
        double rSum = 0;
        for (d = 0; d < nDimensions; ++d) {
            rSum += (rMu[d]*rMu[d]);
        }
        rSum = Math.sqrt(rSum);
        
        double[] gamma = new double[nDimensions];
        double[] t1 = new double[nDimensions];
        for (d = 0; d < nDimensions; ++d) {
            if (rSum > 1.e-17) {
               gamma[d] = Math.min(1., (etaMu[d]/rSum));
            }
            t1[d] = t1Numer[d]/t1Denom[d];
        }
        
        double[] geoMedian2 = new double[nDimensions];
        for (d = 0; d < nDimensions; ++d) {
            //geoMedian2[d] = ((1. - gamma[d])*t1[d]) + (gamma[d]*geoMedian[d]);
            geoMedian2[d] = ((1. - (etaMu[d]/rSum))*t1[d]) + (gamma[d]*geoMedian[d]);
        }
        
        // does T(X) == geoMedian?
        boolean c1 = true;
        for (d = 0; d < nDimensions; ++d) {
            if (Math.abs(geoMedian2[d] - geoMedian[d]) > eps) {
                c1 = false;
                break;
            }
        }
        
        // rSqSum <= etaMu[i]
        boolean c2 = true;
        for (d = 0; d < nDimensions; ++d) {
            if (rSum <= (etaMu[d] + eps)) {
                c2 = false;
                break;
            }
        }
                
        checks[0] = c1;
        checks[1] = c2;
                
        System.arraycopy(geoMedian2, 0, geoMedian, 0, geoMedian.length);
       
        return f(geoMedian);
    }

    /**
     *
     @return
     */
    @Override
    public int getNDimensions() {
        return nDimensions;
    }

    /**
     *
     @return
     */
    @Override
    public double[] getObs() {
        return obs;
    }
    
    /**
     *
     @return
     */
    @Override
    public double getFDEps() {
        return fDEps;
    }
}
