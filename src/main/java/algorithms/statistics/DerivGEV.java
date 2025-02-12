package algorithms.statistics;

import algorithms.misc.MiscMath0;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class DerivGEV {

    /*<pre>
     *                          (   (      ( x-mu))-(1/k))
     *                          (-1*(1 + k*(-----))      )
     *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
     * y = y_const * ----- * exp                           * (1 + k*(-----))
     *               sigma                                   (      (sigma))
     *
     * mu is  the location parameter
     * sigma is the scale parameter and is > 0
     * k is the shape parameter
     *
     * Let z = (1 + k*( (x-mu)/sigma )
     *
     *  y = y_const * (1/sigma) * exp(-1*(z^-(1/k))) * (z^(-1-(1/k)))
     *
     * Components needed in the derivatives:
     *
     *   then deriv of z w.r.t x is
     *   dzdx = k/sigma
     *
     *   dzdk = (x-mu)/sigma
     *
     *   dzdsigma =  -1 * k * (x-mu) * (sigma^-2)
     *
     *   dzdmu = -1*k/sigma

     *   deriv of -z^(-1/k) w.r.t. x is
     *      (1/k) * z^(-1 - (1/k)) * dzdx
     *
     *   deriv of -z^(-1/k) w.r.t. k
     *       use pattern: f(x) = u(x) ^(v(x))
                 ln ( f(x) ) = ln ( u(x) ^(v(x)) )

                 df(x)dx            du(x)dx
                 -------  =  v(x) * -------  +  dv(x)dx * ln(u(x))
                  f(x)               u(x)

                                           du(x)dx
                 df(x)dx = f(x) * ( v(x) * -------  +  dv(x)dx * ln(u(x)) )
                                            u(x)

                 u(k) = -z
                 v(k) = (-1/k)

     *           df(k)dk = -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     *
     *       *Because the ln( negative number) is not defined (though can be approximated with taylor series)
     *        need a different method of creating the derivative of a number with a power that is a function of x.
     *        This derivative is only needed with the derivative with respect to k.
     *
              One solution is to calculate the GEV with slightly different values of k near the given k to get delta GEV/deltaK.
     *
     *   deriv of -z^(-1/k) w.r.t. sigma is
     *       (1/k) * z^(-1 - (1/k)) * dzdsigma
     *
     *   deriv of -z^(-1/k) w.r.t. mu is
     *       (1/k) * z^(-1 - (1/k)) * dzdmu
     *
     *
     *   Let f1 = the first exponential in y
     *          = exp(-1*(z^(-1/k)))
     *
     *   df1dx     = f1 * (1/k) * z^(-1 - (1/k)) * dzdx
     *
     *
     *   df1dk     = f1 * -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     *
     *
     *   df1dsigma = f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
     *
     *
     *   df1dmu    = f1 * (1/k) * z^(-1 - (1/k)) * dzdmu
     *
     *
     *
     *   Let f2 = z^(-1-(1/k))
     *
     *   df2dx = (-1-(1/k)) * z^(-2-(1/k)) * dzdx
     *
     *
     *   df2dk: u(k) = z and v(k) = (-1-(1/k))
     *      using
                                        du(k)dk
                 df2dk = f2 * ( v(k) * -------   +   dv(k)dk * ln(u(k)) )
                                         u(k)
                       = f2 * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
     *
     *   df2dsigma = (-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
     *
     *   df2dmu = (-1-(1/k)) * z^(-2-(1/k)) * dzdmu
     *
     *
     *   Then putting it all together:
     *
     *          yconst
     *   yfit = ------ * f1 * f2
     *           sigma
     *
     *           yconst
     *   dydx =  ------ * ( f1 * df2dx + f2 * df1dx )
     *           sigma
     *
     *           yconst
     *   dydk =  ------ * ( f1 * df2dk + f2 * df1dk )
     *           sigma
     *
     *   dydsigma:
     *        needs to use chain rule once more
     *
     *        f0 = (yconst/sigma)
     *        df0dsigma = -(yconst/sigma^2)
     *
     *        f = f0 * f1 * f2
     *
     *        dydsigma = (  df0dsigma * f1 * f2 ) + ( df1dsigma * f0 * f2 ) + (df2dsigma * f0 * f1 )
     *
     *              yconst
     *   dydmu    = ------ * ( f1 * df2dmu + f2 * df1dmu )
     *              sigma
     * </pre>
     first implemented in project
     https://github.com/nking/two-point-correlation
         w/ Copyright Climb With Your Feet
         and using The MIT License (MIT)
       then moved to this shared library project which has the same copyright
       and license.
     */

    /**
     *
     */


    protected final static Logger log = Logger.getLogger(DerivGEV.class.getName());

    /**
     * calculate the derivative of the GEV w.r.t. x
     *
     * the runtime complexity is O(1), but uses 4 transcendental functions.
     *
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double derivWRTX(double yConst, double mu, double sigma, double k, double x) {

        double z = 1. + k *( (x-mu)/sigma );
        double a, f1, f2;

        double dzdx = k/sigma;
        
        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // For z^(g/h)
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            a = -1. * neg1Pow * (double) Math.pow(invNegZ, (-1./k));
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)) );
        } else {
            a = -1.*Math.pow(z, (-1./k));
            f2 = Math.pow(z, (-1. - (1./k)) );
        }

        f1 = Math.exp( a );
        
        double df2dx, df1dx;

        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            df2dx = neg1Pow * (-1. - (1./k)) * Math.pow(invNegZ, (-2. - (1./k)) ) * dzdx;
            df1dx = f1 * neg1Pow * (1./k) * Math.pow(invNegZ, (-1. - (1./k))) * dzdx;
        } else {
            df2dx = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdx;
            df1dx = f1 * (1./k) * Math.pow(z, (-1. - (1./k))) * dzdx;
        } 

        double dydx = (yConst/sigma) * ( f1 * df2dx + f2 * df1dx );

        return dydx;
    }

    /**
     * calculate the derivative of the GEV w.r.t. k
     *
     the runtime complexity is O(1), but uses 5 transcendental functions.     *
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double derivWRTK(double yConst, double mu, double sigma, double k, double x) {

        double z = 1. + k *( (x-mu)/sigma );
        double a, f1, f2;
        
        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // For z^(g/h)
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)) );
        } else {
            a = -1.*Math.pow(z, (-1./k));
            f2 = Math.pow(z, (-1. - (1./k)) );
        }

        f1 = Math.exp( a );

        // df1dk     = f1 * -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )

        // df2dk = f2 * ( (-1-(1/k)) * dzdk/z + (1/k^2) * ln(z) )

        // any value of z will have trouble with ln(-z) or ln(z) because the built in logarithm doesn't
        //    use a taylor series approximation for negative values plus handling for the number of cycles
        // so df1dk and df2dk are approximated with very small deltas

        double deltaK = 0.0001;

        double k_1 = k + deltaK;
        double z_1 = 1. + k_1 *( (x-mu)/sigma );
        
        boolean z_1IsNegative = (z < 0);
        double a_1, f2_1;
        
        if (z_1IsNegative) {
            double invNegZ = -1.0*(1.0/z_1);
            double neg1Pow = -1.0; // TODO:  revisit this
            a_1 = -1. * neg1Pow * Math.pow(invNegZ, (-1./k_1));
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)) );
            f2_1 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k_1)) );
        } else {
            a_1 = -1. * Math.pow(z_1, (-1./k_1));
            f2_1 = Math.pow(z_1, (-1. - (1./k_1)) );
        }

        f1 = Math.exp( a );

        double df1dk = (Math.exp( a_1 ) - f1)/deltaK;

        double df2dk = (f2_1 - f2)/deltaK;

        // to compare to the derivative:
        /*double dzdk = (x-mu)/sigma;
        if (zIsNegative) {
            double compare_df1dk = f1 * a * ( (1.f/k)*(dzdk/z) + (1.f/(k*k))*Math.log( z ) );
            log.info( String.format("  df1dk   estimate=%4.6f  deriv=%4.6f ", df1dk, compare_df1dk));
        } else {
            double compare_df2dk = f2 * (  (-1.f - (1.f/k))*(dzdk/z) + (1.f/(k*k))*Math.log(z) );
            log.info( String.format("  df2dk   estimate=%4.6f  deriv=%4.6f   (z=%4.6f, k=%4.6f)", df2dk, compare_df2dk, z, k));
        }*/

        double dydk = (yConst/sigma) * ( f1 * df2dk + f2 * df1dk );

        return dydk;
    }

    /**
     * calculate d/dk of GEV using the difference between GEVs given minor changes in k
     *
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    static double estimateDerivUsingDeltaK(double mu, double sigma, double k, double x) {

        double deltaK = 0.0001f*k;

        double d0 = GeneralizedExtremeValue.generateYGEV(x, mu, sigma, (k - deltaK));

        double d1 = GeneralizedExtremeValue.generateYGEV(x, mu, sigma, k);

        double d2 = GeneralizedExtremeValue.generateYGEV(x, mu, sigma, (k + deltaK));

        double d = estimateDerivUsingDelta(d0, d1, d2, deltaK);

        return d;
    }

    /**
     * calculate the derivative giving the 3 values which are separated by delta.
     * The first, d0, was computed with param - delta.
     * The second, d1, was computed with param.
     * The third, d2, was computed with param + delta.
     @param d0
     @param d1
     @param d2
     @param delta
     @return
     */
    protected static double estimateDerivUsingDelta(double d0, double d1, double d2, double delta) {

        double delta0 = d1 - d0;

        double delta1 = d2 - d1;

        double d = (delta0 + delta1)/2.;

        return (d/delta);
    }

    /**
     * calculate the derivative of the GEV w.r.t. sigma
     *
     the runtime complexity is O(1), but uses 5 transcendental functions.     *
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double derivWRTSigma(double yConst, double mu, double sigma, double k, double x) {

        double z = 1. + k *( (x-mu)/sigma );
        double a, f1, f2, df2dsigma, df1dsigma;

        //double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        //double f2 = Math.pow(z, (-1. - (1./k)) );

        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // 
        // For z^(g/h)
        //
        // For the continuous real exponentiation operator, negative base isn't allowed
        // For the discrete real exponentiation operator,
        //    fractional exponents with odd denominators are allowed.
        //    (-z)^(g/h) = ((-z)^g)^(1/h) = ((-1)^g)(z^(g/h))
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        double dzdsigma = -1. * k * (x-mu) * Math.pow(sigma, -2.);
        
        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
            f1 = Math.exp( a );
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)));
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
            df2dsigma = neg1Pow * (-1. - (1./k)) * Math.pow(invNegZ, (-2. - (1./k)) ) * dzdsigma;
            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dsigma =  f1 * neg1Pow * (1./k) * Math.pow(invNegZ, -1. - (1./k)) * dzdsigma;
        } else {
            a = -1. * Math.pow(z, (-1./k));
            f1 = Math.exp( a );
            f2 = Math.pow(z, (-1. - (1.f/k)));
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
            df2dsigma = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdsigma;
            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dsigma =  f1 * (1./k) * Math.pow(z, -1. - (1./k)) * dzdsigma;
        }
               
        double f0 = (yConst/sigma);

        double df0dsigma = -1.f*yConst/(sigma*sigma);

        double dydSigma = (  df0dsigma * f1 * f2 ) + ( df1dsigma * f0 * f2 ) + (df2dsigma * f0 * f1 );

        if (Double.isNaN(dydSigma)) {
            Double d = estimateDerivUsingDeltaSigma(mu, sigma, k, x);
            return (d != null) ? d : 0;
        }

        return dydSigma;
    }

    /**
     * calculate d/dsigma of GEV using the difference between GEVs given minor changes in k
     *
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    static double estimateDerivUsingDeltaSigma(double mu, double sigma, double k, double x) {

        double delta = 0.0001f*sigma;

        double d0 = GeneralizedExtremeValue.generateYGEV(x, mu, (sigma - delta), k);

        double d1 = GeneralizedExtremeValue.generateYGEV(x, mu, sigma, k);

        double d2 = GeneralizedExtremeValue.generateYGEV(x, mu, (sigma + delta), k);

        return estimateDerivUsingDelta(d0, d1, d2, delta);
    }

    /**
     * calculate the derivative of the GEV w.r.t. mu
     *
     the runtime complexity is O(1), but uses 5 transcendental functions.     *
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double derivWRTMu(double yConst, double mu, double sigma, double k, double x) {

        double z = 1. + k *( (x-mu)/sigma );

        //double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        //double f2 = Math.pow(z, (-1. - (1./k)) );

        double a, f1, f2, df1dmu, df2dmu;
        
        double dzdmu = -1. * k/sigma;
        
        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // For z^(g/h)
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
            f1 = Math.exp( a );
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)) );
            
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdmu
            df2dmu = neg1Pow * (-1. - (1./k)) * Math.pow(invNegZ, (-2. - (1./k)) ) * dzdmu;

            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dmu =  f1 * neg1Pow * (1./k) * Math.pow(invNegZ, -1. - (1./k)) * dzdmu;
        } else {
            a = -1.*Math.pow(z, (-1./k));
            f1 = Math.exp( a );
            f2 = Math.pow(z, (-1. - (1./k)) );
          
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdmu
            df2dmu = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdmu;

            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dmu =  f1 * (1./k) * Math.pow(z, -1. - (1./k)) * dzdmu;
        }

        double dydmu = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu );

        if (Double.isNaN(dydmu)) {
            Double d = estimateDerivUsingDeltaMu(mu, sigma, k, x);
            return (d != null) ? d : 0;
        }

        return dydmu;
    }

    /**
     * estimate d/dmu of GEV using the difference between GEVs given minor changes in k
     *
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    static double estimateDerivUsingDeltaMu(double mu, double sigma, double k, double x) {

        double deltaMu = 0.0001f*mu;

        double d0 = GeneralizedExtremeValue.generateYGEV(x, (mu - deltaMu), sigma, k);

        double d1 = GeneralizedExtremeValue.generateYGEV(x, mu, sigma, k);

        double d2 = GeneralizedExtremeValue.generateYGEV(x, (mu + deltaMu), sigma, k);

        return estimateDerivUsingDelta(d0, d1, d2, deltaMu);
    }

    /**
     * for given mu, k, sigma and the data x, normalizedY and normalizedYErr,
     * calculate the changes in k, sigma, and mu which would reduce the chi sq sum.
     *
     * Note that internally, the first and second derivatives of the curve GEV(k, sigma, mu)
     * are used to calculate what the smallest step in k, sigma, or mu would be in
     * order to create a significant change in the GEV curve.
     *
     * The suggested change for k is calculated from d/dk modified by the preconditioner
     * at the point right after the model peak.
     *
     * Note that the suggested changes might be applied by the NonQuadraticConguteSolver as
     * a fraction of the suggested change.
     *
     * runtime cost is
     *
     @param vars array of current values of mu, sigma, and k
     @param varsMin array of minimum allowed values of mu, sigma, and k
     @param varsMax array of maximum allowed values of mu, sigma, and k
     @param chiSqSumReturn item [0] is the best current chiSqSum for given vars.
     *    item [1] is to return the value of the best chiSqSum here to the invoker
     *    if r was set to non-zero values.
     @param x normalized x values of histogram to fit
     @param r the array to pass back values to be added to vars found by this method.
     *     the items are ordered so that
     *          index=0 holds delta mu, index=1 holds delta sigma,
     *     index=2 holds delta k.  the values may be zero if no step was found to
     *     improve the chi square sum.
     @param normalizedY
     @param normalizedYErr
     @param idx0 start of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     @param idx1 stop of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     */
    public static void derivsThatMinimizeChiSqSum(double[] vars, double[] varsMin, double[] varsMax, double[] chiSqSumReturn,
        double[] x, double[] normalizedY, double[] normalizedYErr, double[] r, int idx0, int idx1) {

        Arrays.fill(r, 0);

        // to determine the suggested step for k, sigma, or mu,
        //   we look at the first derivatives of the point in the model GEV right after the maximum
        //   and then modify those derivatives using a preconditioner from ICU0 matrix.
        // that suggested step is tried as plus and minus of current variable and the one
        //   with the smallest resulting chisqsum from the curve produced by the change is the
        //   result returned in derivs for that variable.

        for (int i = idx0; i < (idx1 + 1); i++) {

            double mu = vars[0] + r[0];
            double sigma = vars[1] + r[1];
            double k = vars[2] + r[2];

            double[] yGEV = GeneralizedExtremeValue.genCurve(x, vars[0], vars[1], vars[2]);
            int yMaxIdx = MiscMath0.findYMaxIndex(yGEV);
            if (yMaxIdx < 0) {
                yMaxIdx = 0;
            }
            double yMax = yGEV[yMaxIdx];
            for (int ii = 0; ii < yGEV.length; ii++){
                yGEV[ii] /= yMax;
            }
            double yConst = 1.f/yMax;

            int xIdx = yMaxIdx + 1;
            if (xIdx > x.length - 1) {
                xIdx = x.length - 1;
            }

            double xPoint = x[xIdx];

            double bestChiSqSum = chiSqSumReturn[0];

            double[] yGEVPlus = null;
            double[] yGEVMinus = null;
            double rModified = 0;

            switch(i) {
                // calculate a step size that would affect a change in GEV by using the 1st and 2nd partial derivatives
                case 0: {
                    // k
                    rModified = calculatePreconditionerModifiedResidualK(yConst, mu, sigma, k, xPoint);
                    break;
                }
                case 1: {
                    // sigma
                    rModified = calculatePreconditionerModifiedResidualSigma(yConst, mu, sigma, k, xPoint);
                    break;
                }
                case 2: {
                    // mu
                    rModified = calculatePreconditionerModifiedResidualMu(yConst, mu, sigma, k, xPoint);
                    break;
                }
            }

            switch(i) {
                case 0: {
                    // k
                    if (rModified != 0) {
                        if (((double)(k + rModified) >= varsMin[i]) && ((double)(k + rModified) <= varsMax[i])) {
                            yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, mu, sigma, (k + rModified));
                        }
                        if (((k - rModified) >= varsMin[i]) && ((k - rModified) <= varsMax[i])) {
                            yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, mu, sigma, (k - rModified));
                        }
                    }
                    break;
                }
                case 1: {
                    // sigma
                    if (rModified != 0) {
                        if (((sigma + rModified) >= varsMin[i]) && ((sigma + rModified) <= varsMax[i])) {
                            yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, mu, (sigma + rModified), k);
                        }
                        if (((sigma - rModified) >= varsMin[i]) && ((sigma - rModified) <= varsMax[i])) {
                            yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, mu, (sigma - rModified), k);
                        }
                    }
                    break;
                }
                case 2: {
                    // mu
                    if (rModified != 0) {
                        if (((mu + rModified) >= varsMin[i]) && ((mu + rModified) <= varsMax[i])) {
                            yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, (mu + rModified), sigma, k);
                        }
                        if (((mu - rModified) >= varsMin[i]) && ((mu - rModified) <= varsMax[i])) {
                            yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, (mu - rModified), sigma, k);
                        }
                    }
                    break;
                }
            }

            double chiSqSumPlus = (yGEVPlus != null) ? chiSqSum(yGEVPlus, normalizedY, normalizedYErr) : Double.POSITIVE_INFINITY;
            double chiSqSumMinus = (yGEVMinus != null) ? chiSqSum(yGEVMinus, normalizedY, normalizedYErr) : Double.POSITIVE_INFINITY;
            boolean bestIsPlus = (chiSqSumPlus < chiSqSumMinus);

            if (bestIsPlus && (chiSqSumPlus < bestChiSqSum)) {
                r[i] = (double) rModified;
                bestChiSqSum = chiSqSumPlus;
                chiSqSumReturn[1] = bestChiSqSum;
            } else if (chiSqSumMinus < bestChiSqSum) {
                r[i] = (double) (-1.f * rModified);
                bestChiSqSum = chiSqSumMinus;
                chiSqSumReturn[1] = bestChiSqSum;
            }

            log.finest(" derivs[" + i + "]=" + r[i]);
        }
    }

    /**
      calculate Preconditioner Modified Residual with respect to k as
     resid = dydk/d2ydkdk;
     @param yConst constant factor used in dydk = (yConst/sigma) * ( f1 * df2dk + f2 * df1dk )
     @param mu location parameter of GEV distribution
     @param sigma scale parameter of GEV distribution
     @param k shape parameter of GEV distribution
     @param x the quantile of the distribution
     @return precondition modified residual with respect to k: dydk/d2ydkdk
     */
    public static double calculatePreconditionerModifiedResidualK(double yConst, double mu, double sigma, double k, double x) {

        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        //
        // k component to residuals = d(1,1) * (&#8706;/&#8706;k)
        //       where d(1,1) is 1./(&#8706;<sup>2</sup>/&#8706;k&#8706;k), that is 1./(&#8706;<sup>2</sup>/&#8706;k&#8706;k)

        double dydk = DerivGEV.derivWRTK(yConst, mu, sigma, k, x);

        Double d2ydkdk = estimateDY2DKDK(yConst, mu, sigma, k, x, dydk);

        double resid = dydk/d2ydkdk;

        return resid;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;k&#8706;k empirically.
     *
     * The method uses the tested DerivGEV.derivWRTK() for dfdk and plugs in different k's
     * to estimate d2fdkdk
     *
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double estimateDY2DKDK(double yConst, double mu, double sigma, double k, double x) {

        double dydk = DerivGEV.derivWRTK(yConst, mu, sigma, k, x);

        double d = estimateDY2DKDK(yConst, mu, sigma, k, x, dydk);

        return d;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k empirically.  
     * the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydk
     @return
     */
    public static double estimateDY2DKDK(double yConst, double mu, double sigma, double k, double x, double dydk) {

        double factor = 0.0001f;

        double delta = k*factor;

        double dydk_0 = DerivGEV.derivWRTK(yConst, mu, sigma, (k - delta), x);

        double dydk_2 = DerivGEV.derivWRTK(yConst, mu, sigma, (k + delta), x);

        double d = estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);

        return d;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;k&#8706;sigma, &#8706;<sup>2</sup>/&#8706;k&#8706;sigma empirically.  
     * the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydk
     @return
     */
    public static double estimateDY2DKDSigma(double yConst, double mu, double sigma, double k, double x, double dydk) {

        // &#8706;<sup>2</sup>/&#8706;k&#8706;sigma

        double factor = 0.0001f;

        double delta = (sigma*factor);

        double dydk_0 = DerivGEV.derivWRTK(yConst, mu, (sigma - delta), k, x);

        double dydk_2 = DerivGEV.derivWRTK(yConst, mu, (sigma + delta), k, x);

        double d = estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);

        return d;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;k&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;k&#8706;mu empirically.  
     * the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydk
     @return
     */
    public static double estimateDY2DKDMu(double yConst, double mu, double sigma, double k, double x, double dydk) {

        // &#8706;<sup>2</sup>/&#8706;k&#8706;mu

        double factor = 0.0001f;

        double delta = mu*factor;

        double dydk_0 = DerivGEV.derivWRTK(yConst, (mu - delta), sigma, k, x);

        double dydk_2 = DerivGEV.derivWRTK(yConst, (mu + delta), sigma, k, x);

        double d = estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);

        return d;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma, that is &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydsigma
     @return
     */
    public static double estimateDY2DSigmaDSigma(double yConst, double mu, double sigma, double k, double x, double dydsigma) {

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma = estimate as (dyds_2 - dyds)/ds

        double factor = 0.0001f;

        double delta = sigma*factor;

        double dyds_0 = DerivGEV.derivWRTSigma(yConst, mu, (sigma - delta), k, x);

        double dyds_2 = DerivGEV.derivWRTSigma(yConst, mu, (sigma + delta), k, x);

        double d = estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);

        return d;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;mu&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;mu&#8706;mu empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydmu
     @return
     */
    public static double estimateDY2DMuDMu(double yConst, double mu, double sigma, double k, double x, double dydmu) {

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;mu

        double factor = 0.0001f;

        double delta = mu*factor;

        double dydm_0 = DerivGEV.derivWRTMu(yConst, (mu - delta), sigma, k, x);

        double dydm_2 = DerivGEV.derivWRTMu(yConst, (mu + delta), sigma, k, x);

        double d = estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);

        return d;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;mu&#8706;k, that is &#8706;<sup>2</sup>/&#8706;mu&#8706;k empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydmu
     @return
     */
    public static double estimateDY2DMuDK(double yConst, double mu, double sigma, double k, double x, double dydmu) {

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;k

        double factor = 0.0001f;

        double delta = k*factor;

        double dydm_0 = DerivGEV.derivWRTMu(yConst, mu, sigma, (k - delta), x);

        double dydm_2 = DerivGEV.derivWRTMu(yConst, mu, sigma, (k + delta), x);

        double d = estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);

        return d;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma, that is &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydmu
     @return
     */
    public static double estimateDY2DMuDSigma(double yConst, double mu, double sigma, double k, double x, double dydmu) {

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma

        double factor = 0.0001f;

        double delta = sigma*factor;

        double dydm_0 = DerivGEV.derivWRTMu(yConst, mu, (sigma - delta), k, x);

        double dydm_2 = DerivGEV.derivWRTMu(yConst, mu, (sigma + delta), k, x);

        double d = estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);

        return d;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;sigma&#8706;k, that is &#8706;<sup>2</sup>/&#8706;sigma&#8706;k empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydsigma
     @return
     */
    public static double estimateDY2DSigmaDK(double yConst, double mu, double sigma, double k, double x, double dydsigma) {

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;k = estimate as (dy2_ds_dk - dyds)/dk

        double factor = 0.0001f;

        double delta = (k*factor);

        double dydk_0 = DerivGEV.derivWRTSigma(yConst, mu, sigma, (k - delta), x);

        double dydk_2 = DerivGEV.derivWRTSigma(yConst, mu, sigma, (k + delta), x);

        double d = estimateDerivUsingDelta(dydk_0, dydsigma, dydk_2, delta);

        return d;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @param dydsigma
     @return
     */
    public static double estimateDY2DSigmaDMu(double yConst, double mu, double sigma, double k, double x, double dydsigma) {

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu = estimate as (dy2_ds_dk - dyds)/dk

        double factor = 0.0001f;

        double delta = (sigma*factor);

        double dyds_0 = DerivGEV.derivWRTSigma(yConst, (mu - delta), sigma, k, x);

        double dyds_2 = DerivGEV.derivWRTSigma(yConst, (mu + delta), sigma, k, x);

        double d = estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);

        return d;
    }

    /**
     calculate Preconditioner Modified Residual with respect to sigma
     @param yConst constant factor used in dydk = (yConst/sigma) * ( f1 * df2dsigma + f2 * df1dsigma )
     @param mu location parameter of GEV distribution
     @param sigma scale parameter of GEV distribution
     @param k shape parameter of GEV distribution
     @param x the quantile of the distribution
     @return precondition modified residual with respect to sigma
     */
    public static double calculatePreconditionerModifiedResidualSigma(
        double yConst, double mu, double sigma, double k, double x) {

        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        //
        // sigma component to residuals = d(2,2) * (&#8706;/&#8706;sigma)
        //       where d(2,2) is ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )

        // &#8706;/&#8706;sigma
        double dyds = DerivGEV.derivWRTSigma(yConst, mu, sigma, k, x);

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma
        double d2ydsds = estimateDY2DSigmaDSigma(yConst, mu, sigma, k, x, dyds);


        // (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k) = estimate as (dydsigma_1 - dydsigma)/dsigma
        double d2ydsdk = estimateDY2DSigmaDK(yConst, mu, sigma, k, x, dyds);

        // &#8706;/&#8706;k
        double dydk = DerivGEV.derivWRTK(yConst, mu, sigma, k, x);

        // &#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k
        double d2ydkdk = estimateDY2DKDK(yConst, mu, sigma, k, x, dydk);

        // &#8706;<sup>2</sup>/&#8706;k&#8706;sigma
        double d2ydkds = estimateDY2DKDSigma(yConst, mu, sigma, k, x, dydk);

        double modification = d2ydsds - (d2ydsdk / d2ydkdk) * d2ydkds;

        double resid = dyds / modification;

        return resid;
    }

    /**
     calculate Preconditioner Modified Residual with respect to mu
     @param yConst constant factor used in dydk = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu )
     @param mu location parameter of GEV distribution
     @param sigma scale parameter of GEV distribution
     @param k shape parameter of GEV distribution
     @param x the quantile of the distribution
     @return precondition modified residual with respect to mu
     */
    public static double calculatePreconditionerModifiedResidualMu(double yConst, double mu, double sigma, double k, double x) {

        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        //
        // sigma component to residuals = d(3,3) * (&#8706;/&#8706;mu)
        /*
           where  d(3,3) is 1./(
                     ( (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu) )
                     -
                     (
                        &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma
                        *
                        ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                        *
                        &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu
                     )
                 )

                Let pt1 = (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu)

                Let pt2_1 = ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )

                Let pt2 = ( &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma * pt2_1 * &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu )

            d(3,3) = 1./( pt1 - pt2)
        */

        // &#8706;/&#8706;mu
        double dydmu = DerivGEV.derivWRTMu(yConst, mu, sigma, k, x);

        // &#8706;/&#8706;sigma
        double dyds = DerivGEV.derivWRTSigma(yConst, mu, sigma, k, x);

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma
        double d2ydsds = estimateDY2DSigmaDSigma(yConst, mu, sigma, k, x, dyds);

        if (d2ydsds == 0) {
            return 0;
        }


        // &#8706;<sup>2</sup>/&#8706;mu&#8706;mu
        double d2ydmdm = estimateDY2DMuDMu(yConst, mu, sigma, k, x, dydmu);

        // &#8706;/&#8706;k
        double dydk = DerivGEV.derivWRTK(yConst, mu, sigma, k, x);

        double pt1;

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;k
        double d2ydmdk = estimateDY2DMuDK(yConst, mu, sigma, k, x, dydmu);

        // &#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k
        double d2ydkdk = estimateDY2DKDK(yConst, mu, sigma, k, x, dydk);

        // &#8706;<sup>2</sup>/&#8706;k&#8706;mu
        double d2ydkdm = estimateDY2DKDMu(yConst, mu, sigma, k, x, dydk);

        pt1 = d2ydmdm - (d2ydmdk*d2ydkdm / d2ydkdk);

        // (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k) = estimate as (dydsigma_1 - dydsigma)/dsigma
        double d2ydsdk = estimateDY2DSigmaDK(yConst, mu, sigma, k, x, dyds);

        double d2ydmds = estimateDY2DMuDSigma(yConst, mu, sigma, k, x, dydmu);

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu
        double d2ydsdm = estimateDY2DSigmaDMu(yConst, mu, sigma, k, x, dyds);

        // d(3,3) = 1./( (pt1) - ( &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma * pt2_1 * &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu ))
        //
        // pt2_1 = ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
        // Let pt2 = ( &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma * pt2_1 * &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu )

        double pt2, pt2_1;

        if (d2ydmds == 0) {
            pt2 = 0;
            pt2_1 = 0;
        } else {

            // &#8706;<sup>2</sup>/&#8706;k&#8706;sigma
            double d2ydkds = estimateDY2DKDSigma(yConst, mu, sigma, k, x, dydk);

            pt2_1 = d2ydsds - (d2ydsdk * d2ydkds / d2ydkdk);

            pt2_1 = 1./pt2_1;

            pt2 = d2ydmds * pt2_1 * d2ydsdm;
        }

        double modification = pt1 - pt2;

        double resid = dydmu/modification;

        return resid;
    }

    /**
     *
     @param mu
     @param sigma
     @param k
     @param x
     @param normalizedY
     @param normalizedYErr
     @return
     */
    public static Double chiSqSum(final double mu, final double sigma, final double k, final double[] x,
        final double[] normalizedY, final double[] normalizedYErr) {

        double[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, mu, sigma, k);
        if (yGEV == null) {
            return null;
        }
        return chiSqSum(yGEV, normalizedY, normalizedYErr);
    }

    /**
     * compute the chi square sum for the normalized curves, that is, do not include the factors for
     * yScale of the original unnormalized data.
     *
     @param normalizedYGEV
     @param normalizedY
     @param normalizedYErr
     @return
     */
    public static double chiSqSum(double[] normalizedYGEV, double[] normalizedY, double[] normalizedYErr) {
        return chiSqSum(normalizedYGEV, normalizedY, normalizedYErr, 0, normalizedY.length);
    }


    /**
     * compute the chi square sum for the normalized curves for only part of the indexes of the curve.
     * The calc does not include the factors for
     * yScale of the original unnormalized data.
     *
     @param normalizedYGEV
     @param normalizedY
     @param normalizedYErr
     @param startIdx
     @param stopIdx the last index, exclusive
     @return
     */
    public static double chiSqSum(double[] normalizedYGEV, double[] normalizedY, double[] normalizedYErr,
        int startIdx, int stopIdx) {
        double chiSqSum = 0.0f;
        for (int ii = startIdx; ii < stopIdx; ii++) {
            double z = (normalizedYGEV[ii] - normalizedY[ii]);
            chiSqSum += z*z*normalizedYErr[ii];
        }
        return chiSqSum;
    }

    /**
     * 
     * retrieve the second deriv of GEV &#8706;<sup>2</sup>/&#8706;k&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;k&#8706;mu
     *
     * Note:  the method hasn't been tested yet and will be compared to the results from method estimateDY2DKDMu.
     * 
     * 
     @param yConst
     @param mu
     @param k
     @param sigma
     @param x
     @return
     */
    public static double secondDerivKDerivMu(double yConst, double mu, double sigma, double k, double x) {

        double z = 1. + k *( (x-mu)/sigma );

        //double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        //double f2 = Math.pow(z, (-1. - (1./k)) );

        double a, f1, f2, df1dmu, df2dmu;
        
        double dzdmu = -1. * k/sigma;
        
        boolean zIsNegative = (z < 0);
        // When z is negative, need to use alternative methods for exponentiation:
        // For z^(g/h)
        // For the complex exponentiation operator, that is complex bases, 
        //    the results are not continuous and are infinite.
        //    the principal value is
        //        (-1)^(g/h) * ((1./z)^(g/h))

        if (zIsNegative) {
            double invNegZ = -1.0*(1.0/z);
            double neg1Pow = -1.0; // TODO:  revisit this
            a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
            f1 = Math.exp( a );
            f2 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)) );
            
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdmu
            df2dmu = neg1Pow * (-1. - (1./k)) * Math.pow(invNegZ, (-2. - (1./k)) ) * dzdmu;

            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dmu =  f1 * neg1Pow * (1./k) * Math.pow(invNegZ, -1. - (1./k)) * dzdmu;
        } else {
            a = -1.*Math.pow(z, (-1./k));
            f1 = Math.exp( a );
            f2 = Math.pow(z, (-1. - (1./k)) );
          
            //(-1-(1/k)) * z^(-2-(1/k)) * dzdmu
            df2dmu = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdmu;

            //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
            df1dmu =  f1 * (1./k) * Math.pow(z, -1. - (1./k)) * dzdmu;
        }
        
        double deltaK = 0.0001f;
        double deltaMu = 0.0001f;
        double k_1 = k + deltaK;
        double mu_1 = mu + deltaMu;
        double z_1 = 1. + k_1 *( (x-mu)/sigma );
        double z_1_1 = 1. + k_1 *( (x-mu_1)/sigma );

        double a_1, f2_1;
        if (z_1 < 0) {
            double invNegZ = -1.0*(1.0/z_1);
            double neg1Pow = -1.0; // TODO:  revisit this
            a_1 = -1. * neg1Pow * Math.pow(invNegZ, (-1./k_1));
            f2_1 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k_1)) );
        } else {
            a_1 = -1. * Math.pow(z_1, (-1./k_1));
            f2_1 = Math.pow(z_1, (-1. - (1./k_1)) );
        }
        
        double a_1_1, f2_1_1;
        if (z_1_1 < 0) {
            double invNegZ = -1.0*(1.0/z_1_1);
            double neg1Pow = -1.0; // TODO:  revisit this
            a_1_1 = -1. * neg1Pow * Math.pow(invNegZ, (-1.f/k_1));
            f2_1_1 = neg1Pow * Math.pow(invNegZ, (-1. - (1./k_1)) );
        } else {
            a_1_1 = -1. * Math.pow(z_1_1, (-1.f/k_1));
            f2_1_1 = Math.pow(z_1_1, (-1. - (1./k_1)) );
        }
        
        double f1_1 = Math.exp( a_1 );
        double df1dk = (Math.exp( a_1 ) - f1)/deltaK;
        double df1dk_1 = (Math.exp( a_1_1 ) - f1_1)/deltaK;
        double df2dk = (f2_1 - f2)/deltaK;
        double df2dk_1 = (f2_1_1 - f2_1)/deltaK;


        //pt1_0 = ∂/∂mu( df2dk )  the derivative involves potential log(negative number) so work around:
        double pt1_0 = (df2dk_1 - df2dk)/deltaMu;

        // pt2_0 = ∂/∂mu( df1dk ) the derivative involves potential log(negative number) so work around:
        double pt2_0 = (df1dk_1 - df1dk)/deltaMu;


        double pt2 = f2 * pt2_0 + df1dk * df2dmu;

        double pt1 = (f1 * pt1_0) + (df2dk * df1dmu);

        double d2fdkdmu = (yConst/sigma) * (pt1 + pt2);

        return d2fdkdmu;
    }

    /*
     
     Can see that calculating the formula for the partial derivative &#8706;<sup>2</sup>/&#8706;k&#8706;mu
     is time consuming and error prone so will only implement the "estimate" methods hereafter.


     df1dk     = f1 * (-1*z^(-1/k)) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     df2dk     = f2 * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
     dzdk = (x-mu)/sigma

     &#8706;<sup>2</sup>/&#8706;k&#8706;mu
         = ∂/∂mu ( (yconst/sigma) * ( f1 * df2dk + f2 * df1dk )  )
         = (yconst/sigma) * ∂/∂mu( f1 * df2dk )  + (yconst/sigma) * ∂/∂mu(f2 * df1dk )
         = (yconst/sigma) *  pt1                 + (yconst/sigma) * pt2

         pt1 = ∂/∂mu( f1 * df2dk )
             = f1 * ∂/∂mu( df2dk ) + df2dk * ∂/∂mu(f1)
             = f1 * pt1_0          + df2dk * df1dmu

             pt1_0 = ∂/∂mu( df2dk )
                   = ∂/∂mu( f2 * ( (-1-(1/k)) * dzdk * (1/z)  +  (1/k^2) * ln(z) ) )

                   for a product of more than 2 factors, multiply deriv of one times all
                   other factors then add next deriv times all other factors...

                   = ∂/∂mu( f2 * (-1-(1/k)) * dzdk * (1/z) )  +  ∂/∂mu( f2 * (1/k^2) * ln(z) )

                   = pt1_0_0                                  + pt1_0_1

                 pt1_0_0 = ∂/∂mu( f2 ) * (-1-(1/k)) * dzdk * (1/z)  +  0
                             + ∂/∂mu(dzdk) * f2 * (-1-(1/k)) * (1/z)
                             + ∂/∂mu((1/z)) * f2 * (-1-(1/k)) * dzdk
                         = df2dmu * (-1-(1/k)) * dzdk * (1/z)
                             + pt1_0_0_0 * f2 * (-1-(1/k)) * (1/z)
                             + pt1_0_0_1 * f2 * (-1-(1/k)) * dzdk

                     pt1_0_0_0 = ∂/∂mu(dzdk)
                               = ∂/∂mu( (x-mu)/sigma )
                               = -mu/sigma

                     pt1_0_0_1 = ∂/∂mu(1/z)

                                    z = (1 + k*( (x-mu)/sigma ) = (sigma + k*(x-mu))/sigma
                                    1/z = sigma/(sigma + k*(x-mu))

                               = ∂/∂mu(sigma/(sigma + k*(x-mu)))
                               = ∂/∂mu( sigma * (sigma + k*(x-mu))^-1 )
                               = -1 * sigma * (sigma + k*(x-mu))^-2  * (-k)
                               = k * sigma * (sigma + k*(x-mu))^-2

                 pt1_0_1 = ∂/∂mu( f2 * (1/k^2) * ln(z) )
                         = df2dmu * (1/k^2) * ln(z)
                            + 0
                            + ∂/∂mu( ln(z) ) * f2 * (1/k^2)
                         = df2dmu * (1/k^2) * ln(z)
                            + pt1_0_1_0 * f2 * (1/k^2)

                     pt1_0_1_0 = ∂/∂mu( ln(z) )   where logarithm( u(mu) ) = (1/u(mu)) * d*u(mu)/dmu
                               = ∂/∂mu( ln ( (sigma + k*(x-mu))/sigma )
                               = ( sigma / (sigma + k*(x-mu)) ) * (-k)
                               = (-k * sigma) / ( k * (sigma + k*(x-mu)) )

         pt2 = ∂/∂mu(f2 * df1dk )
             = f2 * ∂/∂mu( df1dk ) + df1dk * df2dmu
             = f2 * pt2_0 + df1dk * df2dmu

             pt2_0 = ∂/∂mu( df1dk )
                   = ∂/∂mu( f1 * (-1*z^(-1/k)) * ( -1*(-1/k) * dzdk * (1/z)  +  (1/k^2) * ln( -z ) ) )
                   = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z) )  +  ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k^2) * ln( -z ) )
                   = pt2_0_0   +   pt2_0_1

                 pt2_0_0 = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z) )
                         = df1dmu * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z)
                             + ∂/∂mu(-1*z^(-1/k)) * f1 * (1/k) * dzdk * (1/z)
                             + 0
                             + ∂/∂mu(dzdk) * f1 * (-1*z^(-1/k)) * (1/k) * (1/z)
                             + ∂/∂mu(1/z) * f1 * (-1*z^(-1/k)) * (1/k) * dzdk
                         = df1dmu * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z)
                             + pt2_0_0_0 * f1 * (1/k) * dzdk * (1/z)
                             + pt2_0_0_1 * f1 * (-1*z^(-1/k)) * (1/k) * (1/z)
                             + pt2_0_0_2 * f1 * (-1*z^(-1/k)) * (1/k) * dzdk

                     pt2_0_0_0 = ∂/∂mu(-1*z^(-1/k))
                               = (1/k) * z^(-1 - (1/k)) * dzdmu  <== from class level comments

                     pt2_0_0_1 = ∂/∂mu(dzdk)
                               = ∂/∂mu( (x-mu)/sigma )
                               = -1/sigm

                     pt2_0_0_2 = ∂/∂mu(1/z)
                               = pt1_0_0_1

                 pt2_0_1 = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k^2) * ln( -z ) )
                         = df1dmu * (-1*z^(-1/k)) * (1/k^2) * ln( -z )
                            + ∂/∂mu(-1*z^(-1/k)) * f1 * (1/k^2) * ln( -z )
                            + 0
                            + ∂/∂mu( ln( -z ) ) * f1 * (-1*z^(-1/k)) * (1/k^2)
                         = df1dmu * (-1*z^(-1/k)) * (1/k^2) * ln( -z )
                            + pt2_0_0_0 * f1 * (1/k^2) * ln( -z )
                            + 0
                            + pt2_0_1_0 * f1 * (-1*z^(-1/k)) * (1/k^2)




                     pt2_0_1_0 = ∂/∂mu( ln ( -1*(sigma + k*(x-mu))/sigma )
                               = (k*sigma)/(sigma + k*(x-mu))


*/
}
