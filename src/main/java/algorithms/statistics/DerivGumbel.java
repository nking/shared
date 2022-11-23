package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;

import java.util.logging.Logger;

public class DerivGumbel {

    /*<pre>

      y = (yconst/sigma) * exp( -(x-mu)/sigma ) * exp(-exp( -(x-mu)/sigma ))

     * mu is  the location parameter
     * sigma is the scale parameter and is > 0
     * k is the shape parameter = 0

     Let z = (x-mu)/sigma
           = (x-mu)/sigma = x/sigma - mu/sigma

     y = (yconst/sigma) * exp(-z) * exp(-exp(-z))

   Components needed in the derivatives:

   then derivs of z  w.r.t. each variable are:
       dz/dx = 1/sigma

       dz/dsigma = -(x-mu)/sigma^2

       dz/dmu = -1/sigma

     y = (1/sigma) * exp(-z) * exp(-exp(-z))
       = (1/sigma) * f1 * exp(-f1)
       = (1/sigma) * f1 * f2

     *   Let f1 = the first exponential in y
     *          = exp(-z)

     *   df1/dx     = f1 * (-1) * dz/dx
     *
     *   df1/dsigma = f1 * (-1) * dz/dsigma
     *
     *   df1/dmu    = f1 * (-1) * dz/dmu
     *
         Let f2 = the second exponential in y
     *          = exp(-f1)

     *   df2/dx = f2 * (-1) * df1/dx

         df2/dsigma = f2 * (-1) * df1/dsigma

         df2/dmu = f2 * (-1) * df1/dmu

     *   Then putting it all together:
     *
     *           yconst
     *   yfit = ------ * f1 * f2
     *           sigma
     *
     *            yconst
     *   dy/dx = ------ * ( f1 * df2/dx + f2 * df1/dx )
     *            sigma
     *
     *             yconst
     *   dy/dmu =  ------ * ( f1 * df2/dmu + f2 * df1/dmu )
     *             sigma
     *
     *   dydsigma:
     *        needs to use chain rule once more
     *
     *        f0 = (yconst/sigma)
     *        df0/dsigma = -(yconst/sigma^2)
     *
     *        f = f0 * f1 * f2
     *
     *        dy/dsigma = (  df0/dsigma * f1 * f2 ) + ( df1/dsigma * f0 * f2 ) + (df2/dsigma * f0 * f1 )
     *
     * </pre>
     */

    protected static Logger log = Logger.getLogger(DerivGumbel.class.getName());

    protected static double calcDZDX(double sigma) {
        return 1./sigma;
    }
    protected static double calcDZDSigma(double mu, double sigma, double x) {
        return -(x-mu)/(sigma*sigma);
    }
    protected static double calcDZDMu(double sigma) {
        return -1./sigma;
    }
    protected static double calcZ(double mu, double sigma, double x) {
        return (x-mu)/sigma;
    }

    /**
     * f1 = exp(-z)
     * @param mu
     * @param sigma
     * @param x
     * @return
     */
    protected static double calcF1(double mu, double sigma, double x) {
        return Math.exp(-1. * calcZ(mu, sigma, x));
    }
    /**
     * f1 = exp(-z)
     * @param z
     * @return
     */
    protected static double calcF1(double z) {
        return Math.exp(-z);
    }
    protected static double calcDF1DX(double mu, double sigma, double x) {
        //f1 * (-1) * dz/dx
        double f1 = calcF1(mu, sigma, x);
        return calcDF1DX(f1, sigma);
    }
    protected static double calcDF1DX(double f1, double sigma) {
        //f1 * (-1) * dz/dx
        double dzdx = calcDZDX(sigma);
        return -1. * f1 * dzdx;
    }
    protected static double calcDF1DSigma(double mu, double sigma, double x) {
        //f1 * (-1) * dz/dsigma
        double f1 = calcF1(mu, sigma, x);
        return calcDF1DSigma(f1, mu, sigma, x);
    }
    protected static double calcDF1DSigma(double f1, double mu, double sigma, double x) {
        //f1 * (-1) * dz/dsigma
        double dzdsigma = calcDZDSigma(mu, sigma, x);
        return -1. * f1 * dzdsigma;
    }
    protected static double calcDF1DMu(double mu, double sigma, double x) {
        //f1 * (-1) * dz/dmu
        double f1 = calcF1(mu, sigma, x);
        return calcDF1DMu(f1, sigma);
    }
    protected static double calcDF1DMu(double f1, double sigma) {
        //f1 * (-1) * dz/dmu
        double dzdmu = calcDZDMu(sigma);
        return -1. * f1 * dzdmu;
    }

    protected static double calcF2(double f1) {
        //exp(-f1)
        return Math.exp(-f1);
    }
    protected static double calcDF2DX(double f1, double f2, double sigma) {
        //f2 * (-1) * df1/dx
        double df1dx = calcDF1DX(f1, sigma);
        return -f2 * df1dx;
    }
    protected static double calcDF2DSigma(double f1, double f2, double mu, double sigma, double x) {
        //f2 * (-1) * df1/dsigma
        double df1dsigma = calcDF1DSigma(f1, mu, sigma, x);
        return -f2 * df1dsigma;
    }
    protected static double calcDF2DMu(double f1, double f2, double sigma) {
        //f2 * (-1) * df1/dmu
        double df1dmu = calcDF1DMu(f1, sigma);
        return -f2 * df1dmu;
    }

    public static double derivWRTX(double yConst, double z, double f1, double f2, double sigma) {

        /*
         *            yconst
         *   dy/dx = ------ * ( f1 * df2/dx + f2 * df1/dx )
         *            sigma
         */
        double df2dx = calcDF2DX(f1, f2, sigma);
        double df1dx = calcDF1DX(f1, sigma);

        double dfdx = (yConst/sigma) * (f1 * df2dx + f2 * df1dx);

        return dfdx;
    }

    public static double derivWRTMu(double yConst, double mu, double sigma, double x) {

        /*
         *             yConst
         *   dy/dmu =  ------ * ( f1 * df2/dmu + f2 * df1/dmu )
         *             sigma
         */
        double z = calcZ(mu, sigma, x);
        double f1 = calcF1(z);
        double f2 = calcF2(f1);

        return derivWRTMu(z, f1, f2, yConst, sigma);
    }
    public static double derivWRTMu(double z, double f1, double f2, double yConst, double sigma) {

        /*
         *             yConst
         *   dy/dmu =  ------ * ( f1 * df2/dmu + f2 * df1/dmu )
         *             sigma
         */
        double df2dmu = calcDF2DMu(f1, f2, sigma);
        double df1dmu = calcDF1DMu(f1, sigma);

        double dfdmu = (yConst/sigma) * (f1 * df2dmu + f2 * df1dmu);

        return dfdmu;
    }

    public static double derivWRTSigma(double yConst, double mu, double sigma, double x) {

        /*
         f0 = (yconst/sigma)
         df0/dsigma = -(yconst/sigma^2)

         f = f0 * f1 * f2

        dy/dsigma = (  df0/dsigma * f1 * f2 ) + ( df1/dsigma * f0 * f2 ) + (df2/dsigma * f0 * f1 )
        */
        double z = calcZ(mu, sigma, x);
        double f1 = calcF1(z);
        double f2 = calcF2(f1);

        return derivWRTSigma(z, f1, f2, yConst, mu, sigma, x);
    }

    public static double derivWRTSigma(double z, double f1, double f2, double yConst, double mu, double sigma, double x) {

        /*
         f0 = (yconst/sigma)
         df0/dsigma = -(yconst/sigma^2)

         f = f0 * f1 * f2

        dy/dsigma = (  df0/dsigma * f1 * f2 ) + ( df1/dsigma * f0 * f2 ) + (df2/dsigma * f0 * f1 )
        */
        double f0 = yConst/sigma;
        double df0dsigma = -f0/sigma;
        double df2dsigma = calcDF2DSigma(f1, f2, mu, sigma, x);
        double df1dsigma = calcDF1DSigma(f1, mu, sigma, x);

        double dfdsigma = (df0dsigma * f1 * f2) + (df1dsigma * f0 * f2) + (df2dsigma * f0 * f1);

        return dfdsigma;
    }

    public static double calcYFit(double yConst, double f1, double f2, double sigma) {
        /*
         *           yconst
         *   yfit = ------ * f1 * f2
         *           sigma
         */
        return (yConst/sigma) * f1 * f2;
    }

    //===============

    /**
     * calculate the derivative giving the 3 values which are separated by delta.
     * The first, d0, was computed with param - delta.
     * The second, d1, was computed with param.
     * The third, d2, was computed with param + delta.
     * @param d0
     * @param d1
     * @param d2
     * @param delta
     * @return
     */
    protected static Double estimateDerivUsingDelta(Double d0, Double d1, Double d2, double delta) {

        if (d0 != null && d1 != null && d2 != null) {

            double delta0 = d1.doubleValue() - d0.doubleValue();
            double delta1 = d2.doubleValue() - d1.doubleValue();

            double d = (delta0 + delta1)/2.;

            return (d/delta);

        } else if (d1 != null && d2 != null) {

            double d = d2.doubleValue() - d1.doubleValue();

            return (d/delta);

        } else if (d0 != null && d1 != null) {

            double d = d1.doubleValue() - d0.doubleValue();

            return (d/delta);

        } else if (d0 != null && d2 != null) {

            double d = d2.doubleValue() - d0.doubleValue();

            return (d/delta);

        } else {

            return null;
        }
    }

    /**
     * calculate d/dsigma of Gumbel using the finite difference method
     *
     * @param mu
     * @param sigma
     * @param x
     * @return
     */
    static Double estimateDerivUsingDeltaSigma(double mu, double sigma, double x, double factor) {

        double delta = factor*sigma;
        Double d0 = Gumbel.pdf(x, mu, (sigma - delta));
        Double d1 = Gumbel.pdf(x, mu, sigma);
        Double d2 = Gumbel.pdf(x, mu, (sigma + delta));
        return estimateDerivUsingDelta(d0, d1, d2, delta);
    }

    /**
     * estimate d/dmu of GEV using the difference between GEVs given minor changes in k
     *
     * @param mu
     * @param sigma
     * @param x
     * @return
     */
    static Double estimateDerivUsingDeltaMu(double mu, double sigma, double x, double factor) {

        double deltaMu = factor*mu;
        Double d0 = Gumbel.pdf(x, (mu - deltaMu), sigma);
        Double d1 = Gumbel.pdf(x, mu, sigma);
        Double d2 = Gumbel.pdf(x, (mu + deltaMu), sigma);

        return estimateDerivUsingDelta(d0, d1, d2, deltaMu);
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma, that is &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same mu, sigma, and x.
     * 
     * @param mu
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    public static double estimateDY2DSigmaDSigma(double yConst, double mu, double sigma, double x, double dydsigma, double factor) {

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma = estimate as (dyds_2 - dyds)/ds

        double delta = sigma*factor;

        Double dyds_0 = DerivGumbel.derivWRTSigma(yConst, mu, (sigma - delta), x);
        Double dyds_2 = DerivGumbel.derivWRTSigma(yConst, mu, (sigma + delta), x);

        Double d = estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);

        return (d != null) ? d.doubleValue() : 0;
    }

    /**
     * estimate &#8706;<sup>2</sup>/&#8706;mu&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;mu&#8706;mu empirically.  
     * the method accepts dydsigma as a given to allow easier
     * reuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.

     * @param mu
     * @param sigma
     * @param x
     * @return
     */
    public static double estimateDY2DMuDMu(double yConst, double mu, double sigma, double x, double dydmu, double factor) {

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;mu

        double delta = mu*factor;

        Double dydm_0 = DerivGumbel.derivWRTMu(yConst, (mu - delta), sigma, x);
        Double dydm_2 = DerivGumbel.derivWRTMu(yConst, (mu + delta), sigma, x);

        Double d = estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);

        return (d != null) ? d.doubleValue() : 0;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma, that is &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     *
     * @param mu
     * @param sigma
     * @param x
     * @return
     */
    public static double estimateDY2DMuDSigma(double yConst, double mu, double sigma, double x, double dydmu, double factor) {

        // &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma

        double delta = sigma*factor;

        Double dydm_0 = DerivGumbel.derivWRTMu(yConst, mu, (sigma - delta), x);
        Double dydm_2 = DerivGumbel.derivWRTMu(yConst, mu, (sigma + delta), x);

        Double d = estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);

        return (d != null) ? d.doubleValue() : 0;
    }

    /**
     * 
     * estimate &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu, that is &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu empirically.  
     * the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * 
     * @param mu
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    public static double estimateDY2DSigmaDMu(double yConst, double mu, double sigma, double x, double dydsigma, double factor) {

        // &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu = estimate as (dy2_ds_dk - dyds)/dk

        double delta = (sigma*factor);

        Double dyds_0 = DerivGumbel.derivWRTSigma(yConst, (mu - delta), sigma, x);
        Double dyds_2 = DerivGumbel.derivWRTSigma(yConst, (mu + delta), sigma, x);

        Double d = estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);

        return (d != null) ? d.doubleValue() : 0;
    }

    public static Double sumSquaredDifferences(final double mu, final double sigma, final double[] x,
        final double[] normalizedY) {

        int n = x.length;

        if (normalizedY.length != n) {
            throw new IllegalArgumentException("x.length must equal normalizedY.length");
        }

        double[] yGEV = Gumbel.generateCurve(x, mu, sigma);

        double yMax = MiscMath0.findMax(yGEV);
        MatrixUtil.multiply(yGEV, 1./yMax);

        double sum = 0;
        double diff;
        for (int i = 0; i < n; ++i) {
            diff = yGEV[i] - normalizedY[i];
            sum += (diff * diff);
        }
        return sum;
    }

}
