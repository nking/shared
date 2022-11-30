package algorithms.statistics;

import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
  <pre>
  Generate curves following the Generalized Extreme Value probability density
  function.
 The GEV is a.k.a. the Fisher-Tippett distribution (though some restrict the
 Fisher-Tippett to the Gumbel distribution).
 
                           (   (      ( x-mu))-(1/k))
                           (-1*(1 + k*(-----))      )
                  1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
  y = y_const * ----- * exp                           * (1 + k*(-----))
                sigma                                   (      (sigma))
 
  mu is  the location parameter
  sigma is the scale parameter and is .gt. 0
  k is the shape parameter
 
 
  if k != 0,
      1 + (k * (x-mu)/sigma) .gt. 0
 
 
  Let z = (x-mu)/sigma
 
  Extreme Value Type I (Gumbel):
                      1
      y = y_const * ----- * exp( -z -exp(-z))
                    sigma
      
      sigma .gt. 0
      k = 0
 
  Extreme Value Type II (Frechet):
                      k     (sigma)^(k+1)      (   (sigma)^k)
      y = y_const * ----- * (-----)       * exp( - (-----)  )
                    sigma   (  x  )            (   (  x  )  )
 
      k .gt. 0, (1/k) .gt. 0
      sigma .gt.  0
      x .gt. 0
      if (1 + k*(x - mu)/sigma) .lte. 0, y = 0
 
  Extreme Value Type III (Weibull):
                      k     (x - mu)^(k+1)      (   (x - mu)^k)
      y = y_const * ----- * (------)       * exp( - (------)  )
                    sigma   (sigma )            (   (sigma )  )
 
                      k
        = y_const * ----- * (z)^(k+1) * exp( -1*(z)^k )
                    sigma
 
      k .lt. 0, (1/k) .lt. 0
      sigma .gt. 0
      x .gt. 0
       if -(1 + k*(x - mu)/sigma) .geq. 0, y = 1
 
 </pre>
  first implemented in project
     https://github.com/nking/two-point-correlation
     w/ Copyright Climb With Your Feet
     and using The MIT License (MIT)
   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class GeneralizedExtremeValue {

    protected final double[] x;
    protected final double[] y;

    protected final double[] dx;
    protected final double[] dy;

    /**
     * calculate a rough estimate of GEV distribution parameters for the given x.
     * Note that the methods works best for sample sizes > 50.
     * <pre>
     * references:
     *     "Estimation of the Generalized Extreme-Value Distribution by the Method of Probability Weighted Moments"
     *    Hosking, Wallis, and Wood 1984
     * </pre>
     * @param x ordered statistic of an observed GEV distribution where X is ordered from min to max.
     * @return
     */
    public static double[] fitUsingMethodOfProbabilityWeightedMoments(double[] x) {

        int n = x.length;

        /* estimate b0, b1, b2 using eqn (4) and plotting position for p_j
        p[j] = (j - 0.35)/n
        eqn (4) : b_r[p[j]] = (1/n) * sum over j=1 to n ( p[j]^r * x[j] )

        below, j is the jth ranked (from smallest to largest) datum.

        b0 = (1/n) * sum over j=1 to n ( 1 * x[j] )
        b1 = (1/n) * sum over j=1 to n ( p[j] * x[j] )
        b2 = (1/n) * sum over j=1 to n ( p[j]^2 * x[j] )

        Alternatively, Landwehr 1979 as referenced in
        IMPROVING PROBABILITY-WEIGHTED MOMENT METHODS FOR THE GENERALIZED EXTREME VALUE DISTRIBUTION
        Diebolt et al 2008
             b_r = (1/n) * sum over j=1 to n ( multiply over l=1 to r( (j-l)/(n-l) ) * X_j,n )
        */

        double[] b = new double[3];
        int i;
        //int j;
        double t;

        double pj;
        for (i = 0; i < n; ++i) {
            //j = n - i; // for the statistic where X is ordered from max to min
            // pj = (j - 0.35)/n; // for the statistic where X is ordered from max to min
            pj = (i - 0.35)/n;
            t = x[i];
            b[0] += t;
            t *= pj;
            b[1] += t;
            t *= pj;
            b[2] += t;
        }
        b[0] /= (double)n;
        b[1] /= (double)n;
        b[2] /= (double)n;

        //System.out.printf("b=%s\n", FormatArray.toString(b, "%.3e"));

        /* eqn (14) for shape estimator
        c = ((2*b1 - b0)/(3*b2-b0)) - math.log(2)/math.log(3)
        kEst = 7.859 * c + 2.9554 * c^2
        */
        double c =  ((2.*b[1] - b[0])/(3.*b[2] - b[0])) - (Math.log(2)/Math.log(3));
        double shapeEst = 7.859*c + 2.9554*c*c;

        /* eqn (15) for location and scale estimators
        scale: alphaEst = (kEst*(2*b1 - b0)/( gamma(1+kEst) * (1 - 2^-kEst)))
        location: epsEst = b0 + (alphaEst*(gamma(1+kEst) - 1)/kEst)
        */
        double gam = Gamma.lanczosGamma9(1. + shapeEst);
        double denom = gam * (1. - Math.pow(2, -shapeEst));
        double scaleEst = shapeEst*(2.*b[1] - b[0])/denom;
        double locationEst = b[0] + (scaleEst*(gam - 1.)/shapeEst);

        // TODO: consider replacing with Castillo and Hadi (1994)
        // Castillo, E., and A. Hadi. (1994).
        //     Parameter and Quantile Estimation for the Generalized Extreme-Value Distribution. Environmetrics 5, 417–432.

        // TODO: consider shape estimator of Dekkers 1989 or Resnick and Starica 1997
        //shapeEst = M1 + 1 − 0.5* (1-((M1)^2)*(M2^(−1)))^(−1),
        // where Mj ≡ (1/k)* sum_over_i=1 to k*( (ln(X_{i:n}) −ln(X_{(k+1):n} )^j, j=1;2.

        return new double[]{locationEst, scaleEst, shapeEst};
    }

    public GeneralizedExtremeValue(double[] xPoints, double[] yPoints, double[] dXPoints, double[] dYPoints) {
        this.x = xPoints;
        this.y = yPoints;
        this.dx = dXPoints;
        this.dy = dYPoints;
    }

    /**
     * generate a normalized GEV
     * @param parameters mu, sigma, and k
     where mu is  the location parameter, sigma is the scale parameter and is > 0,
     and k is the shape parameter.
     * @param yConst normalization factor
     * @return
     */
    public double[] generateNormalizedCurve(double[] parameters, double yConst) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold mu, sigma, and k");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2], yConst);
    }

    /**
     * generate a GEV curve
     * @param parameters mu, sigma, and k
    where mu is  the location parameter, sigma is the scale parameter and is > 0,
    and k is the shape parameter.
     * @return
     */
    public double[] generateNormalizedCurve(double[] parameters) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold mu, sigma, and k");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2]);
    }

    public double[] generateNormalizedCurve(double mu, double sigma, double k) {

        double[] yGEV = generateCurve(x, mu, sigma, k);

        double yConst = determineYConstant(yGEV, mu);

        return generateNormalizedCurve(mu, sigma, k, yConst);
    }

    /**
     * generate a GEV curve w/ parameters k, sigma, and mu and then multiply it
     * by yConst so that the highest peak has value yConst
     *
     * @param k shape parameter
     * @param sigma scale parameter, > 0
     * @param mu location parameter
     * @param yConst
     * @return
     */
    public double[] generateNormalizedCurve(double mu, double sigma, double k, double yConst) {

        double[] yGEV = generateCurve(x, mu, sigma, k);

        double yMax = MiscMath0.findMax(yGEV);

        for (int i = 0; i < yGEV.length; i++) {
            yGEV[i] *= yConst/yMax;
        }
        return yGEV;
    }

    /* <pre>
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
     *
     * Components needed in the derivatives:
     *
     *   Let z = (1 + k*( (x-mu)/sigma )
     *
     *   f1 = exp( -1. * ( z^(-1./k) ) )
     *
     *   f2 = z^(-1. - (1/k))
     *
     *   then y = (yconst/sigma) * f1 * f2
     *</pre>
     */
    public static Double generateYGEV(double xPoint, double mu, double sigma, double k) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Double.isInfinite(1./k)) {
            return generateYEVTypeI(xPoint, sigma, mu);
        }        
        // if k > 0, use Frechet which is Type II
        // if k < 0, use Weibull which is Type III.  Not using k < 0 in this project

        double z = 1. + k * ((xPoint - mu)/sigma);
        double a,b;

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

        if (zIsNegative) {
            double invNegZ = -1.0f*(1.0f/z);
            double neg1Pow = -1.0f; // TODO:  revisit this
            a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
            b = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)));
        } else {
            a = -1. * Math.pow(z, (-1./k));
            b = Math.pow(z, (-1. - (1./k)));
        }
        
        if (Double.isNaN(a) || Double.isNaN(b)) {
            // or return 0?
            throw new IllegalStateException("cannot be NaN");
        } else {
            return ((1./sigma) * Math.exp(a) * b);
        }
    }

    public static Double generateYEVTypeI(double xPoint, double mu, double sigma) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        double z = (xPoint - mu)/sigma;

        double a = (-1.*z - Math.exp(-1.0f*z));

        return ((1./sigma) * Math.exp(a));
    }

    public double[] generateCurve(double[] x1, double mu, double sigma, double k) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma cannot be null");
        }
        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Double.isInfinite(1./k)) {
            return generateEVTypeICurve(x1, mu, sigma);
        }
        
        double[] yGEV = new double[x1.length];

        for (int i = 0; i < x1.length; i++) {

            double z = 1. + k*((x1[i] - mu)/sigma);
            double a,b;

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
            
            if (zIsNegative) {
                double invNegZ = -1.0f*(1.0f/z);
                double neg1Pow = -1.0f; // TODO:  revisit this
                a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
                b = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)));
            } else {
                a = -1. * Math.pow(z, (-1./k));
                b = Math.pow(z, (-1. - (1./k)));
            }

            if (Double.isNaN(a) || Double.isNaN(b)) {
                yGEV[i] = 0;
            } else {
                double t = ((1./sigma) * Math.exp(a) * b);
                if (t < 0 || Double.isNaN(t)) {
                    yGEV[i] = 0;
                } else {
                    yGEV[i] = t;
                }
            }
        }

        return yGEV;
    }

    public static double[] generateEVTypeICurve(double[] x1, double mu, double sigma) {
        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }
        double[] yGEV = new double[x1.length];
        double z;
        double a;
        for (int i = 0; i < x1.length; i++) {
            z = (x1[i] - mu)/sigma;
            a = -(z + Math.exp(-1.*z));
            yGEV[i] = (1./sigma) * Math.exp(a);
        }
        return yGEV;
    }

    public double determineYConstant(double[] yGEV, double mu) {

        int index = MiscMath0.findYMaxIndex(y);
        if (index == -1) {
            return 0;
        }

        return determineYConstant(yGEV, mu, index);
    }

    public double determineYConstant(double[] yGEV, double mu, int index) {

        double yConst = y[index]/yGEV[index];

        return yConst;
    }

    /**
     * calculate sum of the square of the differences between the curves.
     *
     * @param y1 the y data array
     * @param yGEV the generated GEV y model array
     * @return the mean error of the fit
     */
    public static double sumOfSquaredDiff(double[] y1, double[] yGEV) {
        double chiSum = 0.;
        double d;
        for (int i = 0; i < yGEV.length; i++) {
            d = yGEV[i] - y1[i];
            chiSum += (d * d);
        }
        return chiSum;
    }

    /**
     * calculate sum of ( square(y1 - yGEV)/yGEV )
     *
     * @param y1 the y data array
     * @param yGEV the generated GEV y model array
     * @return the mean error of the fit
     */
    public static double chisq(double[] y1, double[] yGEV) {

        double chiSum = 0.;
        double d;
        for (int i = 0; i < yGEV.length; i++) {
            d = (y1[i] - yGEV[i]);
            chiSum += (d*d/yGEV[i]);
        }

        return chiSum;
    }

    /**
     * NOT READY FOR USE.  I dropped some of the terms needed,soneedtoreplace or redeisgn the method.
     <pre>
     return the error in fitting the GEV curve

                                 | dy_fit |^2            | dy_fit |^2            | dy_fit|^2            |dy_fit|^2
       (err_y_fit)^2 =  (err_x)^2|--------|   + (err_k)^2|--------|   + (err_s)^2|-------|   + (err_m)^2|------|
                                 |   dx   |              |   dk   |              |  ds   |              |  dm  |
     </pre>
     * @return
     */
    public static double calculateFittingErrorSquared(GEVYFit yFit, double xPoint) {

        if (yFit == null) {
            return Double.POSITIVE_INFINITY;
        }

        int yMaxIdx = MiscMath0.findYMaxIndex(yFit.getYFit());

        if (yMaxIdx == -1) {
            return Double.POSITIVE_INFINITY;
        }

        double xPeak = yFit.getX()[yMaxIdx];

        // since we don't have an error in the parameters, we'll make a rough
        //   guess with the parameters in bestFit compared to prior bestFit.
        //  the error in the parameters, k, sigma, and mu are assumed to be smaller
        //  than those deltas in order for the best fit to have a better fit.
        //  we can use the full delta to overestimate the error, or assume that
        //  the error has to be smaller than half of the delta otherwise the prior Fit
        //  value would have been kept.  very very rough approx for parameter errors...

        double kDelta = yFit.getKResolution()/2.;
        double sigmaDelta = yFit.getSigmaResolution()/2.;
        double muDelta = yFit.getMuSolutionResolution()/2.;

        double yConst = yFit.getYScale();
        double mu = yFit.getMu();
        double k = yFit.getK();
        double sigma = yFit.getSigma();

        double xDelta = yFit.getX()[1] - yFit.getX()[0];

        double dydx = DerivGEV.derivWRTX(yConst, mu, sigma, k, xPoint);

        double dydk = DerivGEV.derivWRTK(yConst, mu, sigma, k, xPoint);

        double dyds = DerivGEV.derivWRTSigma(yConst, mu, sigma, k, xPoint);

        double dydm = DerivGEV.derivWRTMu(yConst, mu, sigma, k, xPoint);

        // since x is always given as a number, exclude it from propagation
        double err = /*Math.pow(xDelta*dydx, 2)*/ + Math.pow(kDelta*dydk, 2) + Math.pow(sigmaDelta*dyds, 2) + Math.pow(muDelta*dydm, 2);

//System.out.println("error in x alone: " + Math.sqrt( Math.pow(xDelta*dydx, 2) )/yFit.getYScale());

        return err/(yFit.getYScale() * yFit.getYScale());
    }

    /**
      <pre>
      calculate the error in an area / height calculation for y=0 to y > yLimitFraction where y is
      yLimitFraction to to the right of the peak.  This is useful for determining errors for things
      like FWHM for example.
     
      For FWHM we have sum of f = sum(X_i*Y_i)_(i < yLimit)/ Y_i

          err^2 = xError^2*(Y_i/Y_i) = xError^2

          it reduces to the sum of the errors in x.  no pde's...
     </pre>
     
     * @param yFit
     * @param yMaxFactor
     * @return
     */
    public static double calculateWidthFittingError(GEVYFit yFit, double yMaxFactor) {

        if (yFit == null) {
            return Double.POSITIVE_INFINITY;
        }

        int yPeakIdx = MiscMath0.findYMaxIndex(yFit.getOriginalScaleX());
        double yLimit = yMaxFactor * yFit.getOriginalScaleYFit()[yPeakIdx];
        int yLimitIdx = -1;
        for (int i = 0; i < yFit.getOriginalScaleX().length; i++) {
            if (i > yPeakIdx) {
                if (yFit.getOriginalScaleYFit()[i] > yLimit) {
                    yLimitIdx = i;
                } else {
                    break;
                }
            } else {
                yLimitIdx = i;
            }
        }

        // xError[i] should be formally calculated, but the approximation that it can be determined no
        //   better than  the bin center +- binwidth/2  is a minimum error.
        //   a safe addition to that (added in quadrature) would be an error in x derived from chi square
        //   but that is not done here

        double xDelta = yFit.getX()[1] - yFit.getX()[0];
        double xErrorSq = (xDelta*xDelta/4.);

        double sum = 0.0f;
        for (int i = 0; i <= yLimitIdx; i++) {
            sum += xErrorSq;
        }

        sum = Math.sqrt(sum);

        return sum;
    }

    /**
     * calculate the inverse CDF of the GEV, that is, a random variate x given the
     * probability alpha.
     * <pre>
     *     reference is from the book "Extreme Value Distributions, Theory and Applications"
     *     by Kotz and Nadarajah 2000, Section 2.2 Macleod (1989):
     *     when shape != 0:
     *         x = location + scale*(1 - exp(-shape * y))/shape
     *     when shape= 0
     *         x = location + scale * y
     *     where y = -log(-log(p));  0 < p < 1
     * </pre>
     * @param p random variate drawn from U(0,1) where U is the uniform distribution.
     *          NOTE that p must be greater than 0 and less 1, so consider
     *          using a number such as eps=1-11 up to 1e-16 for an offset from 0 or 1.
     *          In other words, p drawn from U(1e-16, 1. - 1e-16).
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @param shape parameter of the distribution function
     * @return
     */
    public static double inverseCdf(double p, double location, double scale, double shape) throws NoSuchAlgorithmException {

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        return inverseCdf(p, location, scale, scale, rand);
    }

    /**
     * calculate the inverse CDF of the GEV, that is, a random variate x given the
     * probability alpha.
     * <pre>
     *     reference is from the book "Extreme Value Distributions, Theory and Applications"
     *     by Kotz and Nadarajah 2000, Section 2.2 Macleod (1989):
     *     when shape != 0:
     *         x = location + scale*(1 - exp(-shape * y))/shape
     *     when shape= 0
     *         x = location + scale * y
     *     where y = -log(-log(p));  0 < p < 1
     * </pre>
     * @param p random variate drawn from U(0,1) where U is the uniform distribution.
     *    NOTE that p must be greater than 0 and less 1, so consider
     *    using a number such as eps=1-11 up to 1e-16 for an offset from 0 or 1.
     *    In other words, p drawn from U(1e-16, 1. - 1e-16).
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @param shape parameter of the distribution function
     * @return
     */
    public static double inverseCdf(double p, double location, double scale, double shape, SecureRandom rand) {
        if (p <= 0) {
            throw new IllegalArgumentException("p must be > 0 and < 1");
        }
        if (p >= 1) {
            throw new IllegalArgumentException("p must be > 0 and < 1");
        }
        double y = -Math.log(-Math.log(p));
        double x;
        if (shape != 0) {
            x = location + scale * (1. - Math.exp(-shape * y)) / shape;
        } else {
            x = location + scale * y;
        }
        return x;
    }

    /**
     * calculate the inverse CDF of the GEV, that is, a random variate x given the
     * probability alpha.
     * <pre>
     *     reference is from the book "Extreme Value Distributions, Theory and Applications"
     *     by Kotz and Nadarajah 2000, Section 2.2 Macleod (1989):
     *     when shape != 0:
     *         x = location + scale*(1 - exp(-shape * y))/shape
     *     when shape= 0
     *         x = location + scale * y
     *     where y = -log(-log(p));  0 < p < 1
     * </pre>
     * @param p array of random variates drawn from U(0,1) where U is the uniform distribution.
     *    NOTE that p elements must be greater than 0 and less 1, so consider
     *    using a number such as eps=1-11 up to 1e-16 for an offset from 0 or 1.
     *    In other words, p elements drawn from U(1e-16, 1. - 1e-16).
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @param shape parameter of the distribution function
     * @return
     */
    public static double[] inverseCdf(double[] p, double location, double scale, double shape, SecureRandom rand) {
        int n = p.length;
        int i;
        double y;
        double[] x = new double[n];
        if (shape != 0) {
            for (i = 0; i < n; ++i) {
                y = -Math.log(-Math.log(p[i]));
                x[i] = location + scale * (1. - Math.exp(-shape * y)) / shape;
            }
        } else {
            for (i = 0; i < n; ++i) {
                y = -Math.log(-Math.log(p[i]));
                x[i] = location + scale * y;
            }
        }
        return x;
    }

    /**
     * randomly draw from the inverse CDF of the GEV to return points that follow a GEV distibution.
     * <pre>
     *     reference is from the book "Extreme Value Distributions, Theory and Applications"
     *     by Kotz and Nadarajah 2000, Section 2.2 Macleod (1989):
     *     when shape != 0:
     *         x = location + scale*(1 - exp(-shape * y))/shape
     *     when shape= 0
     *         x = location + scale * y
     *     where y = -log(-log(p));  0 < p < 1
     * </pre>
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @param shape parameter of the distribution function
     * @return array of randomly drawn numbers from the GEV
     */
    public static double[] sampleRandomlyFrom(double location, double scale, double shape,
                                              int nDraws) throws NoSuchAlgorithmException {

        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        return sampleRandomlyFrom(location, scale, shape, nDraws, rand);
    }

    /**
     * randomly draw from the inverse CDF of the GEV to return points that follow a GEV distibution.
     * <pre>
     *     reference is from the book "Extreme Value Distributions, Theory and Applications"
     *     by Kotz and Nadarajah 2000, Section 2.2 Macleod (1989):
     *     when shape != 0:
     *         x = location + scale*(1 - exp(-shape * y))/shape
     *     when shape= 0
     *         x = location + scale * y
     *     where y = -log(-log(p));  0 < p < 1
     * </pre>
     * @param location parameter of the distribution function
     * @param scale parameter of the distribution function
     * @param shape parameter of the distribution function
     * @return array of randomly drawn numbers from the GEV
     */
    public static double[] sampleRandomlyFrom(double location, double scale, double shape,
                                              int nDraws, SecureRandom rand) {
        double eps = 1e-16;
        double[] p = new double[nDraws];
        for (int i = 0; i < nDraws; ++i) {
            p[i] = rand.nextDouble();
            if (p[i] == 0) {
                p[i] = eps;
            } else if (p[i] == 1.) {
                p[i] -= eps;
            }
        }

        return inverseCdf(p, location, scale, shape, rand);
    }

    public static double[] generateNormalizedCurve(double[] x1, double mu, double sigma, double k) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        double[] yGEV = genCurve(x1, mu, sigma, k);

        double yMax = MiscMath0.findMax(yGEV);

        for (int i = 0; i < yGEV.length; i++) {
            yGEV[i] /= yMax;
        }

        return yGEV;
    }

    public static double[] _genCurve(double[] x1, double mu, double sigma, double k) {
        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }
        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Double.isInfinite(1./k)) {
            return generateEVTypeICurve(x1, mu, sigma);
        }
        double t;
        double[] yGEV = new double[x1.length];
        for (int i = 0; i < x1.length; i++) {
            t = 1. + k * ((x1[i] - mu) / sigma);
            t = Math.pow(t, -1./k);
            yGEV[i] = (1./sigma) * Math.pow(t, k + 1.) * Math.exp(-t);
        }
        return yGEV;
    }

    public static double[] genCurve(double[] x1, double mu, double sigma, double k) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }
        // if k =~ 0, use Gumbel which is Type I
        if (k == 0 || Double.isInfinite(1./k)) {
            return generateEVTypeICurve(x1, mu, sigma);
        }

        double[] yGEV = new double[x1.length];

        for (int i = 0; i < x1.length; i++) {

            double z = 1. + k * ((x1[i] - mu)/sigma);
            double a,b;

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

            if (zIsNegative) {
                double invNegZ = -1.0f*(1.0f/z);
                double neg1Pow = -1.0f; // TODO:  revisit this
                a = -1. * neg1Pow * Math.pow(invNegZ, (-1./k));
                b = neg1Pow * Math.pow(invNegZ, (-1. - (1./k)));
            } else {
                a = -1. * Math.pow(z, (-1./k));
                b = Math.pow(z, (-1. - (1./k)));
            }

            if (Double.isNaN(a) || Double.isNaN(b)) {
                yGEV[i] = 0;
            } else {
                double t = ((1./sigma) * Math.exp(a) * b);
                if (t < 0 || Double.isNaN(t)) {
                    yGEV[i] = 0;
                } else {
                    yGEV[i] = t;
                }
            }
        }

        return yGEV;
    }
}
