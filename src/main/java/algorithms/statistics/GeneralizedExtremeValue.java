package algorithms.statistics;

import algorithms.misc.MiscMath0;

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
  sigma is the scale parameter and is > 0
  k is the shape parameter
 
 
  if k != 0,
      1 + (k * (x-mu)/sigma) > 0
 
 
  Let z = (x-mu)/sigma
 
  Extreme Value Type I (Gumbel):
                      1
      y = y_const * ----- * exp( -z -exp(-z))
                    sigma
      
      sigma > 0
      k = 0
 
  Extreme Value Type II (Frechet):
                      k     (sigma)^(k+1)      (   (sigma)^k)
      y = y_const * ----- * (-----)       * exp( - (-----)  )
                    sigma   (  x  )            (   (  x  )  )
 
      k  > 0, (1/k)>0
      sigma > 0
      x > 0
      if (1 + k*(x - mu)/sigma) <= 0, y = 0
 
  Extreme Value Type III (Weibull):
                      k     (x - mu)^(k+1)      (   (x - mu)^k)
      y = y_const * ----- * (------)       * exp( - (------)  )
                    sigma   (sigma )            (   (sigma )  )
 
                      k
        = y_const * ----- * (z)^(k+1) * exp( -1*(z)^k )
                    sigma
 
      k < 0, (1/k)<0
      sigma > 0
      x > 0
       if -(1 + k*(x - mu)/sigma) >= 0, y = 1
 
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
     * Errors in fitting the curve can be calculated via chain rule of derivatives of the fit:
     *
     *
     * @param y1 the y data array
     * @param yGEV the generated GEV y model array
     * @return the mean error of the fit
     */
    public static double calculateChiSq(double[] y1, double[] yGEV) {

        double yNorm = MiscMath0.findMax(y1);

        double chiSum = 0.;

        for (int i = 0; i < yGEV.length; i++) {

            double z = (yGEV[i] - y1[i])/yNorm;
            z *= z;

            chiSum += z;
        }

        return chiSum;
    }

    /**
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
