package algorithms.misc;

/**
 *
 * @author nichole
 */
public class CDFStandardNormal {
    
    /**
     * 
     * <pre>
     * Best Fit One Parameter Logistic Equation from:
     * "A logistic approximation to the cumulative normal distribution"
     * 2009, Bowling, Khasawneh, Kaewkuekool, and Cho,
     * JIEM, 2009 – 2(1): 114-127
     * https://www.jiem.org/index.php/jiem/article/viewFile/60/27/
     * <pre>
     * @param p can be .gte. -4.5 and .lte. 4.5
     * @return 
     */
    static double _approxBowling(double p) {
        // eqn 10:
        double c = -(0.07056*p*p*p + 1.5976*p); 
        double oneParamLogistic= 1./(1. + Math.exp(c));
        return oneParamLogistic;
    }
    
    
    /**
     an approximation of the CDF of the standard normal gaussian distribution
     by author Haim Short.
     The algorithm is useful for reverse look-ups, that is, given p, estimate 
     the variate X.
     The algorithm is from:
     "Simple Approximations for the Inverse Cumulative Function,
     the Density Function, and the Loss Integral of the Normal Distribution"
    1982, Shore, Appl Statis, 31, No 2, 108-114
    
     NOTE that to compare values to the NIST tables at 
     https://www.itl.nist.gov/div898/handbook/eda/section3/eda3671.htm
      for Area under the Normal Curve from 0 to X (which is p = 0.5 to 1.0),
     one must subtract 0.5 from the search p to use their table.
     example: to find p=0.75, look for value 0.25 in their table.
     the X label on the left hand side is 0.6 and the column label at the top
     is between 0.07 and 0.08, so the X is approx 0.68.
     
    The Shore algorithm in this method produces this for the same example:
    <pre>
    >>> p=0.75
    >>> x1=-5.531*( math.pow((1.-p)/p, 0.1193) - 1.); x1
    0.679421164952841
    *</pre>
     * @param p probability
     * @return the variate, quantile, x.  the range is approximately -3.1 to +3.1.
    */
    public static double approxInverseShort(double p) {
        double z = (p >= 0.5) ? ((1.-p)/p) : (p/(1.-p));
        double x = -5.531 * (Math.pow(z, 0.1193) - 1.);
        if (p <0.5) {
            x*=-1;
        }
        return x;
    }
}
