package algorithms.correlation;

/**
 *
 * implementation of Chaudhuri & Hu 2019.
 * 
 * TODO: consider implementing Brownian Distance Covariance.
 * TODO: add notes for Hilbert-Schmidt independence measure (HSIC) - Lasso.
 * 
 * @author nichole
 */
public class Distance {
    
    /**
     * calculates the distance covariance between univariate vectors x and y as
     * "a weighted  distance between the joint characteristic function and 
     * the product of marginal distributions; 
     * it is 0 if and only if two random vectors  and  are independent. 
     * This measure can detect the presence of a dependence structure when the 
     * sample size is large enough."
     * 
     * This algorithm is an implementation/port of the Matlab code from
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     *      * 
     * Runtime is O(n * lg_2(n)) where n is the number of points in x which is
     * the same as the number in y.
     * 
     * NOTE: redundant points are possible in the rankings as "ties" are handled
     * in the algorithm.  This is one advantage over the similar
     * algorithm of Huo and Szekely (2016).
     * 
     * @param x
     * @param y
     * @return 
     */
    public static double[][] correlation(double[] x, double[] y) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
}
