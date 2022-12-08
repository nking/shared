package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class ModelFit {
    // precision matrix is inverse of the covariance matrix
    // sInv is precision
    // s is covariance
    /**
     * the feature matrix, a.k.a. the design matrix generated from the training vector x.
     */
    protected double[][] phiX;
    /**
     * inverse of the covariance matrix.  this is S^-1 in Bishop's PRML chapters 1 and 3.
     */
    protected double[][] precision;
    /**
     * mean of the regression relationship, roughly similar to a slope term of yTrain/xTrain
     */
    protected double[] mean;
    /**
     * a measure of the joint variablity of pairs of two random elements from vector x.
     * this is S in Bishop's PRML, chapters 1, and 3.
     */
    protected double[][] cov;
    /**
     * alpha is the precision of the prior. the precision is the inverse of the variance.
     */
    protected double alpha;
    /**
     * beta is the precision of the likelihood.  the precision is the inverse of the variance
     */
    protected double beta;

    /**
     *
     @return
     */
    public double[][] getPhiX() {
        return phiX;
    }

    /**
     *
     @param phiX
     */
    public void setPhiX(double[][] phiX) {
        this.phiX = phiX;
    }

    /**
     *
     @return
     */
    public double[][] getPrecision() {
        return precision;
    }

    /**
     *
     @param precision
     */
    public void setPrecision(double[][] precision) {
        this.precision = precision;
    }

    /**
     *
     @return
     */
    public double[] getMean() {
        return mean;
    }

    /**
     *
     @param mean
     */
    public void setMean(double[] mean) {
        this.mean = mean;
    }

    /**
     *
     @return
     */
    public double[][] getCov() {
        return cov;
    }

    /**
     *
     @param cov
     */
    public void setCov(double[][] cov) {
        this.cov = cov;
    }

    /**
     *
     @return
     */
    public double getAlpha() {
        return alpha;
    }

    /**
     *
     @param alpha
     */
    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    /**
     *
     @return
     */
    public double getBeta() {
        return beta;
    }

    /**
     *
     @param beta
     */
    public void setBeta(double beta) {
        this.beta = beta;
    }
}
