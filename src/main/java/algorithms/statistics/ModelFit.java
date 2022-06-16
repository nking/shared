package algorithms.statistics;

public class ModelFit {
    // sInv is precision
    // s is covariance
    protected double[][] phiX;
    protected double[][] precision;
    /**
     * mean of the regression relationship, roughly similar to a slope term of yTrain/xTrain
     */
    protected double[] mean;
    protected double[][] cov;
    protected double alpha;
    protected double beta;

    public double[][] getPhiX() {
        return phiX;
    }

    public void setPhiX(double[][] phiX) {
        this.phiX = phiX;
    }

    public double[][] getPrecision() {
        return precision;
    }

    public void setPrecision(double[][] precision) {
        this.precision = precision;
    }

    public double[] getMean() {
        return mean;
    }

    public void setMean(double[] mean) {
        this.mean = mean;
    }

    public double[][] getCov() {
        return cov;
    }

    public void setCov(double[][] cov) {
        this.cov = cov;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }
}
